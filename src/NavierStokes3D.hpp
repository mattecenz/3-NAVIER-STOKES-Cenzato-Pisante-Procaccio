#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Stokes problem.
class NavierStokes
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    ForcingTerm()
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim - 1; ++i)
        values[i] = 0.0;

      values[dim - 1] = -g;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      if (component == dim - 1)
        return -g;
      else
        return 0.0;
    }

  protected:
    const double g = 0.0;
  };

  // Dirichlet boundary conditions.
  class FunctionG : public Function<dim>
  {
  public:
    // Constructor.
    FunctionG() : Function<dim>(dim + 1)
    {
    }

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/, const unsigned int component = 0) const override
    {
      if (component == 3)
        return 0.;
      else
        return 0.;
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 0.;
      values[1] = 0.;
      values[2] = 0.;
      values[3] = 0.;
    }
  };

  // Neumann boundary conditions.
  class FunctionH : public Function<dim>
  {
  public:
    // Constructor.
    FunctionH()
    {
    }

    // Evaluation.
    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/ = 0) const override
    {
      return 0.;
    }
    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 0.;
      values[1] = 0.;
      values[2] = 0.;
    }
  };

  // Function for the initial condition.
  class FunctionU0 : public Function<dim>
  {
  public:
    FunctionU0()
    {
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int component = 0) const
    {
      if (component == 0)
        return 0.05;
      else
        return 0.;
    }
    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 0.05;
      values[1] = 0.;
      values[2] = 0.;
    }
  };

  // Function for inlet velocity. This actually returns an object with four
  // components (one for each velocity component, and one for the pressure), but
  // then only the first three are really used (see the component mask when
  // applying boundary conditions at the end of assembly). If we only return
  // three components, however, we may get an error message due to this function
  // being incompatible with the finite element space.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
        : Function<dim>(dim + 1)
    {
    }

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 1.;

      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int component = 0) const override
    {
      if (component == 0)
        return 1.;
      else
        return 0.0;
    }

  protected:
    const double alpha = 1.0;
  };

  // Since we're working with block matrices, we need to make our own
  // preconditioner class. A preconditioner class can be any class that exposes
  // a vmult method that applies the inverse of the preconditioner.

  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::Vector &dst,
          const TrilinosWrappers::MPI::Vector &src) const
    {
      dst = src;
    }

  protected:
  };
  class PreconditionBlockIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };

  // Block-diagonal preconditioner.
  class PreconditionBlockDiagonal
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
          solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      SolverControl solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);
    }

    

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;
  };

  // Block-triangular preconditioner.
  class PreconditionBlockTriangular
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass = &pressure_mass_;
      B = &B_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
          solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      tmp.reinit(src.block(1));
      B->vmult(tmp, dst.block(0));
      tmp.sadd(-1.0, src.block(1));

      SolverControl solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
          solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               tmp,
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp;
  };

  class PreconditionSIMPLE
  {
  public:

  void exportmatrix(TrilinosWrappers::SparseMatrix *A, std::string outputFileName)
{
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open())
    {
        int rows = A->m(), cols = A->n();
        // Write dimensions to the first row
        outputFile << rows << " " << cols << std::endl;

        // Write matrix data
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                outputFile << std::setw(8) << std::fixed << std::setprecision(4) << (*A)(i, j) << " ";
            }
            outputFile << std::endl;
        }
        std::cout << "Computed matrix has been written to " << outputFileName << std::endl;

        // Close the file
        outputFile.close();
    }
    else
    {
        std::cerr << "Error opening file for writing." << std::endl;
    }
}

     void exportvector(TrilinosWrappers::MPI::Vector A, std::string outputFileName)
{
    // Write the matrix to the file
    std::ofstream outputFile(outputFileName);
    if (outputFile.is_open())
    {
        int rows = A.size();
        // Write dimensions to the first row
        outputFile << rows << " " << 1 << std::endl;

        // Write matrix data
        for (int i = 0; i < rows; ++i)
        {
            outputFile << std::setw(8) << std::fixed << std::setprecision(4) << A(i) << " "; 
            outputFile << std::endl;
        }
        std::cout << "Computed matrix has been written to " << outputFileName << std::endl;

        // Close the file
        outputFile.close();
    }
    else
    {
        std::cerr << "Error opening file for writing." << std::endl;
    }
}
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               TrilinosWrappers::SparseMatrix &B_)
    {
      B = &B_;
      F = &F_;
      preconditionerF.initialize(F_); 
      TrilinosWrappers::MPI::Vector vett;
      TrilinosWrappers::MPI::Vector tmp;
      for(unsigned int i=0;i<F->m();i++)
        vett(i)= 1.0;
      preconditionerF.vmult(tmp,vett);
      exportvector(tmp,"vett.txt");


      FullMatrix<double> D(F->m(), F->n());
      D.operator=(0);
      for (size_t i = 0; i < F->m(); i++)
      {
        D.set(i, i, - 1 / F->diag_element(i));
      }

      FullMatrix<double> M(B->m(), D.n());
      FullMatrix<double> M2(B->m(), B->n());
      FullMatrix<double> B_full(B->m(), B->n());
      FullMatrix<double> B_fullt(B->n(), B->m());
      B_full.copy_from(*B);

      std::cout << "prima\n";
      B_full.mmult(M, D);
      B_fullt.copy_transposed(B_full);
      M.mmult(M2, B_fullt);

      SparsityPattern sp;
      sp.copy_from(M2);
      S->reinit(sp);
      // Copy values from FullMatrix to SparseMatrix
      for (unsigned int i = 0; i < M2.m(); ++i)
      {
        for (unsigned int j = 0; j < M2.n(); ++j)
        {
          if (M2(i, j) != 0.0)
          {
            S->set(i, j, M2(i, j));            
          }       
        }
        std::cout<<std::endl;
      }

      S->compress(VectorOperation::add);

      exportmatrix(S,"output_S.txt");
      exportmatrix(B,"output_B.txt");
      
      preconditionerS.initialize(*S);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SparsityPattern sp(F->m(), F->n(), 1);
      sp.compress();
      TrilinosWrappers::SparseMatrix D;
      D.reinit(sp);
      for (size_t i = 0; i < F->m(); i++)
      {
        D.set(i, i, - F->diag_element(i));
      }

      SolverControl solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres_velocity(
          solver_control_velocity);

      solver_gmres_velocity.solve(*F,
                                  dst.block(0),
                                  src.block(0),
                                  PreconditionIdentity());
    

      B->vmult(dst.block(1), dst.block(0));
      dst.block(1).sadd(-1.0, src.block(1));
      std::cout << "2\n";

      tmp.reinit(src.block(1));
      tmp2.reinit(src.block(0));
      SolverControl solver_control2(1000,
                                    1e-2 * dst.block(1).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres2(
          solver_control2);
      solver_gmres2.solve(*S,
                          tmp,
                          dst.block(1),
                          PreconditionIdentity());
      std::cout << "3\n";

      D.vmult(tmp, dst.block(0));
      dst.block(1).sadd(1 / alpha, dst.block(1));
      std::cout << "4\n";

      B->transpose();
      B->vmult(tmp2, dst.block(1));
      tmp2.sadd(-1, tmp);
      B->transpose();
      std::cout << "5\n";

      solver_gmres_velocity.solve(D,
                                  dst.block(0),
                                  tmp2,
                                  PreconditionIdentity());
      std::cout << "fine\n";
    }

  protected:
    TrilinosWrappers::SparseMatrix *B;
    const TrilinosWrappers::SparseMatrix *F;
    TrilinosWrappers::SparseMatrix *S;

    TrilinosWrappers::PreconditionILU preconditionerF;
    TrilinosWrappers::PreconditionILU preconditionerS;

    mutable TrilinosWrappers::MPI::Vector tmp;
    mutable TrilinosWrappers::MPI::Vector tmp2;
    const double alpha = 0.5;
  };

	class MyPreconditionSIMPLE
	{
		public:
			void
    	initialize(const TrilinosWrappers::SparseMatrix &F_,
     	          const TrilinosWrappers::SparseMatrix &B_,
								const TrilinosWrappers::SparseMatrix &B_t)
    	{
				F=&F_;
				B=&B_;
				B_T=&B_t;
				//Construct the inverse of the diagonal
				TrilinosWrappers::MPI::Vector diagonal;
				diagonal=0.0;
				for(size_t i=0;i<F_.m();++i){
					diagonal[i]=1./F_.diag_element(i);
					//Save it also in the sparse matrix
					D_inv.set(i,i,1./F_.diag_element(i));
				}
				
				//Create S_tilde
				B_.mmult(S_tilde, B_t, diagonal);

			}
    	void
    	vmult(TrilinosWrappers::MPI::BlockVector &dst,
          	const TrilinosWrappers::MPI::BlockVector &src) const
    	{
				const unsigned int maxiter=100000;
				const double tol=1e-2*src.l2_norm();
      	SolverControl solver_control(maxiter,tol);
      	
				SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
      	
				//Store in temporaries the results
				TrilinosWrappers::MPI::Vector y_u=src.block(0);
				TrilinosWrappers::MPI::Vector y_p=src.block(1);

				TrilinosWrappers::MPI::Vector temp_1=src.block(1);
				
				std::cout<<"Helo"<<std::endl;

				solver.solve(*F, y_u, src.block(0), PreconditionIdentity());
				
				std::cout<<"Helo 1 "<<solver_control.last_step()<<std::endl;

				B->vmult(temp_1,y_u);
				temp_1-=src.block(1);

				std::cout<<"Helo 2"<<std::endl;
				
				solver.solve(S_tilde, y_p, temp_1, PreconditionIdentity()); 
				
				std::cout<<"Helo 3 "<<solver_control.last_step()<<std::endl;

				dst.block(1)=y_p;
				dst.block(1)*=1./alpha;
				temp_1=dst.block(1);

				B_T->vmult(temp_1,dst.block(1));
				//Cannot be same vector
				D_inv.vmult(dst.block(0),temp_1);
				dst.block(0)-=y_u;
				dst.block(0)*=-1.;

				std::cout<<"Helo 4"<<std::endl;
			
			}
		
		protected:

			const double alpha=0.5;

			const TrilinosWrappers::SparseMatrix* F;
			const TrilinosWrappers::SparseMatrix* B_T;
			const TrilinosWrappers::SparseMatrix* B;
			TrilinosWrappers::SparseMatrix S_tilde;
			TrilinosWrappers::SparseMatrix D_inv;

	};

  // Constructor.
  NavierStokes(const std::string &mesh_file_name_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_,
               const double &T_,
               const double &deltat_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), mesh_file_name(mesh_file_name_), degree_velocity(degree_velocity_), degree_pressure(degree_pressure_), T(T_), deltat(deltat_), mesh(MPI_COMM_WORLD)
  {}

  // Setup system.
  void
  setup();

  // Solve system.
  void
  solve();

  // Compute lift and drag.
  void 
  compute_forces();

protected:

  // Drag and lift.
  double drag;
  double lift;

  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble();

  // Solve the problem for one time step.
  void
  solve_time_step();

  // Output results.
  void
  output(const unsigned int &time_step) const;

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu = 1;

  const double rho = 1.;

  // Outlet pressure [Pa].
  const double p_out = 10;

  // Forcing term.
  ForcingTerm forcing_term;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;
  // Final time.
  const double T;
  // TIme step.
  const double deltat;

  // g(x).
  FunctionG function_g;

  // h(x).
  FunctionH function_h;

  // Initial condition.
  FunctionU0 u_0;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;
};

#endif
