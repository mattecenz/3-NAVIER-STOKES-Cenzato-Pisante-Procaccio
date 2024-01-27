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
  static constexpr unsigned int dim = 2;

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

    virtual void
    vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
    {
      values[0] = 0.;
      values[1] = 0.;
      values[2] = 0.;
    }

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/) const override
    {
      return 0.;
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

    virtual double
    value(const Point<dim> & /*p*/, const unsigned int /*component*/) const override
    {
      return 0.;
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
    value(const Point<dim> & /*p*/, const unsigned int component) const
    {
      if (component == 0)
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
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      values[0] = 4 * u_m * p[1] * (0.41 - p[1]) / (0.41 * 0.41) * std::sin(M_PI * get_time());
      values[1] = 0.;
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 0)
        return 4 * u_m * p[1] * (0.41 - p[1]) / (0.41 * 0.41) * std::sin(M_PI * get_time());
      else
        return 0;
    }

  protected:
    double u_m = 1.5;
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


<<<<<<< HEAD
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
=======
  class MyPreconditionSIMPLE
>>>>>>> 56ba5ab9a59d2798f07d08c2839ab57130f30c09
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &B_t)
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;

      // Construct the inverse of the diagonal
      TrilinosWrappers::MPI::Vector diagonal;
      diagonal = 0.0;
      for (size_t i = 0; i < F_.m(); ++i)
      {
        diagonal[i] = 1. / F_.diag_element(i);
        // Save it also in the sparse matrix
        D_inv.set(i, i, 1. / F_.diag_element(i));
      }

      // Create S_tilde
      B_.mmult(S_tilde, B_t, diagonal);

      // Initialize the preconditioners
      preconditioner_F.initialize(*F);
      preconditioner_S.initialize(S_tilde);
    }
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      const unsigned int maxiter = 10000;
      const double tol = 1e-2;
      SolverControl solver_F(maxiter, tol * src.block(0).l2_norm());

      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_gmres(solver_F);

      // Store in temporaries the results
      TrilinosWrappers::MPI::Vector y_u = src.block(0);
      TrilinosWrappers::MPI::Vector y_p = src.block(1);

      TrilinosWrappers::MPI::Vector temp_1 = src.block(1);

      solver_gmres.solve(*F, y_u, src.block(0), preconditioner_F);

      B->vmult(temp_1, y_u);
      temp_1 -= src.block(1);

      SolverControl solver_S(maxiter, tol * temp_1.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector> solver_cg(solver_S);
      solver_cg.solve(S_tilde, y_p, temp_1, preconditioner_S);

      dst.block(1) = y_p;
      dst.block(1) *= 1. / alpha;
      temp_1.reinit(dst.block(0));

      B_T->vmult(temp_1, dst.block(1));
      // Cannot be same vector
      D_inv.vmult(dst.block(0), temp_1);
      dst.block(0) -= y_u;
      dst.block(0) *= -1.;
    }

  protected:
    const double alpha = 0.5;

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix S_tilde;
    TrilinosWrappers::SparseMatrix D_inv;
    TrilinosWrappers::PreconditionILU preconditioner_F;
    TrilinosWrappers::PreconditionILU preconditioner_S;
  };

<<<<<<< HEAD
  class PreconditionaSIMPLE
=======
  class PreconditionASIMPLE
>>>>>>> 56ba5ab9a59d2798f07d08c2839ab57130f30c09
  {
  public:
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &B_,
<<<<<<< HEAD
               const TrilinosWrappers::SparseMatrix &B_t,
               const TrilinosWrappers::MPI::BlockVector &sol_owned)
=======
               const TrilinosWrappers::SparseMatrix &B_t)
>>>>>>> 56ba5ab9a59d2798f07d08c2839ab57130f30c09
    {
      F = &F_;
      B = &B_;
      B_T = &B_t;

<<<<<<< HEAD
      diag_D_inv.reinit(sol_owned.block(0));
      diag_D.reinit(sol_owned.block(0));

      for (unsigned int i : diag_D.locally_owned_elements())
      {
        double temp = F->diag_element(i);
        diag_D[i] = -temp;
        diag_D_inv[i] = 1.0 / temp;
      }

      B->mmult(S, *B_T, diag_D_inv);

      preconditionerF.initialize(*F);
      preconditionerS.initialize(S);
    }
=======
      // Construct the inverse of the diagonal
      TrilinosWrappers::MPI::Vector diagonal;
      diagonal = 0.0;
      for (size_t i = 0; i < F_.m(); ++i)
      {
        diagonal[i] = 1. / F_.diag_element(i);
        // Save it also in the sparse matrix
        D.set(i, i, F_.diag_element(i));
      }

      // Create S_tilde
      B_.mmult(S_tilde, B_t, diagonal);

      preconditionerS.initialize(S_tilde);
      preconditionerF.initialize(F_);
    }

>>>>>>> 56ba5ab9a59d2798f07d08c2839ab57130f30c09
    void
    vmult(TrilinosWrappers::MPI::BlockVector &dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
<<<<<<< HEAD

      tmp.reinit(src.block(1));
      preconditionerF.vmult(dst.block(0), src.block(0));
      dst.block(1) = src.block(1);
      B->vmult(dst.block(1), dst.block(0));
      dst.block(1).sadd(-1.0, src.block(1));
      tmp = dst.block(1);
      preconditionerS.vmult(dst.block(1), tmp);
      dst.block(0).scale(diag_D);
      dst.block(1) *= 1.0 / alpha;
      B_T->vmult_add(dst.block(0), dst.block(1));
      dst.block(0).scale(diag_D_inv);
    }

  protected:
    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix S;

    TrilinosWrappers::PreconditionILU preconditionerF;
    TrilinosWrappers::PreconditionILU preconditionerS;

    TrilinosWrappers::MPI::Vector diag_D;
    TrilinosWrappers::MPI::Vector diag_D_inv;
    mutable TrilinosWrappers::MPI::Vector tmp;
    mutable TrilinosWrappers::MPI::Vector tmp2;
    const double alpha = 0.5;
=======
      const unsigned int maxiter = 100000;
      const double tol = 1e-6;
      SolverControl solver_control(maxiter, tol);

      SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

      // Store temporary results
      TrilinosWrappers::MPI::Vector y_u = src.block(0);
      TrilinosWrappers::MPI::Vector y_p = src.block(1);

      TrilinosWrappers::MPI::Vector temp_2 = src.block(0);

      // 1-st step
      solver.solve(*F, y_u, src.block(0), preconditionerF);

      // 2-nd step
      B->vmult(y_p, y_u);
      y_p.sadd(-1,src.block(1));

      // 3-rd step
      solver.solve(S_tilde, dst.block(1), y_p, preconditionerS);

      // 4-th step
      D.vmult(temp_2, y_u);
      dst.block(1) *= -1./alpha;

      // 5-th step
      B_T->vmult(y_u, dst.block(1));
      temp_2 -= y_u;

      // 6-th step
      solver.solve(D, dst.block(0), temp_2, PreconditionIdentity());

    }

  protected:
    const double alpha = 1.;

    const TrilinosWrappers::SparseMatrix *F;
    const TrilinosWrappers::SparseMatrix *B_T;
    const TrilinosWrappers::SparseMatrix *B;
    TrilinosWrappers::SparseMatrix S_tilde;
    TrilinosWrappers::SparseMatrix D;

    TrilinosWrappers::PreconditionILU preconditionerF;
    TrilinosWrappers::PreconditionILU preconditionerS;
>>>>>>> 56ba5ab9a59d2798f07d08c2839ab57130f30c09
  };

  // Constructor.
  NavierStokes(const std::string &mesh_file_name_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_,
               const double &T_,
               const double &deltat_)
      : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)), mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)), pcout(std::cout, mpi_rank == 0), T(T_), mesh_file_name(mesh_file_name_), degree_velocity(degree_velocity_), degree_pressure(degree_pressure_), deltat(deltat_), mesh(MPI_COMM_WORLD)
  {
  }

  // Setup system.
  void
  setup();

  // Solve system.
  void
  solve();

protected:
  // Assemble system. We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble(const double &time);

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
  const double nu = 1e-3;

  // Forcing term.
  ForcingTerm forcing_term;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh file name.
  const std::string mesh_file_name;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

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

  // Quadrature formula used on boundary lines.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_boundary;

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
