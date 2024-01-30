#include "NavierStokes2D.hpp"

void NavierStokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_boundary = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per boundary cell = " << quadrature_boundary->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0] = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1] = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::none;
        else // other combinations
          coupling[c][d] = DoFTools::always;
      }
    }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
    {
      for (unsigned int d = 0; d < dim + 1; ++d)
      {
        if (c == dim && d == dim) // pressure-pressure term
          coupling[c][d] = DoFTools::always;
        else // other combinations
          coupling[c][d] = DoFTools::none;
      }
    }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
        block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    pcout << "  Initializing the matrices" << std::endl;
    system_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);

    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);
  }
}

// https://www.dealii.org/current/doxygen/deal.II/code_gallery_time_dependent_navier_stokes.html
void NavierStokes::assemble(const double &time)
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_boundary = quadrature_boundary->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_boundary_values(*fe,
                                       *quadrature_boundary,
                                       update_values | update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs = 0.0;
  pressure_mass = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Store the current velocity value in a tensor
  std::vector<Tensor<1, dim>> current_velocity_values(n_q);
  // std::vector<Tensor<2,dim>> current_velocity_gradients(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    cell_matrix = 0.0;
    cell_rhs = 0.0;
    cell_pressure_mass_matrix = 0.0;

    // Retrieve the current solution values.
    fe_values[velocity].get_function_values(solution, current_velocity_values);

    for (unsigned int q = 0; q < n_q; ++q)
    {
      Vector<double> forcing_term_loc(dim);
      forcing_term.vector_value(fe_values.quadrature_point(q),
                                forcing_term_loc);
      Tensor<1, dim> forcing_term_tensor;
      for (unsigned int d = 0; d < dim; ++d)
        forcing_term_tensor[d] = forcing_term_loc[d];

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          // Viscosity term.
          cell_matrix(i, j) +=
              nu *
              scalar_product(fe_values[velocity].gradient(i, q),
                             fe_values[velocity].gradient(j, q)) *
              fe_values.JxW(q);

          // Time derivative discretization.
          cell_matrix(i, j) += fe_values[velocity].value(i, q) *
                               fe_values[velocity].value(j, q) /
                               deltat * fe_values.JxW(q);

          // Convective term.
          cell_matrix(i, j) += current_velocity_values[q] *
                               fe_values[velocity].gradient(j, q) *
                               fe_values[velocity].value(i, q) *
                               fe_values.JxW(q);

          // Pressure term in the momentum equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                               fe_values[pressure].value(j, q) *
                               fe_values.JxW(q);

          // Pressure term in the continuity equation.
          cell_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                               fe_values[pressure].value(i, q) *
                               fe_values.JxW(q);

          // Pressure mass matrix.
          cell_pressure_mass_matrix(i, j) +=
              fe_values[pressure].value(i, q) *
              fe_values[pressure].value(j, q) / nu * fe_values.JxW(q);
        }

        // Forcing term.
        cell_rhs(i) += scalar_product(forcing_term_tensor,
                                      fe_values[velocity].value(i, q)) *
                       fe_values.JxW(q);

        // Time derivative discretization on the right hand side
        cell_rhs(i) += scalar_product(current_velocity_values[q],
                                      fe_values[velocity].value(i, q)) /
                       deltat * fe_values.JxW(q);
      }
    }

    // Boundary integral for Neumann BCs.
    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        // 1 is the inlet velocity and 3 is the outlet
        if (cell->face(f)->at_boundary() &&
            (cell->face(f)->boundary_id() != 1 && cell->face(f)->boundary_id() != 3))
        {
          fe_boundary_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_boundary; ++q)
          {
            Vector<double> neumann_loc(dim);
            function_h.vector_value(fe_boundary_values.quadrature_point(q),
                                    neumann_loc);
            Tensor<1, dim> neumann_loc_tensor;
            for (unsigned int d = 0; d < dim; ++d)
              neumann_loc_tensor[d] = neumann_loc[d];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              cell_rhs(i) +=
                  scalar_product(neumann_loc_tensor,
                                 fe_boundary_values[velocity].value(i, q)) *
                  fe_boundary_values.JxW(q);
            }
          }
        }
      }
    }

    cell->get_dof_indices(dof_indices);

    system_matrix.add(dof_indices, cell_matrix);
    system_rhs.add(dof_indices, cell_rhs);
    pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
  }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;

    // We interpolate first the inlet velocity condition alone, then the wall
    // condition alone, so that the latter "win" over the former where the two
    // boundaries touch.
    inlet_velocity.set_time(time);
    boundary_functions[1] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, false}));

    boundary_functions.clear();
    for (unsigned int i = 1; i < 7; ++i)
    {
      if (i != 1 && i != 3)
      {
        boundary_functions[i] = &function_g;
      }
    }
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                                 {true, true, false}));

    MatrixTools::apply_boundary_values(
        boundary_values, system_matrix, solution, system_rhs, false);
  }
  // pcout<<system_matrix<<std::endl;
}
void NavierStokes::solve_time_step()
{
  pcout << "===============================================" << std::endl;

  const unsigned int maxiter = 100000;
  const double tol = 1e-4 /**system_rhs.l2_norm()*/;

  SolverControl solver_control(maxiter, tol, true);
  // solver_control.enable_history_data();

  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  // PreconditionBlockIdentity preconditioner;

  // PreconditionSIMPLE preconditioner;

  PreconditionaSIMPLE preconditioner;

  pcout << " Assemblying the preconditioner... " << std::endl;

  dealii::Timer timerprec;
  timerprec.restart();

  preconditioner.initialize(
      system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1), solution_owned);

  /*preconditioner.initialize(
      system_matrix.block(0, 0), system_matrix.block(1, 0), system_matrix.block(0, 1));*/

  timerprec.stop();
  pcout << "Time taken to initialize preconditioner: " << timerprec.wall_time() << " seconds" << std::endl;

  time_prec.push_back(timerprec.wall_time());

  // pcout << "done" << std::endl;
  pcout << "===============================================" << std::endl;

  pcout << "Solving the linear system with expected maxiter: " << maxiter;
  pcout << " and tollerance: " << tol << std::endl;

  dealii::Timer timersys;
  timersys.restart();

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);

  timersys.stop();
  pcout << "Time taken to solve Navier Stokes problem: " << timersys.wall_time() << " seconds" << std::endl;

  time_solve.push_back(timersys.wall_time());

  pcout << "Result:  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

void NavierStokes::output(const unsigned int &time_step) const
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> names = {"velocity",
                                    "velocity",
                                    "pressure"};

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::string output_file_name = "output-stokes-2D";
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      time_step,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}

//
void NavierStokes::solve()
{

  pcout << "===============================================" << std::endl;

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    VectorTools::interpolate(dof_handler, u_0, solution_owned);
    solution = solution_owned;

    // Output the initial solution.
    output(0);
    pcout << "===============================================" << std::endl;
  }

  unsigned int time_step = 0;
  double time = 0;

  while (time < T)
  {
    time += deltat;
    ++time_step;

    pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
          << time << ":" << std::flush;

    assemble(time);
    solve_time_step();

    compute_forces();
    output(time_step);
  }
}

void NavierStokes::compute_forces()
{
  pcout << "===============================================" << std::endl;
  pcout << "Computing forces: " << std::endl;

  const unsigned int n_q = quadrature->size();
  const unsigned int n_q_face = quadrature_boundary->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_boundary,
                                   update_values | update_normal_vectors |
                                       update_JxW_values);

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  drag = 0.0;
  lift = 0.0;

  std::vector<Tensor<1, dim>> current_velocity_values(n_q);
  std::vector<double> current_pressure_values(n_q);
  std::vector<Tensor<2, dim>> current_velocity_gradients(n_q);

  double local_lift = 0.0;
  double local_drag = 0.0;

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    fe_values.reinit(cell);

    fe_values[velocity].get_function_values(solution, current_velocity_values);
    fe_values[pressure].get_function_values(solution, current_pressure_values);
    fe_values[velocity].get_function_gradients(solution, current_velocity_gradients);

    if (cell->at_boundary())
    {
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary() &&
            (cell->face(f)->boundary_id() == 5 ||
             cell->face(f)->boundary_id() == 6))
        {
          fe_face_values.reinit(cell, f);

          for (unsigned int q = 0; q < n_q_face; ++q)
          {
            // Get the values
            const double nx = fe_face_values.normal_vector(q)[0];
            const double ny = fe_face_values.normal_vector(q)[1];

            // Construct the tensor
            Tensor<1, dim> tangent;
            tangent[0] = ny;
            tangent[1] = -nx;

            local_drag += (rho * nu * fe_face_values.normal_vector(q) * current_velocity_gradients[q] * tangent * ny -
                           current_pressure_values[q] * nx) *
                          fe_face_values.JxW(q);

            local_lift += (rho * nu * fe_face_values.normal_vector(q) * current_velocity_gradients[q] * tangent * nx +
                           current_pressure_values[q] * ny) *
                          fe_face_values.JxW(q);
          }
        }
      }
    }
  }
  drag = Utilities::MPI::sum(local_drag, MPI_COMM_WORLD);
  lift = Utilities::MPI::sum(local_lift, MPI_COMM_WORLD);
  pcout << "Drag :\t " << drag << " Lift :\t " << lift << std::endl;
  // The mean velocity is defined as 2U(0,H/2,t)/3
  // This is in the case 2D-2 unsteady
  double mean_v = inlet_velocity.getMeanVelocity();
  vec_drag.push_back(drag);
  vec_lift.push_back(lift);
  vec_drag_coeff.push_back((2. * drag) / (mean_v * mean_v * rho * M_PI * 0.1));
  vec_lift_coeff.push_back((2. * lift) / (mean_v * mean_v * rho * M_PI * 0.1));

  pcout
      << "Coeff:\t " << (2. * drag) / (mean_v * mean_v * rho * M_PI * 0.1)
      << " Coeff:\t " << (2. * lift) / (mean_v * mean_v * rho * M_PI * 0.1) << std::endl;

  pcout << "===============================================" << std::endl;
}
