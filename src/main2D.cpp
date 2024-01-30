#include "NavierStokes2D.hpp"

// Main function.
int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string mesh_file_name = argc > 1 ? argv[1] : "../mesh/cilinder_2D_fine.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

  const double T = 8;
  const double deltat = 0.005;

  dealii::Timer timer;
  // Start the timer
  timer.restart();

  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat);

  problem.setup();
  problem.solve();

  // Stop the timer
  timer.stop();

  // Output the elapsed time
  if(rank == 0)
    std::cout << "Time taken to solve ENTIRE Navier Stokes problem: " << timer.wall_time() << " seconds" << std::endl;

  if (rank == 0)
  {
    const std::string output_filename = "forces_results_2D_2case.csv";
    std::ofstream outputFile(output_filename);

    if (!outputFile.is_open())
    {
      std::cerr << "Error opening output file" << std::endl;
      return -1;
    }
    outputFile << "Iteration, Drag, Lift, Coeff Drag, CoeffLift, time prec, time solve" << std::endl;

    for (size_t ite = 0; ite < problem.vec_drag.size(); ite++)
    {
      outputFile << ite * deltat << ", " << problem.vec_drag[ite] << ", " << problem.vec_lift_coeff[ite] << ", " 
                << problem.vec_drag_coeff[ite] << ", " << problem.vec_lift_coeff[ite] << ", "
                << problem.time_prec[ite] << ", " << problem.time_solve[ite]
                << std::endl;
    }
    outputFile.close();
  }

  return 0;
}
