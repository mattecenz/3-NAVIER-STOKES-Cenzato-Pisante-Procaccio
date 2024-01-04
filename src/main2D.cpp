#include "NavierStokes2D.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  const std::string  mesh_file_name  = argc>1?argv[1]:"../mesh/cilinder_2D_coarse.msh";
  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;

	const double T      = 2;
	const double deltat = 0.1;

  NavierStokes problem(mesh_file_name, degree_velocity, degree_pressure, T, deltat);

  problem.setup();
  problem.solve();

  return 0;
}
