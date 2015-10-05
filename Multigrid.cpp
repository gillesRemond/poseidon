#include <boost/numeric/ublas/vector.hpp>

#include <iostream>

#include "constants.h"
#include "Mesh.h"
#include "Multigrid.h"
#include "Equation.h"

using namespace boost::numeric::ublas;


Multigrid::Multigrid( void )
{
}

Multigrid::Multigrid( Parameters parameters )
{
	int nbCells = parameters.nbCells[0] * parameters.nbCells[1];

	// look how many grids are necessary
	m_nbGrids = 1;
	if ( nbCells > 1500 )
		m_nbGrids = ceil( log( nbCells / 1500.0 ) / log( 4.0 ) ) + 1;

	std::cout << "Le nombre de grilles utilisees est de " << m_nbGrids << std::endl;
	
	// create all the grids
	m_grids.resize( m_nbGrids );
	parameters.nbCells[0] = parameters.nbCells[0] / pow( 2.0, m_nbGrids - 1 );
	parameters.nbCells[1] = parameters.nbCells[1] / pow( 2.0, m_nbGrids - 1 );
	
	for ( int gridIndex ( 0 ); gridIndex < m_nbGrids; ++gridIndex )
	{
		// Generate the mesh
		Mesh* mesh = new Mesh( parameters );
		m_grids[ gridIndex ] = mesh;

		// adjust the nb of cells for the next subgrid
		parameters.nbCells[0] *= 2;
		parameters.nbCells[1] *= 2;
	}

	std::cout << "Generation of the meshes." << std::endl;
	stop_chrono();

}


Multigrid::~Multigrid(void)
{
	for ( int gridIndex ( 0 ); gridIndex < m_nbGrids; ++gridIndex )
	{
		// free the memory of the mesh
		delete m_grids[ gridIndex ];
	}
}

void Multigrid::execute( void )
{
	Equation* equation = new Equation();
	boost::numeric::ublas::vector<double> T;
	boost::numeric::ublas::vector<int> prolongator;

	for ( int gridIndex ( 0 ); gridIndex < m_nbGrids; ++gridIndex )
	{
		std::cout << std::endl;
		std::cout << "Run on the grid " << gridIndex << "." << std::endl;

		// Create the linearised equation
		equation -> generateEquation(  m_grids[ gridIndex ] );
		std::cout << "Generation of the equation." << std::endl;
		stop_chrono();

		// Solve the linear equation
		if ( gridIndex == 0 )
			T = equation -> solveGaussElimination();
		else
			T = equation -> solveSOR( 1.8, 0.5 );
		std::cout << "Solving the equation." << std::endl;
		stop_chrono();

		// Update the temperature field
		m_grids[ gridIndex ] -> setTField( T );

		// Prolongate the solution on the next upper grid
		if ( gridIndex < m_nbGrids - 1 )
		{
			int nbCells_x = m_grids[ gridIndex + 1] -> getNbCells( 0 );
			int nbCells_y = m_grids[ gridIndex + 1] -> getNbCells( 1 );
			
			prolongator.resize( nbCells_x * nbCells_y );

			for ( int j ( 0 ); j < nbCells_y; ++j )
			{
				for ( int i ( 0 ); i < nbCells_x; ++i )
				{
					prolongator[ i + j * nbCells_x ] = i / 2 + j / 2 * nbCells_x / 2;
				}
			}

			m_grids[ gridIndex + 1 ] -> prolongateT( m_grids[ gridIndex ], prolongator );

			std::cout << "Prolongate the solution on the next grid." << std::endl;
			stop_chrono();
		}
	}

	// Print the results
	m_grids[ m_nbGrids - 1 ] -> printData();

	std::cout << "Printing the results." << std::endl;
	stop_chrono();

	// free the memory
	delete equation;
}

