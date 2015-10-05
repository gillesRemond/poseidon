#include <fstream>

#include <boost/numeric/ublas/vector.hpp>

#include "constants.h"
#include "Mesh.h"
#include "Node.h"
#include "CVCenter.h"
#include "Cell.h"


using namespace boost::numeric::ublas;


Mesh::Mesh(void)
{
}

Mesh::Mesh( Parameters const& parameters )
{
	m_nbCells.resize( 2 );
	m_nbCells( 0 ) = parameters.nbCells[ 0 ];
	m_nbCells( 1 ) = parameters.nbCells[ 1 ];
	int nbCells = m_nbCells( 0 ) * m_nbCells( 1 );

	// resize the vectors for the points
	m_cvCenters.resize( nbCells );
	m_vertexesV.resize( nbCells + m_nbCells( 0 ) );
	m_vertexesH.resize( nbCells + m_nbCells( 1 ) );
	m_vertexesC.resize( nbCells + m_nbCells( 0 ) + m_nbCells( 1 ) + 1 );

	//////////////////////////////////////////////////////////
	// create all the points

	double cellLength_x = parameters.lengths[ 0 ] / m_nbCells( 0 );
	double cellLength_y = parameters.lengths[ 1 ] / m_nbCells( 1 );

	// create the control volume center
	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			m_cvCenters[ i + j * m_nbCells( 0 ) ] = new CVCenter( i * cellLength_x + cellLength_x / 2, j * cellLength_y + cellLength_y / 2 );
		}
	}

	// create the vertexesV
	for ( int j (0); j < m_nbCells( 1 ) + 1; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			m_vertexesV[ i + j * m_nbCells( 0 ) ] = new Node( i * cellLength_x + cellLength_x / 2, j * cellLength_y );
		}
	}

	// create the vertexesH
	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) + 1; ++i )
		{
			m_vertexesH[ i + j * ( m_nbCells( 0 ) + 1 ) ] = new Node( i * cellLength_x, j * cellLength_y + cellLength_y / 2 );
		}
	}

	// create the vertexesC
	for ( int j (0); j < m_nbCells( 1 ) + 1; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) + 1; ++i )
		{
			m_vertexesC[ i + j * ( m_nbCells( 0 ) + 1 ) ] = new Node( i * cellLength_x, j * cellLength_y );
		}
	}
	
	//////////////////////////////////////////////////////////
	// create all the cells

	m_cells.resize( nbCells );

	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			Cell* cell = new Cell( );
			// P
			cell -> setCenter( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// N
			if ( j < m_nbCells( 1 ) - 1 )
				cell -> setVertex_N( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// S
			if ( j > 0 )
				cell -> setVertex_S( m_cvCenters[ i + ( j - 1 ) * m_nbCells( 0 ) ] );
			// E
			if ( i < m_nbCells( 0 ) - 1 )
				cell -> setVertex_E( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// W
			if ( i > 0 )
				cell -> setVertex_W( m_cvCenters[ ( i - 1 ) + j * m_nbCells( 0 ) ] );
			// n
			cell -> setVertex_n( m_vertexesV[ i + ( j + 1 ) * m_nbCells( 0 ) ] );
			// s
			cell -> setVertex_s( m_vertexesV[ i + j * m_nbCells( 0 ) ] );
			// e
			cell -> setVertex_e( m_vertexesH[ i + 1 + j * ( m_nbCells( 0 ) + 1 ) ] );
			// w
			cell -> setVertex_w( m_vertexesH[ i + j * ( m_nbCells( 0 ) + 1 ) ] );
			// ne
			cell -> setVertex_ne( m_vertexesC[ i + 1 + ( j + 1 ) * ( m_nbCells( 0 ) + 1 ) ] );
			// nw
			cell -> setVertex_nw( m_vertexesC[ i + ( j + 1 ) * ( m_nbCells( 0 ) + 1 ) ] );
			// se
			cell -> setVertex_se( m_vertexesC[ ( i + 1 ) + j * ( m_nbCells( 0 ) + 1 ) ] );
			// sw
			cell -> setVertex_sw( m_vertexesC[ i + j * ( m_nbCells( 0 ) + 1 ) ] );

			// set the boundary types
			vector<int> boundaryTypes ( 4 );

			// north
			if ( j == m_nbCells( 1 ) - 1 )
			{
				boundaryTypes[ NORTH ] = parameters.boundaryTypes[ NORTH ];
			}
			else
			{
				boundaryTypes[ NORTH ] = INNER;
			}

			// south
			if ( j == 0 )
			{
				boundaryTypes[ SOUTH ] = parameters.boundaryTypes[ SOUTH ];
			}
			else
			{
				boundaryTypes[ SOUTH ] = INNER;
			}

			// east
			if ( i == m_nbCells( 0 ) - 1 )
			{
				boundaryTypes[ EAST ] = parameters.boundaryTypes[ EAST ];
			}
			else
			{
				boundaryTypes[ EAST ] = INNER;
			}

			// west
			if ( i == 0 )
			{
				boundaryTypes[ WEST ] = parameters.boundaryTypes[ WEST ];
			}
			else
			{
				boundaryTypes[ WEST ] = INNER;
			}

			cell -> setBoundaryTypes( boundaryTypes );

	
			m_cells[ i + j * m_nbCells( 0 ) ] = cell; 

		} // end of loop in the x direction
	} // end of loop in the y direction
}


Mesh::~Mesh(void)
{
	for ( int i ( 0 ); i < m_cells.size(); ++i )
	{
		delete m_cells[ i ];
	}

	for ( int i ( 0 ); i < m_cvCenters.size(); ++i )
	{
		delete m_cvCenters[ i ];
	}

	for ( int i ( 0 ); i < m_vertexesV.size(); ++i )
	{
		delete m_vertexesV[ i ];
	}

	for ( int i ( 0 ); i < m_vertexesH.size(); ++i )
	{
		delete m_vertexesH[ i ];
	}

	for ( int i ( 0 ); i < m_vertexesC.size(); ++i )
	{
		delete m_vertexesC[ i ];
	}
}

int Mesh::getNbCells( int index )
{
	if ( index >= 0 && index < 2 )
		return m_nbCells( index );
	else
		return 0;
}

Cell* Mesh::getCell( int index )
{
	if ( index >= 0 && index < m_cells.size() )
		return m_cells[ index ];
	else
		return NULL;
}

void Mesh::initialiseUniformField( double T )
{
	int nbCells ( m_cvCenters.size() );

	for ( int i (0); i < nbCells; ++i )
	{
		m_cvCenters( i ) -> setT( T );
	}
}

void Mesh::setTField( vector<double> const& T )
{
	int nbCells ( m_cvCenters.size() );

	for ( int i (0); i < nbCells; ++i )
	{
		m_cvCenters( i ) -> setT( T( i ) );
	}
}

void Mesh::prolongateT( Mesh* lowerMesh, boost::numeric::ublas::vector<int> const& prolongator )
{
	for ( int i (0); i < m_cvCenters.size(); ++i )
	{
		m_cvCenters[ i ] -> setT( lowerMesh -> getT( prolongator[ i ] ) );
	}
}

double Mesh::getT( int index )
{
	if ( index >= 0 && index < m_cvCenters.size() )
		return m_cvCenters[ index ] -> getT();
	else
		return 0;
}

void Mesh::getT( boost::numeric::ublas::vector<double> &T )
{
	int nbCells ( m_cvCenters.size() );

	T.resize( nbCells );

	for ( int i (0); i < nbCells; ++i )
	{
		T( i ) = m_cvCenters( i ) -> getT();
	}
}

void Mesh::generateMesh( Parameters const& parameters )
{
	m_nbCells.resize( 2 );
	m_nbCells( 0 ) = parameters.nbCells[ 0 ];
	m_nbCells( 1 ) = parameters.nbCells[ 1 ];
	int nbCells = m_nbCells( 0 ) * m_nbCells( 1 );

	// resize the vectors for the points
	m_cvCenters.resize( nbCells );
	m_vertexesV.resize( nbCells + m_nbCells( 0 ) );
	m_vertexesH.resize( nbCells + m_nbCells( 1 ) );
	m_vertexesC.resize( nbCells + m_nbCells( 0 ) + m_nbCells( 1 ) + 1 );

	//////////////////////////////////////////////////////////
	// create all the points

	double cellLength_x = parameters.lengths[ 0 ] / m_nbCells( 0 );
	double cellLength_y = parameters.lengths[ 1 ] / m_nbCells( 1 );

	// create the control volume center
	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			m_cvCenters[ i + j * m_nbCells( 0 ) ] = new CVCenter( i * cellLength_x + cellLength_x / 2, j * cellLength_y + cellLength_y / 2 );
		}
	}

	// create the vertexesV
	for ( int j (0); j < m_nbCells( 1 ) + 1; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			m_vertexesV[ i + j * m_nbCells( 0 ) ] = new Node( i * cellLength_x + cellLength_x / 2, j * cellLength_y );
		}
	}

	// create the vertexesH
	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) + 1; ++i )
		{
			m_vertexesH[ i + j * ( m_nbCells( 0 ) + 1 ) ] = new Node( i * cellLength_x, j * cellLength_y + cellLength_y / 2 );
		}
	}

	// create the vertexesC
	for ( int j (0); j < m_nbCells( 1 ) + 1; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) + 1; ++i )
		{
			m_vertexesC[ i + j * ( m_nbCells( 0 ) + 1 ) ] = new Node( i * cellLength_x, j * cellLength_y );
		}
	}
	
	//////////////////////////////////////////////////////////
	// create all the cells

	m_cells.resize( nbCells );

	for ( int j (0); j < m_nbCells( 1 ); ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ); ++i )
		{
			Cell* cell = new Cell( );
			// P
			cell -> setCenter( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// N
			if ( j < m_nbCells( 1 ) - 1 )
				cell -> setVertex_N( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// S
			if ( j > 0 )
				cell -> setVertex_S( m_cvCenters[ i + ( j - 1 ) * m_nbCells( 0 ) ] );
			// E
			if ( i < m_nbCells( 0 ) - 1 )
				cell -> setVertex_E( m_cvCenters[ i + j * m_nbCells( 0 ) ] );
			// W
			if ( i > 0 )
				cell -> setVertex_W( m_cvCenters[ ( i - 1 ) + j * m_nbCells( 0 ) ] );
			// n
			cell -> setVertex_n( m_vertexesV[ i + ( j + 1 ) * m_nbCells( 0 ) ] );
			// s
			cell -> setVertex_s( m_vertexesV[ i + j * m_nbCells( 0 ) ] );
			// e
			cell -> setVertex_e( m_vertexesH[ i + 1 + j * ( m_nbCells( 0 ) + 1 ) ] );
			// w
			cell -> setVertex_w( m_vertexesH[ i + j * ( m_nbCells( 0 ) + 1 ) ] );
			// ne
			cell -> setVertex_ne( m_vertexesC[ i + 1 + ( j + 1 ) * ( m_nbCells( 0 ) + 1 ) ] );
			// nw
			cell -> setVertex_nw( m_vertexesC[ i + ( j + 1 ) * ( m_nbCells( 0 ) + 1 ) ] );
			// se
			cell -> setVertex_se( m_vertexesC[ ( i + 1 ) + j * ( m_nbCells( 0 ) + 1 ) ] );
			// sw
			cell -> setVertex_sw( m_vertexesC[ i + j * ( m_nbCells( 0 ) + 1 ) ] );

			// set the boundary types
			vector<int> boundaryTypes ( 4 );

			// north
			if ( j == m_nbCells( 1 ) - 1 )
			{
				boundaryTypes[ NORTH ] = parameters.boundaryTypes[ NORTH ];
			}
			else
			{
				boundaryTypes[ NORTH ] = INNER;
			}

			// south
			if ( j == 0 )
			{
				boundaryTypes[ SOUTH ] = parameters.boundaryTypes[ SOUTH ];
			}
			else
			{
				boundaryTypes[ SOUTH ] = INNER;
			}

			// east
			if ( i == m_nbCells( 0 ) - 1 )
			{
				boundaryTypes[ EAST ] = parameters.boundaryTypes[ EAST ];
			}
			else
			{
				boundaryTypes[ EAST ] = INNER;
			}

			// west
			if ( i == 0 )
			{
				boundaryTypes[ WEST ] = parameters.boundaryTypes[ WEST ];
			}
			else
			{
				boundaryTypes[ WEST ] = INNER;
			}

			cell -> setBoundaryTypes( boundaryTypes );

	
			m_cells[ i + j * m_nbCells( 0 ) ] = cell; 

		} // end of loop in the x direction
	} // end of loop in the y direction
}

void Mesh::printData()
{
	int nbCells ( m_cells.size() );

	// Create the output file
	std::ofstream fout;
	fout.open ("output.txt", std::ios::trunc);

	/*output << "index\t";
	output << "x (in m)\t";
	output << "y (in m)\t";
	output << "T (in °C)\n";*/
	
	for ( int i (0); i < m_nbCells( 0 ); ++i )
	{
		//fout << m_cells( i ).getIndex() << "\t";
		fout << i << "\t";
		fout << m_cells( i ) -> getVertex_sw() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( i ) -> getVertex_sw() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( i ) -> getT() << "\n";
	}

	fout << 10 << "\t";
	fout << m_cells( m_nbCells( 0 ) - 1 ) -> getVertex_se() -> getPosition( 0 ) << "\t\t";
	fout <<  m_cells( m_nbCells( 0 ) - 1 ) -> getVertex_se() -> getPosition( 1 ) << "\t\t";
	fout << m_cvCenters( m_nbCells( 0 ) - 1 ) -> getT() << "\n";
	
	for ( int j (1); j < m_nbCells( 1 ) ; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) ; ++i )
		{
			//fout << m_cells( i ).getIndex() << "\t";
			fout << i + j * m_nbCells( 0 ) << "\t";
			fout << m_cells( i + j * m_nbCells( 0 ) ) -> getVertex_sw() -> getPosition( 0 ) << "\t\t";
			fout <<  m_cells( i + j * m_nbCells( 0 ) ) -> getVertex_sw() -> getPosition( 1 ) << "\t\t";
			fout << m_cvCenters( i + j * m_nbCells( 0 ) ) -> getT() << "\n";
		}

		fout << m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) << "\t";
		fout << m_cells( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getVertex_se() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getVertex_se() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getT() << "\n";
	}

	for ( int i (0); i < m_nbCells( 0 ); ++i )
	{
		//fout << m_cells( i ).getIndex() << "\t";
		fout << i << "\t";
		fout << m_cells( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getVertex_nw() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getVertex_nw() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getT() << "\n";
	}

	fout << 10 << "\t";
	fout << m_cells( nbCells - 1 ) -> getVertex_ne() -> getPosition( 0 ) << "\t\t";
	fout <<  m_cells( nbCells - 1 ) -> getVertex_ne() -> getPosition( 1 ) << "\t\t";
	fout << m_cvCenters( nbCells - 1 ) -> getT() << "\n";

	fout.close();
}

void Mesh::printData( std::string name )
{
	int nbCells ( m_cells.size() );

	// Create the output file
	std::ofstream fout;
	fout.open ( name, std::ios::trunc );

	/*output << "index\t";
	output << "x (in m)\t";
	output << "y (in m)\t";
	output << "T (in °C)\n";*/
	
	for ( int i (0); i < m_nbCells( 0 ); ++i )
	{
		//fout << m_cells( i ).getIndex() << "\t";
		fout << i << "\t";
		fout << m_cells( i ) -> getVertex_sw() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( i ) -> getVertex_sw() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( i ) -> getT() << "\n";
	}

	fout << 10 << "\t";
	fout << m_cells( m_nbCells( 0 ) - 1 ) -> getVertex_se() -> getPosition( 0 ) << "\t\t";
	fout <<  m_cells( m_nbCells( 0 ) - 1 ) -> getVertex_se() -> getPosition( 1 ) << "\t\t";
	fout << m_cvCenters( m_nbCells( 0 ) - 1 ) -> getT() << "\n";
	
	for ( int j (1); j < m_nbCells( 1 ) ; ++j )
	{
		for ( int i (0); i < m_nbCells( 0 ) ; ++i )
		{
			//fout << m_cells( i ).getIndex() << "\t";
			fout << i + j * m_nbCells( 0 ) << "\t";
			fout << m_cells( i + j * m_nbCells( 0 ) ) -> getVertex_sw() -> getPosition( 0 ) << "\t\t";
			fout <<  m_cells( i + j * m_nbCells( 0 ) ) -> getVertex_sw() -> getPosition( 1 ) << "\t\t";
			fout << m_cvCenters( i + j * m_nbCells( 0 ) ) -> getT() << "\n";
		}

		fout << m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) << "\t";
		fout << m_cells( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getVertex_se() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getVertex_se() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( m_nbCells( 0 ) - 1 + j * m_nbCells( 0 ) ) -> getT() << "\n";
	}

	for ( int i (0); i < m_nbCells( 0 ); ++i )
	{
		//fout << m_cells( i ).getIndex() << "\t";
		fout << i << "\t";
		fout << m_cells( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getVertex_nw() -> getPosition( 0 ) << "\t\t";
		fout <<  m_cells( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getVertex_nw() -> getPosition( 1 ) << "\t\t";
		fout << m_cvCenters( i + (  m_nbCells( 1 ) - 1 ) * m_nbCells( 0 ) ) -> getT() << "\n";
	}

	fout << 10 << "\t";
	fout << m_cells( nbCells - 1 ) -> getVertex_ne() -> getPosition( 0 ) << "\t\t";
	fout <<  m_cells( nbCells - 1 ) -> getVertex_ne() -> getPosition( 1 ) << "\t\t";
	fout << m_cvCenters( nbCells - 1 ) -> getT() << "\n";

	fout.close();
}