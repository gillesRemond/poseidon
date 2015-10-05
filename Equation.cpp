// Equation.cpp
// Class that creates the discretized equation and solves it

#include "constants.h"
#include "Mesh.h"
#include "Equation.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <cmath>

using namespace boost::numeric::ublas;


Equation::Equation()
{}

Equation::Equation( Parameters parameters )
{
	int nbCells_x ( parameters.nbCells[ 0 ] );
	int nbCells_y ( parameters.nbCells[ 1 ] );

	int nbCells ( nbCells_x * nbCells_y );

	m_D0  = zero_vector<double> ( nbCells );
	m_D1  = zero_vector<double> ( nbCells );
	m_D_1 = zero_vector<double> ( nbCells );
	m_D2  = zero_vector<double> ( nbCells );
	m_D_2 = zero_vector<double> ( nbCells );
	m_D3  = zero_vector<double> ( nbCells );
	m_D_3 = zero_vector<double> ( nbCells );
	m_D4  = zero_vector<double> ( nbCells );
	m_D_4 = zero_vector<double> ( nbCells );

	m_D1_index.resize( nbCells );
	m_D_1_index.resize( nbCells );
	m_D2_index.resize( nbCells );
	m_D_2_index.resize( nbCells );
	m_D3_index.resize( nbCells );
	m_D_3_index.resize( nbCells );
	m_D4_index.resize( nbCells );
	m_D_4_index.resize( nbCells );

	for ( int i ( 0 ); i < nbCells; ++i )
	{
		m_D1_index( i ) = 0;
		m_D_1_index( i ) = 0;
		m_D2_index( i ) = 0;
		m_D_2_index( i ) = 0;
		m_D3_index( i ) = 0;
		m_D_3_index( i ) = 0;
		m_D4_index( i ) = 0;
		m_D_4_index( i ) = 0;
	}

	m_B = zero_vector<double> ( nbCells );

	m_T = zero_vector<double> ( nbCells );
}

void Equation::generateEquation( Mesh* mesh )
{
	int nbCells_x ( mesh -> getNbCells( 0 ) );
	int nbCells_y ( mesh -> getNbCells( 1 ) );

	int nbCells = nbCells_x * nbCells_y;

	// initialize all the vectors
	m_D0  = zero_vector<double> ( nbCells );
	m_D1  = zero_vector<double> ( nbCells );
	m_D_1 = zero_vector<double> ( nbCells );
	m_D2  = zero_vector<double> ( nbCells );
	m_D_2 = zero_vector<double> ( nbCells );
	m_D3  = zero_vector<double> ( nbCells );
	m_D_3 = zero_vector<double> ( nbCells );
	m_D4  = zero_vector<double> ( nbCells );
	m_D_4 = zero_vector<double> ( nbCells );

	m_D1_index.resize( nbCells );
	m_D_1_index.resize( nbCells );
	m_D2_index.resize( nbCells );
	m_D_2_index.resize( nbCells );
	m_D3_index.resize( nbCells );
	m_D_3_index.resize( nbCells );
	m_D4_index.resize( nbCells );
	m_D_4_index.resize( nbCells );

	for ( int i ( 0 ); i < nbCells; ++i )
	{
		m_D1_index( i ) = 0;
		m_D_1_index( i ) = 0;
		m_D2_index( i ) = 0;
		m_D_2_index( i ) = 0;
		m_D3_index( i ) = 0;
		m_D_3_index( i ) = 0;
		m_D4_index( i ) = 0;
		m_D_4_index( i ) = 0;
	}

	m_B = zero_vector<double> ( nbCells );

	//m_T = zero_vector<double> ( nbCells );
	mesh -> getT( m_T );

	// create the equation
	for ( int j (0); j < nbCells_y; ++j )
	{
		for ( int i (0); i < nbCells_x; ++i )
		{
			int index = i + j * nbCells_x;

			// we first put the right column indexes for all diagonals stored in the vectors
			if ( index < nbCells - 1 )
				m_D1_index( index ) = index + 1;
			if ( index < nbCells - nbCells_x + 1 )
				m_D2_index( index ) = index + nbCells_x - 1;
			if ( index < nbCells - nbCells_x )
				m_D3_index( index ) = index + nbCells_x;
			if ( index < nbCells - nbCells_x - 1 )
				m_D4_index( index ) = index + nbCells_x + 1;
			if ( index > 0 )
				m_D_1_index( index ) = index - 1;
			if ( index > nbCells_x - 2 )
				m_D_2_index( index ) = index - nbCells_x + 1;
			if ( index > nbCells_x - 1 )
				m_D_3_index( index ) = index - nbCells_x;
			if ( index > nbCells_x )
				m_D_4_index( index ) = index - nbCells_x - 1;
			

			Cell* cell = mesh -> getCell( i + j * nbCells_x );
						
			// we start with the flux going through the northern boundary of the cell
			if ( cell -> getBoundaryTypes( NORTH ) == INNER )
			{
				Cell* cellN = mesh -> getCell( index + nbCells_x );

				double area_n = std::abs( 
					( cell -> getVertex_e() -> getPosition( 0 ) - cellN -> getVertex_w() -> getPosition( 0 ) ) *
					( cellN -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) -
					( cell -> getVertex_e() -> getPosition( 1 ) - cellN -> getVertex_w() -> getPosition( 1 ) ) *
					( cellN -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) ) / 2;

				m_D0[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) +
					    ( cellN -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_e() -> getPosition( 0 ) ) / 4 +
					    ( cell -> getVertex_w() -> getPosition( 0 ) - cellN -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cell -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) +
					    ( cellN -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_e() -> getPosition( 1 ) ) / 4 +
					    ( cell -> getVertex_w() -> getPosition( 1 ) - cellN -> getVertex_w() -> getPosition( 1 ) ) / 4 ) ) / area_n;

				m_D1[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cellN -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_e() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cellN -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_e() -> getPosition( 1 ) ) / 4 ) ) / area_n;

				m_D4[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cellN -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_e() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cellN -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_e() -> getPosition( 1 ) ) / 4 ) ) / area_n;

				m_D_1[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_w() -> getPosition( 0 ) - cellN -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cell -> getVertex_w() -> getPosition( 1 ) - cellN -> getVertex_w() -> getPosition( 1 ) ) / 4 ) ) / area_n;

				m_D2[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_w() -> getPosition( 0 ) - cellN -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cell -> getVertex_w() -> getPosition( 1 ) - cellN -> getVertex_w() -> getPosition( 1 ) ) / 4 ) ) / area_n;

				m_D3[ index ] +=
					( ( cell -> getVertex_nw() -> getPosition( 0 ) - cell -> getVertex_ne() -> getPosition( 0 ) ) * 
					  ( ( cellN -> getVertex_w() -> getPosition( 0 ) - cellN -> getVertex_e() -> getPosition( 0 ) ) +
					    ( cellN -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_e() -> getPosition( 0 ) ) / 4 +
					    ( cell -> getVertex_w() -> getPosition( 0 ) - cellN -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_nw() -> getPosition( 1 ) - cell -> getVertex_ne() -> getPosition( 1 ) ) *
					  ( ( cellN -> getVertex_w() -> getPosition( 1 ) - cellN -> getVertex_e() -> getPosition( 1 ) ) +
					    ( cellN -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_e() -> getPosition( 1 ) ) / 4 +
					    ( cell -> getVertex_w() -> getPosition( 1 ) - cellN -> getVertex_w() -> getPosition( 1 ) ) / 4 ) ) / area_n;
			}

			// we do the flux going through the southern boundary of the cell
			if ( cell -> getBoundaryTypes( SOUTH ) == INNER )
			{
				Cell* cellS = mesh -> getCell( index - nbCells_x );

				double area_s = std::abs( 
					( cellS -> getVertex_e() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) *
					( cell -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_w() -> getPosition( 1 ) ) -
					( cellS -> getVertex_e() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) *
					( cell -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_w() -> getPosition( 0 ) ) ) / 2;

				m_D0[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_w() -> getPosition( 0 ) - cell -> getVertex_e() -> getPosition( 0 ) ) +
					    ( cell -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_e() -> getPosition( 0 ) ) / 4 +
					    ( cellS -> getVertex_w() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_w() -> getPosition( 1 ) - cell -> getVertex_e() -> getPosition( 1 ) ) +
					    ( cell -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_e() -> getPosition( 1 ) ) / 4 +
					    ( cellS -> getVertex_w() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) / 4 )  ) / area_s;

				m_D1[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_e() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_e() -> getPosition( 1 ) ) / 4 )  ) / area_s;

				m_D_2[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_e() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_e() -> getPosition( 1 ) ) / 4 )  ) / area_s;

				m_D_1[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cellS -> getVertex_w() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cellS -> getVertex_w() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) / 4 )  ) / area_s;

				m_D_4[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cellS -> getVertex_w() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cellS -> getVertex_w() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) / 4 )  ) / area_s;

				m_D_3[ index ] +=
					( ( cell -> getVertex_se() -> getPosition( 0 ) - cell -> getVertex_sw() -> getPosition( 0 ) ) * 
					  ( ( cellS -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_w() -> getPosition( 0 ) ) +
					    ( cell -> getVertex_e() -> getPosition( 0 ) - cellS -> getVertex_e() -> getPosition( 0 ) ) / 4 +
					    ( cellS -> getVertex_w() -> getPosition( 0 ) - cell -> getVertex_w() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_se() -> getPosition( 1 ) - cell -> getVertex_sw() -> getPosition( 1 ) ) * 
					  ( ( cellS -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_w() -> getPosition( 1 ) ) +
					    ( cell -> getVertex_e() -> getPosition( 1 ) - cellS -> getVertex_e() -> getPosition( 1 ) ) / 4 +
					    ( cellS -> getVertex_w() -> getPosition( 1 ) - cell -> getVertex_w() -> getPosition( 1 ) ) / 4 )  ) / area_s;
			}

			// we do the flux going through the eastern boundary of the cell
			if ( cell -> getBoundaryTypes( EAST ) == INNER )
			{
				Cell* cellE = mesh -> getCell( index + 1 );

				double area_e = std::abs( 
					( cellE -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) *
					( cellE -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) -
					( cellE -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) *
					( cellE -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) ) / 2;

				m_D0[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) +
					    ( cellE -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) / 4 +
					    ( cell -> getVertex_n() -> getPosition( 0 ) - cellE -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) +
					    ( cellE -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) / 4 +
					    ( cell -> getVertex_n() -> getPosition( 1 ) - cellE -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_e;

				m_D3[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 0 ) - cellE -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 1 ) - cellE -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_e;

				m_D4[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 0 ) - cellE -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 1 ) - cellE -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_e;

				m_D_3[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cellE -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cellE -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) / 4 )  ) / area_e;

				m_D_2[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cellE -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cellE -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) / 4 )  ) / area_e;

				m_D1[ index ] +=
					( ( cell -> getVertex_ne() -> getPosition( 0 ) - cell -> getVertex_se() -> getPosition( 0 ) ) * 
					  ( ( cellE -> getVertex_n() -> getPosition( 0 ) - cellE -> getVertex_s() -> getPosition( 0 ) ) +
					    ( cellE -> getVertex_s() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) / 4 +
					    ( cell -> getVertex_n() -> getPosition( 0 ) - cellE -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_ne() -> getPosition( 1 ) - cell -> getVertex_se() -> getPosition( 1 ) ) * 
					  ( ( cellE -> getVertex_n() -> getPosition( 1 ) - cellE -> getVertex_s() -> getPosition( 1 ) ) +
					    ( cellE -> getVertex_s() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) / 4 +
					    ( cell -> getVertex_n() -> getPosition( 1 ) - cellE -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_e;
			}

			// we do the flux going through the western boundary of the cell
			if ( cell -> getBoundaryTypes( WEST ) == INNER )
			{
				Cell* cellW = mesh -> getCell( index - 1 );

				double area_w = std::abs( 
					( cell -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_n() -> getPosition( 0 ) ) *
					( cell -> getVertex_n() -> getPosition( 1 ) - cellW -> getVertex_s() -> getPosition( 1 ) ) -
					( cell -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_n() -> getPosition( 1 ) ) *
					( cell -> getVertex_n() -> getPosition( 0 ) - cellW -> getVertex_s() -> getPosition( 0 ) ) ) / 2;

				m_D0[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_s() -> getPosition( 0 ) ) +
					    ( cell -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_s() -> getPosition( 0 ) ) / 4 +
					    ( cellW -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_s() -> getPosition( 1 ) ) +
					    ( cell -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_s() -> getPosition( 1 ) ) / 4 +
					    ( cellW -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_w;

				m_D2[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cellW -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cellW -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_w;

				m_D3[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cellW -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cellW -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_w;

				m_D_3[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_s() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_s() -> getPosition( 1 ) ) / 4 )  ) / area_w;

				m_D_4[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_s() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cell -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_s() -> getPosition( 1 ) ) / 4 )  ) / area_w;

				m_D_1[ index ] +=
					( ( cell -> getVertex_sw() -> getPosition( 0 ) - cell -> getVertex_nw() -> getPosition( 0 ) ) * 
					  ( ( cellW -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_n() -> getPosition( 0 ) ) +
					    ( cell -> getVertex_s() -> getPosition( 0 ) - cellW -> getVertex_s() -> getPosition( 0 ) ) / 4 +
					    ( cellW -> getVertex_n() -> getPosition( 0 ) - cell -> getVertex_n() -> getPosition( 0 ) ) / 4 ) +
					  ( cell -> getVertex_sw() -> getPosition( 1 ) - cell -> getVertex_nw() -> getPosition( 1 ) ) * 
					  ( ( cellW -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_n() -> getPosition( 1 ) ) +
					    ( cell -> getVertex_s() -> getPosition( 1 ) - cellW -> getVertex_s() -> getPosition( 1 ) ) / 4 +
					    ( cellW -> getVertex_n() -> getPosition( 1 ) - cell -> getVertex_n() -> getPosition( 1 ) ) / 4 )  ) / area_w;
			}

			if ( i==0 || j==0 )
			{
				m_D0( index ) = 1;
				m_D1( index ) = 0;
				m_D2( index ) = 0;
				m_D3( index ) = 0;
				m_D4( index ) = 0;
				m_D_1( index ) = 0;
				m_D_2( index ) = 0;
				m_D_3( index ) = 0;
				m_D_4( index ) = 0;
				m_B( index ) = 10;
			}
			else if ( i == ( nbCells_x - 1 ) || j == ( nbCells_y - 1) )
			{
				m_D0( index ) = 1;
				m_D1( index ) = 0;
				m_D2( index ) = 0;
				m_D3( index ) = 0;
				m_D4( index ) = 0;
				m_D_1( index ) = 0;
				m_D_2( index ) = 0;
				m_D_3( index ) = 0;
				m_D_4( index ) = 0;
				m_B( index) = 100;
			}
			
		}
	}

	//std::cout << m_D0 << std::endl;
	//std::cout << m_D1 << std::endl;
	//std::cout << m_D_1 << std::endl;
	//std::cout << m_D3 << std::endl;
	//std::cout << m_D_3 << std::endl;
	//std::cout << m_D4 << std::endl;
	//std::cout << m_D2 << std::endl;
	//std::cout << m_D_2 << std::endl;
	//std::cout << m_D_4 << std::endl;
	//std::cout << m_B << std::endl;

	//std::cout << m_D1_index << std::endl;
	//std::cout << m_D2_index << std::endl;
	//std::cout << m_D3_index << std::endl;
	//std::cout << m_D4_index << std::endl;
	//std::cout << m_D_1_index << std::endl;
	//std::cout << m_D_2_index << std::endl;
	//std::cout << m_D_3_index << std::endl;
	//std::cout << m_D_4_index << std::endl;

}


// ajouter une gestion de l'erreur (au cas où la boucle ne converge pas)
vector<double> Equation::solveSOR( double omega, double maxResidual)
{
	int nbCells ( m_T.size() );

	//m_T = m_B;
	vector<double> TPrevious ( m_T );
	vector<double> deltaPrevious ( m_T ); //the difference between T_n and T_n-1
	vector<double> delta ( m_T );
	double residual ( maxResidual + 1 ); // variable estimating the error of the solution during the iterative process of the Gauss-Seidel algorithm
	int nbOfIterations ( 0 ); // variable counting the number of iterations of the Gauss-Seidel algorithm
	double aux ( 0 ); // temporary variable used in the algorithm	

	while ( residual > maxResidual && nbOfIterations < nbCells )
	{
		for ( int i (0); i < nbCells; ++i )
		{
			aux = 0;
			aux += m_D_1( i ) * m_T( m_D_1_index( i ) );
			aux += m_D_2( i ) * m_T( m_D_2_index( i ) );
			aux += m_D_3( i ) * m_T( m_D_3_index( i ) );
			aux += m_D_4( i ) * m_T( m_D_4_index( i ) );
			aux += m_D1( i ) * m_T( m_D1_index( i ) );
			aux += m_D2( i ) * m_T( m_D2_index( i ) );
			aux += m_D3( i ) * m_T( m_D3_index( i ) );
			aux += m_D4( i ) * m_T( m_D4_index( i ) );

			/*std::cout << m_T( i ) << " " << aux << std::endl;*/

			m_T( i ) = omega * ( m_B( i ) - aux ) / m_D0( i ) + ( 1 - omega ) * m_T( i );
		}

		// Calculate the residual. If bigger than maxResidual continue the loop.
		if ( nbOfIterations == 0 ) // For the first iteration, there is no TPrevious and no deltaPrevious -> no residual
		{
			residual = maxResidual+1;
		}
		else if( nbOfIterations == 1 ) // For the second iteration, there is no deltaPrevious -> no residual
		{
			residual = maxResidual+1;
			delta = m_T - TPrevious;
		}
		else
		{
			delta = m_T - TPrevious;

			double lambda1 ( norm_2( delta ) / norm_2( deltaPrevious ) );
			residual = norm_2( delta ) / ( 1 - lambda1 );
		}

		TPrevious = m_T;
		deltaPrevious = delta;
		++nbOfIterations;
	}

	std::cout << std::endl;
	std::cout << "The number of iterations in the SOR algorithm is " << nbOfIterations << std::endl;
	std::cout << "The error is " << residual << std::endl;
	std::cout << std::endl;

	return m_T;
}

vector<double> Equation::solveGaussElimination( )
{
	int nbCells ( m_D0.size() );

	// create the A matrix
	matrix<double> A ( zero_matrix<double> ( nbCells, nbCells ) );

	for ( int i (0); i < nbCells; i++ )
	{
			A( i, m_D1_index( i ) ) = m_D1( i );
			A( i, m_D2_index( i ) ) = m_D2( i );
			A( i, m_D3_index( i ) ) = m_D3( i );
			A( i, m_D4_index( i ) ) = m_D4( i );
			A( i, m_D_4_index( i ) ) = m_D_4( i );
			A( i, m_D_3_index( i ) ) = m_D_3( i );
			A( i, m_D_2_index( i ) ) = m_D_2( i );
			A( i, m_D_1_index( i ) ) = m_D_1( i );
			A( i, i ) = m_D0( i );
	}

	//std::cout << A << std::endl;
	//std::cout << m_B << std::endl;

	// Gauss Elimination
	//Trigonalisation of the A matrix
	for ( int k ( 0 ); k < ( nbCells - 1 ); ++k )
	{
		for ( int l ( k + 1 ); ( l < k + m_D4_index( k ) ) && ( l < nbCells ) ; ++l )
		{
			for ( int p ( k + 1 ); p < nbCells; ++p )
			{
				A( l, p ) -= A( l, k ) / A( k, k ) * A( k, p );
			}

			m_B( l ) -= A( l, k ) / A( k, k ) * m_B( k );

			A( l, k ) = 0;
		}

		/*std::cout << A << std::endl;
		std::cout << m_B << std::endl;*/
	}

	//Solving the linear equation
	for ( int k ( nbCells - 1 ); k >= 0; --k)
	{
		// temporary variable
		double aux ( 0 );

		for ( int p ( k + 1 ); p < nbCells; ++p )
		{
			aux += A( k, p ) * m_T( p );
		}

		m_T( k ) = ( m_B( k ) - aux ) / A( k, k );
	}

	/*std::cout<<prod(A,m_T)<<std::endl;*/

	return m_T;
}