#pragma once

#include "constants.h"
#include "Mesh.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


class Equation
{
private:
	boost::numeric::ublas::vector<double> m_D0;
	boost::numeric::ublas::vector<double> m_D1;
	boost::numeric::ublas::vector<double> m_D_1;
	boost::numeric::ublas::vector<double> m_D2;
	boost::numeric::ublas::vector<double> m_D_2;
	boost::numeric::ublas::vector<double> m_D3;
	boost::numeric::ublas::vector<double> m_D_3;
	boost::numeric::ublas::vector<double> m_D4;
	boost::numeric::ublas::vector<double> m_D_4;
	boost::numeric::ublas::vector<int> m_D1_index; // the jth postion in the column
	boost::numeric::ublas::vector<int> m_D_1_index;
	boost::numeric::ublas::vector<int> m_D2_index;
	boost::numeric::ublas::vector<int> m_D_2_index;
	boost::numeric::ublas::vector<int> m_D3_index;
	boost::numeric::ublas::vector<int> m_D_3_index;
	boost::numeric::ublas::vector<int> m_D4_index;
	boost::numeric::ublas::vector<int> m_D_4_index;
	boost::numeric::ublas::vector<double> m_B;
	boost::numeric::ublas::vector<double> m_T;

public:
	Equation();
	Equation( Parameters parameters );
	const void computeInnerFlux( int );
	const void computeNeumannFlux();
	void generateEquation( Mesh* mesh );
	boost::numeric::ublas::vector<double> solveSOR( double omega, double maxResidual );
	boost::numeric::ublas::vector<double> solveGaussElimination( );

};