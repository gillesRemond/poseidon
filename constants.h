#pragma once

#include <boost/numeric/ublas/vector.hpp>

//using namespace boost::numeric::ublas;

enum{ NORTH, SOUTH, EAST, WEST };
enum{ INNER, DIRICHLET, NEUMANN, ROBIN };

void demarrer_chrono();
void stop_chrono();
void stop_chrono_total();

typedef struct Parameters Parameters;
struct Parameters
{
	boost::numeric::ublas::vector<double> lengths;
	boost::numeric::ublas::vector<int> nbCells;
	boost::numeric::ublas::vector<int> boundaryTypes;
};