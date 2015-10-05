#pragma once

#include <boost/numeric/ublas/vector.hpp>

class Node
{
private :
	boost::numeric::ublas::vector<double> m_position;

public :
	Node();
	Node( boost::numeric::ublas::vector<double> position );
	Node( double x, double y );
	void setPosition(boost::numeric::ublas::vector<double> position);
	boost::numeric::ublas::vector<double> getPosition();
	double getPosition( int coordinate );

};