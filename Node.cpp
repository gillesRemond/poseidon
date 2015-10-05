#include <boost/numeric/ublas/vector.hpp>

#include "Node.h"

using namespace boost::numeric::ublas;

Node::Node(void)
{
}

Node::Node( vector<double> position )
{
	m_position = position;
}

Node::Node( double x, double y )
{
	m_position.resize( 2 );
	m_position[ 0 ] = x;
	m_position[ 1 ] = y;
}

void Node::setPosition( vector<double> position )
{
	m_position = position;
}

vector<double> Node::getPosition()
{
	return m_position;
}

double Node::getPosition( int coordinate )
{
	if ( coordinate >= 0 && coordinate < 2 )
		return m_position[ coordinate ];
	else 
		return 0;
}
