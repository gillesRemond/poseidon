#pragma once

#include "CVCenter.h"
#include "Node.h"

class CVCenter : public Node
{
private :
	int m_index;
	double m_T;

public:
	CVCenter( void );
	CVCenter( double x, double y );
	CVCenter( int index, boost::numeric::ublas::vector<double> position );
	void setIndex( int index );
	void setT( double T );
	int getIndex();
	double getT();
};
