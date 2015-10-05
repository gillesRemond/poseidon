#include "CVCenter.h"


CVCenter::CVCenter(void)
{

}

CVCenter::CVCenter( double x, double y ) : Node( x, y )
{
	m_T = 0;
}

CVCenter::CVCenter( int index, boost::numeric::ublas::vector<double> position ) : Node( position )
{
	m_index = index;
}

void CVCenter::setIndex( int index )
{
	m_index = index;
}

void CVCenter::setT( double T )
{
	m_T = T;
}

int CVCenter::getIndex()
{
	return m_index;
}

double CVCenter::getT()
{
	return m_T;
}
