#include <boost/numeric/ublas/vector.hpp>

#include "Node.h"
#include "CVCenter.h"
#include "Cell.h"
#include "constants.h"

using namespace boost::numeric::ublas;

Cell::Cell(void)
{
}

//Cell::Cell( int const& index, vector< int > const& boundaryTypes, vector< CVCenter* > const& cvCenters,
//		vector< Node* > const& vertexesV, vector< Node* > const& vertexesH, vector< Node* > const& vertexesC )
//{
//	if ( index > 0 && index < cvCenters.size() )
//	{
//		m_center = cvCenters[ index ];
//		// N
//		if ( boundaryTypes[ NORTH ] == INNER )
//			m_vertex_N = cvCenters[ index ];
//
//
//	}
//
//
//
//	if ( i > 1 )
//	m_vertex_N;
//	m_vertex_S;
//	m_vertex_E;
//	m_vertex_W;
//	m_vertex_n;
//	m_vertex_s;
//	m_vertex_e;
//	m_vertex_w;
//	m_vertex_ne;
//	m_vertex_nw;
//	m_vertex_se;
//	m_vertex_sw;
//	boost::numeric::ublas::vector<int> m_boundaryTypes;
//
//
//}



Cell::~Cell(void)
{
}


void Cell::setCenter( CVCenter* cvCenter )
{
	m_center = cvCenter;
}

void Cell::setVertex_N( Node *vertex )
{
	m_vertex_N = vertex;
}

void Cell::setVertex_S( Node *vertex )
{
	m_vertex_S = vertex;
}

void Cell::setVertex_E( Node *vertex )
{
	m_vertex_E = vertex;
}

void Cell::setVertex_W( Node *vertex )
{
	m_vertex_W = vertex;
}

void Cell::setVertex_n( Node *vertex )
{
	m_vertex_n = vertex;
}

void Cell::setVertex_s( Node *vertex )
{
	m_vertex_s = vertex;
}

void Cell::setVertex_e( Node *vertex )
{
	m_vertex_e = vertex;
}

void Cell::setVertex_w( Node *vertex )
{
	m_vertex_w = vertex;
}

void Cell::setVertex_ne( Node *vertex )
{
	m_vertex_ne = vertex;
}

void Cell::setVertex_nw( Node *vertex )
{
	m_vertex_nw = vertex;
}

void Cell::setVertex_se( Node *vertex )
{
	m_vertex_se = vertex;
}

void Cell::setVertex_sw( Node *vertex )
{
	m_vertex_sw = vertex;
}

void Cell::setBoundaryTypes( boost::numeric::ublas::vector<int> boundaryTypes )
{
	if ( boundaryTypes.size() == 4 )
		m_boundaryTypes = boundaryTypes;
}

CVCenter* Cell::getCenter()
{
	return m_center;
}

Node* Cell::getVertex_N()
{
	return m_vertex_N;
}

Node* Cell::getVertex_S()
{
	return m_vertex_S;
}

Node* Cell::getVertex_E()
{
	return m_vertex_E;
}

Node* Cell::getVertex_W()
{
	return m_vertex_W;
}

Node* Cell::getVertex_n()
{
	return m_vertex_n;
}

Node* Cell::getVertex_s()
{
	return m_vertex_s;
}

Node* Cell::getVertex_e()
{
	return m_vertex_e;
}

Node* Cell::getVertex_w()
{
	return m_vertex_w;
}

Node* Cell::getVertex_ne()
{
	return m_vertex_ne;
}

Node* Cell::getVertex_nw()
{
	return m_vertex_nw;
}

Node* Cell::getVertex_se()
{
	return m_vertex_se;
}

Node* Cell::getVertex_sw()
{
	return m_vertex_sw;
}

vector<int> Cell::getBoundaryTypes()
{
	return m_boundaryTypes;
}

int Cell::getBoundaryTypes( int index )
{
	if ( index >= 0 && index < 4 )
		return m_boundaryTypes[ index ];
}
