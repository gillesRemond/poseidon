#pragma once

#include <boost/numeric/ublas/vector.hpp>

#include "Cell.h"
#include "CVCenter.h"
#include "Node.h"

class Cell
{
private:
	CVCenter* m_center;
	Node* m_vertex_N;
	Node* m_vertex_S;
	Node* m_vertex_E;
	Node* m_vertex_W;
	Node* m_vertex_n;
	Node* m_vertex_s;
	Node* m_vertex_e;
	Node* m_vertex_w;
	Node* m_vertex_ne;
	Node* m_vertex_nw;
	Node* m_vertex_se;
	Node* m_vertex_sw;
	boost::numeric::ublas::vector<int> m_boundaryTypes;

public:
	Cell(void);
	/*Cell( int const& index,
		boost::numeric::ublas::vector< int > const& boundaryTypes,
		boost::numeric::ublas::vector< CVCenter* > const& cvCenters,
		boost::numeric::ublas::vector< Node* > const& vertexesV,
		boost::numeric::ublas::vector< Node* > const& vertexesH,
		boost::numeric::ublas::vector< Node* > const& vertexesC );*/
	~Cell(void);
	void setCenter( CVCenter* cvCenter );
	void setVertex_N( Node *vertex );
	void setVertex_S( Node *vertex );
	void setVertex_E( Node *vertex );
	void setVertex_W( Node *vertex );
	void setVertex_n( Node *vertex );
	void setVertex_s( Node *vertex );
	void setVertex_e( Node *vertex );
	void setVertex_w( Node *vertex );
	void setVertex_ne( Node *vertex );
	void setVertex_nw( Node *vertex );
	void setVertex_se( Node *vertex );
	void setVertex_sw( Node *vertex );
	void setBoundaryTypes( boost::numeric::ublas::vector<int> boundaryTypes );
	CVCenter* getCenter();
	Node* getVertex_N();
	Node* getVertex_S();
	Node* getVertex_E();
	Node* getVertex_W();
	Node* getVertex_n();
	Node* getVertex_s();
	Node* getVertex_e();
	Node* getVertex_w();
	Node* getVertex_ne();
	Node* getVertex_nw();
	Node* getVertex_se();
	Node* getVertex_sw();
	boost::numeric::ublas::vector<int> getBoundaryTypes();
	int getBoundaryTypes( int index );
};

