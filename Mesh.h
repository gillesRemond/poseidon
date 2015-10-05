#pragma once

#include <boost/numeric/ublas/vector.hpp>

#include <string>

#include "constants.h"
#include "Node.h"
#include "CVCenter.h"
#include "Cell.h"

class Mesh
{
private:
	boost::numeric::ublas::vector< int > m_nbCells;
	boost::numeric::ublas::vector< CVCenter* > m_cvCenters;
	boost::numeric::ublas::vector< Node* > m_vertexesV;
	boost::numeric::ublas::vector< Node* > m_vertexesH;
	boost::numeric::ublas::vector< Node* > m_vertexesC;
	boost::numeric::ublas::vector< Cell* > m_cells;

public:
	Mesh( void );
	Mesh( Parameters const& parameters );
	~Mesh();
	int getNbCells( int index );
	Cell* getCell( int index );
	void initialiseUniformField( double T );
	void setTField( boost::numeric::ublas::vector<double> const& T );
	void prolongateT( Mesh* lowerMesh, boost::numeric::ublas::vector<int> const& prolongator );
	double getT( int index );
	void getT( boost::numeric::ublas::vector<double> &T );
	void generateMesh( Parameters const& parameters );
	void printData();
	void printData( std::string name );

};

