#pragma once

#include <boost/numeric/ublas/vector.hpp>

#include "constants.h"
#include "Mesh.h"


class Multigrid
{
private:
	int m_nbGrids;
	boost::numeric::ublas::vector<Mesh*> m_grids;

public:
	Multigrid(void);
	Multigrid( Parameters parameters );
	~Multigrid(void);
	void execute( void );
};

