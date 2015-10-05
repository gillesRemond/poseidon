//main.cpp

#include <vector>
#include <iostream>
#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "constants.h"
#include "Multigrid.h"
#include "Mesh.h"
#include "Equation.h"
#include "time.h"

using namespace std;

long int top_chrono;
long int top_chrono_total;
 
void demarrer_chrono() {
        top_chrono = clock();
}
 
void stop_chrono() {
        long int arret_chrono = clock();
        fprintf(stderr, "Le calcul a pris %f secondes.\n",
                (float)(arret_chrono - top_chrono) / CLOCKS_PER_SEC);
		top_chrono = arret_chrono;
}

void stop_chrono_total() {
        long int arret_chrono = clock();
        fprintf(stderr, "Le calcul a pris %f secondes.\n",
                (float)(arret_chrono - top_chrono_total) / CLOCKS_PER_SEC);
		top_chrono = arret_chrono;
}

int main()
{
	Parameters parameters;
	parameters.lengths.resize( 2 );
	parameters.nbCells.resize( 2 );
	parameters.boundaryTypes.resize( 4 );

	// Ask the user for some parameters
	cout << "Welcome in the poseidon software!" << endl << endl;
	cout << "Some parameters for setting up the mesh :" << endl;
	cout << "The length of the rectangular :" << endl;
	cin >> parameters.lengths[0];
	cout << "The height of the rectangular :" << endl;
	cin >> parameters.lengths[1];
	cout << "Number of cells in the x direction :" << endl;
	cin >> parameters.nbCells[0];
	cout << "Number of cells in the y direction :" << endl;
	cin >> parameters.nbCells[1];
	//double omega;
	//cout << "The relaxation parameter omega :" << endl;
	//cin >> omega;

	parameters.boundaryTypes[ 0 ] = DIRICHLET;
	parameters.boundaryTypes[ 1 ] = DIRICHLET;
	parameters.boundaryTypes[ 2 ] = DIRICHLET;
	parameters.boundaryTypes[ 3 ] = DIRICHLET;

	demarrer_chrono();
	top_chrono_total = top_chrono;

	Multigrid* multigrid = new Multigrid( parameters );
	multigrid -> execute();

	cout << "Total time computation." << endl;
	stop_chrono_total();

	// Free the memory
	delete multigrid;

	system("pause");
	return 0;
}