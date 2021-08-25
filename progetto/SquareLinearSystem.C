
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <fstream>
#include "SquareLinearSystem.h"

using namespace std;


int print_precision_mat = 2;
int print_precision_other = 5;


SquareLinearSystem::SquareLinearSystem(unsigned int N_) 
// creo la matrice senza riempirla
// il vettore dei termini noti è posto a 1 di default
{

	N = N_;
	
	b = new double[N];
	b_copy = new double[N];

	x = new double[N];

	matrix = new double[N];
	for (int i = 0; i < N; i++)	matrix[i] = new double[N];

	for (int i = 0; i < N; i++)	{
		b[i] = 1.;
		b_copy[i] = b[i];
	}

}

SquareLinearSystem::SquareLinearSystem(unsigned int N_, int a) 
// creo la matrice riempita con numeri random da 0 ad a
// il vettore dei termini noti è posto a 1 di default
{

	N = N_;
	
	b = new double[N];
	b_copy = new double[N];

	x = new double[N];
	
	matrix = new double[N];
	for (int i = 0; i < N; i++)	matrix[i] = new double[N];

	for (int i = 0; i < N; i++)	{
		b[i] = 1.;
		b_copy[i] = b[i];
	}

	for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				matrix[i][j] = rand() % a;

}


SquareLinearSystem::SquareLinearSystem(double** mat, double* b_)
// copio i puntatori di partenza in quelli dell'oggetto, 
// la matrice e i vettori iniziali saranno quindi spostati dentro l'oggetto
{
	
	b = b_;
	for (int i = 0; i < N; i++)	b_copy[i] = b[i];

	matrix = mat;
	for (int i = 0; i < N; i++)	matrix[i] = mat[i];	


	// per evitare manipolazioni dall'esterno cancello i puntatori esterni
	for (int i = 0; i < N; i++)	delete[] mat[i];
	delete[] mat;
	delete[] b_;

}

SquareLinearSystem::~SquareLinearSystem()
{

	for (int i = 0; i < N; i++)	{
		delete[] matrix[i];
		delete[] triangular_matrix[i];
	}

	delete[] b;
	delete[] b_copy;
	delete[] x;

}

void set_matrix(const double** mat)
// copia gli elementi della matrice di partenza in quelli dell'oggetto
// tutte e due le matrici vengono mantenute singolarmente
{

	for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				matrix[i][j] = mat[i][j];

}

void set_known_terms(const double* b_)
// copia i termini noti singolarmente
// tutti e due i vettori vengono mantenuti
{
	
	for (int i = 0; i < N; i++){
		b[i] = b_[i];
		b_copy = b[i];	
	}

}

void set_soliutions(const double* x_)
// copia le soluzioni singolarmente
// tutti e due i vettori vengono mantenuti
{

	for (int i = 0; i < N; i++)
		x[i] = x_[i];	

}

void print_mat() const
{

	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << " " << mat[i][j];
		}
		cout << endl;
	}

	cout << setprecision(print_precision_other);

}

void print_known_terms() const
{

	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) 
		cout << " " << b[i];
	
	cout << endl;
	cout << setprecision(print_precision_other);

}

void print_solutions() const
{

	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) 
		cout << " " << x[i];
	
	cout << endl;
	cout << setprecision(print_precision_other);

}

double** copy_mat() const
{

	double **copy = new double*[N];
	for (int i = 0; i < N; i++)	copy[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			copy[i][j] = mat[i][j];

	return copy;

}

void gauss_jordan_elimination()
// riduce una qualunque matrice a triangolare con l'algoritmo di Gauss-Jordan
{

	clock_t time = clock();

	triangular_matrix = copy_mat(matrix);
	double* tmp;

	double tmp_double, scale;

	double *column = new double[N];
	int max_index;

	bool exit = false;


	for (int i = 0; i < N-1; i++) {		
		

		// PIVOTING

		for (int k = 0; k < N-i; k++)	column[k] = abs(triangular_matrix[k+i][i]);


		// trovo l'indice del pivot con elemento massimo in modulo e lo metto per primo
		
		max_index = distance(column, max_element(column, &column[N-i-1])) + i;
		

		tmp = triangular_matrix[i];
		triangular_matrix[i] = triangular_matrix[max_index];
		triangular_matrix[max_index] = tmp;

		tmp_double = b[i];
		b[i] = b[max_index];
		b[max_index] = tmp_double;


		// eseguo la sottrazione
		for (int k = i+1; k < N; k++) {

			if (triangular_matrix[i][i]==0)	break;

			scale = triangular_matrix[k][i] / triangular_matrix[i][i];


			b[k] = 1E-15*trunc(1E15*(b[k] - b[i] * scale));

			for (int j = i+1; j < N; j++)
				triangular_matrix[k][j] = 1E-15*trunc(1E15*(triangular_matrix[k][j] - triangular_matrix[i][j] * scale));

			triangular_matrix[k][i] = 0;

		}


		// azzero i contari per ripartire
		exit = false;

	}

	delete[] column;
	delete[] tmp;

	time = clock() - time;
	time_to_triangularize = double(time)/CLOCKS_PER_SEC;

}

bool is_triangular(const double **mat) const
{

	// cerco l'indice del primo elemento non nullo della primariga
	int index;
	for (index = 0; index < N; index++)
		if (mat[0][index] != 0) break;



	for (int i = 0; i < N; i++)
		for (int j = 0; j < index+i && j < N; j ++)
			if (mat[i][j] != 0)
				return false;

	return true;

}

void solve_triangular()
// calcola le soluzioni di un sistema triangolare con la back-substitution
{

	if (!is_triangular(triangular_matrix)) 
		cout << "La matrice non è stata triangolarizzata: le soluzioni non possono essere calcolate" << endl;
	else {

		for (int i = N-1; i >= 0; i--) {
		
			x[i] = b[i];
		
			for (int j = i+1; j < N ; j++)
				x[i] -=  triangular_matrix[i][j] * x[j];
		
			x[i] /= triangular_matrix[i][i];
		
		}

	}

}

void evaluate_error()
{

	double *tmp = new double[N];

	for (int i = 0; i < N; i++) {
		
		tmp[i] = 0.;

		for(int j = 0; j < N; j++)
			tmp[i] += matrix[i][j] * x[j];

		tmp[i] = abs(tmp[i] - b_copy[i]);

	}


	error = *max_element(tmp, &tmp[N-1]);
	delete[] tmp;

}




