


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

using namespace std;


void print_mat(double **mat);
void gauss_jordan_elimination(double **mat);
bool is_triangular(double** mat);
double* solve_triangular (double **mat);
double **copy_mat(double **mat);


int print_precision_mat = 2;
int print_precision_other = 5;

unsigned int N = 33;
double *b;




int main() {


	// definisco la matrice e il vettore dei termini noti

	cout << endl << "Creazione matrici e termini noti: ";

	double **matrix = new double*[N];
	for (int i = 0; i < N; i++)	matrix[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			matrix[i][j] = rand() % 10;

	cout << endl << "matrice iniziale:" << endl;
	print_mat(matrix);
	cout << endl;


	b = new double[N];
	double *b_1 = new double[N];

	for (int i = 0; i < N; i++) {
		b[i] = 1.;
		b_1[i] = b[i];
 	}

	// faccio una copia della matrice per poi calcolare l'errore
	double **matrix_copy = copy_mat(matrix);

	cout << "OK" << endl;


	// riduco la matrice
	cout << "Riduzione matrice principale: ";
	gauss_jordan_elimination(matrix);
	cout << "OK" << endl;

	cout << "matrice ridotta:" << endl;
	print_mat(matrix);
	cout << endl;

	// calcolo le soluzioni
	cout << "Calcolo delle soluzioni: ";
	double *x = solve_triangular(matrix);
	cout << "OK" << endl;


	// calcolo l'errore 
	// il calcolo viene fatto risostituendo il vettore delle soluzioni nella matrice iniziale, 
	// l'errore è quindi quotato come la massima differenza in modulo tra il nuovo vettore dei termini
	// noti calcolaro come matrix * x e quello inziale 

	cout << "Calcolo dell'errore: ";
	double *b_2 = new double[N];


	for (int i = 0; i < N; i++) {
		
		b_2[i] = 0.;

		for(int j = 0; j < N; j++)
			b_2[i] += matrix_copy[i][j] * x[j];

		b_2[i] = abs(b_2[i] - b_1[i]);

	}

	cout << "OK\t\t";

	double error = *max_element(b_2, &b_2[N-1]);
	cout << endl << "   Errore massimo: " << error << endl;
	cout << "soluzioni:" << endl;
	for (int i = 0; i < N; i++)
		cout << "  " << x[i];
	cout << endl;


	cout << "termini noti:" << endl;
	for (int i = 0; i < N; i++)
		cout << "  " << b[i];
	cout << endl;
	

	cout << "Distruzione delle matrici: ";
	for (int i = 0; i < N; i++) delete[] matrix[i];
	for (int i = 0; i < N; i++) delete[] matrix_copy[i];
	delete[] matrix;
	delete[] matrix_copy;
	delete[] x;
	delete[] b;
	delete[] b_1;
	delete[] b_2;
	cout << "OK" << endl << endl; 

	return 0;
}






void print_mat(double **mat) {
	
	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << " " << mat[i][j];
		}
		cout << endl;
	}

	cout << setprecision(print_precision_other);

}

void gauss_jordan_elimination(double **mat) {

	double* tmp;
	double tmp_double, scale, max_el;
	int max_index;
	double *column = new double[N];

	for (int i = 0; i < N-1; i++) {
		
		// PIVOTING


		for (int k = 0; k < N-i; k++)	column[k] = abs(mat[k+i][i]);


		// trovo l'indice del pivot con elemento massimo in modulo e lo metto per primo
		
		//max_index = distance(column, max_element(column, &column[N-i-1])) + i;
		max_el = 0.;
		max_index = 0;

		for (int l = 0; l < N-i; l++) {

			if (abs(column[l]) > max_el) {
				max_el = column[l];
				max_index = l;
			}
		}

		max_index += i;

		

		tmp = mat[i];
		mat[i] = mat[max_index];
		mat[max_index] = tmp;

		tmp_double = b[i];
		b[i] = b[max_index];
		b[max_index] = tmp_double;





		// eseguo la sottrazione
		for (int k = i+1; k < N; k++) {

			if (mat[i][i]==0)	break;

			scale = mat[k][i] / mat[i][i];
			
			b[k] = (b[k] - b[i] * scale);


			for (int j = i+1; j < N; j++) 
				mat[k][j] = (mat[k][j] - mat[i][j] * scale);

			mat[k][i] = 0;

		}

	}

	delete[] column;

}



bool is_triangular(double** mat) {

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




double* solve_triangular (double **mat) {

	// data una matrice triangolare e un vettore di termini noti la funzione restituisce un puntatore 
	// al vettore di lunghezza N contenente le soluzioni trovate con la back-substitution

	if (!is_triangular(mat)) {
		cout << "La matrice non è triangolare: return NULL" << endl;
		return NULL;
	}

	double *x = new double[N];
	
	for (int i = N-1; i >= 0; i--) {
	
		x[i] = b[i];
	
		for (int j = i+1; j < N ; j++)
			x[i] -=  mat[i][j] * x[j];
	
		x[i] /= mat[i][i];
	
	}

	return x;

}


double **copy_mat(double **mat) {

	double **copy = new double*[N];
	for (int i = 0; i < N; i++)	copy[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			copy[i][j] = mat[i][j];

	return copy;

}
