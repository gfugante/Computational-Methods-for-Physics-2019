

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

typedef float Dato_t;

void print_mat(Dato_t **mat);
void gauss_jordan_elimination(Dato_t **mat);
bool is_triangular(Dato_t** mat);
Dato_t* solve_triangular (Dato_t **mat);
Dato_t **copy_mat(Dato_t **mat);


int print_precision_mat = 2;
int print_precision_other = 5;

unsigned int N = 10;

int range_mat = 1E6;
Dato_t range_b = 1.;

Dato_t *b;


bool print = false;


int main() {

	ofstream out("dati/CPU_double_mat1E6_b1.txt");

	cout << "-----------------------------------------" << endl;
	cout << "N" << "\t" << "time" << "\t\t" << "error" << endl;
	cout << "-----------------------------------------" << endl;


	Dato_t **matrix, **matrix_copy;
	Dato_t *b_1, *b_2, *x;
	Dato_t error;
	int incremento = 100;

	clock_t time;
	Dato_t time_total = 0.;


	while(time_total < 60*60) {
		
		time = clock();
		// creo la matrice random e una sua copia per il calcolo dell'errore
		
		matrix = new Dato_t*[N];
	
		for (int i = 0; i < N; i++)	matrix[i] = new Dato_t[N];

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				matrix[i][j] = rand() % range_mat;

		matrix_copy = copy_mat(matrix);


		// creo i vettori dei termini noti necessari e li riempio
		
		b = new Dato_t[N];
		b_1 = new Dato_t[N];
		b_2 = new Dato_t[N];

		for (int i = 0; i < N; i++) {
			b[i] = range_b;
			b_1[i] = b[i];
		}


		// riduco la matrice
		
		gauss_jordan_elimination(matrix);

		// trovo le soluzioni

		x = solve_triangular(matrix);


		// calcolo l'errore
		
		for (int i = 0; i < N; i++) {
		
			b_2[i] = 0.;

			for(int j = 0; j < N; j++)
				b_2[i] += matrix_copy[i][j] * x[j];

			b_2[i] = abs(b_2[i] - b_1[i]);

		}


		error = *max_element(b_2, &b_2[N-1]);

		time = clock() - time;
		time_total += Dato_t(time)/CLOCKS_PER_SEC;

		cout << N << "\t" << scientific << Dato_t(time)/CLOCKS_PER_SEC  << "\t" << error << fixed <<endl;
		 out << N << "\t" << scientific << Dato_t(time)/CLOCKS_PER_SEC  << "\t" << error << fixed <<endl;

		N+=incremento;

	}

	N-=incremento;


	cout << endl << "Distruzione delle matrici: ";
	for (int i = 0; i < N; i++) delete[] matrix[i];
	for (int i = 0; i < N; i++) delete[] matrix_copy[i];
	delete[] matrix;
	delete[] matrix_copy;
	delete[] x;
	delete[] b;
	delete[] b_1;
	delete[] b_2;
	cout << "OK" << endl; 



	cout << endl;
	return 0;
}







void print_mat(Dato_t **mat) {
	
	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << " " << mat[i][j];
		}
		cout << endl;
	}

	cout << setprecision(print_precision_other);

}


void gauss_jordan_elimination(Dato_t **mat) {

	Dato_t* tmp;
	Dato_t tmp_Dato_t, scale, max_el;
	int max_index;
	Dato_t *column = new Dato_t[N];

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

		tmp_Dato_t = b[i];
		b[i] = b[max_index];
		b[max_index] = tmp_Dato_t;





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



bool is_triangular(Dato_t** mat) {

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




Dato_t* solve_triangular (Dato_t **mat) {

	// data una matrice triangolare e un vettore di termini noti la funzione restituisce un puntatore 
	// al vettore di lunghezza N contenente le soluzioni trovate con la back-substitution

	if (!is_triangular(mat)) {
		cout << "La matrice non è triangolare: return NULL" << endl;
		return NULL;
	}

	Dato_t *x = new Dato_t[N];
	
	for (int i = N-1; i >= 0; i--) {
	
		x[i] = b[i];
	
		for (int j = i+1; j < N ; j++)
			x[i] -=  mat[i][j] * x[j];
	
		x[i] /= mat[i][i];
	
	}

	return x;

}


Dato_t **copy_mat(Dato_t **mat) {

	Dato_t **copy = new Dato_t*[N];
	for (int i = 0; i < N; i++)	copy[i] = new Dato_t[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			copy[i][j] = mat[i][j];

	return copy;

}

