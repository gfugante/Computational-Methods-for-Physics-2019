
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;



void print_mat(double **mat, int N);
void gauss_jordan_elimination(double **mat, int N, bool print=false);
bool is_diagonal(double** mat, int N);





void print_mat(double **mat, int N) {

	int print_precision_mat = 2;
	int print_precision_other = 5;
	
	cout << setprecision(print_precision_mat);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << " " << mat[i][j];
		}
		cout << endl;
	}

	cout << setprecision(print_precision_other);

}









void gauss_jordan_elimination(double **mat, int N, bool print) {

	double* tmp;
	double tmp_double;
	bool exit = false;
	double scale;

	for (int i = 0; i < N-1; i++) {
		
		if(print) cout << endl;		
		

		// PIVOTING
		
		double *column = new double[N-i];
		int max_index;

		for (int k = 0; k < N-i; k++)	column[k] = abs(mat[k+i][i]);


		// trovo l'indice del pivot con elemento massimo in modulo e lo metto per primo
		
		// cout << "max index: " << distance(column, max_element(column, &column[N])) << " + " << i << " = "  << distance(column, max_element(column, &column[N-i-1])) + i << endl;
		//max_index = distance(column, max_element(column, &column[N-i-1])) + i;
		//delete[] column;
		

		tmp = mat[i];
		mat[i] = mat[max_index];
		mat[max_index] = tmp;

		//tmp_double = b_1[i];
		//b_1[i] = b_1[max_index];
		//b_1[max_index] = tmp_double;


		if(print) {

			cout << "scambi" << endl << endl;
			cout << "  " << i << " <-> " << max_index << endl;
			// SYNCTHREADS! max_element usa delle if!
	

			// ho trovato il pivot e l'ho messo nella riga giusta per cominciare a fare la sottrazione
			cout << endl;
			print_mat(mat, N);
			cout << endl;
			cout << "sottrazioni" << endl << endl;

		}




		// eseguo la sottrazione
		for (int k = i+1; k < N; k++) {

			if (mat[i][i]==0)	break;

			scale = mat[k][i] / mat[i][i];
			//b_1[k] = 1E-15*trunc(1E15*(b_1[k] - b_1[i] * scale));


			for (int j = N-1; j >= i; j--) {

				if (print)	cout << "  " << mat[k][j]  << " - " << mat[i][j] << " * " << mat[k][i] << " / " << mat[i][i] << " = " << mat[k][j] - mat[i][j] * scale << endl;
				mat[k][j] = 1E-15*trunc(1E15*(mat[k][j] - mat[i][j] * scale));

			}
			cout << endl;
		}


		if (print) {
			cout << endl;
			print_mat(mat, N);
			cout << endl;
			cout << "-------------------------------------------------" << endl;
		}


		// azzero i contari per ripartire
		exit = false;

	}

}









bool is_diagonal(double** mat, int N) {

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
