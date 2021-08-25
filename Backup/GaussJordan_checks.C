


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

using namespace std;


void print_mat(double **mat);
void gauss_jordan_elimination(double **mat);
bool is_diagonal(double** mat);


int print_precision_mat = 2;
int print_precision_other = 5;

unsigned int N = 10;
double *b_1;
double *b_2;

double scale;

bool print = false;


int main() {


	// definisco la matrice e il vettore x

	double **matrix = new double*[N];
	double *x = new double[N];
	b_1 = new double[N];
	b_2 = new double[N];
	
	for (int i = 0; i < N; i++)	matrix[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			matrix[i][j] = rand() % 2;

	for (int i = 0; i < N; i++)	{
		x[i] = 1.;
		b_1[i] = 0.;
		b_2[i] = 0.;
	}

	// moltiplico matrix e x, trovo un vettore soluzione b che dovrò riarrangiare nello stesso modo in cui riarrangio matrix
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
	 		b_1[i] += matrix[i][j]*x[j];
	


	// print matrix
	
	if (print) {

		cout << endl << "-------------------------------------------------" << endl;
		cout << endl << "Matrice iniziale" << endl << endl;
		print_mat(matrix);
		cout << endl << "-------------------------------------------------" << endl;
	
	}

	gauss_jordan_elimination(matrix);

	if (print) {

		cout << endl;
		cout << endl << "-------------------------------------------------" << endl;
		cout << endl << "Matrice ridotta" << endl << endl;
		print_mat(matrix);
		cout << endl << "-------------------------------------------------" << endl;

	}


	if (is_diagonal(matrix)) 	cout << "Matrice ridotta con successo" << endl;
	else 						cout << "L'algoritmo fallisce per N = " << N << endl;
	


	// ri-moltiplico matrix e x e trovo un nuovo vettore di soluzioni, quoto l'errore commesso come il massimo elemento 
	// della differenza tra questo vettore e quello trovato prima di ridurre la matrice
	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++)	 b_2[i] += matrix[i][j]*x[j];
		b_2[i] = 1E-15*trunc(1E15*abs(b_2[i] - b_1[i]));

	}

	double error = *max_element(b_2, &b_2[N]);
	cout << "Errore massimo: " << scientific << error << endl;
	

	for (int i = 0; i < N; i++) delete[] matrix[i];
	delete[] matrix;
	delete[] x;
	delete[] b_1;
	delete[] b_2;

	//cout << endl;
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
	double tmp_double;
	bool exit = false;

	for (int i = 0; i < N-1; i++) {
		
		//cout << endl;		
		

		// PIVOTING
		
		double *column = new double[N-i];
		int max_index;

		for (int k = 0; k < N-i; k++)	column[k] = abs(mat[k+i][i]);


		// trovo l'indice del pivot con elemento massimo in modulo e lo metto per primo
		
		//cout << "max index: " << distance(column, max_element(column, &column[N])) << " + " << i << " = "  << distance(column, max_element(column, &column[N-i-1])) + i << endl;
		max_index = distance(column, max_element(column, &column[N-i-1])) + i;
		delete[] column;
		

		tmp = mat[i];
		mat[i] = mat[max_index];
		mat[max_index] = tmp;

		tmp_double = b_1[i];
		b_1[i] = b_1[max_index];
		b_1[max_index] = tmp_double;



		if (print) {

			cout << "scambi" << endl << endl;
			cout << "  " << i << " <-> " << max_index << endl;
			// SYNCTHREADS! max_element usa delle if!
			
	
			// ho trovato il pivot e l'ho messo nella riga giusta per cominciare a fare la sottrazione
			cout << endl;
			print_mat(mat);
			cout << endl;
			cout << "sottrazioni" << endl << endl;
		
		}




		// eseguo la sottrazione
		for (int k = i+1; k < N; k++) {

			if (mat[i][i]==0)	break;

			scale = mat[k][i] / mat[i][i];
			b_1[k] = 1E-15*trunc(1E15*(b_1[k] - b_1[i] * scale));


			for (int j = i+1; j < N; j++) {

				if (print) cout << "  " << mat[k][j]  << " - " << mat[i][j] << " * " << scale << " = " << mat[k][j] - mat[i][j] * scale << endl;
				mat[k][j] = 1E-15*trunc(1E15*(mat[k][j] - mat[i][j] * scale));

			}
			mat[k][i] = 0;
			if (print) cout << endl;
		}


		if (print) {

			cout << endl;
			print_mat(mat);
			cout << endl;
			cout << "-------------------------------------------------" << endl;

		}
		// azzero i contari per ripartire
		exit = false;

	}

}

bool is_diagonal(double** mat) {

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


double* solve_triangular (double **mat, double *b) {

	// data una matrice triangolare e un vettore di termini noti la funzione restituisce un puntatore 
	// al vettore di lunghezza N contenente le soluzioni trovate con la back-substitution

	if (!is_diagonal(mat)) {
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

