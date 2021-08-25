
#ifndef __SQUARE_LINEAR_SYSTEM__
#define __SQUARE_LINEAR_SYSTEM__

class SquareLinearSystem {

	protected: 
		unsigned int N;							// dimensione
		double **matrix, **triangular_matrix;	// matrice e matrice ridotta
		double *b, *b_copy;						// termini noti
		double *x;								// soluzioni
		double time_to_triangularize;			// tempo per triangolrizzare la matrice (in secondi)
		double error;							// errore massimo nella risoluzione del sistema
	
	public:
		SquareLinearSystem(unsigned int N_);
		SquareLinearSystem(unsigned int N_, int a);
		SquareLinearSystem(double** mat, double* b_);	// attenzione non vengono copiati gli elementi ma vengono acquisiti i puntatori
		~SquareLinearSystem();


		double** get_matrix() 		{return matrix};
		double* get_known_terms() 	{return b};
		double* get_solutions()		{return x};

		void set_matrix(const double** mat);
		void set_known_terms(const double* b_);
		void set_soliutions(const double* x_);

		void print_mat() const;
		void print_known_terms() const;
		void print_solutions() const;

		double** copy_mat() const;

		void gauss_jordan_elimination();
		bool is_triangular(const double **mat) const;
		void solve_triangular();
		void evaluate_error();

};

#endif 

