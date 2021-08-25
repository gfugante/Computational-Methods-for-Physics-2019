#include <iostream> 
#include <cstdlib>
#include <cmath>
#define EPS 1e-6
#include <ctime>
using namespace std;
#define MAX 30
#define THREADS 32
typedef float dato;
//funzione per controllare la coerenza dei risultati cpu-gpu
void check(dato *a, dato *b,int N) {
  dato diff;
	for(int i = 0; i < N; ++i)
		for(int j = 0 ; j < N ; ++j )
			if(fabs(a[i*N+j]/b[i*N+j]-1.)>EPS) {
			  diff = a[i*N+j]-b[i*N+j];
				cerr << "Errore in posizione (" << i << "," << j <<"), a: "
				<< a[i*N+j] << ", check: " << b[i*N+j] << "   " << diff << endl; 
				return;
			}
	cerr << "OK" << endl;
	return;
};

//Versione cpu del prodotto matrice matrice
void matrix_product_cpu(dato *a,dato *b, dato *c,int N){

	dato temp;
	for(int i = 0; i < N; ++i) 
		for(int j = 0; j < N; ++j) {
			temp = 0; 
			for(int k = 0; k < N;++k) 
				temp += a[i*N+k]*b[k*N+j];	
			c[i*N+j] = temp;
		}
	return;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
//Prodotto matriciale: versione ottimizzata. 
//Si considerano matrici quadrate di ordine N  multiplo della dimensione del blocco di threads THREADS (quadrato).
//E' possibile ottimizzare il calcolo di un gruppo di righe per un gruppo di colonne
//specificamente di un blocco di elementi di A di dimensione THREADSxN e un blocco di B di dimensione Nx THREADS. 
//Per ridurre il numero di accessi alla memoria e' conveniente sfruttare la memoria condivisa della scheda grafica.
//La memoria condivisa e' presente in quantita' limitata e quindi diventa necessario ridurre il problema in tanti sottoproblemi.
//Piu' precisamente, il sottoproblema elementare e' costituito da un prodotto di due matrici di dimensione THREADSxTHREADS
//memorizzate sulla memoria veloce (shared) della GPU

__global__ void matrix_product( dato* A, dato* B, dato* C,int width) {



    int ib = blockIdx.y;
    int jb = blockIdx.x;	        
    int it = threadIdx.y;
    int jt = threadIdx.x;
        
    //Indice della prima sottomatrice di A elaborata dal blocco
    // width e' un multiplo intero di THREADS
    // aBegin  include un certo numero ib  di gruppi di blocchi rettangolari THREADSxwidth
    int aBegin = width*THREADS*ib;
    
    //Indice dell'ultima sottomatrice di A elaborata dal blocco
    //
    int aEnd   = aBegin + width - 1;
    
    // numero di colonne tra una sottomatrice e la successiva
    int aStep  = THREADS;
    
    
    //indice della prima sottomatrice di B elaborata dal blocco
    //bBegin include un certo numero jb di blocchi di colonne, blocchi larghi THREADS 
    int bBegin = THREADS*jb;
         
    // numero di elementi tra una sottomatrice e la successiva    
    int bStep  = THREADS * width;
     

    // Csub e' usata come variabile in cui memorizzare il valore dell'elemento di C
    // calcolato dal thread
    // Viene aggiornato ripetutamente nel ciclo for seguente
    dato Csub = 0;    
    
    
    //Iterazione sulle sottomatrici in cui viene suddiviso il calcolo degli elementi del blocco 
    for (int a = aBegin, b = bBegin; a <= aEnd; a += aStep, b += bStep) {

        //Dichiarazione della variabile in cui salvare la sottomatrice di A in esame
	__shared__ dato As[THREADS][THREADS];
        
        __shared__ dato Bs[THREADS][THREADS];

        //Carico gli elementi di ciascuna sottomatrice in memoria condivisa: ogni thread del 
	//blocco carica un elemento! 	      
	As[it][jt] = A[a +  width*it + jt]; 
	Bs[it][jt] = B[b +  width*it + jt];
	
        //Sincronizzo per assicurarmi che ogni thread del blocco abbia caricato gli elementi.
	__syncthreads();

 	//Calcolo i contributi agli elementi di matrice di C dovute alla sottomatrici in esame        
	for (int k = 0; k < THREADS ; ++k ) 
        	Csub += As[it][k]*Bs[k][jt];
	//l'elemento C[it][jt] viene aggiornato in nblocks passaggi, pari al numero di iterazioni
	//del ciclo for

   			
        //Sincronizzo per assicurarmi che il calcolo precedente sia terminato prima di caricare nuove
	//sottomatrici
	__syncthreads();
    }

    //Scrivo i risultati in C. Ogni thread elabora un elemento di C. 
    int c = width*THREADS*ib + THREADS*jb;              
    C[c +  width*it + jt] = Csub;
 
    
};

//crea cornice di zeri
__global__ void completa_matrice(dato *a, int N) {

	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int offset = gridDim.x*blockDim.x*i + j;
	
	if(i >= N || j >= N) 
		a[offset] = 0.;


};

///////////////////////////////////////////////////////////////////////////////////
void matrix(dato *a,dato *b,dato *c,int N) {
	
	dato *a_dev,*b_dev,*c_dev;	
	int nblock = (N+THREADS-1)/THREADS;	
	int width = nblock*THREADS;

	
	cout << "width=" << width << "  dim matrix=" << N << "\n";

	cudaMalloc((void**)&a_dev,width*width*sizeof(dato));
	cudaMalloc((void**)&b_dev,width*width*sizeof(dato));
	cudaMalloc((void**)&c_dev,width*width*sizeof(dato));
	
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
		
	cudaMemcpy2D(a_dev,width*sizeof(dato),a,N*sizeof(dato),N*sizeof(dato),N,cudaMemcpyHostToDevice);
	cudaMemcpy2D(b_dev,width*sizeof(dato),b,N*sizeof(dato),N*sizeof(dato),N,cudaMemcpyHostToDevice);

	dim3 c_threads(THREADS,THREADS);
	dim3 c_blocks(nblock,nblock);
		
	dim3 threads(THREADS,THREADS);
	dim3 blocks(nblock,nblock);
	
	if(N % THREADS != 0) {
		completa_matrice<<<c_blocks,c_threads>>>(a_dev,N);
	        completa_matrice<<<c_blocks,c_threads>>>(b_dev,N);
	}
	
	matrix_product<<<blocks,threads>>>(a_dev,b_dev,c_dev,width);
	
	
	//copio il risultato in c
	cudaMemcpy2D(c,N*sizeof(dato),c_dev,width*sizeof(dato),N*sizeof(dato),N,cudaMemcpyDeviceToHost);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapsed;
	cudaEventElapsedTime(&elapsed,start,stop);
	cout << N << " " << elapsed/1000. << endl;	
	
	//dealloco 
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(a_dev);
	cudaFree(b_dev);
	cudaFree(c_dev);			
};



int main(int argc, char **argv) {

	cudaSetDevice(0);
	const int N = atoi(argv[1]);	

	dato *a_host = new dato[N*N];
	dato *b_host = new dato[N*N];
	dato *c_host = new dato[N*N];
	//	dim matrix=100;
	dato *c_check = new dato[N*N];		
	
	srand48(time(NULL));
	
	
	//riempio matrici a e b
	for(int i = 0; i < N ;++i)
		for(int j = 0; j < N;++j) {
			int index = i*N+j;
			a_host[index] = drand48()*MAX;
			b_host[index] = drand48()*MAX;	
		}
	
	float elapsed_cpu;		
	matrix(a_host,b_host,c_host,N);
	
	clock_t start,stop;	
	start = clock();
	matrix_product_cpu(a_host,b_host,c_check,N);
	stop = clock();
	elapsed_cpu = (float)(stop - start)/CLOCKS_PER_SEC;
	cout << elapsed_cpu << endl;
	
		check(c_host,c_check,N);
	delete[](a_host);
	delete[](b_host);
	delete[](c_host);
	delete[](c_check);
	
	return 0;
};
