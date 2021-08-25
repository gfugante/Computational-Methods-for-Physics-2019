#include <iostream>
using namespace std;
#include <cstdlib>

#define imin(a,b) (a<b?a:b)
const int N = 100000;

const int threadsPerBlock = 512;
//trucco per ottimizzare il numero di blocchi
const int blocksPerGrid = imin(32,(N+threadsPerBlock-1)/threadsPerBlock);


///////////////////////////////////////////////////////////////////////////////////////////
__global__ void scalar_prod( float *a,  float *b,  float *c){

  // la memoria shared e' interna e comune ai threads di un singolo blocco
  __shared__ float cache[threadsPerBlock];

  int tid = threadIdx.x+blockDim.x*blockIdx.x;

  int cacheIndex = threadIdx.x;
  float temp = 0;

  while(tid < N){
    temp +=a[tid]*b[tid];

    //incremento dell'indice della componente che vi4ene sommata
    tid += blockDim.x*gridDim.x;
  };

  //il numero di elementi da valutare puo' essere maggiore del numero di thread * numero di blocchi
  //la somma (temp) corre sui diversi gruppi di blocchi combinando le componenti associate 
  //allo stesso valore di tid (modulo il numero totale di threads)

  //scopo del prodotto scalare e' di restituire un numero; parte di questo numero viene calcolato sulla GPU
  //la somma parziale appena calcolata viene messa nella memoria shared del blocco in corrispondenza del valore
  // di threadIdx.x
  cache[cacheIndex] = temp;

  //risulta necessario attendere che tutti i thread abbiano terminato di valutare le loro somme
  //prima di passare alla somma delle somme parziali
  __syncthreads();

  //si puo' quindi passare a valutare la somma parziale dei risultati dei vari threads di un blocco

int i=blockDim.x/2;
  while(i!=0){
    if (cacheIndex < i)
      cache[cacheIndex]+=cache[cacheIndex+i];


//       c[0]   c[1]  c[2]  c[3]   c[4]  c[5]  c[6] c[7]

      //questa somma combina i valori della meta' superiore del blocco di threads con quelli della 
      //meta' inferiore e li associa a indici della meta' inferiore


    //e' necessario attendere che tutti i threads abbiano effettuato la somma
    __syncthreads();

    //tutti gli elementi rilevanti hanno indici da zero a meta' del blocco
    //nella seconda iterazione si sommano gli elementi del secondo quarto con quelli del primo
    //quarto del blocco, e cosi' via
    i/=2;

    //poiche' i e' una variabile intera, nel caso 1/2 = 0.5 viene restituito 0
  };

  if (cacheIndex ==0)
    //il risultato finale e' la somma parziale di un blocco
    c[blockIdx.x] = cache[0];

}

///////////////////////////////////////////////////////////////////////////////////////////////////


void sp_cpu(float *a, float *b, float &check){

  int i;
  check=0;
  for (i=0; i<N; i++){
    check+=a[i]*b[i];

  }

}





int main( void ){

  float *a, *b, *partial_c;
  float *dev_a, *dev_b, *dev_partial_c;
  float  sp, sp2;

  a = (float *)malloc(N*sizeof(float));
  b = (float *)malloc(N*sizeof(float));
  partial_c = (float *)malloc(blocksPerGrid*sizeof(float));

  cudaEvent_t start1,stop1;
  float time1;

  for( int i=0; i<N; i++){
    a[i] = drand48();
    b[i] = drand48();
  };

  cudaEventCreate(&start1);
  cudaEventCreate(&stop1);
  

  cudaMalloc( (void **)&dev_a, N*sizeof(float)  );
  cudaMalloc( (void **)&dev_b, N*sizeof(float)  );

  //la memoria allocata serve per recuperare le somme parziali associate a
  //ciascun blocco
  cudaMalloc( (void **)&dev_partial_c, blocksPerGrid*sizeof(float)  );


  // all'indirizzo di c viene copiato 
  // quanto si trova sul device a partire da dev_c
  // per una lunghezza pari a sizeof(int)

  cudaEventRecord(start1,0);

  
  cudaMemcpy( dev_a, a, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_b, b, N*sizeof(float), cudaMemcpyHostToDevice);




// il kernel viene invocato, la GPU viene suddivisa in 32 blocchi
// in ciascun blocco lavorano 512 thread
 scalar_prod<<<blocksPerGrid,threadsPerBlock>>>(dev_a, dev_b, dev_partial_c);



// all'indirizzo di partial_c viene copiato 
// quanto si trova sul device a partire da dev_partial_c
// per una lunghezza pari a blocksPerGrid*sizeof(float)

  cudaMemcpy( partial_c, dev_partial_c, blocksPerGrid*sizeof(float), cudaMemcpyDeviceToHost);

  sp=0;
  //infine le somme parziali dei vari blocchi vengono combinate 
   for( int i=0; i<blocksPerGrid; i++){
     sp += partial_c[i];
   };
     cout << "scalar product gpu=" << sp <<"\n";


  cudaFree( dev_a);
  cudaFree( dev_b);
  cudaFree( dev_partial_c);



  cudaEventRecord(stop1,0);
  cudaEventSynchronize(stop1);
  cudaEventElapsedTime(&time1, start1, stop1);
  cout << "time1=" << time1 <<"\n";

  cudaEventDestroy(start1);
  cudaEventDestroy(stop1);


  ///////////////////////////////////////////
  cudaEventCreate(&start1);
  cudaEventCreate(&stop1);
  cudaEventRecord(start1,0);

  //check del prodotto scalare, calcolato sulla CPU
  sp_cpu(a,b,sp2);

  cudaEventRecord(stop1,0);
  cudaEventSynchronize(stop1);
  cudaEventElapsedTime(&time1, start1, stop1);
     cout << "scalar product cpu=" << sp2 <<"\n";

  cout << "time1=" << time1 <<"\n";

  cudaEventDestroy(start1);
  cudaEventDestroy(stop1);



  free(a);
  free(b);
  free(partial_c);


  return 0;
}

