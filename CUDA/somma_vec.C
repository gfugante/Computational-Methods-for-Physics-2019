#include <iostream>
using namespace std;
#include <cstdlib>

#define N 65535

// ogni kernel viene dichiarato con __global__
__global__ void somma_vec( int *a, int *b, int *c){

  int tid = blockIdx.x;

  if(tid < N)
    c[tid]=a[tid]+b[tid];
  
}

int main( void ){

  int a[N],b[N],c[N];
  int *dev_a, *dev_b, *dev_c;

  cudaEvent_t start1,stop1;
  float time1;

  for( int i=0; i<N; i++){
    a[i] = 2*i;
    b[i] = -i;
  };


  // la possibilita' di misurare la durata di un evento viene resa possibile
  // dall'utilizzo di tag, che possono essere associati a due istanti temporali
  // e al relativo utilizzo della GPU
  cudaEventCreate(&start1);
  cudaEventCreate(&stop1);

// la memoria richiesta sul device deve essere allocata:
// dev_c puntatore a intero, sulla scheda viene allocata memoria pari al
// secondo argomento di cudaMalloc associata all'indirizzo di dev_c
// (void **) pointer a pointer
  cudaMalloc( (void **)&dev_a, N*sizeof(int)  );
  cudaMalloc( (void **)&dev_b, N*sizeof(int)  );
  cudaMalloc( (void **)&dev_c, N*sizeof(int)  );


  //da questo istante misuriamo il tempo trascorso sulla GPU
  cudaEventRecord(start1,0);

  //i dati contenuti in a e b vengono copiati dall'HOST al DEVICE
  // e vengono associati a dev_a e a dev_b
  cudaMemcpy( dev_a, a, N*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_b, b, N*sizeof(int), cudaMemcpyHostToDevice);


// il kernel viene invocato  chiedendo che ci siano N blocchi ciascuno 
// contenente un singolo thread
// il risultato viene associato a dev_c
  somma_vec<<<N,1>>>(dev_a, dev_b, dev_c);


// nella variabile c viene copiato dal DEVICE all'HOST
// quanto si trova sul device a partire da dev_c
// per una lunghezza pari a N*sizeof(int)

  cudaMemcpy( c, dev_c, N*sizeof(int), cudaMemcpyDeviceToHost);

  for( int i=0; i<N; i++){
    cout << "c["<< i <<"] = " << a[i] <<"  " << b[i] << " = " <<c[i] <<"\n";
  };


  cudaFree( dev_a);
  cudaFree( dev_b);
  cudaFree( dev_c);

  //il contatore finale viene misurato
  cudaEventRecord(stop1,0);
  //si chiede che la GPU abbia completato le sue operazioni prima di restituire
  // il valore del contatore stop1 
  cudaEventSynchronize(stop1);

  //la differenza tra start1 e stop1 corrisponde al tempo trascorso
  //misurato in millisecondi
  cudaEventElapsedTime(&time1, start1, stop1);
  cout << "time1=" << time1 <<"\n";

  cudaEventDestroy(start1);
  cudaEventDestroy(stop1);

  return 0;
}

