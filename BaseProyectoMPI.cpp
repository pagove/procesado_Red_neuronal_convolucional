 
/****************** BaseProyecto.cpp ********************/ /**
 *
 * @File: BaseProyecto.cpp
 *
 * @Brief Este programa es el programa base para las practicas de ARC
 *        donde intentaremos paralelizar el siguiente codigo.
 *                                             
 * @Author: Ivan Gregoria Bou
 * @Author: Pablo Gomar Velló
 * @date 14/02/2020                               
 * 
 */ /*********************************************************/
 
#include <iostream>
#include <vector> 
#include <math.h>
#include<algorithm>
#include <sys/time.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <mpi.h>

using namespace std;

const int REPS = 10;

int main (int argc, char* args[]){
    
            
    /**
     * Declaracion de varibales.
     */
    string caracter = args[2];
        long n = atoi(args[1]);
    long m = 5; 
    long s = 20;    
    float *a = 0;  
    float *c = new float[s*m*m];
    float *m1 = 0; 
    float aux = 0; 
    double long equacion;
    struct timeval ini, fin;
    double tiempo;
    ofstream fich_out("ResultadosTiempo.txt", ios::app);
    ofstream fich_out1("ResultadosFlops.txt", ios::app);
    
    int rank, size;
    long n_local, resto;
    int *partes = 0;
    int  *desplazamiento = 0;
    int a_parcial;
    
    MPI_Init(&argc, &args);
    
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    n_local = n / size;
    resto = n % size;
    /*
     * =======================================================================        INICIO PROGRAMA                     ========================================================================
     */
    
    
    if (fich_out.is_open() && fich_out1.is_open()){
        
        if(rank == 0){
                        a = new float[n*n];
                        partes = new int[size];
                        desplazamiento = new int[size];
                        m1 = new float[(n/2)*(n/2)];
                        //Creacion de la matriz principal (imagen).
                for(int x = 0; x < n; x++){
                    for(int y = 0; y < n; y++){
                        a[(x*n) + y] = ( ((x*y)+x+y) / (3* pow(n,2)) ) * 1e3;
                    }
                }
                
                
                //Creacion de la bateria de filtros.
                srand(5); // Se inicializa la semilla siempre con el numero 5 no con el tiempo del pc
                for (int k = 0; k < s; k++){
                    for(int i = 0; i < m; i++){
                        for(int j = 0; j < m; j++){
                            c[(k*m*m) + (i*m) + j] = (rand() * pow(-1, j+k))/ RAND_MAX;
                        }
                    }
                }


        
        partes[0] = (n_local+resto+2)*n;
        desplazamiento[0] = 0;
        
        for(int i = 1; i < (size - 1); i++){
                partes[i] = (n_local+4)*n;
                desplazamiento[i] = desplazamiento[i-1] + partes[i-i] - 4 * n;
                }

                partes[size-1] = (n_local+2)*n;
        desplazamiento[size-1] = desplazamiento[size-2] + partes[size -2] -4 * n;
        
                        gettimeofday(&ini, NULL);

                }//fin rank0
        
        /*
        * =======================================================================        INICIO CALCULO PROGRAMA               ========================================================================
        */
        
        if(( rank == 0) || (rank == size -1)){
                a_parcial = 2;
                } else {
                        a_parcial = 4;
                }

                float *a_local = new float[(n_local + a_parcial) * n];
        float *b_local = new float[s*n_local*n];
        float *r_local = new float[s*(n_local/2)*(n/2)];
        float *m1_local = new float[(n_local/2)*(n/2)];

        
        for (int t = 0; t < REPS; t++){
            aux = 0;
            
            for(int i = 0; i < s*n_local*n; i++){
                b_local[i] = 0;
            }
            
            MPI_Bcast(c, s*m*m, MPI_FLOAT, 0, MPI_COMM_WORLD);
            MPI_Scatterv(a, partes, desplazamiento, MPI_FLOAT, a_local, (n_local+a_parcial)*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
            
            if (rank == 0){
                    for(int z = 0; z < s; z++){ // parte izquierda
                        for(int x = 2; x < n_local;  x++){
                            for(int y = 2; y <= (n - 3); y++){
                                b_local[(z*n_local*n) + (x*n) + y] = 0;
                                for(int i = -m/2; i <= m/2; i++){
                                    for(int j = -m/2; j <= m/2; j++){                  
                                        b_local[(z*n_local*n) + (x*n) + y] += a_local[(n* (x + i)) + y + j] * c[(m*m*z) + (m*((m/2) + i)) + ((m/2) + j)];
                                    }
                                }
                                
                            }
                        }
                    }

                        } else if (rank == (size - 1)){ // parte derecha
                                for(int z = 0; z < s; z++){
                        for(int x = 0; x <= (n_local - 3);  x++){
                            for(int y = 2; y <= (n - 3); y++){
                                b_local[(z*n_local*n) + (x*n) + y] = 0;
                                for(int i = -m/2; i <= m/2; i++){
                                    for(int j = -m/2; j <= m/2; j++){                  
                                        b_local[(z*n_local*n) + (x*n) + y] += a_local[(n* (x + i + 2)) + y + j] * c[(m*m*z) + (m*((m/2) + i)) + ((m/2) + j)];
                                    }
                                }
                                
                            }
                        }
                    }
                        } else { // cdntro matriz
                                //Calculo de la convolucion de la imagen.
                    for(int z = 0; z < s; z++){
                        for(int x = 0; x < n_local;  x++){
                            for(int y = 2; y <= (n - 3); y++){
                                b_local[(z*n_local*n) + (x*n) + y] = 0;
                                for(int i = -m/2; i <= m/2; i++){
                                    for(int j = -m/2; j <= m/2; j++){                  
                                        b_local[(z*n_local*n) + (x*n) + y] += a_local[(n* (x + i + 2)) + y + j] * c[(m*m*z) + (m*((m/2) + i)) + ((m/2) + j)];
                                    }
                                }
                                
                            }
                        }
                    }

                        }
            
            
            
            //Aplicacion de la funcion no lineal.
            for(int z = 0; z < s; z++){
                for(int x = 0; x < n_local; x++ ){
                    for(int y = 0; y < n; y++){
                        b_local[(z*n_local*n) + (n*x) + y] = 1.0/(1 + exp(-b_local[(z*n_local*n) + (n*x) + y]));
                    }
                }
            }
            
            //Pooling
            for(int z = 0; z < s; z++){
                for(int x = 0; x < n_local/2; x++){
                    for(int y = 0; y < n/2; y++){
                        r_local[(z*(n_local/2)*(n/2)) + (x*(n/2)) + y] = fmax( fmax(b_local[(z*n_local*n)+ (n*2*x) + 2*y], b_local[(z*n_local*n) + (n*2*x) + ((2*y) +1)]), fmax(b_local[(z*n_local*n) + (((2*x) +1)*n) + (2*y)], (b_local[(z*n_local*n) + (((2*x)+1)*n) + ((2*y)+1)])));
                    }
                }
            }
            
            //promedio.
            for(int x = 0; x < n_local/2; x++){
                for(int y = 0; y < n/2; y++){
                    aux = 0;
                    for(int k = 0; k < s; k++){
                        aux = aux + r_local[(k*(n_local/2)* (n/2)) + ((n/2 )*x) + y];
                    }
                    m1_local[(x*(n_local/2)) + y] = aux/s ;
                }
            }
              
        } // bucle de 10 reps para la medida del tiempo
        
                MPI_Gather(m1_local, (n/2)*(n_local/2), MPI_FLOAT, m1, (n/2)*(n_local/2), MPI_FLOAT, 0, MPI_COMM_WORLD);


                if (rank == 0){
                        gettimeofday(&fin, NULL);
    
    
    /*
    * =======================================================================        SALIDA PROGRAMA                     ========
    */

                if (caracter == "d"){   
                    cout << "Matriz A \t Matriz M" << endl;
                    
                    for(int i = 0; i < 4; i++){
                        cout << fixed << setprecision(4) << a[(i*n) + i] << "\t         ";
                        cout << fixed << setprecision(4) << m1[(i*n/2) + i];
                        cout << endl;
                    }
                    
                } else {
                    equacion = ( (2*m*m*n*n*s) - (16*m*m*n*s) + (32*m*m*s) - (14*n*n*s) + (n*n) ) / 1e6;
                    tiempo = ((fin.tv_sec - ini.tv_sec) + ((fin.tv_usec - ini.tv_usec)/1e6))/10;
                    cout << "Iter: " << n << endl;
                    cout << "Tiempo: " << tiempo << "s" << endl;
                    cout << "Flops:  " << "  " << equacion/tiempo << endl;
                    fich_out << n << "  " << tiempo << endl;
                    fich_out1 << n << "  " << equacion/tiempo << endl;
                }          
                }
                    
        } else {
                cout << " No se ha podido abrir el archivo txt " << endl;
        }

        MPI_Finalize();

    return 0;
}


