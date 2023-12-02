/****************** BaseProyectoOpenMP.cpp ********************/ /**
 *
 * @File: BaseProyectoOpenMP.cpp
 *
 * @Brief Este programa es la paralelizacion mediante OpenMP del programa
 *        BaseProyecto.
 *                                             
 * @Author: Ivan Gregoria Bou
 * @Author: Pablo Gomar Velló
 * @date 09/03/2020                               
 * 
 */ /*********************************************************/
 
#include <iostream>
#include <vector> 
#include <math.h>
#include <algorithm>
#include <sys/time.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <omp.h>


using namespace std;

const int REPS = 10;

int main (int argc, char* args[]){
    
    /**
    * float *A_mem = new float[]
    */
            
    /**
     * Declaracion de varibales.
     */
    long n = atoi(args[1]);
    long m = 5; 
    long s = 20;    
    float *a = new float[n*n];
    float *c = new float[s*m*m];
    float *b = new float[s*n*n];
    float *r = new float[s*(n/2)*(n/2)];
    float *m1 = new float[(n/2)*(n/2)];
    float aux; 
    double long equacion;
    float suma;
    string caracter = args[2];
    struct timeval ini, fin;
    double tiempo;
    ofstream fich_out("ResultadosTiempo.txt", ios::app);
    ofstream fich_out1("ResultadosFlops.txt", ios::app);
    //int hilos = omp_get_max_threads();
    //int hilos = 24;
      int hilos = atoi(args[3]);
     /*
     * =======================================================================        INICIO PROGRAMA                     ========================================================================
     */
    
    //omp_set_num_threads(hilos);
    cout<< "\nNumero de hilos: " << hilos << endl;
    
    if (fich_out.is_open() && fich_out1.is_open()){
        
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
        
        
        /*
        * =======================================================================        INICIO CALCULO PROGRAMA               ========================================================================
        */
        
        gettimeofday(&ini, NULL);
        
        
        for (int t = 0; t < REPS; t++){
            //Calculo de la convolucion de la imagen.
            aux = 0;
            #pragma omp for schedule(dynamic, 2)
            for(int i = 0; i < s*n*n; i++){
                b[i] = 0;
            }
            
            #pragma omp parallel 
            {
                #pragma omp for schedule (dynamic, 2) reduction(+:suma)
                for(int z = 0; z < s; z++){
                    for(int x = 2; x <= (n - 3);  x++){
                        for(int y = 2; y <= (n - 3); y++){
                            suma = 0.0;
                            for(int i = -m/2; i <= m/2; i++){
                                for(int j = -m/2; j <= m/2; j++){                  
                                    suma += a[(n* (x + i)) + y + j] * c[(m*m*z) + (m*((m/2) + i)) + ((m/2) + j)];
                                }
                            }
                            b[(z*n*n) + (x*n) + y] = suma;
                            
                        }
                    }
                }
            }
            
            #pragma omp parallel shared(b)
            {
                //Aplicacion de la funcion no lineal.
                #pragma omp for schedule(dynamic, 2) 
                for(int z = 0; z < s; z++){
                    for(int x = 0; x < n; x++ ){
                        for(int y = 0; y < n; y++){
                            b[(z*n*n) + (n*x) + y] = 1.0/(1 + exp(-b[(z*n*n) + (n*x) + y]));
                        }
                    }
                }
            }
            
            #pragma omp parallel for shared(r, b) schedule (dynamic, 2)
                //Pooling
                for(int z = 0; z < s; z++){
                    for(int x = 0; x < n/2; x++){
                        for(int y = 0; y < n/2; y++){
                            r[(z*n/2*n/2) + (x*n/2) + y] = fmax( fmax(b[(z*n*n)+ (n*2*x) + 2*y], b[(z*n*n) + (n*2*x) + ((2*y) +1)]), fmax(b[(z*n*n) + (((2*x) +1)*n) + (2*y)], b[(z*n*n) + (((2*x)+1)*n) + ((2*y)+1)]));
                        }
                    }
                }
            
            #pragma omp parallel for shared(r) schedule(dynamic, 2) reduction(+:aux)
                //promedio.
                for(int x = 0; x < n/2; x++){
                    for(int y = 0; y < n/2; y++){
                        aux = 0;
                        for(int k = 0; k < s; k++){
                            aux = aux + r[(k*n/2 * n/2) + (n/2 *x) + y];
                        }
                        m1[(x*n/2) + y] = aux/s ;
                    }
                }
            
        }
            
        gettimeofday(&fin, NULL);
        
        
        /*
        * =======================================================================        SALIDA PROGRAMA                     ========================================================================
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
    
    } else {
        cout << " No se ha podido abrir el archivo txt " << endl;
    }    
    
    return 0;
}

