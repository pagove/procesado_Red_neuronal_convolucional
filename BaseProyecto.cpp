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

using namespace std;


int main (int argc, char* args[]){
    
    /**
     * Declaracion de varibales.
     */
    int n = atoi(args[1]);
    int m = 5; 
    int s = 20; 
    float a[n][n]; 
    float c[s][m][m]; 
    float b[s][n][n];
    float r[s][n/2][n/2];
    float m1[n/2][n/2];
    float aux; 
    string caracter = args[2];
	struct timeval tim;
    double tiniu, tinis, tfinu, tfins;
    
    /*
     * =======================================================================        INICIO PROGRAMA                     ========================================================================
     */
    
    //Creacion de la matriz principal (imagen).
    for(int x = 0; x < n; x++){
        for(int y = 0; y < n; y++){
        
            a[x][y] = ( ((x*y)+x+y) / (3* pow(n,2)) ) * 1e3;
        }
    }
    
    
    //Creacion de la bateria de filtros.
     srand(5); // Se inicializa la semilla siempre con el numero 5 no con el tiempo del pc
     for (int k = 0; k < s; k++){
         for(int i = 0; i < m; i++){
             for(int j = 0; j < m; j++){
                 c[k][i][j] = (rand() * pow(-1, j+k))/ RAND_MAX;
             }
        }
    }
	
	
	/*
     * =======================================================================        INICIO CALCULO PROGRAMA               ========================================================================
     */
    
	gettimeofday(&tim, NULL);
    tiniu = tim.tv_usec;
    tinis = tim.tv_sec;
    
    //Calculo de la convolucion de la imagen.
    for(int z = 0; z < s; z++){
        for(int x = 2; x < (n - 3);  x++){
            for(int y = 2; y < (n - 3); y++){
            
                for(int i = -m/2; i <= m/2; i++){
                    for(int j = -m/2; j <= m/2; j++){
                        b[z][x][y] = a[x+i][y+j] * c[z][ (m/2) + i][ (m/2) + j]; // Como hace la multiplicacion;                         
                    }
                }
                
            }
        }
    }
    
    //Aplicacion de la funcion no lineal.
    for(int z = 0; z < s; z++){
        for(int x = 0; x < n; x++ ){
            for(int y = 0; y < n; y++){
                b[z][x][y] = 1.0/(1 + exp(-b[z][x][y]));
            }
        }
    }
    
    //Pooling
    for(int z = 0; z < s; z++){
        for(int x = 0; x < n/2; x++){
            for(int y = 0; y < n/2; y++){
            
                r[z][x][y] = fmax( fmax(b[z][2*x][2*y], b[z][2*x][(2*y) +1]), fmax(b[z][(2*x) +1][2*y], b[z][(2*x)+1][(2*y)+1]));
            }
        }
    }
    
    //promedio.
    for(int x = 0; x < n/2; x++){
        for(int y = 0; y < n/2; y++){
            aux = 0;
            for(int k = 0; k < s; k++){
                aux = aux + r[k][x][y];
            }
            
            m1[x][y] = (1.0/s) * aux;
        }
    }
    gettimeofday(&tim, NULL);
    tfinu = tim.tv_usec;
    tfins = tim.tv_sec;
    
    /*
     * =======================================================================        SALIDA PROGRAMA                     ========================================================================
     */

    if (caracter == "d"){   
        cout << "Matriz A \t Matriz M" << endl;
        
        for(int i = 0; i < 4; i++){
            cout << fixed << setprecision(4) << a[i][i] << "\t         ";
            cout << fixed << setprecision(4) << m1[i][i];
            cout << endl;
        }
        
    } else {
        cout << "Tiempo: " << ((tfins-tinis)*1e6 + tfinu-tiniu )/1e6 << "s" << endl;
    }
    
    
    return 0;
}
