#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char const *argv[]){


    char **Consultas;

    Consultas = (char **) malloc(sizeof (char *)*3816);
    for (int i = 0; i < 106*6*3*2; i++)
        Consultas[i] = (char *) malloc(sizeof (char)*500);

    int contador=0;
    for (int i = 0; i < 106; ++i)
    {
        if (i<9){
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/left/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/left/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/left/ring_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/right/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/right/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/00%d/right/ring_%d.bmp", i+1,j+1);
                contador++;
            }
        } if (i>=9 && i < 99){
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/left/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/left/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/left/ring_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/right/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/right/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/0%d/right/ring_%d.bmp", i+1,j+1);
                contador++;
            }
        }
        if (i>=99){
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/left/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/left/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/left/ring_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/right/index_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/right/middle_%d.bmp", i+1,j+1);
                contador++;
            }        
            for (int j = 0; j < 6; ++j){
                sprintf(Consultas[contador], "/home/rhernandez/SDUMLA/roi/%d/right/ring_%d.bmp", i+1,j+1);
                contador++;
            }
        }

	//printf ("%s \n", Consultas[contador-1]);
    }


    omp_set_num_threads(12);
   #pragma omp parallel
    {
        int tid = omp_get_thread_num(); //ID del thread
        int procs = omp_get_num_threads(); //Nro. total de threads
        int i=0,a=0,ii=0;
        for (int t = tid+1; t <= 106; t+=procs){
            for (int j = 1; j <= 2; ++j){
                for (int k = 1; k <= 3; ++k){
                    for (int l = 1; l <= 6; ++l){
                        i= (t-1)*36+(j-1)*18+(k-1)*6+(l-1);
                        ii= (int)(i/36)+1;

                        for (int q = 1; q <= 106; ++q){
                            for (int w = 1; w <= 2; ++w){
                                for (int e = 1; e <= 3; ++e){
                                    for (int r = 1; r <= 6; ++r){
                                        a=(q-1)*36+(w-1)*18+(e-1)*6+(r-1);
                                        char ejecutar[500];
                                        sprintf(ejecutar, "touch Resultados_persona_%d.txt && a=$(/home/rhernandez/CPM_Daisy/CPM %s %s) && echo %d,%d,%d,%d,%d,%d,%d,%d,$a >> Resultados_persona_%d.txt", ii, Consultas[i], Consultas[a], ii,j,k,l,q,w,e,r,ii);
                                        system(ejecutar);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}
