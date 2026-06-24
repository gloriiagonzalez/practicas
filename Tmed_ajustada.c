#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define X0 5000 // ini cell population
#define tstop 50.0
#define tchange 50.0
#define rpd 95000 // reactions per dot
#define D 0.80
#define D2 0.50
#define D3 0.90
#define cbir 1.1
#define sigma 5.0
#define NP 20 // numero de poblaciones

#define T0 0.6
#define alpha 5.0

int main(){

    int mu;
    int n;
    int trans = 0;

    double t;
    double tau;
    double a0;
    double X[NP];
    double cb[NP];
    double cd[NP];
    double T[NP];
    double r1, r2;
    double a[2*NP];
    double paso = 1.0/NP;
    double tol;
    double Xtot;
    double sum;
    double Tmed;

    n = 0;
    t = 0;

    X[0] = X0;
    T[0] = paso;

    for (int i = 1; i < NP; i++){
        X[i] = X0;
        T[i] = T[i-1] + paso;
    }

    //ajuste nacimiento a la tolerancia 
    for (int i = 0; i < NP; i++) {
        cb[i] = cbir* pow(T[i], alpha) / (pow(T0, alpha) + pow(T[i], alpha));
    }
 


    for (int i = 0; i < NP; i++){
        if (D > T[i]) {
            cd[i] = 1;
        } else {
            cd[i] = exp(sigma*(D - T[i]));
        }
    }

    srand(time(NULL));

    FILE *f = fopen("migracion.txt", "w");
    FILE *g = fopen("threshold_medio.txt", "w");

    while (t < tstop){

        for (int i = 0; i < NP; i++){
            a[i] = X[i]*cb[i];        // birth
            a[i+NP] = X[i]*cd[i]; // death
        }

        a0 = 0.0;
        for (int i = 0; i < 2*NP; i++) {
            a0 += a[i];
        }

        if (a0 <= 0.0) break;
        
        r1 = (double)rand()/RAND_MAX;
        r2 = (double)rand()/RAND_MAX;

        tau = -log(r1)/a0;
        t = t + tau;

        tol = a[0];
        mu = 1;

        while(r2*a0 > tol){
            tol += a[mu];
            mu++;
        }

        if (mu < NP+1){
            int dest = rand() % NP;
            X[dest] = X[dest] + 1;
        } else {
            X[mu-1-NP] = X[mu-1-NP] - 1;
        }

        // transiciones
        if (t > tchange && trans == 0){
            trans = 1;
            for (int i = 0; i < NP; i++){
                if (D2 > T[i]) {
                    cd[i] = 1;
                } else {
                    cd[i] = exp(sigma*(D2 - T[i]));
                }
            }
        }
        else if (t > 2.0*tchange && trans == 1){
            trans = 2;
            for (int i = 0; i < NP; i++){
                if (D3 > T[i]) {
                    cd[i] = 1;
                } else {
                    cd[i] = exp(sigma*(D3 - T[i]));
                }
            }
        }

        // Threshold medio
        Xtot = 0.0;
        sum = 0.0;

        for (int i = 0; i < NP; i++){
            Xtot += X[i];
            sum += X[i]*T[i];
        }

        Tmed = sum/Xtot;

        if (n % rpd == 0){
            fprintf(g, "%.8f %.8f\n", t, Tmed);

            fprintf(f, "%.8f", t);
            for (int i = 0; i < NP; i++)
                fprintf(f, " %.0f", X[i]);

            fprintf(f, "\n");
        }

        n++;
    }

    fclose(f);
    fclose(g);

    return 0;
}