#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define X0      7500
#define tstop   40.0
#define tchange 20.0
#define rpd     18000
#define D       0.60
#define D2      0.80
#define D3      0.80
#define cb      0.95
#define sigma   3.0
#define NP      20

int main(){

    int mu;
    int n;
    int trans = 0;

    double t;
    double tau;
    double a0;
    double X[NP];
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

    for (int i = 0; i < NP; i++){
        if (D > T[i]) {
            cd[i] = 1;
        } else {
            cd[i] = exp(sigma*(D - T[i]));
        }
    }

    srand(time(NULL));

    FILE *f = fopen("m_trans3.txt", "w");
    FILE *g = fopen("t_trans3.txt", "w");

    while (t < tstop){

        for (int i = 0; i < NP; i++){
            a[i]    = X[i] * cb;       // birth
            a[i+NP] = X[i] * cd[i];   // death
        }

        a0 = 0.0;
        for (int i = 0; i < 2*NP; i++){
            a0 += a[i];
        }

        if (a0 <= 0.0) break;

        r1 = (double)rand() / RAND_MAX;
        r2 = (double)rand() / RAND_MAX;

        /* avoid log(0) */
        if (r1 == 0.0) r1 = 1e-15;

        tau = -log(r1) / a0;
        t  += tau;

        /* --- reaction selection (0-based, bounded) --- */
        tol = a[0];
        mu  = 0;
        while (mu < 2*NP - 1 && r2*a0 > tol){
            mu++;
            tol += a[mu];
        }

        if (mu < NP){
            /* birth: move newborn to a neighbouring subpopulation */
            int i      = mu;
            int offset = (rand() % 3) - 1;
            int dest   = (i + offset + NP) % NP;
            X[dest] += 1;
        } else {
            /* death: guard against negative population */
            int victim = mu - NP;
            if (X[victim] > 0)
                X[victim] -= 1;
        }

        /* --- transition updates --- */
        if (t > tchange && trans == 0){
            trans = 1;
            for (int i = 0; i < NP; i++){
                cd[i] = (D2 > T[i]) ? 1.0 : exp(sigma*(D2 - T[i]));
            }
        } else if (t > 2.0*tchange && trans == 1){
            trans = 2;
            for (int i = 0; i < NP; i++){
                cd[i] = (D3 > T[i]) ? 1.0 : exp(sigma*(D3 - T[i]));
            }
        }

        /* --- mean threshold --- */
        Xtot = 0.0;
        sum  = 0.0;
        for (int i = 0; i < NP; i++){
            Xtot += X[i];
            sum  += X[i] * T[i];
        }

        if (Xtot <= 0.0) break;   /* total extinction */
        Tmed = sum / Xtot;

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
