/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */
#define MIN_AMB_RES 3         /* min number of ambiguities for ILS-AR */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x, y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

#define PAR_EL 1

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D) {
    int i, j, k, info = 0;
    double a, *A = mat(n, n);

    memcpy(A, Q, sizeof(double) * n * n);
    for (i = n - 1; i >= 0; i--) {
        if ((D[i] = A[i + i * n]) <= 0.0) {
            info = -1;
            break;
        }
        a = sqrt(D[i]);
        for (j = 0; j <= i; j++) L[i + j * n] = A[i + j * n] / a;
        for (j = 0; j <= i - 1; j++) for (k = 0; k <= j; k++) A[j + k * n] -= L[i + k * n] * L[i + j * n];
        for (j = 0; j <= i; j++) L[i + j * n] /= L[i + i * n];
    }
    tracemat(3,D,n,1,10,3);
    free(A);
    if (info) fprintf(stderr, "%s : LD factorization error\n", __FILE__);
    return info;
}

/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j) {
    int k, mu;

    if ((mu = (int) ROUND(L[i + j * n])) != 0) {
        for (k = i; k < n; k++) L[k + n * j] -= (double) mu * L[k + i * n];
        for (k = 0; k < n; k++) Z[k + n * j] -= (double) mu * Z[k + i * n];
    }
}

/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z) {
    int k;
    double eta, lam, a0, a1;

    eta = D[j] / del;
    lam = D[j + 1] * L[j + 1 + j * n] / del;
    D[j] = eta * D[j + 1];
    D[j + 1] = del;
    for (k = 0; k <= j - 1; k++) {
        a0 = L[j + k * n];
        a1 = L[j + 1 + k * n];
        L[j + k * n] = -L[j + 1 + j * n] * a0 + a1;
        L[j + 1 + k * n] = eta * a0 + lam * a1;
    }
    L[j + 1 + j * n] = lam;
    for (k = j + 2; k < n; k++) SWAP(L[k + j * n], L[k + (j + 1) * n]);
    for (k = 0; k < n; k++) SWAP(Z[k + j * n], Z[k + (j + 1) * n]);
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z) {
    int i, j, k;
    double del;

    j = n - 2;
    k = n - 2;
    while (j >= 0) {
        if (j <= k) for (i = j + 1; i < n; i++) gauss(n, L, Z, i, j);
        del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
        //tracemat(1,D,n,1,10,3);
        if (del + 1E-6 < D[j + 1]) { /* compared considering numerical error */
            perm(n, L, D, j, del, Z);
            k = j;
            j = n - 2;
        } else j--;
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------
* args   : n      I  number of float parameters
*          m      I  number of fixed solution
           L,D    I  transformed covariance matrix
           zs     I  transformed double-diff phase biases
           zn     O  fixed solutions
           s      O  sum of residuals for fixed solutions                    */
extern int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s) {
    int i, j, k, c, nn = 0, imax = 0;
    double newdist, maxdist = 1E99, y;
    double *S = zeros(n, n);
    double *dist = mat(n, 1);
    double *zb = mat(n, 1);
    double *z = mat(n, 1);
    double *step = mat(n, 1);

    k = n - 1;
    dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = ROUND(zb[k]);
    y = zb[k] - z[k];
    step[k] = SGN(y);  /* step towards closest integer */
    for (c = 0; c < LOOPMAX; c++) {
        newdist = dist[k] + y * y / D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */
        if (newdist < maxdist) {
            /* Case 1: move down */
            if (k != 0) {
                dist[--k] = newdist;
                for (i = 0; i <= k; i++)
                    S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
                zb[k] = zs[k] + S[k + k * n];
                z[k] = ROUND(zb[k]); /* next valid integer */
                y = zb[k] - z[k];
                step[k] = SGN(y);
            }
                /* Case 2: store the found candidate and try next valid integer */
            else {
                if (nn < m) {  /* store the first m initial points */
                    if (nn == 0 || newdist > s[imax]) imax = nn;
                    for (i = 0; i < n; i++) zn[i + nn * n] = z[i];
                    s[nn++] = newdist;
                } else {
                    if (newdist < s[imax]) {
                        for (i = 0; i < n; i++) zn[i + imax * n] = z[i];
                        s[imax] = newdist;
                        for (i = imax = 0; i < m; i++) if (s[imax] < s[i]) imax = i;
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0]; /* next valid integer */
                y = zb[0] - z[0];
                step[0] = -step[0] - SGN(step[0]);
            }
        }
            /* Case 3: exit or move up */
        else {
            if (k == n - 1) break;
            else {
                k++;  /* move up */
                z[k] += step[k];  /* next valid integer */
                y = zb[k] - z[k];
                step[k] = -step[k] - SGN(step[k]);
            }
        }
    }
    for (i = 0; i < m - 1; i++) { /* sort by s */
        for (j = i + 1; j < m; j++) {
            if (s[i] < s[j]) continue;
            SWAP(s[i], s[j]);
            for (k = 0; k < n; k++) SWAP(zn[k + i * n], zn[k + j * n]);
        }
    }
    free(S);
    free(dist);
    free(zb);
    free(z);
    free(step);

    if (c >= LOOPMAX) {
        fprintf(stderr, "%s : search loop count overflow\n", __FILE__);
        return -2;
    }
    return 0;
}

/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1) (double-diff phase biases)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s) {
    int info;
    double *L, *D, *Z, *z, *E,*ZL,*W,*ZQ;

    if (n <= 0 || m <= 0) return -1;
    L = zeros(n, n);D = mat(n, 1);Z = eye(n);
    z = mat(n, 1);E = mat(n, m); W = mat(n,n);ZQ=mat(n,n);
    /* LD (lower diaganol) factorization (Q=L'*diag(D)*L) */
    if (!(info = LD(n, Q, L, D))) {
        /* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) */
        reduction(n, L, D, Z);
        matmul("TN", n, 1, n, 1.0, Z, a, 0.0, z); /* z=Z'*a */
        matmul("TN",n,n,n,1.0,Z,Q,0.0,W);
        matmul("NN",n,n,n,1.0,W,Z,0.0,ZQ);

        /* mlambda search
            z = transformed double-diff phase biases
            L,D = transformed covariance matrix */
        if (!(info = search(n, m, L, D, z, E, s))) {  /* returns 0 if no error */
            info = solve("T", Z, E, n, m, F); /* F=Z'\E */
        }
    }
    free(L);free(D);free(Z);free(z);free(E);
    return info;
}

/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z) {
    double *L, *D;
    int i, j, info;

    if (n <= 0) return -1;

    L = zeros(n, n);
    D = mat(n, 1);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            Z[i + j * n] = i == j ? 1.0 : 0.0;
        }
    /* LD factorization */
    if ((info = LD(n, Q, L, D))) {
        free(L);
        free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n, L, D, Z);

    free(L);
    free(D);
    return 0;
}

/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s) {
    double *L, *D;
    int info;

    if (n <= 0 || m <= 0) return -1;

    L = zeros(n, n);
    D = mat(n, 1);

    /* LD factorization */
    if ((info = LD(n, Q, L, D))) {
        free(L);
        free(D);
        return info;
    }
    /* mlambda search */
    info = search(n, m, L, D, a, F, s);

    free(L);
    free(D);
    return info;
}

static void matswap(int s, int n, int i, double *A){
    if(s > n ) return;
    double *m = mat(n,n);
    matcpy(m, A, n, n);
    for(int j = s; j < n; j++){
        A[j + i * n] = m[j + (i+1)*n];
        A[j + (i+1) * n] = m[j + i*n];
    }
    free(m);
}

void maxsort(double *a,int *b,int n){
    double temp1,temp2;
    int maxindex,lastindex;
    temp2 = -1;
    for(int i = 0; i < n; i++){
        if(*(a+i) > temp2){
            temp2 = *(a+i);
            maxindex = i;
        }
    }
    lastindex = maxindex;
    temp1 = temp2;
    *b  = maxindex;
    for(int j = 1; j < n; j++){
        temp1 = -1;
        for(int i = 0; i < n; i++){
            if(*(a+i) > temp1 && *(a+i) <= temp2 && i != lastindex){
                temp1 = *(a+i);
                maxindex = i;
            }
        }
        lastindex = maxindex;
        temp2=temp1;
        *(b+j) = maxindex;
    }
}

void minsort(int *a,int *b,int n){
    int temp1,temp2;
    int minindex,lastindex;
    temp2 = (int)1e8;
    for(int i = 0; i < n; i++){
        if(*(a+i) < temp2){
            temp2 = *(a+i);
            minindex = i;
        }
    }
    lastindex = minindex;
    temp1 = temp2;
    *b  = minindex;
    for(int j = 1; j < n; j++){
        temp1 = (int)1e8;
        for(int i = 0; i < n; i++){
            if(*(a+i) < temp1 && *(a+i) >= temp2 && i != lastindex){
                temp1 = *(a+i);
                minindex = i;
            }
        }
        lastindex = minindex;
        temp2=temp1;
        *(b+j) = minindex;
    }

}

extern int reshape(double *A, int na, int nx, int ny, int *index){
    int row, col;
    double  *B;

    B = mat(nx, ny);

    for (row = 0; row < nx; row++){
        for (col = 0; col < ny; col++) {
            if (row < na && col < na){
                B[row + col * ny] = A[row + col * ny];
            } else if(row < na){
                B[row + col * ny] = A[row + index[col - na] * ny];
            } else if(col < na && ny != 1){
                B[row + col * ny] = A[index[row - na] + col * ny];
            } else if(ny == 1){
                B[row + col * ny] = A[index[row - na]];
            } else{
                B[row + col * ny] = A[index[row - na] + index[col - na] * ny];
            }
        }
    }
    matcpy(A, B, nx, ny);
    free(B);

}
