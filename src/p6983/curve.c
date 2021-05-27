#include "curve.h"

// Curve y² = x³ + x, with CM by √-4, j-invariant 1728
void init_curve(proj *E0) {
  E0->x = fp2_0;
  E0->z = fp2_1;
}

const int dist_trace = 0;
const int dist_norm = 1;

// distorsion map on E0, twisted Edwards form
void ted0_dist(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_mul3(&Q->x, &Pcopy.t, &fp2_i);
    Q->y = Pcopy.z;
    Q->z = Pcopy.y;
    fp2_mul3(&Q->t, &Pcopy.x, &fp2_i);
}

// distorsion map on E0, montgomery form
void mont0_dist(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_neg2(&Q->x, &Pcopy.x);
    Q->z = Pcopy.z;
}

// distorsion map on E0, montgomery form
void montxy0_dist(proj2 *Q, const proj2 *P) {
    proj2 Pcopy = *P;
    fp2_neg2(&Q->x, &Pcopy.x);
    fp2_mul3(&Q->y, &Pcopy.y, &fp2_i);
    Q->z = Pcopy.z;
}

// frobenius map on E0, montgomery form
void ted0_frob(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
}

// frobenius map on the twist of E0, montgomery form
void ted0_frob_twist(point *Q, const point *P) {
    point Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
    fp2_frob2(&Q->t, &Pcopy.t);
    //    fp2_mul2(&Q->x, &ted0_twister_frob);
    //fp2_mul2(&Q->t, &ted0_twister_frob);
}

// frobenius map on E0, montgomery form
void mont0_frob(proj *Q, const proj *P) {
    proj Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->z, &Pcopy.z);
}

// distorsion map on E0, montgomery form
void montxy0_frob(proj2 *Q, const proj2 *P) {
    proj2 Pcopy = *P;
    fp2_frob2(&Q->x, &Pcopy.x);
    fp2_frob2(&Q->y, &Pcopy.y);
    fp2_frob2(&Q->z, &Pcopy.z);
}
