#include "curve.h"

// cube root of unity = -√-3/2 - 1/2
static const fp2 w = {
  { 0xf93d1742150ce0b0, 0x902188c8f3cc0437, 0x4232a1eae7c0f7c,
    0xd9f078b7b9f9fd61, 0x3e52553cd09ac6b2, 0xba425a5b30354c0a, 0xcb2 },
  { 0xf93d1742150ce0b0, 0x902188c8f3cc0437, 0x4232a1eae7c0f7c,
    0xd9f078b7b9f9fd61, 0x3e52553cd09ac6b2, 0xba425a5b30354c0a, 0xcb2 }
};
// fourth root of unity i := (-3)^(p-1)/4
static const fp i = {
  0xab20c114ec65a19d, 0xf33e9dc2c0c309c9, 0x63ce7b6a01d8ac2,
  0x2db4cd44f3721bd5, 0x56cea44edacdb49, 0xebf3e524c4cf83dd, 0x305
};

// Curve y² = x³ + √3 x² + x, with CM by √-3, j-invariant 0
void init_curve(proj *E0) {
  E0->x =  (fp2){ fp_0, i };
  E0->z = fp2_1;
}

const int dist_trace = -1;
const int dist_norm = 1;
const int q_norm = 3;

// distorsion map on E0, twisted Edwards form
// X ↦ X  ( 4Z + 2(Y-Z)(ωi + 1) )   ( 4Y + 2(Y-Z)(ωi + ω) )
// Y ↦ Z  2ω²(Y + Z)                ( 4Z + 2(Y-Z)(ωi - ω) )
// Z ↦ Z  2ω²(Y + Z)                ( 4Y + 2(Y-Z)(ωi + ω) )   
// T ↦ X  ( 4Z + 2(Y-Z)(ωi + 1) )   ( 4Z + 2(Y-Z)(ωi - ω) )
void ted0_dist(point *Q, const point *P) {
  // 2(ωi + ω) = -j(1+i) - (1+i)
  // 2(ωi - ω) = j(1-i) + (1-i)
  // 2(ωi + 1) = -ji + (2-i)
  // 2ω² = j - 1
  fp2 wipw, wimw, wi1, w2;
  fp_neg2(&wi1.im, &i);
  wi1.re = wi1.im; wi1.re.x.c[0] += 2;
  wipw.im = wi1.im; wipw.im.x.c[0] -= 1;
  wipw.re = wipw.im; 
  wimw.im = wi1.im; wimw.im.x.c[0] += 1;
  wimw.re = wimw.im;
  w2.im = fp_1;
  fp_neg2(&w2.re, &fp_1);
  
  fp2 YmZ, YpZ, Y4, Z4, A, B, C, D;
  fp2_sub3(&YmZ, &P->y, &P->z);
  fp2_add3(&YpZ, &P->y, &P->z);
  fp2_add3(&Y4, &P->y, &P->y);
  fp2_add2(&Y4, &Y4);
  fp2_add3(&Z4, &P->z, &P->z);
  fp2_add2(&Z4, &Z4);
  
  fp2_mul3(&A, &YmZ, &wi1);
  fp2_add2(&A, &Z4);
  fp2_mul2(&A, &P->x);
  fp2_mul3(&B, &YpZ, &w2);
  fp2_mul2(&B, &P->z);
  fp2_mul3(&C, &YmZ, &wipw);
  fp2_add2(&C, &Y4);
  fp2_mul3(&D, &YmZ, &wimw);
  fp2_add2(&D, &Z4);

  fp2_mul3(&Q->x, &A, &C);
  fp2_mul3(&Q->y, &B, &D);
  fp2_mul3(&Q->z, &B, &C);
  fp2_mul3(&Q->t, &A, &D);
}

// distorsion map on E0, montgomery form
// (X:Y:Z) ↦ ( ω(X - iωZ) : Y : Z )
void mont0_dist(proj *Q, const proj *P) {
  proj Pcopy = *P;
  Q->x.re = i; Q->x.im = fp_0;
  fp2_mul2(&Q->x, &w);
  fp2_mul2(&Q->x, &Pcopy.z);
  fp2_sub3(&Q->x, &Pcopy.x, &Q->x);
  fp2_mul2(&Q->x, &w);
  Q->z = Pcopy.z;
}

// distorsion map on E0, montgomery form
// just the same
void montxy0_dist(proj2 *Q, const proj2 *P) {
  proj2 Pcopy = *P;
  Q->x.re = i; Q->x.im = fp_0;
  fp2_mul2(&Q->x, &w);
  fp2_mul2(&Q->x, &Pcopy.z);
  fp2_sub3(&Q->x, &Pcopy.x, &Q->x);
  fp2_mul2(&Q->x, &w);
  Q->y = Pcopy.y;
  Q->z = Pcopy.z;
}

// frobenius map on E0, twisted Edwards form
// ( X : Y : Z : T ) ↦ ( iT : Z : Y : iX )^p
void ted0_frob(point *Q, const point *P) {
  point Pcopy = *P;
  Q->t.re = i;
  Q->t.im = fp_0;
  fp2_mul3(&Q->x, &Q->t, &Pcopy.t);
  Q->y = Pcopy.z;
  Q->z = Pcopy.y;
  fp2_mul2(&Q->t, &Pcopy.x);
  fp2_frob2(&Q->x, &Q->x);
  fp2_frob2(&Q->y, &Q->y);
  fp2_frob2(&Q->z, &Q->z);
  fp2_frob2(&Q->t, &Q->t);
}

// frobenius map on the twist of E0, twisted Edwards form
// ( X : Y : Z : T ) ↦ ( -T : Z : Y : -X )^p
void ted0_frob_twist(point *Q, const point *P) {
  point Pcopy = *P;
  Q->x = Pcopy.t;
  Q->y = Pcopy.z;
  Q->z = Pcopy.y;
  Q->t = Pcopy.x;
  fp2_neg1(&Q->x);
  fp2_neg1(&Q->t);
  fp2_frob2(&Q->x, &Q->x);
  fp2_frob2(&Q->y, &Q->y);
  fp2_frob2(&Q->z, &Q->z);
  fp2_frob2(&Q->t, &Q->t);
}

// frobenius map on E0, montgomery form
void mont0_frob(proj *Q, const proj *P) {
  proj Pcopy = *P;
  fp2_frob2(&Q->x, &Pcopy.x);
  fp2_neg1(&Q->x);
  fp2_frob2(&Q->z, &Pcopy.z);
}

// frobenius map on E0, montgomery form
void montxy0_frob(proj2 *Q, const proj2 *P) {
  proj2 Pcopy = *P;
  fp2_frob2(&Q->x, &Pcopy.x);
  fp2_neg1(&Q->x);
  Q->y = (fp2){ i, fp_0 };
  fp2_mul2(&Q->y, &Pcopy.y);
  fp2_frob2(&Q->y, &Q->y);
  fp2_frob2(&Q->z, &Pcopy.z);
}

// sqrt(-q) map in E0, Montgomery form
// (X:Z) ↦ ( -X (X + Z√3)² : Z (X√3 + Z)² )
void mont0_sqrt_minus_q(proj *Q, const proj *P) {
  proj Pcopy = *P;
  Q->x.re = fp_0;
  Q->x.im = i;
  fp2_mul3(&Q->z, &Pcopy.x, &Q->x);
  fp2_add2(&Q->z, &Pcopy.z);
  fp2_sq1(&Q->z);
  fp2_mul2(&Q->z, &Pcopy.z);
  
  fp2_mul2(&Q->x, &Pcopy.z);
  fp2_add2(&Q->x, &Pcopy.x);
  fp2_sq1(&Q->x);
  fp2_mul2(&Q->x, &Pcopy.x);
  fp2_neg1(&Q->x);
}

// sqrt(-q) map in E0, twisted Edwards form
void ted0_sqrt_minus_q(point *Q, const proj *E, const point *P) {
    point A = *P;
    ted0_dist(&A,&A);
    ted_double(&A, E, &A);
    ted_add(&A, E, &A, P);
    ted_neg(&A, &A);
    *Q = A;
}
