#include "fp2.h"

const fp2 fp2_0 = { { 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0 } };

/* Arithmetic modulo X^2 + c */
static const uint64_t c = 19;  // No point in sharing this, hopefull
static const fp minvc = { 0x6f4b7e2298080773, 0x5d6d1db84760b3e3, 0xe6cd4274877b69,
  0x88d03fbf8c8ba21a, 0xec4b731f63861282, 0x91e4773789bb3f93, 0x6f }; // -1/c

// √-c turns out to be a non-square
fp2 fp2_non_residue() {
  return fp2_i;
}

void fp2_mul3(fp2 *x, fp2 const *y, fp2 const *z) {
  fp xsum, xim;
  fp_add3(&xsum, &y->re, &y->im);
  fp_add3(&xim, &z->re, &z->im);
  fp_mul2(&xsum, &xim);
  fp_mul3(&xim, &y->im, &z->im);
  fp_mul3(&x->re, &y->re, &z->re);
  fp_sub3(&x->im, &xsum, &xim);
  fp_sub2(&x->im, &x->re);
  fp_mul3_64(&xim, &xim, c);
  fp_sub2(&x->re, &xim);
}

void fp2_sq2(fp2 *x, fp2 const *y) {
  fp sum, diff, cim;
  fp_mul3_64(&cim, &y->im, c);
  fp_add3(&sum, &y->re, &cim);
  fp_sub3(&diff, &y->re, &y->im);
  fp_mul3(&x->im, &y->re, &y->im);
  fp_mul3_64(&cim, &x->im, c - 1);
  fp_add2(&x->im, &x->im);
  fp_mul3(&x->re, &sum, &diff);
  fp_sub2(&x->re, &cim);
}

void fp2_inv(fp2 *x) {
  fp inorm, im2;
  fp_sq2(&inorm, &x->re);
  fp_sq2(&im2, &x->im);
  fp_mul3_64(&im2, &im2, c);
  fp_add2(&inorm, &im2);
  fp_inv(&inorm);
  fp_mul2(&x->re, &inorm);
  fp_mul2(&x->im, &inorm);
  fp_neg1(&x->im);
}

bool fp2_issquare(const fp2 *x) {
  fp inorm, im2;
  fp_sq2(&inorm, &x->re);
  fp_sq2(&im2, &x->im);
  fp_mul3_64(&im2, &im2, c);
  fp_add2(&inorm, &im2);
  return fp_issquare(&inorm);
}

void fp2_frob2(fp2 *x, const fp2 *y) {
  x->re = y->re;
  fp_neg2(&x->im, &y->im);
}

void fp2_exp(fp2 *res, fp2 const *x, uintbig const *k)
{
    if (fp2_iszero(x)) { *res = *x; return; }
    const fp2 xcopy = *x;
    *res = fp2_1;

    unsigned long i = BITS;
    while (--i && !uintbig_bit(k, i));
    do {
        fp2_sq1(res);
        if (uintbig_bit(k, i)) {
            fp2_mul2(res, &xcopy);
        }
    } while (i--);
}



// dlp of h in basis g, which has order ell, naive implementation
bool fp2_dlp_naive(long *res, const fp2 *h, const fp2 *g, long ell) {
    long logarithm = 0;
    fp2 x = fp2_1;

    for (int i = 0; i < ell; ++i) {
        if (fp2_equal(h,&x)) { *res = logarithm; return true;}
        logarithm++;
        fp2_mul2(&x,g);
    }

    return false;
}

void fp2_sqrt(fp2 *x) {
    if (fp_iszero(&x->im)) {
        fp x_re_copy = x->re;

        if (fp_issquare(&x_re_copy)) {
            fp_sqrt(&x->re);
            return;
        }
        else {
	  fp_mul3(&x->im, &x->re, &minvc);
	  fp_sqrt(&x->im);
	  fp_set(&x->re, 0);
	  return;
        }
    }

    fp sdelta, re, tmp1, inv2, im;

    // sdelta = sqrt(re^2 + c im^2)
    fp_sq2(&sdelta, &x->re);
    fp_sq2(&tmp1, &x->im);
    fp_mul3_64(&tmp1, &tmp1, c);
    fp_add2(&sdelta, &tmp1);

    fp_sqrt(&sdelta);

    fp_set(&inv2,2);
    fp_inv(&inv2);

    fp_add3(&re,&x->re,&sdelta);
    fp_mul2(&re,&inv2);

    if (!fp_issquare(&re)) {
        fp_sub3(&re,&x->re,&sdelta);
        fp_mul2(&re,&inv2);
    }

    fp_sqrt(&re);

    im = re;

    fp_inv(&im);
    fp_mul2(&im,&inv2);
    fp_mul2(&im,&x->im);

    x->re = re;
    x->im = im;
}
