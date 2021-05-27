#include "uintbig.h"
#include <gmp.h>

#define N_WORDS 7
#define N_LIMBS (N_WORDS * 64 / GMP_LIMB_BITS)

const uintbig uintbig_1 = { 1, 0, 0, 0, 0, 0, 0 };

void uintbig_set(uintbig *x, uint64_t y) {
  x->c[0] = y;
  for (int i = 1; i < N_WORDS; i++)
    x->c[i] = 0;
}

bool uintbig_bit(uintbig const *x, uint64_t k) {
  return x->c[k / 64] >> (k % 64) & 1;
}

bool uintbig_add3(uintbig *x, uintbig const *y, uintbig const *z) {
  return mpn_add_n(x->c, y->c, z->c, N_LIMBS);
}
bool uintbig_sub3(uintbig *x, uintbig const *y, uintbig const *z) {
  return mpn_sub_n(x->c, y->c, z->c, N_LIMBS);
}

void uintbig_mul3_64(uintbig *x, uintbig const *y, uint64_t z) {
  mpn_mul_1(x->c, y->c, N_LIMBS, z);
}
uint64_t uintbig_div3_64(uintbig *x, uintbig const *y, uint64_t z) {
  return mpn_divmod_1(x->c, y->c, N_LIMBS, z);
}


void gentobig(uintbig *res, GEN a) {
    // assert(gsigne(a) >= 0);
    pari_sp ltop = avma;
    GEN b;
    res->c[0] = umodi2n(a,32);
    b = shifti(a, -32);
    res->c[0] += (umodi2n(b,32) << 32);
    b = shifti(b, -32);
    for (int i = 1; i < N_WORDS; i++) {
      res->c[i] = umodi2n(b,32);
      b = shifti(b, -32);
      res->c[i] += (umodi2n(b,32) << 32);
      b = shifti(b, -32);
    }
    //assert(isexactzero(b));
    avma = ltop;
}
