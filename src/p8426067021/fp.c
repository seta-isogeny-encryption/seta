#include "fp.h"
#include "rng.h"
#include "constants.h"
#include <gmp.h>
#include <pari/pari.h>
#include <assert.h>

#define FP_LIMBS 7

const fp fp_0 = {0, 0, 0, 0, 0, 0, 0};
const fp fp_1 = {1, 0, 0, 0, 0, 0, 0};
const uintbig p = { 0xf27a2e842a19c161, 0x20431191e798086f, 0x846543d5cf81ef9,
  0xb3e0f16f73f3fac2, 0x7ca4aa79a1358d65, 0x7484b4b6606a9814, 0x1965 };
const uintbig mp = { 0x44ae023bb1f00634, 0x3af1ce6a4ca4abb2, 0xeb01332fdc38d3d5,
  0x7341d22fa6e160eb, 0x8d5963746c74b128, 0x798ec126c5cb80d0, 0x3c6 }; // 2^(7·64) mod p
// 2^4-th roots of unity with the square roots of their inverse
const struct {
  uint64_t limb;
  fp root;
} roots_of_unity[16] = {
  { 0x1, { 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 } },
  { 0x3fc912e708b9bc94, { 0x49628593e22d3f07, 0x3806d750cc070161, 0x5c44fe281c38b5db, 0xe7808dbca8f08a24, 0x83ab7e05471fe289, 0xd16d7dbc1cbb5fe7, 0xedf } },
  { 0x6384ab2a23e3f47e, { 0xf2853f72585b4d0, 0xe42bf3d86edb5fb4, 0x65f2d7c047679390, 0xa8a55c8d92b5e21a, 0x3692c8beae1f7a59, 0x787662aaa6335a38, 0x1df } },
  { 0xf9d209392b46864e, { 0x5978b8e983aebb1c, 0x74cb701c3fe7d994, 0xa8b9824f76a8abc5, 0x5e4e785421783e30, 0x61743bbb8dd8764d, 0xc00f82a4d9456ee8, 0x867 } },
  { 0x47596d6f3db41fc4, { 0x906e5c74b87d4155, 0xc4bb32e8cd0891d4, 0x69bca847bb2012c6, 0xb1d2bb3a5754fc5d, 0x85603c95d027e0bb, 0x66469fe9eaf9730f, 0x138 } },
  { 0xd41aff8b08b27c40, { 0x2b69c91f6be62a62, 0xe121d81f0a2046b6, 0x7d12810a6c06db9f, 0xa2a530576b2430cd, 0xd3126366c87232d, 0xcf15d2ca109f0ce2, 0x76f } },
  { 0x620bd20f719c800c, { 0x1e5f2ef921674521, 0xe61bc57cff8bc932, 0x19a0469dca6e28d5, 0xd02692adc8970fde, 0x33f530578fd9f210, 0x51f95b4f9180fce3, 0x15b4 } },
  { 0xe351da8d04940c91, { 0xcb44678912ebf0f9, 0xaf3541183b44086e, 0x9c98936247870b9e, 0xdc1da753a8367274, 0x944b86f5dc3e3bcc, 0xed73c2f73fd30685, 0xc10 } },
  { 0xf27a2e842a19c160, { 0xab20c114ec65a19d, 0xf33e9dc2c0c309c9, 0x63ce7b6a01d8ac2, 0x2db4cd44f3721bd5, 0x56cea44edacdb49, 0xebf3e524c4cf83dd, 0x305 } },
  { 0xb2b11b9d216004cd, { 0xcf4a476a289a73f5, 0x41a466d1cb57b33, 0x4b79bff5e179c688, 0xe39b6a73b8f8869a, 0x51b1c3df7ae1d6b0, 0x597a5ae381ad065, 0x8f6 } },
  { 0x8ef5835a0635cce3, { 0xf8a8254afed33b13, 0xaae44912ab07c973, 0xcba2bff8b7c09773, 0xac22df978af485d0, 0x9988a99ba2a79119, 0xa79382611676fe51, 0x1194 } },
  { 0xf8a8254afed33b13, { 0x70aec77ed171a012, 0x2c93b7e1956523bc, 0x6ac9172ba5f18843, 0x9d88640716f0f0fc, 0x11dc4ea38bb8c438, 0xb69f07d3ad1a6e98, 0x1056 } },
  { 0xab20c114ec65a19d, { 0x8ef5835a0635cce3, 0x7c4fb9885857136f, 0x61faa2695ae77c86, 0xb1d2bb3a5754fbeb, 0x85603c95d027e0bb, 0x66469fe9eaf9730f, 0x138 } },
  { 0x1e5f2ef921674521, { 0x9ed5d3c87646dc20, 0xb1e88995bde72a59, 0x694f15a04e4d9dcb, 0x72147e3f00503d02, 0x83b96a675ddc852e, 0x154efee9f7f746e3, 0x83b } },
  { 0x906e5c74b87d4155, { 0xb2b11b9d216004cd, 0xe6f8dfe63d4fd142, 0x340a7cb325c8f003, 0x3b1317f13f3f66b, 0x5005d31d6375fa56, 0x98d7a4e3990a86b, 0xde0 } },
  { 0xf2853f72585b4d0, { 0x244e3faae8a5a4b2, 0xad9972aeb7375a2b, 0xde6eed9f18abb25f, 0x55c3508bad50b99e, 0x61306861e1abbeb2, 0x1ce225897837dcc, 0x1249 } },
};

void fp_set(fp *x, uint64_t y) {
  x->x.c[0] = y;
  for (int i = 1; i < FP_LIMBS; i++)
    x->x.c[i] = 0;
}

void fp_cswap(fp *x, fp *y, bool c) {
  uint64_t tmp;
  for (int i = 0; i < FP_LIMBS*c; i++) {
    tmp = y->x.c[i];
    y->x.c[i] = x->x.c[i];
    x->x.c[i] = tmp;
  }
}

void fp_enc(fp *x, uintbig const *y) {
  for (int i = 0; i < FP_LIMBS; i++)
    x->x.c[i] = y->c[i];
}
void fp_dec(uintbig *x, fp const *y) {
  for (int i = 0; i < FP_LIMBS; i++)
    x->c[i] = y->x.c[i];
}

void fp_add2(fp *x, fp const *y) { fp_add3(x, x, y); }
void fp_sub2(fp *x, fp const *y) { fp_sub3(x, x, y); }
void fp_mul2(fp *x, fp const *y) { fp_mul3(x, x, y); }

void fp_add3(fp *x, fp const *y, fp const *z) {
  mp_limb_t carry = mpn_add_n(x->x.c, y->x.c, z->x.c, FP_LIMBS);
  if (carry) {
    mpn_add_n(x->x.c, x->x.c, mp.c, FP_LIMBS);
  } else if (mpn_cmp(x->x.c, p.c, FP_LIMBS) >= 0) {
    mpn_sub_n(x->x.c, x->x.c, p.c, FP_LIMBS);
  }
}

void fp_sub3(fp *x, fp const *y, fp const *z) {
  mp_limb_t borrow = mpn_sub_n(x->x.c, y->x.c, z->x.c, FP_LIMBS);
  if (borrow) {
    mpn_add_n(x->x.c, x->x.c, p.c, FP_LIMBS);
  }
}

void fp_mul3(fp *x, fp const *y, fp const *z) {
  uint64_t tmp[FP_LIMBS*2], thrash[FP_LIMBS+1];
  mpn_mul_n(tmp, y->x.c, z->x.c, FP_LIMBS);
  mpn_tdiv_qr(thrash, x->x.c, 0, tmp, 2 * FP_LIMBS, p.c, FP_LIMBS);
}

void fp_mul3_64(fp *x, fp const *y, uint64_t z) {
  uint64_t tmp[FP_LIMBS + 1], thrash[2];
  mpn_mul(tmp, y->x.c, FP_LIMBS, &z, 1);
  mpn_tdiv_qr(thrash, x->x.c, 0, tmp, FP_LIMBS + 1, p.c, FP_LIMBS);
}

void fp_sq1(fp *x) { fp_sq2(x, x); }
void fp_sq2(fp *x, fp const *y) { fp_mul3(x, y, y); }
void fp_inv(fp *x) {
  mpz_t res, mpzx, mpzp;
  mpz_init(res);
  mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
  mpz_roinit_n(mpzp, p.c, FP_LIMBS);
  mpz_invert(res, mpzx, mpzp);
  int i = 0;
  for (; i < res->_mp_size; ++i) {
    x->x.c[i] = ((uint64_t*)res->_mp_d)[i];
  }
  for (; i < FP_LIMBS; ++i) {
    x->x.c[i] = 0;
  }
  mpz_clear(res);
}
bool fp_issquare(fp *x) {
  mpz_t mpzx, mpzp;
  mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
  mpz_roinit_n(mpzp, p.c, FP_LIMBS);
  int s = mpz_legendre(mpzx, mpzp);
  return s+1;
}

// Tonelli - Shanks
void fp_sqrt(fp *x) {
  mpz_t mpzx, mpzp, mpzsqrt, mpzcof, mpzQ, mpzexp;
  mpz_init(mpzsqrt);
  mpz_init(mpzcof);
  mpz_init(mpzexp);
  
  mpz_roinit_n(mpzx, x->x.c, FP_LIMBS);
  mpz_roinit_n(mpzp, p.c, FP_LIMBS);
  mpz_roinit_n(mpzQ, p_even_cofactor.c, FP_LIMBS);
  mpz_add_ui(mpzexp, mpzQ, 1);
  mpz_div_ui(mpzexp, mpzexp, 2);
  mpz_powm(mpzsqrt, mpzx, mpzexp, mpzp);
  mpz_powm(mpzcof, mpzx, mpzQ, mpzp);

  int i = 0;
  for (; i < mpzsqrt->_mp_size; ++i) {
    x->x.c[i] = ((uint64_t*)mpzsqrt->_mp_d)[i];
  }
  for (; i < FP_LIMBS; ++i) {
    x->x.c[i] = 0;
  }

  // shortcut: there's few enough 2^n roots of unity, that we can just
  // use a table
  for (i = 0; i < 1 << (two_tors_height - 1); i++) {
    if (roots_of_unity[i].limb == ((uint64_t*)mpzcof->_mp_d)[0]) {
      fp_mul2(x, &roots_of_unity[i].root);
      break;
    }
  }
  assert(i < 1 << (two_tors_height - 1));
    
  mpz_clear(mpzsqrt);
  mpz_clear(mpzcof);
  mpz_clear(mpzexp);
}

void fp_random(fp *x) {
  uint64_t thrash;
  randombytes(x->x.c + 0, 8*FP_LIMBS);
  mpn_tdiv_qr(&thrash, x->x.c, 0, x->x.c, FP_LIMBS, p.c, FP_LIMBS);
}
