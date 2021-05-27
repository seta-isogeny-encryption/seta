#include "fp.h"
#include "constants.h"
#include "rng.h"
#include <gmp.h>
#include <pari/pari.h>
#include <assert.h>

#define FP_LIMBS 7

const fp fp_0 = {0, 0, 0, 0, 0, 0, 0};
const fp fp_1 = {1, 0, 0, 0, 0, 0, 0};
const uintbig p = { 0x86a90941ba75a7c1, 0x31828544879e2296, 0x1b685fe43ce3748,
  0x3f212b8be3c80cb, 0xdc28deed6b1eff8, 0xfb987c1cb8e3c599, 0xd3 };
const uintbig mp = { 0x9dc63624671879e, 0x45946ae936d7e095, 0x691a474146c7f6ef,
  0x297e6130c82d43d2, 0x2a48c4b06ac25191, 0x1bf59ddb9d0bd84c, 0x97 }; // 2^(7·64) mod p
// 2^5-th roots of unity with the square roots of their inverse
const struct {
  uint64_t limb;
  fp root;
} roots_of_unity[32] = {
  { 0x1, { 0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 } },
  { 0x613ec2eed657a62f, { 0x4310088acccf0c0f, 0x7cf934dd1832f4db, 0xe18eb36f5fc1f7fb, 0x160fe2de7f8dd830, 0x4b9b170fbbb85aaa, 0xbc220c04e7870ee0, 0x95 } },
  { 0x71a2e11189bca5a6, { 0xa16811b7b02a8851, 0xc697e0b3604b939f, 0x2b35380a09d98866, 0x84f91c8498370540, 0x6dd51297a0eb66de, 0x1f0bf12003fbf9b2, 0x99 } },
  { 0xf017710f15ed4217, { 0x9b12511f143614d8, 0x2b7c8221ae21ad94, 0xb58d2e3d2f417153, 0x995ac8c7c2e8b7f8, 0xe8ce5da060b6093d, 0xcf9d82f98e47bf74, 0x4e } },
  { 0x6ba6d5d613899978, { 0x88cfc28bd145cf9a, 0xa4efb70590753adc, 0xe8737c776cad3de5, 0xbf185a4ffc6d36fd, 0x5c8e7243a43dabec, 0x467b868414a7cf34, 0xb9 } },
  { 0xb2afb2de448a6c21, { 0x91071f72b80f600f, 0xd64e33ff4c48942, 0x6cbd23b1de0e9740, 0x84320a42488d0548, 0x239d93e9c8708a40, 0x56dafa0e1a209468, 0x7c } },
  { 0xb24d69c9631375fb, { 0xadfab22c57e0c1c5, 0xc084c59e7781af9c, 0x845bc2e60d330e55, 0x3af423a229432f45, 0x8705def7580191b1, 0x3b4c37f3be810ade, 0x91 } },
  { 0xce22ab51d5b2d5a9, { 0x16e97bd62721fc19, 0x9a43368d90cb6415, 0x4a1c4b2f826049ff, 0x6c6c6482f0ed3ea7, 0x40d73a92a52cc01c, 0xe94f5fde593d5afd, 0x93 } },
  { 0xc7b1b74064686bcd, { 0xcc9e5e9e7e874b67, 0x647b5d79b76058ae, 0xf355d556b7e8c7ec, 0xc6c39090c5f2061b, 0xb54432abd7b6e08c, 0xce7bc0ca6493ef2a, 0xb5 } },
  { 0x62fd8c664ae84aa4, { 0xd65761dfc6cd9e72, 0xa9392d5248fa0d3, 0xd37a3fbf00afca32, 0xebee8f1be0e4e4de, 0xd20352e96eb58fb3, 0xa09b718b7bb71f32, 0x4c } },
  { 0xcbae27b4fb12d8f3, { 0x3339a3d3cbba52d9, 0xdeade525a731258a, 0x904e2a4b38affc51, 0x1be87aefe9dfe531, 0x6edc39d069c76e0b, 0x184d6e9f17d862e6, 0x79 } },
  { 0x536f656deebb54e8, { 0x57ecaf903c6570cb, 0xe2de43d9e68cd81e, 0x8a4738c24d6ad955, 0xf9fea8fb302efbf1, 0x43dd8adec8c28672, 0xf47fba620fff399e, 0xb3 } },
  { 0xba0aaaa33bee5c5a, { 0xbafae18cbf62cece, 0x4b25bf24f5a36a37, 0x5d2713e30b5925f2, 0x12b800ff7f51ac54, 0x3915af892bf88b36, 0x82a72c288a3e216c, 0xbc } },
  { 0xd8ae57156294e5fc, { 0xcf09236d7c592800, 0x119db480d345da01, 0xcddcbcad431f3585, 0x3caf296b9344cfdf, 0xebfd4d401a76ab24, 0x4a59568c36dcd663, 0x9a } },
  { 0xfdd946b5e92fd827, { 0x23ab7cdb6f8d5d1d, 0x6102cc155079d0e4, 0xda41a374e1f4b1ed, 0xdfcd9cc3db4e9d81, 0x43cd9d68c19e7c24, 0xb0768009ababac74, 0x53 } },
  { 0xe540f78a0a4b1f70, { 0x860b629df247dab3, 0x8493464862572efb, 0x7262b71be5371674, 0x11b2840eca953793, 0x75934f0211fa1f5e, 0x58a07ca54d9446de, 0xbe } },
  { 0x86a90941ba75a7c0, { 0xbef75201560d3bf4, 0xff349ef3ba8ebb92, 0x47f31bc1395a7d63, 0x78e7f71088c58a38, 0xf49617a525d1e83e, 0xd766b0c3cf6fbf9d, 0x88 } },
  { 0x256a4652e41e0192, { 0xe1d416867a01d947, 0xe9dfaae7f00dad2b, 0xe26c82849e4c8d82, 0x548e8f13592972bf, 0xe09980facda8bd56, 0x655732cf6b322e2a, 0x1e } },
  { 0x1506283030b9021b, { 0xb8865defe4c2d218, 0x2a151d91ac334d1b, 0x3daf9e243ebe390, 0x1d146705d0ba0cdd, 0xd5dc905e1f294145, 0x9acb7040c9b26bbc, 0x74 } },
  { 0x96919832a48865aa, { 0xedfcfe6a9c9b3e1, 0x8b3f4d3bdcd398b7, 0x2e276f36aceef1fe, 0xd04430fb774598d, 0x9702b68097ef5850, 0xed4980fb9a3a219f, 0xb5 } },
  { 0x1b02336ba6ec0e49, { 0xd45b9f78576231c6, 0x7bf998d4a13856c9, 0xf4d276c1d6382ad8, 0x338b14d60b3b292e, 0x3535b163820c6ee2, 0x2794b790e4e1db2, 0x2c } },
  { 0xd3f9566375eb3ba0, { 0xc794dcef959d6500, 0xa7dea40deb2903c7, 0x421fd37c86eda793, 0x8c537b84a833b498, 0x2fbe05cd41efc5d, 0x168e7da68558d46e, 0x2b } },
  { 0xd45b9f78576231c6, { 0xd3f9566375eb3ba0, 0x413c4146d28e955d, 0x304edbb598d359bf, 0x57ca7ed0874728ed, 0x87f6d2a2c76c37a2, 0x1e0879411a4eba26, 0x3d } },
  { 0xb8865defe4c2d218, { 0xe218c3810a0a335e, 0xf88414af22535194, 0x258ec719fa13e1df, 0x2868b4714bf63a00, 0xda596d121063054b, 0xf8aed25bfa175d98, 0x70 } },
  { 0xbef75201560d3bf4, { 0x1b02336ba6ec0e49, 0x98f92b2ae859042a, 0x8a7be8bb1997584d, 0xc6c39090c5f20630, 0xb54432abd7b6e08c, 0xce7bc0ca6493ef2a, 0xb5 } },
  { 0x23ab7cdb6f8d5d1d, { 0xac528cbfe756fd8b, 0xf1dfea125a9e16e4, 0x1e92bc47edf0a06, 0x19e64d24300ec009, 0xb9dbcf50d20d7f06, 0x1c8ddac911d4189e, 0x86 } },
  { 0xbafae18cbf62cece, { 0x96919832a48865aa, 0x2db300189e97da69, 0x85ef4ac6d8d219ae, 0xbdd5af110fcf7705, 0xaf160bdd62518dba, 0x9c09ed9c4520c82c, 0xad } },
  { 0x3339a3d3cbba52d9, { 0xee8b8ddb19bd5a22, 0xcfff1c248d093b1e, 0xa2d7fd12b026a92b, 0x68c740f988cf5139, 0x7343e75763fb7ef, 0xab4b02f64bf45ea3, 0xa9 } },
  { 0xcc9e5e9e7e874b67, { 0x1506283030b9021b, 0x27a075b325f35389, 0xf8c23979af5e3b6c, 0xb9d1722bb7290be0, 0x62701c314b310b1d, 0xe1d2742a872bd1d3, 0x2e } },
  { 0xadfab22c57e0c1c5, { 0x392d04e8d7d2cb47, 0xd2b4260525890e01, 0x83459099e83ce77f, 0x2ad42bd6c69531ea, 0xdcc0c779bec5130d, 0x357253d3111b913f, 0x5f } },
  { 0x88cfc28bd145cf9a, { 0x256a4652e41e0192, 0x8b8f60e4efa4499d, 0xb27d1f558e2c5703, 0xd8dc0bcc402aae01, 0xe368eea4b7b18b92, 0x7f6a0616856c7455, 0xc5 } },
  { 0xa16811b7b02a8851, { 0xf03b3b4cabf6e89a, 0x411e73940733376d, 0xc21b4bfc1e7ccb7f, 0xe46093fe23f2bcfc, 0x5d700e6a4fd491b2, 0x123dbcaa74eef5a6, 0xc8 } },
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
