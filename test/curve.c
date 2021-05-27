#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "mont.h"
#include "constants.h"
#include "curve.h"

void test_frob(proj Q, const proj *A) {
  proj Ed, P; proj2 PP, QQ; point Pe, Qe;
  mont_to_ted(&Ed, A, false);

  P = Q;
  xtoxy(&PP, A, &P); QQ = PP;
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;
  
  // π² = -p = ±1 (depending on p mod 4)
  mont0_frob(&P, &P);
  montxy0_frob(&PP, &PP);
  ted0_frob(&Pe, &Pe);
  assert(is_on_curve(&P, A));
  assert(xy_is_on_curve(A, &PP));
  assert(ted_is_on_curve(&Pe, &Ed));

  mont0_frob(&P, &P);
  montxy0_frob(&PP, &PP);
  ted0_frob(&Pe, &Pe);
  if (!curve_order_is_p_plus_one) {
    xyNEG(&PP, &PP);
    ted_neg(&Pe, &Pe);
  }
  assert(mont_equal(&P, &Q));
  assert(xy_equal(&PP, &QQ));
  assert(ted_equal(&Pe, &Pe));
}

void test_frob_twist(proj Q, const proj *A) {
  proj Ed, P; point Pe, Qe;
  mont_to_ted(&Ed, A, true);

  P = Q;
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;
  
  // π² = -p = ±1 (depending on p mod 4)
  mont0_frob(&P, &P);
  ted0_frob_twist(&Pe, &Pe);
  assert(!is_on_curve(&P, A));
  assert(ted_is_on_curve(&Pe, &Ed));

  mont0_frob(&P, &P);
  ted0_frob_twist(&Pe, &Pe);
  if (curve_order_is_p_plus_one) {
    ted_neg(&Pe, &Pe);
  }
  assert(mont_equal(&P, &Q));
  assert(ted_equal(&Pe, &Qe));

  // test coherence with ted0_frob
  //    π
  // E -----> E       : ax² + y² = 1 + d(xy)²
  // ^        |/t
  // |*t      v
  // E* ----> E*      : Bax² + y² = 1 + Bd(xy)² ,   B = t² ,  t^(p-1) = B^(p-1)/2
  //
  //  frob_twist = B^(p-1)/2 ~ frob
  //
  uintbig p12;
  fp2 twister = fp2_non_residue();
  uintbig_sub3(&p12, &p, &uintbig_1);
  uintbig_div3_64(&p12, &p12, 2);
  fp2_exp(&twister, &twister, &p12);
  ted0_frob(&Pe, &Pe);
  fp2_mul2(&Pe.x, &twister);
  fp2_mul2(&Pe.t, &twister);
  ted0_frob_twist(&Qe, &Qe);
  assert(ted_equal(&Pe, &Qe));
}

void test_dist(proj Q, const proj *A, bool twist) {
  proj Ed, P; proj2 PP, QQ, tmp; point Pe, Qe, tmpe;
  mont_to_ted(&Ed, A, twist);

  P = Q;
  if (!twist) {
    xtoxy(&PP, A, &P);
    QQ = PP;
  }
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;
				    
  // ι² - tι + n = 0   (t = dist_trace, n = dist_norm)
  mont0_dist(&P, &P);
  if (!twist) montxy0_dist(&PP, &PP);
  ted0_dist(&Pe, &Pe);
  assert(is_on_curve(&P, A) != twist);
  assert(twist || xy_is_on_curve(A, &PP));
  assert(ted_is_on_curve(&Pe, &Ed));

  if (!twist) {
    xytox(&Q, &PP);
    assert(mont_equal(&Q, &P));
  }
  ted_to_mont_point(&Q, &Pe);
  assert(mont_equal(&Q, &P));

  int i;
  if (dist_trace > 0) {
    i = dist_trace;
    if (!twist) xyNEG(&tmp, &QQ);
    ted_neg(&tmpe, &Qe);
  } else {
    i = -dist_trace;
    tmp = QQ;
    tmpe = Qe;
  }
  for (; i > 0; i--) {
    if (!twist) xyADD(&PP, A, &PP, &tmp);
    ted_add(&Pe, &Ed, &Pe, &tmpe);
  }
  
  if (!twist) montxy0_dist(&PP, &PP);
  ted0_dist(&Pe, &Pe);
  
  if (dist_norm < 0) {
    i = -dist_norm;
    if (!twist) xyNEG(&tmp, &QQ);
    ted_neg(&tmpe, &Qe);
  } else {
    i = dist_norm;
    tmp = QQ;
    tmpe = Qe;
  }
  for (; i > 0; i--) {
    if (!twist) xyADD(&PP, A, &PP, &tmp);
    ted_add(&Pe, &Ed, &Pe, &tmpe);
  }
  assert(twist || xy_is_zero(&PP));
  assert(ted_iszero(&Pe));
}

void test_anticomm(proj P, const proj *A, bool twist) {
  proj Ed; proj2 PP, QQ, tmp; point Pe, Qe, tmpe;
  mont_to_ted(&Ed, A, twist);

  if (!twist) {
    xtoxy(&PP, A, &P);
    QQ = PP;
  }
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;

  // ιπ = π(t - ι)   (t = dist_trace)
  if (!twist) montxy0_frob(&PP, &PP);
  twist ? ted0_frob_twist(&Pe, &Pe) : ted0_frob(&Pe, &Pe);
  mont0_dist(&P, &P);
  if (!twist) montxy0_dist(&PP, &PP);
  ted0_dist(&Pe, &Pe);
  assert(is_on_curve(&P, A) != twist);
  assert(twist || xy_is_on_curve(A, &PP));
  assert(ted_is_on_curve(&Pe, &Ed));

  tmp = QQ;
  tmpe = Qe;
  if (!twist) montxy0_dist(&QQ, &QQ);
  ted0_dist(&Qe, &Qe);
  xyNEG(&QQ, &QQ);
  ted_neg(&Qe, &Qe);

  int i;
  if (dist_trace < 0) {
    i = -dist_trace;
    if (!twist) xyNEG(&tmp, &tmp);
    ted_neg(&tmpe, &tmpe);
  } else {
    i = dist_trace;
  }
  for (; i > 0; i--) {
    if (!twist) xyADD(&QQ, A, &QQ, &tmp);
    ted_add(&Qe, &Ed, &Qe, &tmpe);
  }

  if (!twist) montxy0_frob(&QQ, &QQ);
  twist ? ted0_frob_twist(&Qe, &Qe) : ted0_frob(&Qe, &Qe);

  assert(twist || xy_equal(&PP, &QQ));
  assert(ted_equal(&Pe, &Qe));
}

void test_comm(proj Q, const proj *A, int k, bool twist) {
  proj Ed, P; proj2 PP, QQ; point Pe, Qe;
  mont_to_ted(&Ed, A, twist);

  P = Q;
  if (!twist) {
    xtoxy(&PP, A, &P);
    QQ = PP;
  }
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;
  
  uintbig mul;
  uintbig_set(&mul, k);

  mont0_dist(&P, &P);
  if (!twist) montxy0_dist(&PP, &PP);
  ted0_dist(&Pe, &Pe);
  mont0_frob(&P, &P);
  if (!twist) montxy0_frob(&PP, &PP);
  twist ? ted0_frob_twist(&Pe, &Pe) : ted0_frob(&Pe, &Pe);
  xMUL(&P, A, &P, &mul);
  if (!twist) xyMUL(&PP, A, &PP, &mul);
  ted_mul(&Pe, &Pe, &Ed, &mul);
  assert(is_on_curve(&P, A) != twist);
  assert(twist || xy_is_on_curve(A, &PP));
  assert(ted_is_on_curve(&Pe, &Ed));

  xMUL(&Q, A, &Q, &mul);
  if (!twist) xyMUL(&QQ, A, &QQ, &mul);
  ted_mul(&Qe, &Qe, &Ed, &mul);
  mont0_dist(&Q, &Q);
  if (!twist) montxy0_dist(&QQ, &QQ);
  ted0_dist(&Qe, &Qe);
  mont0_frob(&Q, &Q);
  if (!twist) montxy0_frob(&QQ, &QQ);
  twist ? ted0_frob_twist(&Qe, &Qe) : ted0_frob(&Qe, &Qe);
  assert(mont_equal(&P, &Q));
  assert(twist || xy_equal(&PP, &QQ));
  assert(ted_equal(&Pe, &Pe));
}

void test_sqrt_minus_q(proj Q, const proj *A, bool twist) {
  proj Ed, P; point Pe, Qe;
  mont_to_ted(&Ed, A, twist);

  P = Q;
  mont_to_ted_point(&Pe, A, &P); Qe = Pe;

  uintbig q;
  uintbig_set(&q, q_norm);
  
  // π² = -q
  mont0_sqrt_minus_q(&P, &P);
  ted0_sqrt_minus_q(&Pe, &Ed, &Pe);
  assert(is_on_curve(&P, A) != twist);
  assert(ted_is_on_curve(&Pe, &Ed) == ted_is_on_curve(&Qe, &Ed));

  mont0_sqrt_minus_q(&P, &P);
  ted0_sqrt_minus_q(&Pe, &Ed, &Pe);
  ted_neg(&Pe, &Pe);
  xMUL(&Q, A, &Q, &q);
  ted_mul(&Qe, &Qe, &Ed, &q);
  assert(mont_equal(&P, &Q));
  assert(ted_equal(&Pe, &Qe));

  // Check coherence with dist     (is this really necessary?)
  // 2dist = -(√-q - dist_trace)   (the other choice would be √-q + dist_trace)
  point tmpe = Qe = Pe;
  ted0_sqrt_minus_q(&Qe, &Ed, &Qe);
  int i;
  if (dist_trace < 0) {
    i = -dist_trace;
  } else {
    i = dist_trace;
    ted_neg(&tmpe, &tmpe);
  }
  for (; i > 0; i--)
    ted_add(&Qe, &Ed, &Qe, &tmpe);
  ted_neg(&Qe, &Qe); // -(√-q - t)P
  
  ted0_dist(&Pe, &Pe);
  ted_double(&Pe, &Ed, &Pe); // 2*dist(P)
  assert(ted_equal(&Pe, &Qe));
}

int main() {
  srand48(1);

  proj A; init_curve(&A);
  proj P;
  
  for (int i = 0; i < 50; i++) {
    do {
      fp2_random(&P.x); fp2_random(&P.z);
    } while (!is_on_curve(&P, &A));
    
    test_frob(P, &A);
    test_dist(P, &A, false);
    test_anticomm(P, &A, false);
    test_comm(P, &A, i+2, false);
    test_sqrt_minus_q(P, &A, false);

    do {
      fp2_random(&P.x); fp2_random(&P.z);
    } while (is_on_curve(&P, &A));

    test_frob_twist(P, &A);
    test_dist(P, &A, true);
    test_anticomm(P, &A, true);
    test_comm(P, &A, i+2, true);
    test_sqrt_minus_q(P, &A, true);
  }

  printf("    \033[1;32mAll tests passed\033[0m\n");
  exit(0);
}
