
#ifndef SETA_H
#define SETA_H

#include <pari/pari.h>
#include "idiso.h"

#define C_LEN 27
#define T_LEN 4

// a point with a part on the curve and a part on the twist
typedef struct hybrid_point {
    proj E;
    long twist_mult[T_LEN];
    long curve_mult[C_LEN];
    proj twist_pt;
    proj curve_pt;
} hybrid_point;

typedef struct public_param {
  hybrid_point N1_basis[3];
  hybrid_point N2_basis[3];

  GEN N1_twist;
  GEN N1_curve;
  GEN N1_primary_twist; // primary factorisation
  GEN N1_primary_curve;

  GEN N2_twist;
  GEN N2_curve;
  GEN N2_primary_twist;
  GEN N2_primary_curve; 

  GEN D; 
  GEN d; 
  GEN fac_det; 
} public_param;

typedef struct public_key {
  hybrid_point N1_basis[3];
  hybrid_point N2_basis[3];
} public_key;

typedef struct ciphertext {
  hybrid_point N2_basis[3];
} ciphertext;


typedef struct secret_key {
    uintbig ker_psi_plus_1;
    uintbig ker_psi_plus_2;
    uintbig ker_psi_minus_1;
    uintbig ker_psi_minus_2;
    uintbig ker_psi_dual_plus_1;
    uintbig ker_psi_dual_plus_2;
    uintbig ker_psi_dual_minus_1;
    uintbig ker_psi_dual_minus_2;
} secret_key;


extern const long N1_len;
extern const long N1_fact[];
extern const long N1_mult[];

extern const long N2_len;
extern const long N2_fact[];
extern const long N2_mult[];



// Compute basis of the N1 and N2 torsion
void seta_setup(public_param *param);

void seta_genkey(public_key *pk, secret_key *sk, const public_param *param);

// Only for testing
void seta_genkey_invalid_for_testing(public_key *pk, secret_key *sk, const public_param *param);

void seta_eval(ciphertext *ctxt, const public_key *pk, const uintbig *x);

void seta_inverse(odd_isogeny *phi_m, const ciphertext *ctxt, const public_key *pk, const secret_key *sk, const public_param *param);


#endif
