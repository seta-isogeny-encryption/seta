// p = 2·6334792777^12 - 1
//
// p-1 = 2^6 * 3^3 * 5 * 7 * 13 * 17 * 53 * 109 * 199 * 271 * 277 * 433 * 461 * 547 * 3373 * 10181 * 16273 * 104173 * 600109 * 2229307 * 3254137 * 3290561 * 5790487 * 6109057 * 9918889 * 12011149 * 20650789 * 24685489 * 37015501
// p+1 = 2 * 5939^12 * 1066643^12


#include "constants.h"

const long class_mod_4 = 1;
const long two_tors_height = 6;

const long security_level = 128;

const long signing_length=1000;
const long signing_length_two_tors_height_step = 31;
const long last_step_length = 10;


const char* p_str = "8352547421173864960624333703667331271395143445202337668463192809880783896176785654136882064526006534223584195029346241";
const char* all_the_torsion_str = "69765048422958181897960319077352793214907650293450379167416508328959226975464456756427967097183686868327009326439047901292946921759400717445370636114701166739708048309054296279586015687227479047950787096565629295598124341383191860830080";

const uintbig p_plus_odd_cofactor = { 2, 0, 0, 0, 0, 0, 0 };
const uintbig p_minus_odd_cofactor = { 0x40, 0, 0, 0, 0, 0, 0 };
const uintbig p_even_cofactor = { 0x5a1aa42506e9d69f, 0x20c60a15121e788a,
  0x2c06da17f90f38dd, 0xe00fc84ae2f8f203, 0x64370a37bb5ac7bf,
  0x4fee61f072e38f16, 0x3 };

#define M_LEN 28
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  { 3, 5, 7, 13, 17, 53, 109, 199, 271, 277, 433, 461, 547, 3373, 10181,
    16273, 104173, 600109, 2229307, 3254137, 3290561, 5790487, 6109057,
    9918889, 12011149, 20650789, 24685489, 37015501 };
const long p_minus_mult[M_LEN] =
  { 3, 1, 1,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,    1,     1,
        1,      1,      1,       1,       1,       1,       1,       1,
          1,        1,        1,        1,        1 };

#define P_LEN 2
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] = { 5939, 1066643 };
const long p_plus_mult[P_LEN] = {   12,      12 };


const long N1_mult_minus[M_LEN] =
  { 0, 0, 0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,    0,     0,
        0,      0,      0,       0,       0,       0,       0,       0,
          0,        0,        0,        0,        0 };
const long N1_mult_plus[P_LEN] =
  {    0,      12 };

const long N2_mult_minus[M_LEN] =
  { 3, 1, 1,  1,  1,  1,   0,   1,   1,   1,   1,   1,   1,    1,     1,
        1,      1,      1,       1,       1,       1,       1,       1,
          1,        1,        1,        1,        0 };
const long N2_mult_plus[P_LEN] =
  {   12,       0 };

/* // the multiplicities to take to obtain log2(p) bits of torsion (for commitment) */
/* const long p_minus_mult_com[M_LEN] = */
/*   { 0,  1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, */
/*        1,    1,    1,    1,    1 }; */
/* const long p_plus_mult_com[P_LEN] = */
/*   { 0, 2,  1,  1,  1,   1,   1,   1,   1,    1,    1,    1 }; */


/* // the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge) */
/* const long p_minus_mult_cha[M_LEN] = */
/*   { 53,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, */
/*        0,    0,    0,    0,    0 }; */
/* const long p_plus_mult_cha[P_LEN] = */
/*   { 21, 0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0 }; */
