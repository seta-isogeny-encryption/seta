#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "uintbig.h"

// The class of p mod 4, and its consequence for the order of curves
extern const long class_mod_4;
#define curve_order_is_p_plus_one (class_mod_4 == 3)

extern const long two_tors_height;
// The cofactor of 2^two_tors_height in p±1
extern const uintbig p_even_cofactor;
extern const long security_level;

// the signing isogeny has degree 2^signing_length
extern const long signing_length;
// we have the equality signin_length = two_tors_height * (signing_length_two_tors_height_step -1 ) + last_step_length
extern const long signing_length_two_tors_height_step;
extern const long last_step_length;

extern const char* p_str;
extern const char* all_the_torsion_str;

// The useful odd factors in p-1
extern const long p_minus_len;
extern const long p_minus_fact[];
extern const long p_minus_mult[];
// The cofactor of the useful odd torsion in p-1
extern const uintbig p_minus_odd_cofactor;

// The useful odd factors in p+1
extern const long p_plus_len;
extern const long p_plus_fact[];
extern const long p_plus_mult[];
// The cofactor of the useful odd torsion in p+1
extern const uintbig p_plus_odd_cofactor;

// Same as above, but along the curve/twist dichotomy
#define on_curve_len (curve_order_is_p_plus_one ? p_plus_len : p_minus_len)
#define on_curve_fact (curve_order_is_p_plus_one ? p_plus_fact : p_minus_fact)
#define on_curve_mult (curve_order_is_p_plus_one ? p_plus_mult : p_minus_mult)
#define on_curve_odd_cofactor (curve_order_is_p_plus_one ? p_plus_odd_cofactor : p_minus_odd_cofactor)
#define on_twist_len (!curve_order_is_p_plus_one ? p_plus_len : p_minus_len)
#define on_twist_fact (!curve_order_is_p_plus_one ? p_plus_fact : p_minus_fact)
#define on_twist_mult (!curve_order_is_p_plus_one ? p_plus_mult : p_minus_mult)
#define on_twist_odd_cofactor (!curve_order_is_p_plus_one ? p_plus_odd_cofactor : p_minus_odd_cofactor)

// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
extern const long p_minus_mult_com[];
extern const long p_plus_mult_com[];

// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
extern const long p_minus_mult_cha[];
extern const long p_plus_mult_cha[];

// Séta constants
extern const long N1_mult_minus[];
extern const long N1_mult_plus[];
extern const long N2_mult_minus[];
extern const long N2_mult_plus[];

// Same as above, but along the curve/twist dichotomy
#define N1_mult_curve (curve_order_is_p_plus_one ? N1_mult_plus : N1_mult_minus)
#define N1_mult_twist (!curve_order_is_p_plus_one ? N1_mult_plus : N1_mult_minus)
#define N2_mult_curve (curve_order_is_p_plus_one ? N2_mult_plus : N2_mult_minus)
#define N2_mult_twist (!curve_order_is_p_plus_one ? N2_mult_plus : N2_mult_minus)

#endif
