// p = 2·8426067021^12 - 1
//
// p-1 = 2^5 * 5 * 7 * 13 * 17 * 19 * 23 * 41 * 73 * 313 * 1009 * 2857 * 3733 * 5519 * 6961 * 53113 * 499957 * 763369 * 2101657 * 2616791 * 7045009 * 11959093 * 17499277 * 20157451 * 30298249 * 33475999 * 39617833 * 45932333
// p+1 = 2 * 3^24 * 43^12 * 257^12 * 84719^12


#include "constants.h"

const long class_mod_4 = 1;
const long two_tors_height = 5;

const long security_level = 128;

const long signing_length=1000;
const long signing_length_two_tors_height_step = 31;
const long last_step_length = 10;


const char* p_str = "256170376103179994761914482856266627791915646373582727175785233336597312364178751161184887957805876215573849680566141281";
const char* all_the_torsion_str = "65623261592844692118321438425339452115173228005679839727238382594155983305126436954282908454249767432625394800791080317067422349066481244932552342489741556853813851868814195191266951532939668125495137254071411275260925216193796110052320960";

const uintbig p_plus_odd_cofactor = { 6, 0, 0, 0, 0, 0, 0 };
const uintbig p_minus_odd_cofactor = { 0x20, 0, 0, 0, 0, 0, 0 };
const uintbig p_even_cofactor = { 0x7f93d1742150ce0b, 0xc902188c8f3cc043,
  0x104232a1eae7c0f7, 0x2d9f078b7b9f9fd6, 0xa3e52553cd09ac6b,
  0x2ba425a5b30354c0, 0xcb };


#define M_LEN 27
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  { 5, 7, 13, 17, 19, 23, 41, 73, 313, 1009, 2857, 3733, 5519, 6961, 53113, 499957,
    763369, 2101657, 2616791, 7045009, 11959093, 17499277, 20157451, 30298249, 33475999,
    39617833, 45932333 };
const long p_minus_mult[M_LEN] =
  { 1, 1,  1,  1,  1,  1,  1,  1,   1,    1,    1,    1,    1,    1,     1,      1,
         1,       1,       1,       1,        1,        1,        1,        1,        1,
     1,        1 };

#define P_LEN 4
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] = {  3, 43, 257, 84719 };
const long p_plus_mult[P_LEN] = { 23, 12,  12,    12 };

// const long N1_mult_minus[M_LEN] =
//   { 0, 0,  0,  0,  0,  0,  0,  0,   0,    0,    0,    0,    0,    0,     0,      0,
//          0,       0,       0,       0,        0,        0,        0,        0,        0,
//      0,        0 };
// const long N1_mult_plus[P_LEN] =
//   {  0, 12,   0,    10 };

// const long N2_mult_minus[M_LEN] =
//   { 1, 1,  1,  1,  1,  1,  1,  1,   1,    1,    1,    1,    1,    1,     1,      1,
//          1,       1,       1,       1,        1,        1,        1,        1,        1,
//      1,        0 };
// const long N2_mult_plus[P_LEN] =
//   { 22,  0,  12,     0 };

const long N1_mult_minus[M_LEN] =
  { 0, 0,  0,  0,  0,  0,  0,  0,   0,    0,    0,    0,    0,    0,     0,      0,
         0,       0,       0,       0,        0,        0,        0,        0,        0,
     0,        0 };
const long N1_mult_plus[P_LEN] =
  {  0, 12,   0,    11 };

const long N2_mult_minus[M_LEN] =
  { 1, 1,  1,  1,  1,  1,  0,  1,   1,    1,    1,    1,    1,    1,     1,      1,
         1,       1,       1,       1,        1,        1,        1,        0,        1,
     1,        1 };
const long N2_mult_plus[P_LEN] =
  { 21,  0,  12,     0 };


// // SMALLER VALUES FOR TESTING
// const long N1_mult_minus[M_LEN] =
//   { 0, 0,  0,  0,  0,  0,  0,  0,   0,    0,    0,    0,    0,    0,     0,      0,
//          0,       0,       0,       0,        0,        0,        0,        0,        0,
//      0,        0 };
// const long N1_mult_plus[P_LEN] =
//   {  0, 11,   0,    0 };

// const long N2_mult_minus[M_LEN] =
//   { 1, 1,  1,  1,  1,  1,  1,  1,   1,    1,    1,    1,    1,    1,     1,      0,
//          0,       0,       0,       0,        0,        0,        0,        0,        0,
//      0,        0 };
// const long N2_mult_plus[P_LEN] =
//   { 23,  0,  12,     0 };










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
