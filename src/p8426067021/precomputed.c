
#include <assert.h>
#include "precomputed.h"
#include "curve.h"
#include "uintbig.h"
#include "precomputed_generated.h"

struct precomp_struct global_setup;

void init_precomputations() {
  init_curve(&global_setup.E0);
  init_precomputations_generated();
  
   	long var = fetch_var(); 
     GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC); 
     GEN p = strtoi(p_str); 
  /*   GEN B = alg_hilbert(nf, stoi(-1), negi(p), var, 0); */

    long q = q_norm;
    long c = -1;
    GEN multtable = mkvec4(
    mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
    mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q,0)),
    mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(p),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(p),gen_0,gen_0)),
    mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q,0), mkcol4(gen_0,p,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q),p)),gen_0,gen_0,gen_0))
    );

    GEN B = alg_csa_table(nf, multtable, var,0);

     GEN B_1 = mkcol4s(1,0,0,0); 
     GEN B_i = mkcol4s(0,1,0,0); 
     GEN B_j = mkcol4s(0,0,1,0); 
     GEN B_ji = mkcol4s(0,0,0,1); 


  /*   GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2 */
  /*   GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2 */

     global_setup.p = p; 
     global_setup.q = q;
     global_setup.B = B; // the quaternion algebra 
     global_setup.qf = mkmat4(mkcol4s(1,0,0,0), 
                              mkcol4s(0,q,0,0), 
                              mkcol4(gen_0,gen_0,p,gen_0), 
                              mkcol4(gen_0,gen_0,gen_0,gmul(p,stoi(q)))); // quadratic form defined by the reduced norm 

     // global_setup.torsion_fm = Z_factor_limit(strtoi( 
     //    all_the_torsion_str 
     //    ), 70000000); 

     global_setup.gen_p_plus_fact = cgetg(3, t_MAT); 
     gel(global_setup.gen_p_plus_fact,1) = cgetg(p_plus_len+1, t_COL); 
     gel(global_setup.gen_p_plus_fact,2) = cgetg(p_plus_len+1, t_COL); 
     global_setup.gen_p_minus_fact = cgetg(3, t_MAT); 
     gel(global_setup.gen_p_minus_fact,1) = cgetg(p_minus_len+1, t_COL); 
     gel(global_setup.gen_p_minus_fact,2) = cgetg(p_minus_len+1, t_COL); 

     global_setup.gen_p_plus_primary = cgetg(p_plus_len+1, t_VEC); 
     global_setup.gen_p_minus_primary = cgetg(p_minus_len+1, t_VEC); 

     for (int i = 0; i < p_plus_len; ++i) { 
         gel(gel(global_setup.gen_p_plus_fact,1),i+1) = stoi(p_plus_fact[i]); 
         gel(gel(global_setup.gen_p_plus_fact,2),i+1) = stoi(p_plus_mult[i]); 
         gel(global_setup.gen_p_plus_primary,i+1) = powuu(p_plus_fact[i],p_plus_mult[i]); 
     } 

     for (int i = 0; i < p_minus_len; ++i) { 
         gel(gel(global_setup.gen_p_minus_fact,1),i+1) = stoi(p_minus_fact[i]); 
         gel(gel(global_setup.gen_p_minus_fact,2),i+1) = stoi(p_minus_mult[i]); 
         gel(global_setup.gen_p_minus_primary,i+1) = powuu(p_minus_fact[i],p_minus_mult[i]); 
     } 

     global_setup.gen_odd_torsion = gmul(ZV_prod(global_setup.gen_p_plus_primary), 
     	ZV_prod(global_setup.gen_p_minus_primary)); 

	GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
    GEN B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    GEN B3 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));
    GEN B4 = mkcol4(gen_0,gdiv(gen_1,stoi(q)),gen_0,gneg(gdiv(stoi(c),stoi(q))));

    global_setup.O0 = alglathnf(B,mkmat4(B1,B2,B3,B4), gen_0); // the cannonical maximal order 
  /*   global_setup.O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0); */


     global_setup.one = B_1; 
     global_setup.i = B_i; 
     global_setup.j = B_j; 
     global_setup.ji = B_ji; 

     global_setup.O0_b1 = B1; 
     global_setup.O0_b2 = B2; 
     global_setup.O0_b3 = B3; 
     global_setup.O0_b4 = B4; 
     global_setup.O0_to_standard = mkmat4(B1, B2, B3, B4); 
     global_setup.standard_to_O0 = RgM_inv(global_setup.O0_to_standard); 

     global_setup.smooth_famat_for_klpt =Z_factor_limit(global_setup.gen_odd_torsion, 7045009+1);
         // output(global_setup.smooth_famat_for_klpt);

  /*   // output(alg_O0_to_standard(mkcol4s(1,0,0,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,1,0,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,0,1,0))); */
  /*   // output(alg_O0_to_standard(mkcol4s(0,0,0,1))); */

}


long ell_to_index(long ell, bool *twist) {
	switch (ell) {
		case 3: *twist = true; return 0;
		case 43: *twist = true; return 1;
		case 257: *twist = true; return 2;
		case 84719: *twist = true; return 3;
		case 5: *twist = false; return 0;
		case 7: *twist = false; return 1;
		case 13: *twist = false; return 2;
		case 17: *twist = false; return 3;
		case 19: *twist = false; return 4;
		case 23: *twist = false; return 5;
		case 41: *twist = false; return 6;
		case 73: *twist = false; return 7;
		case 313: *twist = false; return 8;
		case 1009: *twist = false; return 9;
		case 2857: *twist = false; return 10;
		case 3733: *twist = false; return 11;
		case 5519: *twist = false; return 12;
		case 6961: *twist = false; return 13;
		case 53113: *twist = false; return 14;
		case 499957: *twist = false; return 15;
		case 763369: *twist = false; return 16;
		case 2101657: *twist = false; return 17;
		case 2616791: *twist = false; return 18;
		case 7045009: *twist = false; return 19;
		case 11959093: *twist = false; return 20;
		case 17499277: *twist = false; return 21;
		case 20157451: *twist = false; return 22;
		case 30298249: *twist = false; return 23;
		case 33475999: *twist = false; return 24;
		case 39617833: *twist = false; return 25;
		case 45932333: *twist = false; return 26;
	}
	assert(0);
	return 0;
}

long ell_to_e(long ell) {
	switch (ell) {
		case 3: return 23;
		case 43: return 12;
		case 257: return 12;
		case 84719: return 12;
		case 2: return 5;
		default: return 1;
	}
}
