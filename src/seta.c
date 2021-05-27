#include <stdint.h>
//#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>
#include <gmp.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "qfsolve_factor.h"
#include "printing.h"


#define FP_LIMBS (4 * 64 / GMP_LIMB_BITS)

#include "seta.h"
// #include "isomorphism.h"
// #include "uintbig.h"





static isog_degree mont_order(const proj *P, const proj *A, bool twist) {
  
    const long *fact, *mult;
    long len;
      
    if(!twist) {
        fact = on_curve_fact; mult = on_curve_mult;
        len = on_curve_len;
    }
    else {
        fact = on_twist_fact; mult = on_twist_mult;
        len = on_twist_len;
    }

  proj tmp;
  uintbig cof;
  isog_degree deg = degree_co((isog_degree){ 0 }, mult, len);
  for (int j = 0; j < len; j++) {
    degree_unset(&deg, j);
    degree_to_uint(&cof, deg, fact, len);
    xMUL(&tmp, A, P, &cof);
    uintbig_set(&cof, fact[j]);
    uint8_t v = 0;
    for ( ; !mont_iszero(&tmp); v++) {
      xMUL(&tmp, A, &tmp, &cof);
    }
    degree_set(&deg, j, v);
  }
  return deg;
}


static GEN alg_standard_to_O0(GEN elt) {
    return RgM_RgC_mul(global_setup.standard_to_O0, elt);
}


// static void fp2_random_replayable(fp2 *x) {
//   uint64_t thrash;
//   x->re.x.c[0] = random_Fl(0xffffffffffffffff);
//   x->re.x.c[1] = random_Fl(0xffffffffffffffff);
//   x->re.x.c[2] = random_Fl(0xffffffffffffffff);
//   x->re.x.c[3] = random_Fl(0xffffffffffffffff);
//   mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->re.x.c, 0, (mp_srcptr)x->re.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
//   x->im.x.c[0] = random_Fl(0xffffffffffffffff);
//   x->im.x.c[1] = random_Fl(0xffffffffffffffff);
//   x->im.x.c[2] = random_Fl(0xffffffffffffffff);
//   x->im.x.c[3] = random_Fl(0xffffffffffffffff);
//   mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->im.x.c, 0, (mp_srcptr)x->im.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
// }


static void random_point(proj *P, proj const *A, long ell, long e, bool twist) {
    uintbig cofactor;
    if ((!twist && curve_order_is_p_plus_one) || (twist && !curve_order_is_p_plus_one)) 
        uintbig_add3(&cofactor, &p, &uintbig_1);
    else { uintbig_sub3(&cofactor, &p, &uintbig_1); }
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    for (int i = 0; i < e; ++i) {
        uintbig_div3_64(&cofactor, &cofactor, ell); 
    }
    proj Z;

    while (1) {
        fp2_random(&P->x); P->z = fp2_1;
        if (twist == is_on_curve(P, A)) continue;
        xMUL(P, A, P, &cofactor);
        Z = *P;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&Z, A, &Z, &ell_big);
        }
        if (!fp2_iszero(&Z.z)) { 
            // a final test
            xMUL(&Z, A, &Z, &ell_big);
            assert(fp2_iszero(&Z.z));
            return;
        }
    }
}


static void random_basis(proj *P1, proj *P2, point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
    point P1_mul, P2_mul, tmp;
    uintbig ell_big;
    uintbig_set(&ell_big, ell);
    fp2 weil;

    proj E;
    mont_to_ted(&E, A, twist);

    random_point(P1, A, ell, e, twist);

    mont_to_ted_point(P1_ted, A, P1);
    P1_mul = *P1_ted;
    for (int i = 0; i < e-1; ++i) {
        ted_mul(&P1_mul, &P1_mul, &E, &ell_big);
    }

    assert(ted_is_on_curve(&P1_mul,&E));
    assert(!ted_iszero(&P1_mul));
    ted_mul(&tmp, &P1_mul, &E, &ell_big);
    assert(ted_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        mont_to_ted_point(P2_ted, A, P2);
        P2_mul = *P2_ted;
        for (int i = 0; i < e-1; ++i) {
            ted_mul(&P2_mul, &P2_mul, &E, &ell_big);
        }
        assert(ted_is_on_curve(&P2_mul,&E));
        assert(!ted_iszero(&P2_mul));
        ted_mul(&tmp, &P2_mul, &E, &ell_big);
        assert(ted_iszero(&tmp));

        ted_weil(&weil, &E, &P1_mul, &P2_mul, &ell_big);
        fp2_sub2(&weil, &fp2_1);
    } while (fp2_iszero(&weil));

}


void copy_hybrid_point(hybrid_point *P, const hybrid_point *Q){
    for (long i = 0; i < on_twist_len; ++i) {
        P->twist_mult[i] = Q->twist_mult[i];
    }

    for (long i = 0; i < on_curve_len; ++i) {
        P->curve_mult[i] = Q->curve_mult[i];
    }

    P->E = Q->E;
    P->twist_pt = Q->twist_pt;
    P->curve_pt = Q->curve_pt;
}


void endo_to_ker(uintbig *a, uintbig *b, const GEN endo, long len, const long N_fact[], const long N_mult[], const GEN N_primary) {

    GEN ker_1 = zerovec(len), ker_2 = zerovec(len);
    GEN m1,m2,m3,m4;
    uintbig_set(a, 0);
    uintbig_set(b, 0);

    for (int i = 0; i < len; ++i) {
        long ell = N_fact[i], e = N_mult[i];
        action_from_elle(&m1, &m2, &m3, &m4, ell, e);
        GEN ker_ell = endo_to_kernel_action(endo, m1, m2, m3, m4, ell, e);

        gel(ker_1,i+1) = gel(ker_ell,1);
        gel(ker_2,i+1) = gel(ker_ell,2);

    }

    gentobig(a, ZV_chinese(ker_1, N_primary, NULL));
    gentobig(b, ZV_chinese(ker_2, N_primary, NULL));
}




void seta_genkey(public_key *pk, secret_key *sk, const public_param *param) {

    const hybrid_point *N1_basis = param->N1_basis;
    const hybrid_point *N2_basis = param->N2_basis;
    GEN N1 = gmul(param->N1_twist, param->N1_curve), N2 = gmul(param->N2_twist, param->N2_curve);

    GEN B = global_setup.B;
    GEN O0 = global_setup.O0;
    proj E0 = global_setup.E0;
    clock_t t;



    t = tic();

    GEN qf = mkmat4(mkcol4(gmul(global_setup.p,stoi(global_setup.q)),gen_0,gen_0,gen_0),
                    mkcol4(gen_0,global_setup.p,gen_0,gen_0),
                    mkcol4(gen_0,gen_0,stoi(global_setup.q),gen_0),
                    mkcol4(gen_0,gen_0,gen_0,gneg(param->D)));
    GEN sol = qfsolve_factor(qf, param->fac_det);  

    TOC(t,"quadratic system");
    // output(gmul(param->D,gmul(gsqr(stoi(global_setup.q)),gsqr(global_setup.p))));
    // output(famat_prod(param->fac_det));
    // output(gneg(ZM_det(qf)));
    assert(gcmp(famat_prod(param->fac_det),gneg(ZM_det(qf))) == 0);
    // output(qf);

        // find an order containing (a*i+b*j+c*k)/u, where sol = (a,b,c,u)
    // printf("sol = ");
    // output(sol);
    GEN theta_u = mkcol4(gen_0, gel(sol,3), gel(sol,2), gel(sol,1));
    GEN u = gabs(gel(sol,4),0);

    // output(algnorm(B, theta_u,0));
    // output(u);
    // output(param->D);
    assert(gcmp(algnorm(B, theta_u,0),gmul(gsqr(u),param->D)) == 0);
    assert(alglatcontains(B, O0, theta_u, NULL));

    GEN I = lideal_create(B, O0, theta_u, u); // I = O0*theta*u + O0*u

    // small correction if N(I) < u
    theta_u = gdiv(gmul(theta_u, lideal_norm(I)),u);
    u = lideal_norm(I);


    GEN OL = alglatlefttransporter(B, lideal_lattice(I), lideal_lattice(I));
    GEN OR = alglatrighttransporter(B, lideal_lattice(I), lideal_lattice(I));


    GEN index;
    int subset = alglatsubset(B, OL, O0, &index);
    assert((subset) && (gcmp(index, gen_1) == 0));

    assert(alglatcontains(B, OR, algdivr(B, theta_u,  alg_scalar(B, u)), NULL));


    GEN one = mkcol4s(1,0,0,0);
    GEN O0_ideal = lideal_create(B, O0, one, gen_1);
    printf("KLPT\n");
    GEN J = klpt_general_power(O0_ideal, I, gen_2);


    GEN NJ = lideal_norm(J);
    printf("len %f\n", dbllog2r(itor(NJ,10)));
    GEN v = NJ;


    GEN alpha = lideal_isom(I, J); // I*alpha = J
    assert(alpha);
    //klpt_special_smooth(I, famat_sqr(famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact)));


    GEN eta_v = algmul(B, alg_conj(B,alpha) , algmul(B,theta_u,alpha));

    assert(alglatcontains(B, O0, eta_v, NULL));

    GEN content, N1eta_plus_d_v = algadd(B, gmul(eta_v, N1), alg_scalar(B, gmul(param->d,v)));
    alg_primitive(&content, B, O0, N1eta_plus_d_v);
    printf("content of eta_v+d "); output(content);
    //printf("norm eta_v+d "); output(algnorm(B,N1eta_plus_d_v,0));
    //printf("N2^2 "); output(gmul(gmul(N2,N2),gmul(v,v)));
    


    t = tic();

    two_walk_long phi;
    GEN L;
    special_isogeny phi_L;

    init_trivial_two_walk_long(&phi);

    printf("ideal_to_isogeny_O0_two_long...\n");

    ideal_to_isogeny_O0_two_long(&phi, &L, &phi_L, J, false);

    TOC(t,"ideal_to_isogeny_O0_two_long");



    proj zip[12];
    zip[0] = N1_basis[0].curve_pt; 	zip[1] = N1_basis[1].curve_pt; 	zip[2] = N1_basis[2].curve_pt;
    zip[3] = N2_basis[0].curve_pt; 	zip[4] = N2_basis[1].curve_pt; 	zip[5] = N2_basis[2].curve_pt;
    zip[6] = N1_basis[0].twist_pt; 	zip[7] = N1_basis[1].twist_pt; 	zip[8] = N1_basis[2].twist_pt;
    zip[9] = N2_basis[0].twist_pt; 	zip[10] = N2_basis[1].twist_pt; zip[11] = N2_basis[2].twist_pt;


    proj E_theta = E0;

    eval_walk_long_mult(&phi, &E_theta, zip, 12);


    for (int i = 0; i < 3; ++i) {       
        copy_hybrid_point(pk->N1_basis+i, N1_basis+i);
        copy_hybrid_point(pk->N2_basis+i, N2_basis+i);
        pk->N1_basis[i].E = E_theta;
        pk->N2_basis[i].E = E_theta;
    }

    pk->N1_basis[0].curve_pt = zip[0]; pk->N1_basis[1].curve_pt = zip[1]; pk->N1_basis[2].curve_pt = zip[2]; 
    pk->N2_basis[0].curve_pt = zip[3]; pk->N2_basis[1].curve_pt = zip[4]; pk->N2_basis[2].curve_pt = zip[5]; 
    pk->N1_basis[0].twist_pt = zip[6]; pk->N1_basis[1].twist_pt = zip[7]; pk->N1_basis[2].twist_pt = zip[8]; 
    pk->N2_basis[0].twist_pt = zip[9]; pk->N2_basis[1].twist_pt = zip[10]; pk->N2_basis[2].twist_pt = zip[11]; 



    #ifndef NDEBUG

    isog_degree ord1,ord2; 

    for (int j = 0; j < 3; ++j)
    {
	    degree_one(&ord1);
	    for (long i = 0; i < on_twist_len; ++i) {
	        degree_set(&ord1, i, pk->N1_basis[0].twist_mult[i]);
	    }
	    ord2 = mont_order(&pk->N1_basis[j].twist_pt, &E_theta, true);
	    assert(ord1.val == ord2.val);

	    degree_one(&ord1);
	    for (long i = 0; i < on_curve_len; ++i) {
	        degree_set(&ord1, i, pk->N1_basis[0].curve_mult[i]);
	    }
	    ord2 = mont_order(&pk->N1_basis[j].curve_pt, &E_theta, false);
	    assert(ord1.val == ord2.val);
    }
    #endif




    // compute the action of eta

    GEN endo_v = alg_standard_to_O0(eta_v);



    // compute the kernel of psi (coefficients of linear combination of phi_m(P) and phi_m(Q))
    GEN endo_eta_mod_N2 = gmul(endo_v,ginvmod(v, N2));
    //printf("endo_eta_mod_N2 "); output(endo_eta_mod_N2);



    GEN psi_mod_N2 = gmul(endo_eta_mod_N2,N1);
    gel(psi_mod_N2,1) = gadd(param->d,gel(psi_mod_N2,1));


    GEN psi_dual_mod_N2 = gneg(gmul(endo_eta_mod_N2,N1));
    gel(psi_dual_mod_N2,1) = gadd(param->d,gel(psi_dual_mod_N2,1));

    if(curve_order_is_p_plus_one) {
        endo_to_ker(&sk->ker_psi_plus_1, &sk->ker_psi_plus_2, psi_mod_N2, on_curve_len, on_curve_fact, N2_basis[0].curve_mult, param->N2_primary_curve);
        endo_to_ker(&sk->ker_psi_minus_1, &sk->ker_psi_minus_2, psi_mod_N2, on_twist_len, on_twist_fact, N2_basis[0].twist_mult, param->N2_primary_twist);
    
        endo_to_ker(&sk->ker_psi_dual_plus_1, &sk->ker_psi_dual_plus_2, psi_dual_mod_N2, on_curve_len, on_curve_fact, N2_basis[0].curve_mult, param->N2_primary_curve);
        endo_to_ker(&sk->ker_psi_dual_minus_1, &sk->ker_psi_dual_minus_2, psi_dual_mod_N2, on_twist_len, on_twist_fact, N2_basis[0].twist_mult, param->N2_primary_twist);
    } 
    else {
        endo_to_ker(&sk->ker_psi_minus_1, &sk->ker_psi_minus_2, psi_mod_N2, on_curve_len, on_curve_fact, N2_basis[0].curve_mult, param->N2_primary_curve);
        endo_to_ker(&sk->ker_psi_plus_1, &sk->ker_psi_plus_2, psi_mod_N2, on_twist_len, on_twist_fact, N2_basis[0].twist_mult, param->N2_primary_twist);
    
        endo_to_ker(&sk->ker_psi_dual_minus_1, &sk->ker_psi_dual_minus_2, psi_dual_mod_N2, on_curve_len, on_curve_fact, N2_basis[0].curve_mult, param->N2_primary_curve);
        endo_to_ker(&sk->ker_psi_dual_plus_1, &sk->ker_psi_dual_plus_2, psi_dual_mod_N2, on_twist_len, on_twist_fact, N2_basis[0].twist_mult, param->N2_primary_twist);
    }


}





void seta_genkey_invalid_for_testing(public_key *pk, secret_key *sk, const public_param *param) {

    proj E_theta;

    // define E_theta as the target of a short random walk
    two_walk phi_two;
    random_two_walk(&phi_two);
    dual_walk(&phi_two);
    E_theta = phi_two.A;

    // find a (random) basis of the torsion of E_theta

    proj basis_twist[on_twist_len][3], basis_curve[on_curve_len][3];
    point ted_0 = {fp2_0, fp2_1, fp2_1, fp2_0};
    point basis_twist_ted[on_twist_len][3], basis_curve_ted[on_curve_len][3];

    proj E_ted_twist;
    proj E_ted_curve;
    mont_to_ted(&E_ted_twist, &E_theta, true);
    mont_to_ted(&E_ted_curve, &E_theta, false);


    point pt;

    

    for (int i = 0; i < on_curve_len; i++) {
        long ell = on_curve_fact[i], e = param->N1_basis[0].curve_mult[i];
        // printf("%ld^%ld\n",ell,e);

        if (e == 0) {
          basis_curve_ted[i][0] = ted_0; basis_curve_ted[i][1] = ted_0; basis_curve_ted[i][2] = ted_0;
          basis_curve[i][0].x = fp2_1; basis_curve[i][1].x = fp2_1; basis_curve[i][2].x = fp2_1;
          basis_curve[i][0].z = fp2_0; basis_curve[i][1].z = fp2_0; basis_curve[i][2].z = fp2_0;
          continue;
        }

        random_basis(&basis_curve[i][0], &basis_curve[i][1], &basis_curve_ted[i][0], &basis_curve_ted[i][1], &E_theta, ell, e, false);
        ted_add(&basis_curve_ted[i][2], &E_ted_curve, &basis_curve_ted[i][0], &basis_curve_ted[i][1]);
        ted_to_mont_point(&basis_curve[i][2], &basis_curve_ted[i][2]);
    }


    for (int i = 0; i < on_twist_len; i++) {
        long ell = on_twist_fact[i], e = param->N1_basis[0].twist_mult[i];

        if (e == 0) {
          basis_twist_ted[i][0] = ted_0; basis_twist_ted[i][1] = ted_0; basis_twist_ted[i][2] = ted_0;
          basis_twist[i][0].x = fp2_1; basis_twist[i][1].x = fp2_1; basis_twist[i][2].x = fp2_1;
          basis_twist[i][0].z = fp2_0; basis_twist[i][1].z = fp2_0; basis_twist[i][2].z = fp2_0;
          continue;
        }

        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_theta, ell, e, true);
        ted_add(&basis_twist_ted[i][2], &E_ted_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);
    }


    for (int i = 0; i < 3; ++i) {       
        copy_hybrid_point(pk->N1_basis+i, param->N1_basis+i);
        pk->N1_basis[i].E = E_theta;

        ted_add(&pt, &E_ted_curve, &basis_curve_ted[0][i], &basis_curve_ted[1][i]);
        for (int j = 2; j < on_curve_len; j++) {
            ted_add(&pt, &E_ted_curve, &pt, &basis_curve_ted[j][i]);
        }
        ted_to_mont_point(&pk->N1_basis[i].curve_pt, &pt);

        ted_add(&pt, &E_ted_twist, &basis_twist_ted[0][i], &basis_twist_ted[1][i]);
        for (int j = 2; j < on_twist_len; j++) {
            ted_add(&pt, &E_ted_twist, &pt, &basis_twist_ted[j][i]);
        }
        ted_to_mont_point(&pk->N1_basis[i].twist_pt, &pt);
    }









    

    for (int i = 0; i < on_curve_len; i++) {
        long ell = on_curve_fact[i], e = param->N2_basis[0].curve_mult[i];
        // printf("%ld^%ld\n",ell,e);

        if (e == 0) {
          basis_curve_ted[i][0] = ted_0; basis_curve_ted[i][1] = ted_0; basis_curve_ted[i][2] = ted_0;
          basis_curve[i][0].x = fp2_1; basis_curve[i][1].x = fp2_1; basis_curve[i][2].x = fp2_1;
          basis_curve[i][0].z = fp2_0; basis_curve[i][1].z = fp2_0; basis_curve[i][2].z = fp2_0;
          continue;
        }

        random_basis(&basis_curve[i][0], &basis_curve[i][1], &basis_curve_ted[i][0], &basis_curve_ted[i][1], &E_theta, ell, e, false);
        ted_add(&basis_curve_ted[i][2], &E_ted_curve, &basis_curve_ted[i][0], &basis_curve_ted[i][1]);
        ted_to_mont_point(&basis_curve[i][2], &basis_curve_ted[i][2]);
    }


    for (int i = 0; i < on_twist_len; i++) {
        long ell = on_twist_fact[i], e = param->N2_basis[0].twist_mult[i];

        if (e == 0) {
          basis_twist_ted[i][0] = ted_0; basis_twist_ted[i][1] = ted_0; basis_twist_ted[i][2] = ted_0;
          basis_twist[i][0].x = fp2_1; basis_twist[i][1].x = fp2_1; basis_twist[i][2].x = fp2_1;
          basis_twist[i][0].z = fp2_0; basis_twist[i][1].z = fp2_0; basis_twist[i][2].z = fp2_0;
          continue;
        }

        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_theta, ell, e, true);
        ted_add(&basis_twist_ted[i][2], &E_ted_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);
    }


    for (int i = 0; i < 3; ++i) {       
        copy_hybrid_point(pk->N2_basis+i, param->N2_basis+i);
        pk->N2_basis[i].E = E_theta;

        ted_add(&pt, &E_ted_curve, &basis_curve_ted[0][i], &basis_curve_ted[1][i]);
        for (int j = 2; j < on_curve_len; j++) {
            ted_add(&pt, &E_ted_curve, &pt, &basis_curve_ted[j][i]);
        }
        ted_to_mont_point(&pk->N2_basis[i].curve_pt, &pt);

        ted_add(&pt, &E_ted_twist, &basis_twist_ted[0][i], &basis_twist_ted[1][i]);
        for (int j = 2; j < on_twist_len; j++) {
            ted_add(&pt, &E_ted_twist, &pt, &basis_twist_ted[j][i]);
        }
        ted_to_mont_point(&pk->N2_basis[i].twist_pt, &pt);
    }







    // Randomly define some pseudo-action of theta on the torsion

    long N_LIMBS = 7;
    randombytes((sk->ker_psi_plus_1).c + 0, 8*N_LIMBS);
    uintbig_set(&sk->ker_psi_plus_2,1);
    randombytes((sk->ker_psi_minus_1).c + 0, 8*N_LIMBS);
    uintbig_set(&sk->ker_psi_minus_2,1);
    randombytes((sk->ker_psi_dual_plus_1).c + 0, 8*N_LIMBS);
    uintbig_set(&sk->ker_psi_dual_plus_2,1);
    randombytes((sk->ker_psi_dual_minus_1).c + 0, 8*N_LIMBS);
    uintbig_set(&sk->ker_psi_dual_minus_2,1);


}

void seta_eval(ciphertext *ctxt, const public_key *pk, const uintbig *x) {
    pari_sp ltop = avma;

    proj E_theta = pk->N1_basis[0].E;

    odd_isogeny phi_m;
    degree_one(&phi_m.deg_plus);
    degree_one(&phi_m.deg_minus);
    // phi_m.kernel_plus.x = fp2_1;
    // phi_m.kernel_plus.z = fp2_0;

    if(curve_order_is_p_plus_one) {
        for (long i = 0; i < p_minus_len; ++i) {
            degree_set(&phi_m.deg_minus, i, pk->N1_basis[0].twist_mult[i]);
        }

        for (long i = 0; i < p_plus_len; ++i) {
            degree_set(&phi_m.deg_plus, i, pk->N1_basis[0].curve_mult[i]);
        }

        xBIDIM(&phi_m.kernel_minus, &E_theta, &pk->N1_basis[0].twist_pt, &uintbig_1, &pk->N1_basis[1].twist_pt, x, &pk->N1_basis[2].twist_pt);
        xBIDIM(&phi_m.kernel_plus, &E_theta, &pk->N1_basis[0].curve_pt, &uintbig_1, &pk->N1_basis[1].curve_pt, x, &pk->N1_basis[2].curve_pt);
    } 
    else {
        for (long i = 0; i < p_minus_len; ++i) {
            degree_set(&phi_m.deg_minus, i, pk->N1_basis[0].curve_mult[i]);
        }

        for (long i = 0; i < p_plus_len; ++i) {
            degree_set(&phi_m.deg_plus, i, pk->N1_basis[0].twist_mult[i]);
        }

        xBIDIM(&phi_m.kernel_minus, &E_theta, &pk->N1_basis[0].curve_pt, &uintbig_1, &pk->N1_basis[1].curve_pt, x, &pk->N1_basis[2].curve_pt);
        xBIDIM(&phi_m.kernel_plus, &E_theta, &pk->N1_basis[0].twist_pt, &uintbig_1, &pk->N1_basis[1].twist_pt, x, &pk->N1_basis[2].twist_pt);
    }


    #ifndef NDEBUG
    isog_degree ord; 
    ord = mont_order(&phi_m.kernel_minus, &E_theta, curve_order_is_p_plus_one);
    assert(ord.val == phi_m.deg_minus.val);
    ord = mont_order(&phi_m.kernel_plus, &E_theta, !curve_order_is_p_plus_one);
    assert(ord.val == phi_m.deg_plus.val);
    #endif





    copy_hybrid_point(ctxt->N2_basis+0, &pk->N2_basis[0]);
    copy_hybrid_point(ctxt->N2_basis+1, &pk->N2_basis[1]);
    copy_hybrid_point(ctxt->N2_basis+2, &pk->N2_basis[2]);

    proj zip[6];
    zip[0] = ctxt->N2_basis[0].twist_pt;
    zip[1] = ctxt->N2_basis[1].twist_pt;
    zip[2] = ctxt->N2_basis[2].twist_pt;
    zip[3] = ctxt->N2_basis[0].curve_pt;
    zip[4] = ctxt->N2_basis[1].curve_pt;
    zip[5] = ctxt->N2_basis[2].curve_pt;

    proj E_m = E_theta;
    eval_mult(&E_m, &phi_m, zip, 6);

    ctxt->N2_basis[0].E = E_m;
    ctxt->N2_basis[1].E = E_m;
    ctxt->N2_basis[2].E = E_m;

    ctxt->N2_basis[0].twist_pt = zip[0];
    ctxt->N2_basis[1].twist_pt = zip[1];
    ctxt->N2_basis[2].twist_pt = zip[2];
    ctxt->N2_basis[0].curve_pt = zip[3];
    ctxt->N2_basis[1].curve_pt = zip[4];
    ctxt->N2_basis[2].curve_pt = zip[5];



    #ifndef NDEBUG
    isog_degree ord1, ord2, ord3, ord4;
    degree_one(&ord1);
    degree_one(&ord2);
    for (long i = 0; i < on_curve_len; ++i) {
      degree_set(&ord1, i, pk->N2_basis[0].curve_mult[i]);
    }
    for (long i = 0; i < on_twist_len; ++i) {
      degree_set(&ord2, i, pk->N2_basis[0].twist_mult[i]);
    }
    for (int j = 0; j < 3; ++j)
    {
        assert(is_on_curve(&ctxt->N2_basis[j].curve_pt, &E_m));
	    ord3 = mont_order(&ctxt->N2_basis[j].curve_pt, &E_m, false);
    	assert(ord1.val == ord3.val);
        assert(!is_on_curve(&ctxt->N2_basis[j].twist_pt, &E_m));
	    ord4 = mont_order(&ctxt->N2_basis[j].twist_pt, &E_m, true);
    	assert(ord2.val == ord4.val);
    }
    #endif

    avma = ltop;

}


// len = p_minus_len
// psi_basis_ted = uptosign_psi_basis_minus_ted
void seta_ker_psi_minus_d(uintbig *x, uintbig *y, long len, const long N1_fact[], const long N1_mult[],
    const point basis_ted[][3], const point psi_basis_ted[], const proj *E_ted, const public_param *param) {

    pari_sp ltop = avma;

    GEN N1 = gen_1;
    GEN N1_primary = zerovec(len);
    for (long i = 0; i < len; ++i) {
        gel(N1_primary, i+1) = powuu(N1_fact[i],N1_mult[i]);
        N1 = gmul(N1,gel(N1_primary, i+1));
    }


    GEN psi_matrix_N1 = zerovec(len);
    GEN sol_1 = zerovec(len);
    GEN sol_2 = zerovec(len);


    // Compute the action of psi on the N1-torsion and find the N1-kernel of psi - d, i.e., of psi_m

    for (int i = 0; i < len; i++) {
        long ell = N1_fact[i], e = N1_mult[i];
        if (e == 0) {
          gel(psi_matrix_N1, i+1) =  mkmat2(mkcol2(gen_1,gen_0),mkcol2(gen_0,gen_1));
          continue;
        }


        GEN gelle = powuu(ell,e);

        point P;

        uintbig z;

        gentobig(&z,diviiexact(N1,gelle));
        //assert(uintbig_equal(&z,&z_bis));


        GEN A,B,C,D;
        ted_mul(&P, psi_basis_ted + 0, E_ted, &z);
        long found;
        found = ted_bidim_log(&A, &C, E_ted, &P, basis_ted[i] + 0, basis_ted[i]+1, ell, e);
        assert(found);

        ted_mul(&P, psi_basis_ted + 1, E_ted, &z);
        found = ted_bidim_log(&B, &D, E_ted, &P, basis_ted[i] + 0, basis_ted[i]+1, ell, e);
        assert(found);

        gel(psi_matrix_N1, i+1) = gmod(gmul(mkmat2(mkcol2(A,C),mkcol2(B,D)),  ginvmod(diviiexact(N1,gelle), gelle) ),gelle);
        //output(psi_matrix_N1);

        GEN t = gtrace(gel(psi_matrix_N1, i+1));
        if (!isexactzero(gmod(gsub(t,gadd(param->d,param->d)), gelle))) {
            //printf("FLIP!\n");
            gel(psi_matrix_N1, i+1) = gneg(gel(psi_matrix_N1, i+1));
        }
        t = gtrace(gel(psi_matrix_N1, i+1));
        assert(isexactzero(gmod(gsub(t,gadd(param->d,param->d)), gelle)));

        GEN phi_eta_matrix = gsub(gel(psi_matrix_N1, i+1), mkmat2(mkcol2(param->d,gen_0),mkcol2(gen_0,param->d)));

        GEN ker = matkermod(phi_eta_matrix, gelle, NULL);


        GEN sol, sol_reduced;



        int dim_ker = lg(ker)-1;
        printf("dim_ker %d\n", dim_ker);
        assert(dim_ker > 0);
        for (int i = 1; i <= dim_ker; ++i){
          sol = gel(ker,i);
          sol_reduced = gmodgs(sol,ell);
          if (!isexactzero(sol_reduced)) break;
          assert(i < dim_ker);
        }

        gel(sol_1,i+1) = gel(sol, 1);
        gel(sol_2,i+1) = gel(sol, 2);


    }

    gentobig(x, ZV_chinese(sol_1, N1_primary, NULL));
    gentobig(y, ZV_chinese(sol_2, N1_primary, NULL));

    avma = ltop;
}




void seta_inverse(odd_isogeny *phi_m, const ciphertext *ctxt, const public_key *pk, const secret_key *sk, const public_param *param) {

    pari_sp ltop = avma;

    clock_t t;
    proj E_m = ctxt->N2_basis[0].E;
    // proj E_theta = pk->N1_basis[0].E;




    //Reconstitute psi from the secret key and the N2-basis of the ciphertext

    odd_isogeny psi1 = trivial_odd_isogeny(), psi2 = trivial_odd_isogeny();

    if(curve_order_is_p_plus_one) {
        for (long i = 0; i < p_plus_len; ++i) {
          degree_set(&psi1.deg_plus, i, ctxt->N2_basis[0].curve_mult[i]);
          degree_set(&psi2.deg_plus, i, ctxt->N2_basis[0].curve_mult[i]);
        }
        for (long i = 0; i < p_minus_len; ++i) {
          degree_set(&psi1.deg_minus, i, ctxt->N2_basis[0].twist_mult[i]);
          degree_set(&psi2.deg_minus, i, ctxt->N2_basis[0].twist_mult[i]);
        }
    
        xBIDIM(&psi1.kernel_plus, &E_m, &ctxt->N2_basis[0].curve_pt, &sk->ker_psi_plus_1, &ctxt->N2_basis[1].curve_pt, &sk->ker_psi_plus_2, &ctxt->N2_basis[2].curve_pt);
        xBIDIM(&psi1.kernel_minus, &E_m, &ctxt->N2_basis[0].twist_pt, &sk->ker_psi_minus_1, &ctxt->N2_basis[1].twist_pt, &sk->ker_psi_minus_2, &ctxt->N2_basis[2].twist_pt);
    
        xBIDIM(&psi2.kernel_plus, &E_m, &ctxt->N2_basis[0].curve_pt, &sk->ker_psi_dual_plus_1, &ctxt->N2_basis[1].curve_pt, &sk->ker_psi_dual_plus_2, &ctxt->N2_basis[2].curve_pt);
        xBIDIM(&psi2.kernel_minus, &E_m, &ctxt->N2_basis[0].twist_pt, &sk->ker_psi_dual_minus_1, &ctxt->N2_basis[1].twist_pt, &sk->ker_psi_dual_minus_2, &ctxt->N2_basis[2].twist_pt);
    }
    else {
        for (long i = 0; i < p_plus_len; ++i) {
          degree_set(&psi1.deg_plus, i, ctxt->N2_basis[0].twist_mult[i]);
          degree_set(&psi2.deg_plus, i, ctxt->N2_basis[0].twist_mult[i]);
        }
        for (long i = 0; i < p_minus_len; ++i) {
          degree_set(&psi1.deg_minus, i, ctxt->N2_basis[0].curve_mult[i]);
          degree_set(&psi2.deg_minus, i, ctxt->N2_basis[0].curve_mult[i]);
        }
    
        xBIDIM(&psi1.kernel_plus, &E_m, &ctxt->N2_basis[0].twist_pt, &sk->ker_psi_plus_1, &ctxt->N2_basis[1].twist_pt, &sk->ker_psi_plus_2, &ctxt->N2_basis[2].twist_pt);
        xBIDIM(&psi1.kernel_minus, &E_m, &ctxt->N2_basis[0].curve_pt, &sk->ker_psi_minus_1, &ctxt->N2_basis[1].curve_pt, &sk->ker_psi_minus_2, &ctxt->N2_basis[2].curve_pt);
    
        xBIDIM(&psi2.kernel_plus, &E_m, &ctxt->N2_basis[0].twist_pt, &sk->ker_psi_dual_plus_1, &ctxt->N2_basis[1].twist_pt, &sk->ker_psi_dual_plus_2, &ctxt->N2_basis[2].twist_pt);
        xBIDIM(&psi2.kernel_minus, &E_m, &ctxt->N2_basis[0].curve_pt, &sk->ker_psi_dual_minus_1, &ctxt->N2_basis[1].curve_pt, &sk->ker_psi_dual_minus_2, &ctxt->N2_basis[2].curve_pt);
    }

    #ifndef NDEBUG
    isog_degree ord;
    ord = mont_order(&psi1.kernel_plus, &E_m, !curve_order_is_p_plus_one);
    assert(ord.val == psi1.deg_plus.val);
    ord = mont_order(&psi1.kernel_minus, &E_m, curve_order_is_p_plus_one);
    assert(ord.val == psi1.deg_minus.val);

    ord = mont_order(&psi2.kernel_plus, &E_m, !curve_order_is_p_plus_one);
    assert(ord.val == psi2.deg_plus.val);
    ord = mont_order(&psi2.kernel_minus, &E_m, curve_order_is_p_plus_one);
    assert(ord.val == psi2.deg_minus.val);
    #endif

    special_isogeny psi = trivial_special_isogeny();


    psi.source = E_m;
    psi.target = E_m;

    psi.phi2_set = false;
    psi.phi2_dual_set = true;

    psi.phi1 = psi1;
    psi.phi2_dual = psi2;

    // #ifndef NDEBUG
    // proj D = psi.source;
    // eval_special(&D, &psi, &ctxt->N2_basis[0].curve_pt);
    // #endif




    // generate N1-basis

    proj basis_twist[on_twist_len][3], basis_curve[on_curve_len][3];
    point ted_0 = {fp2_0, fp2_1, fp2_1, fp2_0};
    point basis_twist_ted[on_twist_len][3], basis_curve_ted[on_curve_len][3];


    proj E_ted_twist;
    proj E_ted_curve;
    mont_to_ted(&E_ted_twist, &E_m, true);
    mont_to_ted(&E_ted_curve, &E_m, false);
    t = tic();

    // twist part
    for (int i = 0; i < on_twist_len; i++) {
        long ell = on_twist_fact[i], e = pk->N1_basis[0].twist_mult[i];

        if (e == 0) {
          basis_twist_ted[i][0] = ted_0; basis_twist_ted[i][1] = ted_0; basis_twist_ted[i][2] = ted_0;
          basis_twist[i][0].x = fp2_1; basis_twist[i][1].x = fp2_1; basis_twist[i][2].x = fp2_1;
          basis_twist[i][0].z = fp2_0; basis_twist[i][1].z = fp2_0; basis_twist[i][2].z = fp2_0;
          continue;
        }


        // GEN gelle = powuu(ell,e);
        // GEN inv2 = Fp_inv(gen_2, gelle);

        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_m, ell, e, true);

        ted_add(&basis_twist_ted[i][2], &E_ted_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);


        #ifndef NDEBUG
        point test;
        proj test2;
        uintbig ell_big;
        gentobig(&ell_big, powuu(ell,e));
        ted_mul(&test, &basis_twist_ted[i][0], &E_ted_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_mul(&test, &basis_twist_ted[i][1], &E_ted_twist, &ell_big);
        assert(ted_iszero(&test));

        ted_to_mont_point(&test2, &basis_twist_ted[i][0]);
        assert(!is_on_curve(&test2, &E_m));
        xMUL(&test2, &E_m, &test2, &ell_big);
        assert(mont_iszero(&test2));

        ted_to_mont_point(&test2, &basis_twist_ted[i][1]);
        assert(!is_on_curve(&test2, &E_m));
        xMUL(&test2, &E_m, &test2, &ell_big);
        assert(mont_iszero(&test2));


        ted_mul(&test, &basis_twist_ted[i][2], &E_ted_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_to_mont_point(&test2, &basis_twist_ted[i][2]);
        assert(!is_on_curve(&test2, &E_m));
        xMUL(&test2, &E_m, &test2, &ell_big);
        assert(mont_iszero(&test2));
        #endif



    }

    // plus part
    for (int i = 0; i < on_curve_len; i++) {
        long ell = on_curve_fact[i], e = pk->N1_basis[0].curve_mult[i];

        if (e == 0) {
          basis_curve_ted[i][0] = ted_0; basis_curve_ted[i][1] = ted_0; basis_curve_ted[i][2] = ted_0;
          basis_curve[i][0].x = fp2_1; basis_curve[i][1].x = fp2_1; basis_curve[i][2].x = fp2_1;
          basis_curve[i][0].z = fp2_0; basis_curve[i][1].z = fp2_0; basis_curve[i][2].z = fp2_0;
          continue;
        }

        random_basis(&basis_curve[i][0], &basis_curve[i][1], &basis_curve_ted[i][0], &basis_curve_ted[i][1], &E_m, ell, e, false);

        ted_add(&basis_curve_ted[i][2], &E_ted_curve, &basis_curve_ted[i][0], &basis_curve_ted[i][1]);
        ted_to_mont_point(&basis_curve[i][2], &basis_curve_ted[i][2]);
    }



    TOC(t,"basis");


    proj basis_twist_sum[3], basis_curve_sum[3];
    point basis_twist_ted_sum[3], basis_curve_ted_sum[3];


    ted_add(&basis_twist_ted_sum[0], &E_ted_twist, &basis_twist_ted[0][0], &basis_twist_ted[1][0]);
    ted_add(&basis_twist_ted_sum[1], &E_ted_twist, &basis_twist_ted[0][1], &basis_twist_ted[1][1]);
    for (int i = 2; i < on_twist_len; i++) {
        ted_add(&basis_twist_ted_sum[0], &E_ted_twist, &basis_twist_ted_sum[0], &basis_twist_ted[i][0]);
        ted_add(&basis_twist_ted_sum[1], &E_ted_twist, &basis_twist_ted_sum[1], &basis_twist_ted[i][1]);
    }
    ted_add(&basis_twist_ted_sum[2], &E_ted_twist, &basis_twist_ted_sum[0], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[0], &basis_twist_ted_sum[0]);
    ted_to_mont_point(&basis_twist_sum[1], &basis_twist_ted_sum[1]);
    
    ted_to_mont_point(&basis_twist_sum[2], &basis_twist_ted_sum[2]);


    ted_add(&basis_curve_ted_sum[0], &E_ted_curve, &basis_curve_ted[0][0], &basis_curve_ted[1][0]);
    ted_add(&basis_curve_ted_sum[1], &E_ted_curve, &basis_curve_ted[0][1], &basis_curve_ted[1][1]);
    for (int i = 2; i < on_curve_len; i++) {
        ted_add(&basis_curve_ted_sum[0], &E_ted_curve, &basis_curve_ted_sum[0], &basis_curve_ted[i][0]);
        ted_add(&basis_curve_ted_sum[1], &E_ted_curve, &basis_curve_ted_sum[1], &basis_curve_ted[i][1]);
    }
    ted_add(&basis_curve_ted_sum[2], &E_ted_curve, &basis_curve_ted_sum[0], &basis_curve_ted_sum[1]);
    ted_to_mont_point(&basis_curve_sum[0], &basis_curve_ted_sum[0]);
    ted_to_mont_point(&basis_curve_sum[1], &basis_curve_ted_sum[1]);
    
    ted_to_mont_point(&basis_curve_sum[2], &basis_curve_ted_sum[2]);




    // evaluate psi on the N1-basis, and find the images in the twisted Edwards model, up to sign

    proj psi_basis_twist[3];
    proj psi_basis_curve[3];
    proj C = psi.source;
    t = tic();
    // TODO: do the following three at once

    proj zip[6];
    zip[0] = basis_twist_sum[0];
    zip[1] = basis_twist_sum[1];
    zip[2] = basis_twist_sum[2];
    zip[3] = basis_curve_sum[0];
    zip[4] = basis_curve_sum[1];
    zip[5] = basis_curve_sum[2];

    // psi_basis_twist[0] = eval_special(&C, &psi, basis_twist_sum + 0);
    // psi_basis_twist[1] = eval_special(&C, &psi, basis_twist_sum + 1);
    // psi_basis_twist[2] = eval_special(&C, &psi, basis_twist_sum + 2);
    // psi_basis_curve[0] = eval_special(&C, &psi, basis_curve_sum + 0);
    // psi_basis_curve[1] = eval_special(&C, &psi, basis_curve_sum + 1);
    // psi_basis_curve[2] = eval_special(&C, &psi, basis_curve_sum + 2);

    eval_special_mult(&C, &psi, zip, 6);

    psi_basis_twist[0] = zip[0];
    psi_basis_twist[1] = zip[1];
    psi_basis_twist[2] = zip[2];
    psi_basis_curve[0] = zip[3];
    psi_basis_curve[1] = zip[4];
    psi_basis_curve[2] = zip[5];

    TOC(t,"eval_special");

    proj j1,j2;
    
    #ifndef NDEBUG
    jinv256(&j1, &C);
    jinv256(&j2, &psi.source);
    // if (!mont_equal(&j1,&j2)) fprintf(stderr,"ERROR in seta_inverse: j-invariants of psi.source and the target of psi do not match\n");
    assert(mont_equal(&j1,&j2));
    #endif


    // minus part
    point uptosign_psi_basis_twist_ted[3];
    mont_to_ted_point(&uptosign_psi_basis_twist_ted[0], &E_m, &psi_basis_twist[0]);
    mont_to_ted_point(&uptosign_psi_basis_twist_ted[1], &E_m, &psi_basis_twist[1]);
    //mont_to_ted_point(&uptosign_psi_basis_twist_ted[2], &E_m, &psi_basis_twist[2]);

    ted_add(&uptosign_psi_basis_twist_ted[2], &E_ted_twist, &uptosign_psi_basis_twist_ted[0], &uptosign_psi_basis_twist_ted[1]);
    proj test_pt;
    ted_to_mont_point(&test_pt, &uptosign_psi_basis_twist_ted[2]);
    if (!mont_equal(&test_pt, psi_basis_twist+2)) {
      ted_neg(&uptosign_psi_basis_twist_ted[1], &uptosign_psi_basis_twist_ted[1]);
      ted_add(&uptosign_psi_basis_twist_ted[2], &E_ted_twist, &uptosign_psi_basis_twist_ted[0], &uptosign_psi_basis_twist_ted[1]);
    }
    ted_to_mont_point(&test_pt, &uptosign_psi_basis_twist_ted[2]);
    assert(mont_equal(&test_pt, psi_basis_twist+2));

    // plus part
    point uptosign_psi_basis_curve_ted[3];
    mont_to_ted_point(&uptosign_psi_basis_curve_ted[0], &E_m, &psi_basis_curve[0]);
    mont_to_ted_point(&uptosign_psi_basis_curve_ted[1], &E_m, &psi_basis_curve[1]);
    //mont_to_ted_point(&uptosign_psi_basis_curve_ted[2], &E_m, &psi_basis_curve[2]);

    ted_add(&uptosign_psi_basis_curve_ted[2], &E_ted_curve, &uptosign_psi_basis_curve_ted[0], &uptosign_psi_basis_curve_ted[1]);

    ted_to_mont_point(&test_pt, &uptosign_psi_basis_curve_ted[2]);
    if (!mont_equal(&test_pt, psi_basis_curve+2)) {
      ted_neg(&uptosign_psi_basis_curve_ted[1], &uptosign_psi_basis_curve_ted[1]);
      ted_add(&uptosign_psi_basis_curve_ted[2], &E_ted_curve, &uptosign_psi_basis_curve_ted[0], &uptosign_psi_basis_curve_ted[1]);
    }
    ted_to_mont_point(&test_pt, &uptosign_psi_basis_curve_ted[2]);
    assert(mont_equal(&test_pt, psi_basis_curve+2));








    odd_isogeny phi_m_dual;
    degree_one(&phi_m_dual.deg_plus);
    degree_one(&phi_m_dual.deg_minus);

    if (curve_order_is_p_plus_one) {
        for (long i = 0; i < p_minus_len; ++i) {
          degree_set(&phi_m_dual.deg_minus, i, pk->N1_basis[0].twist_mult[i]);
        }
        for (long i = 0; i < p_plus_len; ++i) {
          degree_set(&phi_m_dual.deg_plus, i, pk->N1_basis[0].curve_mult[i]);
        }
    }
    else {
        for (long i = 0; i < p_minus_len; ++i) {
          degree_set(&phi_m_dual.deg_minus, i, pk->N1_basis[0].curve_mult[i]);
        }
        for (long i = 0; i < p_plus_len; ++i) {
          degree_set(&phi_m_dual.deg_plus, i, pk->N1_basis[0].twist_mult[i]);
        }
    }


    uintbig x, y;
    proj phi_m_source = E_m;

    for (int i = 0; i < 4; ++i) {
        seta_ker_psi_minus_d(&x, &y, on_twist_len, on_twist_fact, pk->N1_basis[0].twist_mult, basis_twist_ted, uptosign_psi_basis_twist_ted, &E_ted_twist, param);
        xBIDIM(&phi_m_dual.kernel_minus, &E_m, &basis_twist_sum[0], &x, &basis_twist_sum[1], &y, &basis_twist_sum[2]);

        seta_ker_psi_minus_d(&x, &y, on_curve_len, on_curve_fact, pk->N1_basis[0].curve_mult, basis_curve_ted, uptosign_psi_basis_curve_ted, &E_ted_curve, param);
        xBIDIM(&phi_m_dual.kernel_plus, &E_m, &basis_curve_sum[0], &x, &basis_curve_sum[1], &y, &basis_curve_sum[2]);

        if(!curve_order_is_p_plus_one){ // swap
            proj tmp = phi_m_dual.kernel_plus;
            phi_m_dual.kernel_plus = phi_m_dual.kernel_minus;
            phi_m_dual.kernel_minus = tmp;
        }

        // proj zip[1];
        // zip[0] = ctxt->N2_basis[0].curve_pt;
        *phi_m = phi_m_dual;
        dual(&phi_m_source, phi_m);

        jinv256(&j1, &phi_m_source);
        jinv256(&j2, &pk->N1_basis[0].E);
        

        if (!mont_equal(&j1,&j2)) {
            printf("%dThis should probably not happen\n", i);
            for (int j = 0; j < 3; ++j)
            {
                if (i % 2 == 0){
                    ted_neg(&uptosign_psi_basis_twist_ted[j], &uptosign_psi_basis_twist_ted[j]);
                }
                else {
                    ted_neg(&uptosign_psi_basis_curve_ted[j], &uptosign_psi_basis_curve_ted[j]);
                }
            }
        }
        else { break; }

    }

    assert(mont_equal(&j1,&j2));

    avma = ltop;
}









void seta_setup(public_param *param) {

    pari_sp ltop = avma;

    assert(C_LEN == on_curve_len);
    assert(T_LEN == on_twist_len);
    // minus part, i.e., twist
    // long N1_fact[19] =
    //   {  3, 43, 103, 109, 199, 227, 419, 491, 569, 631, 677, 857, 859, 883,
    //     1019, 1171, 1879, 2713, 4283 };

    // plus part, i.e., non-twist
    // long N2_fact[12] =
    //   {  5, 7, 11, 31, 83, 107, 137, 751, 827, 3691, 4019, 6983 };
        
    // long N1_mult[19] =
    //     { 0,  0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    //        0,    0,    0,    0,    0 };
    // long N2_mult[12] =
    //     { 21, 1,  1,  1,  0,   1,   0,   1,   1,    0,    1,    0 };
       


    // long N1_mult[19] =
    //     { 21,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    //        0,    0,    0,    0,    0 };
    // long N2_mult[12] =
    //     { 21, 2,  1,  1,  1,   1,   1,   1,   1,    0,    0,    0 };

        
    // long N1_mult[19] =
    //     { 16,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    //        0,    0,    0,    0,    0 };
    // long N2_mult[12] =
    //     { 21, 2,  1,  1,  1,   0,   1,   1,   1,    0,    1,    0 };


    proj E0 = global_setup.E0;

    hybrid_point *N1_basis = param->N1_basis;
    hybrid_point *N2_basis = param->N2_basis;


    N1_basis[0].E = E0;
    N1_basis[1].E = E0;
    N1_basis[2].E = E0;
    N2_basis[0].E = E0;
    N2_basis[1].E = E0;
    N2_basis[2].E = E0;

    GEN coeff_1_N1 = zerovec(on_twist_len), 
        coeff_1_N2 = zerovec(on_twist_len), 
        coeff_0 = zerovec(on_twist_len);
    for (long i = 0; i < on_twist_len; ++i) {
      gel(coeff_1_N1, i+1) = powuu(on_twist_fact[i],on_twist_mult[i] - N1_mult_twist[i]);
      gel(coeff_1_N2, i+1) = powuu(on_twist_fact[i],on_twist_mult[i] - N2_mult_twist[i]);
    }
    N1_basis[0].twist_pt = coeff_to_E0(mkvec2(coeff_1_N1,coeff_0), true);
    N1_basis[1].twist_pt = coeff_to_E0(mkvec2(coeff_0,coeff_1_N1), true);
    N1_basis[2].twist_pt = coeff_to_E0(mkvec2(coeff_1_N1,coeff_1_N1), true);

    N2_basis[0].twist_pt = coeff_to_E0(mkvec2(coeff_1_N2,coeff_0), true);
    N2_basis[1].twist_pt = coeff_to_E0(mkvec2(coeff_0,coeff_1_N2), true);
    N2_basis[2].twist_pt = coeff_to_E0(mkvec2(coeff_1_N2,coeff_1_N2), true);




    coeff_1_N1 = zerovec(on_curve_len);
    coeff_1_N2 = zerovec(on_curve_len);
    coeff_0 = zerovec(on_curve_len);
    for (long i = 0; i < on_curve_len; ++i) {
      gel(coeff_1_N1, i+1) = powuu(on_curve_fact[i],on_curve_mult[i] - N1_mult_curve[i]);
      gel(coeff_1_N2, i+1) = powuu(on_curve_fact[i],on_curve_mult[i] - N2_mult_curve[i]);
    }

    N1_basis[0].curve_pt = coeff_to_E0(mkvec2(coeff_1_N1,coeff_0), false);
    N1_basis[1].curve_pt = coeff_to_E0(mkvec2(coeff_0,coeff_1_N1), false);
    N1_basis[2].curve_pt = coeff_to_E0(mkvec2(coeff_1_N1,coeff_1_N1), false);

    N2_basis[0].curve_pt = coeff_to_E0(mkvec2(coeff_1_N2,coeff_0), false);
    N2_basis[1].curve_pt = coeff_to_E0(mkvec2(coeff_0,coeff_1_N2), false);
    N2_basis[2].curve_pt = coeff_to_E0(mkvec2(coeff_1_N2,coeff_1_N2), false);


    for (long i = 0; i < on_twist_len; ++i) {
      N1_basis[0].twist_mult[i] = N1_mult_twist[i];
      N1_basis[1].twist_mult[i] = N1_mult_twist[i];
      N1_basis[2].twist_mult[i] = N1_mult_twist[i];
      N2_basis[0].twist_mult[i] = N2_mult_twist[i];
      N2_basis[1].twist_mult[i] = N2_mult_twist[i];
      N2_basis[2].twist_mult[i] = N2_mult_twist[i];
    }
    for (long i = 0; i < on_curve_len; ++i) {
      N1_basis[0].curve_mult[i] = N1_mult_curve[i];
      N1_basis[1].curve_mult[i] = N1_mult_curve[i];
      N1_basis[2].curve_mult[i] = N1_mult_curve[i];
      N2_basis[0].curve_mult[i] = N2_mult_curve[i];
      N2_basis[1].curve_mult[i] = N2_mult_curve[i];
      N2_basis[2].curve_mult[i] = N2_mult_curve[i];
    }


    #ifndef NDEBUG

    isog_degree ord1,ord2; 
    
    degree_one(&ord1);
    for (long i = 0; i < on_twist_len; ++i) {
        degree_set(&ord1, i, N1_basis[0].twist_mult[i]);
    }
    ord2 = mont_order(&N1_basis[0].twist_pt, &E0, true);

    assert(ord1.val == ord2.val);

    degree_one(&ord1);
    for (long i = 0; i < on_curve_len; ++i) {
        degree_set(&ord1, i, N1_basis[0].curve_mult[i]);
    }
    ord2 = mont_order(&N1_basis[0].curve_pt, &E0, false);
    assert(ord1.val == ord2.val);
    #endif


    clock_t t;


    param->N1_twist = gen_1;
    param->N1_curve = gen_1;
    param->N1_primary_twist = zerovec(on_twist_len);
    param->N1_primary_curve = zerovec(on_curve_len);
    param->N2_twist = gen_1;
    param->N2_curve = gen_1;
    param->N2_primary_twist = zerovec(on_twist_len);
    param->N2_primary_curve = zerovec(on_curve_len);

    // GEN N1_minus_cof = gen_1, N1_plus_cof = gen_1;
    // GEN N2_minus_cof = gen_1, N2_plus_cof = gen_1;

    for (long i = 0; i < on_twist_len; ++i) {
        gel(param->N1_primary_twist, i+1) = powuu(on_twist_fact[i],N1_basis[0].twist_mult[i]);
        param->N1_twist = gmul(param->N1_twist,gel(param->N1_primary_twist, i+1));
        gel(param->N2_primary_twist, i+1) = powuu(on_twist_fact[i],N2_basis[0].twist_mult[i]);
        param->N2_twist = gmul(param->N2_twist,gel(param->N2_primary_twist, i+1));
    }
    for (long i = 0; i < on_curve_len; ++i) {
        gel(param->N1_primary_curve, i+1) = powuu(on_curve_fact[i],N1_basis[0].curve_mult[i]);
        param->N1_curve= gmul(param->N1_curve,gel(param->N1_primary_curve, i+1));
        gel(param->N2_primary_curve, i+1) = powuu(on_curve_fact[i],N2_basis[0].curve_mult[i]);
        param->N2_curve = gmul(param->N2_curve,gel(param->N2_primary_curve, i+1));
    }




    GEN N1 = gmul(param->N1_twist, param->N1_curve), N2 = gmul(param->N2_twist, param->N2_curve);

    param->d = gmod(N2, gmul(N1,N1));
    param->D = diviiexact(gsub(gmul(N2,N2), gmul(param->d,param->d)),gmul(N1,N1));
    // output(param->d);
    // output(param->D);

    assert(Fp_issquare(gmul(param->D,stoi(global_setup.q)), global_setup.p));

    // while ((gcmp(ggcd(param->d,N2),gen_1) != 0) || !Fp_issquare(gmul(param->D,stoi(global_setup.q)), global_setup.p)) {
    //     // printf("gcd(d,N2) = "); output(ggcd(param->d,N2));
    //     param->d = gadd(param->d, gmul(N1,N1));
    //     param->D = diviiexact(gsub(gmul(N2,N2), gmul(param->d,param->d)),gmul(N1,N1));
    //     assert(gcmp(param->d,N2) < 0); // didn't find d such that D is a square mod p and (N2,d) = 1
    // }

    // assert(Fp_issquare(gmul(param->D,stoi(global_setup.q)), global_setup.p));

    for (int i = 0; i < on_twist_len; ++i) {
        if (N1_basis[0].twist_mult[i] > 0) {
            long ell = on_twist_fact[i];
            assert(!Fp_issquare(gneg(param->D), stoi(ell)));
            assert(!isexactzero(gmod(gneg(param->D), stoi(ell))));
        }
    }
    for (int i = 0; i < on_curve_len; ++i) {
        if (N1_basis[0].curve_mult[i] > 0) {
            long ell = on_curve_fact[i];
            assert(!Fp_issquare(gneg(param->D), stoi(ell)));
            assert(!isexactzero(gmod(gneg(param->D), stoi(ell))));
        }
    }


    t = tic();
    // output(param->D);


    // GEN D1 = gsub(N2, param->d);
    // GEN D2 = gadd(N2, param->d);
    // D1 = diviiexact(param->D, ggcd(param->D,D1));
    // D2 = diviiexact(param->D, D1);
    // output(D1);
    // output(D2);


    printf("\\\\ FACTORISATION OF D\n"); 
    param->fac_det = absZ_factor(param->D);
    // printf("WARNING: param->fac_det is hardcoded for testing purpuses\n");
    // param->fac_det = mkmat2(mkcoln(8,strtoi("1571"),strtoi("1281383"),strtoi("6827983"),strtoi("8400253"),strtoi("85316653726537"),strtoi("549074433904903"),strtoi("844031771502222229118375835277"),strtoi("74915729981907191583653570582561")), 
    //     mkcoln(8,gen_1,gen_1,gen_1,gen_1,gen_1,gen_1,gen_1,gen_1));
    

    output(param->fac_det);
    printf("\n\n"); 



    param->fac_det = famat_mul(param->fac_det, famat_mul(to_famat(global_setup.p, gen_2),to_famat(stoi(global_setup.q), gen_2))); 
    TOC(t,"factor");

    gerepileall(ltop, 11, &param->N1_twist, &param->N1_curve, &param->N1_primary_twist, &param->N1_primary_curve,
                &param->N2_twist, &param->N2_curve, &param->N2_primary_twist, &param->N2_primary_curve,
                &param->D, &param->d, &param->fac_det);







}






