#define _XOPEN_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <pari/pari.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>


#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"

#include "curve.h"
#include "mont.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"

/*
static char* fp2_hash(fp2 x) {
    return pari_sprintf("h%lu", (x.re.x.c[0]+3*x.re.x.c[1]+5*x.re.x.c[2]+7*x.re.x.c[3]
                    +11*x.im.x.c[0]+13*x.im.x.c[1]+17*x.im.x.c[2]+23*x.im.x.c[3]) % 100003);
}
static fp2 fp2_ratio(fp2 *x, fp2 *y) {
    fp2 tmp;
    tmp = *y;
    assert(!fp2_iszero(&tmp));
    fp2_inv(&tmp);
    fp2_mul2(&tmp, x);
    //printf("%lld %lld %lld %lld\n", tmp.re.x.c[0], tmp.re.x.c[1], tmp.im.x.c[0], tmp.im.x.c[1]);
    return tmp;
}
static char* proj_hash(proj x) {
    return fp2_hash(fp2_ratio(&x.x,&x.z));
}



static GEN norm0(GEN x) {
    return algnorm(global_setup.B, x,0);
}
*/



static void mont0_id_plus_dist(proj *Q, const proj *P, const proj *A) {
    proj Pcopy = *P;
    fp2 x,x2,x3;
    x = Pcopy.x;
    fp2_sq2(&x2, &Pcopy.x);
    fp2_mul3(&x3, &x2, &Pcopy.x);

    fp2_mul3(&x, &x, &Pcopy.z);
    fp2_mul3(&x, &x, &Pcopy.z); // x*z^2
    fp2_mul3(&x2, &x2, &Pcopy.z); // x^2*z

    fp2_mul3(&Q->z, &x2, &A->z);
    fp2_add2(&Q->z, &Q->z); // multiplication by 2

    fp2_mul2(&x, &A->z);
    fp2_mul2(&x2, &A->x);
    fp2_mul2(&x3, &A->z);
    fp2_add3(&Q->x, &x, &x2);
    fp2_add2(&Q->x, &x3);
    fp2_mul2(&Q->x, &fp2_i);
    fp2_neg1(&Q->x);
}




static point ted_action(GEN coeff, const point *P, const proj *E, bool twist) {
    point A, B, C, D;

    A = B = C = D = *P;

    ted0_dist(&B,&B);

    if (!twist) ted0_frob(&C,&C);
    else ted0_frob_twist(&C,&C);

    D = B;
    if (!twist) ted0_frob(&D,&D);
    else ted0_frob_twist(&D,&D);

    uintbig x;

    gentobig(&x, gel(coeff,1));
    ted_mul(&A, &A, E, &x);
    gentobig(&x, gel(coeff,2));
    ted_mul(&B, &B, E, &x);
    gentobig(&x, gel(coeff,3));
    ted_mul(&C, &C, E, &x);
    gentobig(&x, gel(coeff,4));
    ted_mul(&D, &D, E, &x);

    point Q = A;
    ted_add(&Q, E, &Q, &B);
    ted_add(&Q, E, &Q, &C);
    ted_add(&Q, E, &Q, &D);

    return Q;
}















static bool is_in_kernel_of_endo_two(GEN endo, GEN v, proj pt) {

    proj P0 = torsion_basis_two[0], P1 = torsion_basis_two[1], P2 = torsion_basis_two[2];
    uintbig x,y;

    proj K;

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&K, &global_setup.E0, &P0, &x, &P1, &y, &P2);

    if (!(mont_equal(&K,&pt))) return false;


    proj iP0, iP1, iP2;



    mont0_dist(&iP0, &P0);
    mont0_dist(&iP1, &P1);
    mont0_dist(&iP2, &P2);


    proj P0_iP0, P1_iP1, P2_iP2;

    mont0_id_plus_dist(&P0_iP0, &P0, &global_setup.E0);
    mont0_id_plus_dist(&P1_iP1, &P1, &global_setup.E0);
    mont0_id_plus_dist(&P2_iP2, &P2, &global_setup.E0);


    proj a_ib_0, a_ib_1, a_ib_2, a_ib_K;

    gentobig(&x, gel(endo,1));
    gentobig(&y, gel(endo,2));
    xBIDIM(&a_ib_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&a_ib_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&a_ib_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&a_ib_K, &global_setup.E0, &a_ib_0, &x, &a_ib_1, &y, &a_ib_2);


    proj c_id_0, c_id_1, c_id_2, c_id_K;

    gentobig(&x, gel(endo,3));
    gentobig(&y, gel(endo,4));
    xBIDIM(&c_id_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&c_id_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&c_id_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&c_id_K, &global_setup.E0, &c_id_0, &x, &c_id_1, &y, &c_id_2);

    proj jc_jid_K;
    mont0_frob(&jc_jid_K, &c_id_K);


    return mont_equal(&jc_jid_K, &a_ib_K);


}




static bool is_in_kernel_of_endo_ell(GEN endo, GEN v, proj *pt, long ell) {
    long e = ell_to_e(ell);
    bool twist;
    //unsigned long index = ell_to_index(ell, &twist);
    ell_to_index(ell, &twist);
    uintbig x,y;

    GEN gelle = powuu(ell,e);

    GEN cof = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);
    cof = gdiv(cof, gelle);

    uintbig cof_big;

    gentobig(&cof_big,cof);

    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;

    proj P0 = basis[0], P1 = basis[1], P2 = basis[2];
    xMUL(&P0, &global_setup.E0, &P0, &cof_big);
    xMUL(&P1, &global_setup.E0, &P1, &cof_big);
    xMUL(&P2, &global_setup.E0, &P2, &cof_big);


    proj ker_alt;
    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&ker_alt, &global_setup.E0, &P0, &x, &P1, &y, &P2);

    assert(mont_equal(&ker_alt, pt));


    // GEN m_id, m_i, m_j, m_ji;
    // action_from_elle(&m_id, &m_i, &m_j, &m_ji, ell, e);


    proj iP0, iP1, iP2; 

    mont0_dist(&iP0, &P0);
    mont0_dist(&iP1, &P1);
    mont0_dist(&iP2, &P2);


    proj P0_iP0, P1_iP1, P2_iP2;


    mont0_id_plus_dist(&P0_iP0, &P0 , &global_setup.E0);
    mont0_id_plus_dist(&P1_iP1, &P1 , &global_setup.E0);
    mont0_id_plus_dist(&P2_iP2, &P2 , &global_setup.E0);


    proj a_ib_0, a_ib_1, a_ib_2, a_ib_K;


    //endo = gmod(algmul(global_setup.B, endo, mkcol4s(0,1,0,0)), gelle);


    gentobig(&x, gel(endo,1));
    gentobig(&y, gel(endo,2));
    xBIDIM(&a_ib_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&a_ib_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&a_ib_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&a_ib_K, &global_setup.E0, &a_ib_0, &x, &a_ib_1, &y, &a_ib_2);

    proj c_id_0, c_id_1, c_id_2, c_id_K;

    gentobig(&x, gel(endo,3));
    gentobig(&y, gel(endo,4));
    xBIDIM(&c_id_0, &global_setup.E0, &P0, &x, &iP0, &y, &P0_iP0);
    xBIDIM(&c_id_1, &global_setup.E0, &P1, &x, &iP1, &y, &P1_iP1);
    xBIDIM(&c_id_2, &global_setup.E0, &P2, &x, &iP2, &y, &P2_iP2);

    gentobig(&x, gel(v,1));
    gentobig(&y, gel(v,2));
    xBIDIM(&c_id_K, &global_setup.E0, &c_id_0, &x, &c_id_1, &y, &c_id_2);

    proj jc_jid_K;
    mont0_frob(&jc_jid_K, &c_id_K);

    //printf("compare: %s %s\n", proj_hash(a_ib_K), proj_hash(jc_jid_K));

    return (mont_equal(&jc_jid_K, &a_ib_K));

}















static bool is_in_kernel_of_endo_odd(GEN endo, GEN v, proj *pt, bool twist) {

    long len = (twist) ? p_minus_len : p_plus_len;
    const long *fact = (twist) ? p_minus_fact : p_plus_fact;
    const long *mult = (twist) ? p_minus_mult : p_plus_mult;


    for (int i = 0; i < len; ++i) {
        long ell = fact[i];
        long e = mult[i];

        GEN gelle = powuu(ell,e);
        GEN cof = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);
        cof = gdiv(cof, gelle);
        uintbig cof_big;

        gentobig(&cof_big,cof);

        proj P;
        xMUL(&P, &global_setup.E0, pt, &cof_big);

        if (!(is_in_kernel_of_endo_ell(gmod(endo,gelle), gmod(v,gelle), &P, ell))) {
            printf("is_in_kernel_of_endo_odd failed at ell = %ld\n", ell);
            return false;
        }
    }

    // the above method does not detect distorsions. The one below should be finer

    proj E;
    mont_to_ted(&E, &global_setup.E0, twist);
    point K, K1, K2, K_endo;

    uintbig x,y;
    proj K_mont, K_mont_bis;


    const point *basis_ted = (twist) ? torsion_basis_twist_ted_sum : torsion_basis_ted_sum;
    const proj *basis = (twist) ? torsion_basis_twist_sum : torsion_basis_sum;
    GEN order = (twist) ? famat_prod(global_setup.gen_p_minus_fact) : famat_prod(global_setup.gen_p_plus_fact);

    K1 = basis_ted[0];
    K2 = basis_ted[1];
    gentobig(&x, gel(v,1));
    ted_mul(&K1, &K1, &E, &x);
    gentobig(&y, gel(v,2));
    ted_mul(&K2, &K2, &E, &y);

    ted_add(&K, &E, &K1, &K2);
    endo = gmod(endo,order);
    K_endo = ted_action(endo, &K, &E, twist);
    ted_to_mont_point(&K_mont, &K);
    xBIDIM(&K_mont_bis, &global_setup.E0, &basis[0], &x, &basis[1], &y, &basis[2]);

    assert(mont_equal(&K_mont, &K_mont_bis));

    return ted_iszero(&K_endo);


}




int test_rootn() {
    float accumulated_time = 0.;
    int repetitions = 5;
    clock_t t;

    unsigned int length_n = 10;

    GEN A = global_setup.B;
    GEN p = global_setup.p;
    GEN order = global_setup.O0;

    // GEN s = global_setup.gen_odd_torsion;
    GEN n, rho;

    do {
        n = randomprime(powiu(gen_2, length_n));
        rho = Fp_sqrt(n, p);
    } while (!rho);

    GEN sol = qfsolve(mkmat4(mkcol4(p,gen_0,gen_0,gen_0),
        mkcol4(gen_0,p,gen_0,gen_0),
        mkcol4(gen_0,gen_0,gen_1,gen_0),
        mkcol4(gen_0,gen_0,gen_0,gneg(n))));

    GEN eta = mkcol4(gen_0, gel(sol,3), gel(sol,1), gel(sol,2));
    GEN s = gel(sol,4);

    //GEN x = algmul(A,eta,eta);
    //x = diviiexact(gel(x,1),gmul(s,s));

    GEN I = lideal_create(A, order, eta, s);

    GEN O = alglatrighttransporter(A, lideal_lattice(I), lideal_lattice(I));

    GEN theta = algdivr(A,eta,alg_scalar(A,s));

    assert(alglatcontains(A, order, eta, NULL));
    assert(alglatcontains(A, O, theta, NULL));

    // output(p);
    // output(sol);
    // output(eta);
    // output(n);
    // output(x);

    // GEN k;
    // do {
    //     k = randomi();
    // }

    return 1;
}






int test_diamond() {
    float accumulated_time = 0.;
    int repetitions = 5;
    clock_t t;
    GEN A = global_setup.B;
    GEN order = global_setup.O0;
    GEN fm = famat_mul(global_setup.gen_p_plus_fact, global_setup.gen_p_minus_fact);

    for (int i = 0; i < repetitions; i++) {
        t = tic();

        long e1 = 33;
        long e2 = 33;
        GEN Nodd = gsqr(global_setup.gen_odd_torsion);
        GEN Ntwo = powiu(gen_2, e1+e2);


        // STEP 1: Find an endomorphism gamma of norm Nodd*Ntwo = T^2 * 2^e1 * 2^e2,
        //         where T = global_setup.gen_odd_torsion

        GEN gamma = NULL;
        while (!gamma) {
            // parity option is 1 so the 2-walk is not backtracking
            gamma = norm_equation_special(global_setup.p, gmul(Nodd,Ntwo), 1, true);
        }
        gamma = gtrans(gamma);
        GEN n;
        gamma = alg_primitive(&n, A, order, gamma);

        GEN N = algnorm(A,gamma,0);
        Nodd = ggcd(Nodd, N);

        assert(gcmp(N, gmul(Nodd,Ntwo)) == 0);
        assert(Z_lval(N, 2) == e1+e2);

        // long m = remove_1_i(A,order,&gamma);
        // assert(m < e1);
        // e1 -= m;

        GEN gamma_conj = alg_conj(A, gamma);
        
        // m = remove_1_i(A,order,&gamma_conj);
        // assert(m < e2);
        // e2 -= m;
        // gamma = alg_conj(A, gamma_conj);
        // printf("e1 = %ld; e2 = %ld\n", e1, e2);

        // GEN gamma_i = algmul(A, gamma, mkcol4s(0,1,0,0));
        // GEN gamma_conj_i = algmul(A, gamma_conj, mkcol4s(0,1,0,0));


        // STEP 2: Compute the kernel of the degree T and 2^e1 initial isogenies of 
        //         gamma, and the degree T and 2^e2 initial isogenies of gamma_conj
        //         The corresponding isogenies are psi_1 (T), phi_1 (2^e1), 
        //         psi_2 (T), and phi_2 (2^e2), 

        GEN H1_odd = lideal_create(A, order, gamma, ggcd(global_setup.gen_odd_torsion, Nodd));
        GEN H1_two = lideal_create(A, order, gamma, powuu(2, e1));

        GEN H2_odd = lideal_create(A, order, gamma_conj, gdiv(Nodd, lideal_norm(H1_odd)));
        GEN H2_two = lideal_create(A, order, gamma_conj, powuu(2, e2));

        odd_isogeny psi_1 = ideal_to_isogeny_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));
        odd_isogeny psi_2 = ideal_to_isogeny_O0_T(H2_odd, famat_Z_gcd(fm,lideal_norm(H2_odd)));

        proj psi_1_source = global_setup.E0;
        proj psi_2_source = global_setup.E0;
        two_walk phi_1 = ideal_to_isogeny_O0_two(H1_two);
        two_walk phi_2 = ideal_to_isogeny_O0_two(H2_two);



        // STEP 3: Check the isogeny correctness of the kernels from previous step

        GEN coeff_1 = ideal_to_kernel_O0_T(H1_odd, famat_Z_gcd(fm,lideal_norm(H1_odd)));
        GEN v_plus_1 = torsion_crt_compose(gel(coeff_1,1), false);
        GEN v_minus_1 = torsion_crt_compose(gel(coeff_1,2), true);

        assert(is_in_kernel_of_endo_odd(gamma, v_plus_1, &psi_1.kernel_plus, false));
        assert(is_in_kernel_of_endo_odd(gamma, v_minus_1, &psi_1.kernel_minus, true));

        GEN coeff_2 = ideal_to_kernel_O0_T(H2_odd, famat_Z_gcd(fm,lideal_norm(H2_odd)));
        GEN v_plus_2 = torsion_crt_compose(gel(coeff_2,1), false);
        GEN v_minus_2 = torsion_crt_compose(gel(coeff_2,2), true);

        assert(is_in_kernel_of_endo_odd(gamma_conj, v_plus_2, &psi_2.kernel_plus, false));
        assert(is_in_kernel_of_endo_odd(gamma_conj, v_minus_2, &psi_2.kernel_minus, true));

        GEN v1 = ideal_to_kernel_O0_ell(H1_two, 2);
        GEN endo1 = gmod(gamma, powuu(2,33));
        assert(is_in_kernel_of_endo_two(endo1, v1, phi_1.ker));

        GEN v2 = ideal_to_kernel_O0_ell(H2_two, 2);
        GEN endo2 = gmod(gamma_conj, powuu(2,33));
        assert(is_in_kernel_of_endo_two(endo2, v2, phi_2.ker));


        // STEP 3: Compute psi_1*phi_1, phi_1*psi_1, psi_2*phi_2, phi_2*psi_2

        two_walk phi_1_pushed = push_two_walk_through_odd_isogeny(&phi_1, &psi_1, &psi_1_source);
        two_walk phi_2_pushed = push_two_walk_through_odd_isogeny(&phi_2, &psi_2, &psi_2_source);
        
        proj psi_1_pushed_source = psi_1_source;
        proj psi_2_pushed_source = psi_2_source;
        odd_isogeny psi_1_pushed = push_odd_isogeny_through_two_walk(&psi_1, &psi_1_pushed_source, &phi_1);
        odd_isogeny psi_2_pushed = push_odd_isogeny_through_two_walk(&psi_2, &psi_2_pushed_source, &phi_2);

        proj A1, A2, pt1 = phi_1_pushed.ker, pt2 = phi_2_pushed.ker;
        isomorphism isom1, isom2;
        two_walk phi_1_pushed_adjusted;
        two_walk phi_2_pushed_adjusted;
        eval_walk_isom(&isom1, &phi_1_pushed_adjusted, &A1, &pt1, &phi_1_pushed, &pt1);
        eval_walk_isom(&isom2, &phi_2_pushed_adjusted, &A2, &pt2, &phi_2_pushed, &pt2);

        assert(fp2_iszero(&pt1.z) && !fp2_iszero(&pt1.x));
        assert(fp2_iszero(&pt2.z) && !fp2_iszero(&pt2.x));

        pt1 = psi_1_pushed.kernel_plus;
        pt2 = psi_2_pushed.kernel_plus;

        proj B1 = psi_1_pushed_source;
        eval(&B1, & psi_1_pushed, &pt1);
        proj B2 = psi_2_pushed_source;
        eval(&B2, & psi_2_pushed, &pt2);

        assert(fp2_iszero(&pt1.z) && !fp2_iszero(&pt1.x));
        assert(fp2_iszero(&pt2.z) && !fp2_iszero(&pt2.x));

        // STEP 3: Check that psi_1*phi_1, phi_1*psi_1, psi_2*phi_2, phi_2*psi_2
        //         all have the same target curve

        proj jA1,jA2;
        jinv256(&jA1, &A1);
        jinv256(&jA2, &A2);

        proj jB1,jB2;
        jinv256(&jB1, &B1);
        jinv256(&jB2, &B2);

        assert(mont_equal(&jA1,&jB1));
        assert(mont_equal(&jA2,&jB2));
        assert(mont_equal(&jA1,&jA2));

        accumulated_time += toc(t);
    }

    printf("average time\t [%f ms]\n",  (accumulated_time / repetitions));
    return 1;
}


// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    test_rootn();

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);
}



