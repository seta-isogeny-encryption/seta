#define _XOPEN_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <math.h>
#include <assert.h>

#include "ideal.h"
#include "toolbox.h"
#include "klpt.h"
#include "precomputed.h"


struct klpt_test_setup_t {
    GEN p; // the prime
    GEN B; // the quaternion algebra
    GEN O0; // the cannonical maximal order
    GEN torsion_fm; // factorisation matrix of the available torsion
    long q;
    long bits;
};

struct klpt_test_setup_t test_setup;




GEN norm0(GEN x) {
    return algnorm(test_setup.B, x,0);
}

/*
static GEN alg_O0_to_standard(GEN elt) {
    return RgM_RgC_mul(test_setup.O0_to_standard, elt);
}

static GEN alg_standard_to_O0(GEN elt) {
    return RgM_RgC_mul(test_setup.standard_to_O0, elt);
}
*/


int test_lideal1() {
    GEN B_O0 = test_setup.O0;
    GEN B = test_setup.B;

    GEN gen = mkcol4s(2,7,3,-5);
    GEN lideal = lideal_create(B, B_O0, gen, NULL);

    //output(lideal);
    output(norm0(gen));
    output(lideal_norm(lideal));


    gel(lideal,3) = NULL;
    GEN z = lideal_generator(lideal);

    output(norm0(z));
    GEN lideal2 = lideal_create(B, B_O0, z, lideal_norm(lideal));

    //output(lideal2);
    gel(lideal2,2) = gen_0;
    output(lideal_norm(lideal2));

    return 0;
}



int test_lideal2() {
    GEN O0 = test_setup.O0;
    GEN B = test_setup.B;

    GEN gen = lattice_random(B, O0, stoi(10000)); // = mkcol4s(2,7,0,0);

    output(algnorm(B, gen, 0));

    GEN lideal = lideal_create(B, O0, gen, NULL);

    printf("LLL\n");
    //output(lideal_lll(lideal));

    //printf("Short\n");
    output(norm0(gel(lideal_lll(lideal),1)));
    output(norm0(gel(lideal_lll(lideal),2)));
    output(norm0(gel(lideal_lll(lideal),3)));
    output(norm0(gel(lideal_lll(lideal),4)));

    //lideal_equiv_prime(lideal);


    GEN I1 = lideal_random_2e(B, O0, 256);

    output(lideal_norm(I1));

    GEN I2 = lideal_equiv_prime(I1,NULL);

    output(lideal_norm(I2));


    lideal_isom(I1, I2);

    output(algmultable(B));

    return 0;
}

int test_norm_eq() {
    float accumulated_time_ms = 0.;
    int repetitions = 10000;
    clock_t t;

    //GEN fm = famat_sqr(test_setup.torsion_fm);

    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        long ctr = 1;
        GEN sol = NULL, target;
        //GEN N = randomi(mpshift(gen_1, 128));
        //GEN p_div_N = gadd(truedivii(test_setup.p, N), gen_1);
        //GEN fm_1;

        t = tic();
        while (!sol) {
            //fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 2));

            //target = gmul(N,famat_prod(fm_1));
            target = randomi(gmul(test_setup.p, stoi(ctr/20 + 1)));
            sol = norm_equation_special(test_setup.p, target, 0, false);
            ctr++;
        }
        //output(sol);
        accumulated_time_ms += toc(t);

        GEN sum1 = gadd(gsqr(gel(sol,1)), gsqr(gel(sol,2)));
        GEN sum2 = gmul(test_setup.p,gadd(gsqr(gel(sol,3)), gsqr(gel(sol,4))));
        GEN sum = gadd(sum1,sum2);
        if (gcmp(sum,target) != 0) {
            printf("error in norm_equation_special: incorrect result!\n");
            break;
        }

        avma = av;
    }
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 0;
}

int test_famat_rand() {
    float accumulated_time_ms = 0.;
    int repetitions = 100;
    clock_t t;

    GEN fm = famat_sqr(test_setup.torsion_fm);

    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN target;
        GEN N = randomi(mpshift(gen_1, 128));
        GEN B = gmul(gadd(truedivii(test_setup.p, N), gen_1),stoi(3));
        GEN fm_1;

        t = tic();
        for (int i = 0; i < 100; ++i) {
            //fm_1 = famat_random(fm, gmulgs(p_div_N,(ctr/10) + 2));
            fm_1 = famat_random(fm, B);
            target = gmul(N,famat_prod(fm_1));
        }
        //output(sol);
        accumulated_time_ms += toc(t);

        avma = av;
    }
    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    return 0;
}

int test_famat() {
    GEN fm = mkmat2(mkcol4s(2,3,5,7), mkcol4s(4,2,1,2));
    GEN fm2;
    output(famat_pop(fm,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));
    output(famat_pop(fm2,&fm2));



    fm = mkmat2(mkcol4s(2,3,5,7), mkcol4s(30,2,3,2));
    output(fm);

    for (int i = 0; i < 20; ++i) {
        fm2 = famat_random(fm,stoi(10));
        output(fm2);
    }

    return 0;
}


int test_klpt() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 16;

    GEN torsion_fm = test_setup.torsion_fm;
    GEN B = test_setup.B;
    GEN O0 = test_setup.O0;

    pari_sp av = avma;
    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(B, O0, test_setup.bits);

        // assert(lideal_cyclic(&I));
        // output(lideal_lattice(I));

        GEN fm = famat_sqr(torsion_fm);
        gel(gel(fm,2),1) = gen_0; // remove power of 2

        clock_t t = tic();
        // printf("klpt_special_smooth call\n");
        GEN J = klpt_special_smooth(I, fm);
        // printf("klpt_special_smooth done\n");
        accumulated_time_ms += toc(t);
        //accumulated_time_ms += toc(t);

        GEN NJ = lideal_norm(J);
        accumulated_bitlength += dbllog2r(itor(NJ,10));
        //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));

        // check norm

        int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
        if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

        // check isomorphism

        GEN alpha = lideal_isom(I, J); // I*alpha = J
        if (!alpha) { printf("output of klpt is not isomorphic to input\n"); break; }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));
    return 0;
}

int test_equiv_nearprime() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 100;
    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(test_setup.B, test_setup.O0, 130);
        GEN fm = famat_sqr(test_setup.torsion_fm);
        gel(gel(fm,2),1) = gen_0; // remove power of 2

        clock_t t = tic();
        GEN J = lideal_equiv_nearprime(I,fm,0);
        // GEN J = lideal_equiv_prime_except(I,NULL,NULL);
        accumulated_time_ms += toc(t);

        //GEN NJ = lideal_norm(J);
        //if (!ispseudoprime(NJ,0)) { printf("output of lideal_equiv_prime_except is not of prime norm\n"); break; }

        // check isomorphism

        GEN alpha = lideal_isom(I, J); // I*alpha = J
        if (!alpha) { printf("output of lideal_equiv_* is not isomorphic to input\n"); break; }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));

    return 0;
}

int test_klpt2e() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 16;
    pari_sp av = avma;

    for (int i = 0; i < repetitions; ++i) {
        GEN I = lideal_random_2e(test_setup.B, test_setup.O0, 64);
        GEN fm = famat_sqr(test_setup.torsion_fm);
        gel(gel(fm,2),1) = gen_0; // remove power of 2

        clock_t t = tic();
        GEN J = klpt_special_smooth_small_2e_input(I, fm);
        accumulated_time_ms += toc(t);
        //accumulated_time_ms += toc(t);

        if (J) {
            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));
            //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));

            // check norm

            int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
            if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

            // check isomorphism

            GEN alpha = lideal_isom(I, J); // I*alpha = J
            if (!alpha) { printf("output of klpt is not isomorphic to input\n"); }
        }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / repetitions));

    return 0;
}

int test_klpt_general() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 15;
    int win = 0;


    GEN p = test_setup.p;
    GEN B = test_setup.B;
    GEN O0 = test_setup.O0;
    long q = 1;



    pari_sp av = avma;
    for (int i = 0; i < repetitions; ++i) {
        unsigned int length_NI = 64;

        GEN NI = NULL;

        do {
            NI = randomprime(powiu(gen_2, length_NI));
        } while (Fp_issquare(gen_2,NI));

        GEN alpha = NULL;

        unsigned int margin = 6;
        while (!alpha) {
            alpha = norm_equation_special_q(p, q, gmul(NI,randomi(powiu(gen_2, 256-length_NI + margin))));
            ++margin;
        }


        GEN I = lideal_create(B, O0, gtrans(alpha), NI);

        GEN K = lideal_random_2e(B, O0, 130);
        K = lideal_equiv_prime(K,NULL);
        clock_t t = tic();
        GEN J = klpt_general_power(I, K, gen_2);
        accumulated_time_ms += toc(t);


        if (J) {
            ++win;

            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));

            //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));


            // check norm

            //int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
            //if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

            //output(Z_factor_limit(NJ, 3));

            // check isomorphism

            GEN I1 = lideal_inter(I,K);
            GEN I2 = lideal_inter(I,J);

            GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
            if (!alpha) { printf("output of klpt is not isomorphic to input!\n"); break; }
        }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

    return 0;
}

int test_klpt_power_simple() {

    float accumulated_time_ms = 0., accumulated_bitlength = 0.;
    int repetitions = 15;
    int win = 0;

    // GEN torsion_fm = test_setup.torsion_fm;
    // GEN p = test_setup.p;
    GEN B = test_setup.B;
    GEN O0 = test_setup.O0;
    // long q = test_setup.q;
    GEN B_1 = mkcol4s(1,0,0,0);

    pari_sp av = avma;
    for (int i = 0; i < repetitions; ++i) {

        GEN I = lideal_create(B, O0, B_1, gen_1);

        GEN K = lideal_random_2e(B, O0, test_setup.bits);
        K = lideal_equiv_prime(K,NULL);
        clock_t t = tic();
        printf("klpt_general_power call\n");
        GEN J = klpt_general_power(I, K, gen_2);
        printf("klpt_general_power done\n");
        accumulated_time_ms += toc(t);


        if (J) {
            ++win;

            GEN NJ = lideal_norm(J);
            accumulated_bitlength += dbllog2r(itor(NJ,10));

            //printf("%ld bits\n", (long)dbllog2r(itor(NJ,10)));


            // check norm

            //int smooth_norm = (gcmp(famat_prod(famat_Z_gcd(fm, NJ)), NJ) == 0);
            //if (!smooth_norm) { printf("output of klpt does not have a valid norm\n"); break; }

            //output(Z_factor_limit(NJ, 3));

            // check isomorphism

            GEN I1 = lideal_inter(I,K);
            GEN I2 = lideal_inter(I,J);

            GEN alpha = lideal_isom(I1, I2); // I1*alpha = I2
            if (!alpha) { printf("output of klpt is not isomorphic to input!\n"); break; }
        }

        avma = av;
    }

    printf("average time\t [%f ms]\n",  (accumulated_time_ms / repetitions));
    printf("average length\t %d bits\n", (int) (accumulated_bitlength / win));

    return 0;
}



// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(80000000, 1<<18);

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    // long var = fetch_var();
    // GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);

    // GEN a = stoi(-1),
    //     p = strtoi("73743043621499797449074820543863456997944695372324032511999999999999999999999"),
    //     b = negi(p);

    // GEN B = alg_hilbert(nf, a, b, var, 0);

    // GEN B_1 = mkcol4s(1,0,0,0);
    // GEN B_i = mkcol4s(0,1,0,0);
    // GEN B_j = mkcol4s(0,0,1,0);
    // GEN B_ji = mkcol4s(0,0,0,1);
    // //GEN B_ij = mkcol4s(0,0,0,-1);

    // GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
    // GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2

    // GEN B_O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0);

    // test_setup.p = p;
    // test_setup.B = B; // the quaternion algebra
    // test_setup.qf = mkmat4(mkcol4s(1,0,0,0),
    //                          mkcol4s(0,1,0,0),
    //                          mkcol4(gen_0,gen_0,p,gen_0),
    //                          mkcol4(gen_0,gen_0,gen_0,p)); // quadratic form defined by the reduced norm


    // test_setup.O0 = B_O0; // the cannonical maximal order
    // test_setup.one = B_1;
    // test_setup.i = B_i;
    // test_setup.j = B_j;
    // test_setup.ji = B_ji;

    // test_setup.O0_b1 = B_1;
    // test_setup.O0_b2 = B_i;
    // test_setup.O0_b3 = B_1k_2;
    // test_setup.O0_b4 = B_ij_2;
    // test_setup.O0_to_standard = mkmat4(B_1, B_i, B_1k_2, B_ij_2);
    // test_setup.standard_to_O0 = RgM_inv(test_setup.O0_to_standard);

    long var = fetch_var();
    GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);

    GEN B_1 = mkcol4s(1,0,0,0);
    GEN B_i = mkcol4s(0,1,0,0);



    test_setup.torsion_fm = Z_factor_limit(strtoi(
        "197530174297949459837634878151545563369632855190375548677707409417459236752253845947265965991865263091519488000000000000000000000"
        ), 30000);

    test_setup.torsion_fm = famat_sqr(test_setup.torsion_fm);


    long c,q;
    GEN multtable, B1, B2, B3, B4;






    /*************** p = 1 mod 8 and j = 0 ********************/
    test_setup.bits = 397;
    test_setup.p = strtoi("256170376103179994761914482856266627791915646373582727175785233336597312364178751161184887957805876215573849680566141281");
    test_setup.q = 3;
    c = 1;
    test_setup.B = alg_hilbert(nf, stoi(-test_setup.q), negi(test_setup.p), var, 0);

    q = test_setup.q;
    multtable = mkvec4(
        mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
        mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q,0)),
        mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(test_setup.p),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(test_setup.p),gen_0,gen_0)),
        mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q,0), mkcol4(gen_0,test_setup.p,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q),test_setup.p)),gen_0,gen_0,gen_0))

        );

    test_setup.B = alg_csa_table(nf, multtable, var,0);

    B1 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    B2 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));
    B3 = mkcol4(gen_0,gdiv(gen_1,stoi(test_setup.q)),gen_0,gneg(gdiv(stoi(c),stoi(test_setup.q))));
    B4 = mkcol4(gen_0,gen_0,gen_0,gneg(gen_1));

    // GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
    // GEN B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    // GEN B3 = mkcol4(gen_0,gdiv(gen_2,stoi(11)),gen_0,gneg(gdiv(gen_1,stoi(11))));
    // GEN B4 = mkcol4(gen_0,gneg(gdiv(gen_1,stoi(11))),ghalf,(gdiv(gen_1,stoi(22))));
    test_setup.O0 = alglathnf(test_setup.B,mkmat4(B1,B2,B3,B4), gen_0); // the cannonical maximal order



    // test_klpt2e();
    // test_klpt();
    test_klpt_power_simple();





    /*************** p = 1 mod 8 ********************/
    test_setup.bits = 112;
    test_setup.p = strtoi("4565346563423453435142634556476321");
    test_setup.q = 11;
    c = 5;
    test_setup.B = alg_hilbert(nf, stoi(-test_setup.q), negi(test_setup.p), var, 0);

    q = test_setup.q;
    multtable = mkvec4(
        mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
        mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q,0)),
        mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(test_setup.p),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(test_setup.p),gen_0,gen_0)),
        mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q,0), mkcol4(gen_0,test_setup.p,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q),test_setup.p)),gen_0,gen_0,gen_0))

        );

    test_setup.B = alg_csa_table(nf, multtable, var,0);

    B1 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    B2 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));
    B3 = mkcol4(gen_0,gdiv(gen_1,stoi(test_setup.q)),gen_0,gneg(gdiv(stoi(c),stoi(test_setup.q))));
    B4 = mkcol4(gen_0,gen_0,gen_0,gneg(gen_1));

    // GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
    // GEN B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    // GEN B3 = mkcol4(gen_0,gdiv(gen_2,stoi(11)),gen_0,gneg(gdiv(gen_1,stoi(11))));
    // GEN B4 = mkcol4(gen_0,gneg(gdiv(gen_1,stoi(11))),ghalf,(gdiv(gen_1,stoi(22))));
    test_setup.O0 = alglathnf(test_setup.B,mkmat4(B1,B2,B3,B4), gen_0); // the cannonical maximal order



    test_klpt2e();
    test_klpt();
    test_klpt_power_simple();





    

    /*************** p = 5 mod 8 ********************/
    test_setup.bits = 112;
    test_setup.p = strtoi("4565346563423453435142634556475949");
    test_setup.q = 2;
    test_setup.B = alg_hilbert(nf, stoi(-test_setup.q), negi(test_setup.p), var, 0);
    GEN B_2_i_ij_4 = mkcol4(ghalf,gneg(gmul(ghalf,ghalf)),gen_0,gneg(gmul(ghalf,ghalf)));
    GEN B_1_i_j_2 = mkcol4(gneg(ghalf),ghalf,ghalf,gen_0);
    test_setup.O0 = alglathnf(test_setup.B,mkmat4(B_1, B_i, B_2_i_ij_4, B_1_i_j_2), gen_0); // the cannonical maximal order



    test_klpt2e();
    test_klpt();
    test_klpt_power_simple();





    /*************** p = 3 mod 4 ********************/
    test_setup.bits = 112;
    test_setup.p = strtoi("4565346563423453435142634556476019");
    test_setup.q = 1;
    test_setup.B = alg_hilbert(nf, stoi(-test_setup.q), negi(test_setup.p), var, 0);
    GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
    GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2
    test_setup.O0 = alglathnf(test_setup.B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0); // the cannonical maximal order



    test_klpt2e();
    test_klpt();
    //test_klpt_general();
    test_klpt_power_simple();




    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);
}
