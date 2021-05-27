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

#include "mont.h"
#include "tedwards.h"
#include "constants.h"
#include "precomputed.h"
#include "curve.h"


// const uintbig cof_p_plus1 = {1ULL,0,0,0};
// const uintbig cof_p_plus2 = {1ULL,0,0,0};



// const char* p_str = "8352547421173864960624333703667331271395143445202337668463192809880783896176785654136882064526006534223584195029346241";
// const char* all_the_torsion_str = "69765048422958181897960319077352793214907650293450379167416508328959226975464456756427967097183686868327009326439047901292946921759400717445370636114701166739708048309054296279586015687227479047950787096565629295598124341383191860830080";



// const uintbig cof_p_plus1 = {517434778561ULL,0,0,0};
// const uintbig cof_p_plus2 = {26602537156291ULL,0,0,0};



// const char* p_str = "73743043621499797449074820543863456997944695372324032511999999999999999999999";
// const char* all_the_torsion_str = "197530174297949459837634878151545563369632855190375548677707409417459236752253845947265965991865263091519488000000000000000000000";

struct quaternion_setup_t {
    GEN p; // the prime
    uintbig p_uint; // the prime
    GEN B; // the quaternion algebra
    GEN qf; // the quaternion algebra
    GEN O0; // the cannonical maximal order
    GEN one;
    GEN i;
    GEN j;
    GEN ji;
    GEN torsion_fm; // factorisation matrix of the available torsion

    GEN O0_b1;
    GEN O0_b2;
    GEN O0_b3;
    GEN O0_b4;
    GEN O0_to_standard;
    GEN standard_to_O0;

    proj E0;
};

struct quaternion_setup_t precomp_setup;



uintbig stobig(long long x) {
    uintbig x_big;
    uintbig_set(&x_big, x);
    return x_big;
}

char* pari_int_code(GEN i) {
    if (is_bigint(i)) return pari_sprintf("strtoi(\"%Ps\")", i);
    else return pari_sprintf("stoi(%PsULL)", i);
}

char* pari_2x2_matrix_code(GEN M) {
    return  pari_sprintf("mkmat2(mkcol2(%s,%s),mkcol2(%s,%s))", pari_int_code(gcoeff(M,1,1)),
            pari_int_code(gcoeff(M,2,1)),
            pari_int_code(gcoeff(M,1,2)),
            pari_int_code(gcoeff(M,2,2)));
}


char* fp_code(const fp *x) {
    return pari_sprintf("{ %luULL, %luULL, %luULL, %luULL, %luULL, %luULL, %luULL }", x->x.c[0], x->x.c[1], x->x.c[2], x->x.c[3] , x->x.c[4], x->x.c[5], x->x.c[6]);
}

char* fp2_code(const fp2 *x) {
    return pari_sprintf("{ %s,\n %s }", fp_code(&x->re), fp_code(&x->im));
}

char* proj_code(const proj *P) {
    if (fp2_iszero(&P->z))
        return pari_sprintf("{ %s,\n %s }", fp2_code(&P->x), fp2_code(&P->z));
    else {
        fp2 x,y,z,tmp = P->z;
        fp2_inv(&tmp);
        fp2_mul3(&x,&P->x,&tmp);
        fp2_mul3(&z,&P->z,&tmp);
        return pari_sprintf("{ %s,\n %s }", fp2_code(&x), fp2_code(&z));
    }
}

char* proj2_code(const proj2 *P) {
    if (fp2_iszero(&P->z))
        return pari_sprintf("{ %s,\n %s,\n %s }", fp2_code(&P->x), fp2_code(&P->y), fp2_code(&P->z));
    else {
        fp2 x,y,z,tmp = P->z;
        fp2_inv(&tmp);
        fp2_mul3(&x,&P->x,&tmp);
        fp2_mul3(&y,&P->y,&tmp);
        fp2_mul3(&z,&P->z,&tmp);
        return pari_sprintf("{ %s,\n %s,\n %s }", fp2_code(&x), fp2_code(&y), fp2_code(&z));
    }
}

char* ted_code(const point *P) {
    return pari_sprintf("{ %s,\n %s,\n %s,\n %s }", fp2_code(&P->x), fp2_code(&P->y), fp2_code(&P->z), fp2_code(&P->t));
}

void print_loop(char* prefix, char* name1, char* name2) {
    printf("%sfor (int i = 0; i < 3; ++i) {\n",prefix);
    printf("%s\tfp_enc( &(&(&%s[i])->x)->re, &(%s[i][0][0]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->x)->im, &(%s[i][0][1]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->z)->re, &(%s[i][1][0]) );\n",prefix,name1,name2);
    printf("%s\tfp_enc( &(&(&%s[i])->z)->im, &(%s[i][1][1]) );\n",prefix,name1,name2);
    printf("%s}\n",prefix);
}

GEN norm0(GEN x) {
    return algnorm(precomp_setup.B, x,0);
}

// // this is broken
// void fp2_random(fp2 *x) {
//   uint64_t thrash;
//   for (int i = 0; i < FP_LIMBS; ++i) {
//     x->re.x.c[i] = random_Fl(0xffffffffffffffff);
//     x->im.x.c[i] = random_Fl(0xffffffffffffffff);
//   }
//   mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->re.x.c, 0, (mp_srcptr)x->re.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
//   mpn_tdiv_qr((mp_ptr)&thrash, (mp_ptr)x->im.x.c, 0, (mp_srcptr)x->im.x.c, FP_LIMBS, (mp_srcptr)p.c, FP_LIMBS);
// }


void random_point(proj *P, proj const *A, long ell, long e, bool twist) {
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

void random_basis(proj *P1, proj *P2, point *P1_ted, point *P2_ted, proj const *A, long ell, long e, bool twist) {
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

void random_basis_2e(proj *P1, proj *P2, proj const *A, long e, bool twist) {
    proj P1_mul, P2_mul, tmp;
    uintbig ell_big;
    long ell = 2;
    uintbig_set(&ell_big, ell);

    random_point(P1, A, ell, e, twist);

    P1_mul = *P1;
    for (int i = 0; i < e-1; ++i) {
        xMUL(&P1_mul, A, &P1_mul, &ell_big);
    }

    assert(is_on_curve(&P1_mul,A));
    assert(!mont_iszero(&P1_mul));
    xMUL(&tmp, A, &P1_mul, &ell_big);
    assert(mont_iszero(&tmp));

    do {
        random_point(P2, A, ell, e, twist);
        P2_mul = *P2;
        for (int i = 0; i < e-1; ++i) {
            xMUL(&P2_mul, A, &P2_mul, &ell_big);
        }
        assert(is_on_curve(&P2_mul,A));
        assert(!mont_iszero(&P2_mul));
        xMUL(&tmp, A, &P2_mul, &ell_big);
        assert(mont_iszero(&tmp));

    } while (mont_equal(&P1_mul,&P2_mul));

}

bool fp2_ord_trivial(long *res, const fp2 *g, long bound) {
    long order = 1;
    fp2 x = *g;

    for (int i = 0; i < bound; ++i) {
        if (fp2_equal(&fp2_1,&x)) { *res = order; return true; }
        order++;
        fp2_mul2(&x,g);
    }

    return false;
}

// log of Q is base P1,P2, a basis of the 2-torsion
// ASSUMES Q is in the 2-torsion
void bidim_log_2(long long *a, long long *b, const proj *Q, const proj *P1, const proj *P2) {
    if (mont_iszero(Q))         { *a = 0; *b = 0; }
    else if (mont_equal(Q,P1))  { *a = 1; *b = 0; }
    else if (mont_equal(Q,P2))  { *a = 0; *b = 1; }
    else                        { *a = 1; *b = 1; }
}

// Far from optimal
// P12 = P1 + P2
bool bidim_log_2e(long long *a, long long *b, const proj *A, const proj *Q, const proj *P1, const proj *P2, const proj *P12, long e) {
    uintbig big_2, x, y, tmp_big;
    uintbig_set(&big_2, 2);

    long long log_1 = 0, log_2 = 0, tmp;
    proj Ps1[e], Ps2[e], Ps12[e], Qs[e], R;

    Ps1[0] = *P1;
    Ps2[0] = *P2;
    Ps12[0] = *P12;
    Qs[0] = *Q;

    for (int i = 1; i < e; ++i) {
        xMUL(&Ps1[i], A, &Ps1[i-1], &big_2);
        xMUL(&Ps2[i], A, &Ps2[i-1], &big_2);
        xMUL(&Ps12[i], A, &Ps12[i-1], &big_2);
        xMUL(&Qs[i], A, &Qs[i-1], &big_2);
    }

    // first bit
    bidim_log_2(&log_1, &log_2, &Qs[e-1], &Ps1[e-1], &Ps2[e-1]);

    // next bits
    for (int i = 1; i < e; ++i) {

        uintbig_set(&x, log_1);
        uintbig_set(&y, log_2);

        tmp = (1ULL << i);
        uintbig_set(&tmp_big, tmp);

        xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
        if (mont_equal(&Qs[e-i-1], &R)) {
            // do nothing
        }
        else {
            uintbig_add3(&x, &x, &tmp_big);
            xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
            if (mont_equal(&Qs[e-i-1], &R)) {
                log_1 += tmp;
            }
            else {
                uintbig_set(&x, log_1);
                uintbig_add3(&y, &y, &tmp_big);
                xBIDIM(&R, A, &Ps1[e-i-1], &x, &Ps2[e-i-1], &y, &Ps12[e-i-1]);
                if (mont_equal(&Qs[e-i-1], &R)) {
                    log_2 += tmp;
                }
                else {
                    log_1 += tmp;
                    log_2 += tmp;
                }
            }
                
        }
    }
    *a = log_1;
    *b = log_2;

    // test
    uintbig_set(&x, *a);
    uintbig_set(&y, *b);
    xBIDIM(&R, A, P1, &x, P2, &y, P12);
    assert(mont_equal(Q, &R));
    return true;
}

GEN action_two_3_4(GEN m_i, GEN m_j, long e) {
    GEN m_ij2,m_1ji2;
    GEN gelle = stoi(1LL<<e);
    GEN gelle1 = stoi(1LL<<(e-1));
    GEN iMi, MM, test1, test2;
    GEN pplus14 = gdiv(gadd(precomp_setup.p,gen_1),stoi(4));
    GEN id = mkmat2(mkcol2s(1,0),mkcol2s(0,1));


    for (int i11 = 0; i11 < 2; ++i11){
        for (int i12 = 0; i12 < 2; ++i12){
            for (int i21 = 0; i21 < 2; ++i21){
                for (int i22 = 0; i22 < 2; ++i22){
                    m_ij2 = gdiv(gadd(m_i, m_j),gen_2);
                    if (i11) gcoeff(m_ij2,1,1) = gadd(gcoeff(m_ij2,1,1),gelle1);
                    if (i12) gcoeff(m_ij2,1,2) = gadd(gcoeff(m_ij2,1,2),gelle1);
                    if (i21) gcoeff(m_ij2,2,1) = gadd(gcoeff(m_ij2,2,1),gelle1);
                    if (i22) gcoeff(m_ij2,2,2) = gadd(gcoeff(m_ij2,2,2),gelle1);
                    printf("B\n");
                    output(gmul(m_i,m_ij2));
                    printf("C\n");
                    output(gmul(gmul(m_i,m_ij2),m_i));
                    printf("D\n");
                    output(gmod(gmul(gmul(m_i,m_ij2),m_i),gelle));
                    printf("E\n");

                    iMi = gmod(gmul(gmul(m_i,m_ij2),m_i),gelle);


                    printf("done\n");

                    test1 = gmod(gsub(m_ij2,m_i),gelle);

                    MM = gmod(gmul(m_ij2,m_ij2),gelle);
                    test2 = gmod(gneg(gmul(id,pplus14)),gelle);



                    if (gequal(iMi,test1) && gequal(MM,test2)) { // a candidate
                        m_1ji2 = gmod(gneg(gmul(m_ij2,m_i)),gelle);



                        printf("\t/* candidate %d%d%d%d */\n", i11,i12,i21,i22);
                        printf("\taction_two_3 = %s;\n", pari_2x2_matrix_code(m_1ji2));
                        printf("\taction_two_4 = %s;\n", pari_2x2_matrix_code(m_ij2));
                    }
                }
            }
        }
    }


    return NULL;
}


void compute_action(GEN *M, void (*endo)(point*, const point*), const point *P1, const point *P2, const proj *E, long ell, long e) {
    point Q;
    GEN a,b,c,d;
    assert(ted_is_on_curve(P1,E));
    endo(&Q,P1);
    // proj E1;
    // printf("E0 = %s\n", proj_code(E));
    // fp2_frob2(&E1.x, &E->x);
    // fp2_frob2(&E1.z, &E->z);
    assert(ted_is_on_curve(&Q,E));

    assert(ted_bidim_log(&a, &c, E, &Q, P1, P2, ell, e));

    endo(&Q,P2);
    assert(ted_bidim_log(&b, &d, E, &Q, P1, P2, ell, e));

    *M = mkmat2(mkcol2(a,c),mkcol2(b,d));
}

void check_action(GEN M, void (*endo)(point*, const point*), void (*mont_endo)(proj*, const proj*), const point *P1, const point *P2, const proj *E, const proj *E_mont, long ell, long e) {
    GEN x,y,u,v;
    uintbig X,Y,U,V;
    point Q, endoQ, endoQ_combination, tmp;

    GEN gelle = powuu(ell,e);

    point P12;
    proj M1,M2,M12,MQ,MendoQ,MendoQ_combination;

    ted_to_mont_point(&M1, P1);
    ted_to_mont_point(&M2, P2);
    ted_add(&P12, E, P1, P2);
    ted_to_mont_point(&M12, &P12);


    for (int i = 0; i < 10; ++i) {    
        x = randomi(gelle);
        y = randomi(gelle);
        gentobig(&X, x);
        gentobig(&Y, y);

        ted_mul(&Q, P1, E, &X);
        ted_mul(&tmp, P2, E, &Y);
        ted_add(&Q, E, &Q, &tmp);

        endo(&endoQ,&Q);

        GEN A,B,C,D;
        A = gcoeff(M,1,1);
        B = gcoeff(M,1,2);
        C = gcoeff(M,2,1);
        D = gcoeff(M,2,2);

        u = gadd(gmul(A,x),gmul(B,y));
        v = gadd(gmul(C,x),gmul(D,y));
        gentobig(&U, u);
        gentobig(&V, v);

        ted_mul(&endoQ_combination, P1, E, &U);
        ted_mul(&tmp, P2, E, &V);
        ted_add(&endoQ_combination, E, &endoQ_combination, &tmp);

        assert(ted_equal(&endoQ,&endoQ_combination));


        // same check in montgomery form
        ted_to_mont_point(&MQ, &Q);
        mont_endo(&MendoQ,&MQ);

        xBIDIM(&MendoQ_combination, E_mont, &M1, &U, &M2, &V, &M12);

        assert(mont_equal(&MendoQ,&MendoQ_combination));
    }
}

// P12 = P1 + P2
void compute_action_2e(GEN *M, void (*endo)(proj*, const proj*), const proj *P1, const proj *P2, const proj *P12, const proj *A, long e) {
    proj Q, R;
    long long a,b,c,d,x,y;
    uintbig biga, bigc;

    endo(&Q,P1); 
    assert(is_on_curve(&Q,A));

    assert(bidim_log_2e(&a, &c, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, a);
    uintbig_set(&bigc, c);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));

    endo(&Q,P2);
    assert(bidim_log_2e(&b, &d, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, b);
    uintbig_set(&bigc, d);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));


    endo(&Q,P12);
    assert(bidim_log_2e(&x, &y, A, &Q, P1, P2, P12, e));

    uintbig_set(&biga, x);
    uintbig_set(&bigc, y);
    xBIDIM(&R, A, P1, &biga, P2, &bigc, P12);
    assert(mont_equal(&Q, &R));

    if (((x != (a+b)%(1LL<<e)) && ((1LL<<e)-x != (a+b)%(1LL<<e))) || ((y != (c+d)%(1LL<<e)) && ((1LL<<e)-y != (c+d)%(1LL<<e)))) {
        b = (1LL<<e)-b;
        d = (1LL<<e)-d;
    }

    assert( ((x == (a+b)%(1LL<<e)) && (y == (c+d)%(1LL<<e)))  ||  (((1LL<<e)-x == (a+b)%(1LL<<e)) && ((1LL<<e)-y == (c+d)%(1LL<<e))));

    *M = mkmat2(mkcol2(stoi(a),stoi(c)),mkcol2(stoi(b),stoi(d)));
}

bool check_action_2e(GEN M, void (*endo)(proj*, const proj*), void (*endoxy)(proj2*, const proj2*), const proj2 *P1xy, const proj2 *P2xy, const proj *E, long e) {
    long long x,y;
    uintbig X,Y,U,V;
    proj2 Qxy,Rxy,exy1,exy2,exy12,eQxy_combination,eQxy;

    proj2 P12xy;

    GEN t;

    for (int i = 0; i < 20; ++i) {  
        x = random_Fl(1ULL<<e);
        y = random_Fl(1ULL<<e);
        uintbig_set(&X, x);
        uintbig_set(&Y, y);

        // printf("(x,y) ");
        // output(mkcol2s(x,y));

        endoxy(&exy1,P1xy);
        endoxy(&exy2,P2xy);
        endoxy(&exy12,&P12xy);


        xyMUL(&Qxy, E, P1xy, &X);
        // printf("xP1xy = %s\n",proj2_code(&Qxy));
        xyMUL(&Rxy, E, P2xy, &Y);
        // output(stoi(y));
        // printf("P2xy = %s\n",proj2_code(P2xy));
        // printf("yP2xy = %s\n",proj2_code(&Rxy));
        // assert(!xy_is_zero(&Rxy));
        xyADD(&Qxy, E, &Qxy, &Rxy);

        endoxy(&eQxy,&Qxy);



        long long a,b;
        proj pt1, pt2, pt12,eQ,eQ_combination,tmp;
        proj2 P12xy;
        xytox(&pt1, P1xy);
        xytox(&pt2, P2xy);
        xyADD(&P12xy, E, P1xy, P2xy);
        xytox(&pt12, &P12xy);
        xytox(&eQ, &eQxy);
        assert(bidim_log_2e(&a, &b, E, &eQ, &pt1, &pt2, &pt12, e));

        // xMUL(&tmp, E, &pt1, &X);
        // printf("xP1 = %s\n",proj_code(&tmp));
        // xMUL(&tmp, E, &pt2, &Y);
        // printf("yP2 = %s\n",proj_code(&tmp));


        t = gmod(RgM_RgC_mul(M, mkcol2s(x,y)), stoi(1ULL<<e));
        //printf("e = %ld\n", e);
        // printf("bidim [%ld,%ld] or [%ld,%ld] \n",(1ULL<<e)-a,(1ULL<<e)-b,a,b);
        // printf("coeff ");
        // output(t);
        // printf("mat ");
        // output(M);




        uintbig_set(&U, itos_or_0(gel(t,1)));
        uintbig_set(&V, itos_or_0(gel(t,2)));

        xyMUL(&eQxy_combination, E, P1xy, &U);
        xyMUL(&Rxy, E, P2xy, &V);
        xyADD(&eQxy_combination, E, &eQxy_combination, &Rxy);
        proj2 neg;
        xyNEG(&neg, &eQxy);


        // printf("eQ = %s\n",proj_code(&eQ));
        // xBIDIM(&eQ_combination, E, &pt1, &U, &pt2, &V, &pt12);
        // printf("eQ_combination1 = %s\n",proj_code(&eQ_combination));
        // uintbig_set(&U, a);
        // uintbig_set(&V, b);
        // xBIDIM(&eQ_combination, E, &pt1, &U, &pt2, &V, &pt12);
        // printf("eQ_combination2 = %s\n",proj_code(&eQ_combination));
        // printf("eQxy = %s\n eQxy_combination = %s\n",proj2_code(&eQxy),proj2_code(&eQxy_combination));

        assert(xy_is_on_curve(E, &eQxy_combination));
        assert(xy_equal(&eQxy,&eQxy_combination) || xy_equal(&neg,&eQxy_combination));

        if (xy_equal(&neg,&eQxy_combination)) return false; // flip sign!
    }
    return true;
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

    long var = fetch_var();
    GEN nf = nfinit(pol_x(fetch_var()),LOWDEFAULTPREC);
    
    GEN m1 = mkmat2(mkcol2s(1,0),mkcol2s(0,1));

    // // p6334792777

    // long q = 19;
    // long c = 3;

    // precomp_setup.E0.x =  (fp2){
    // fp_0,
    // { 0xe887b88259aa7fbe, 0xeb339ca64c8d30ac, 0xebc6554e12df73c9,
    //   0x1b4d274bf1982cea, 0x4bd0eb80201b471c, 0x28e9cb4b6564bef7, 0x1d}
    // };
    // precomp_setup.E0.z = fp2_1;


    // p8426067021

    long q = 3;
    long c = -1;

    // Curve y² = x³ + √3 x² + 1, with CM by √-3, j-invariant 0
    init_curve(&precomp_setup.E0);

    // GEN a = stoi(-1);
    GEN p_pari = strtoi(p_str);
    // GEN b = negi(p_pari);

    // GEN B = alg_hilbert(nf, a, b, var, 0);

    GEN B_1 = mkcol4s(1,0,0,0);
    GEN B_i = mkcol4s(0,1,0,0);
    GEN B_j = mkcol4s(0,0,1,0);
    GEN B_ji = mkcol4s(0,0,0,1);
    // //GEN B_ij = mkcol4s(0,0,0,-1);


    // GEN B_1k_2 = mkcol4(ghalf,gen_0,gen_0,gneg(ghalf)); // (1-ji)/2
    // GEN B_ij_2 = mkcol4(gen_0,ghalf,ghalf,gen_0); // (i+j)/2

    // GEN B_O0 = alglathnf(B,mkmat4(B_1, B_i, B_1k_2, B_ij_2), gen_0);



    GEN multtable = mkvec4(
    mkmat4(mkcol4s(1,0,0,0), mkcol4s(0,1,0,0), mkcol4s(0,0,1,0), mkcol4s(0,0,0,1)),
    mkmat4(mkcol4s(0,1,0,0), mkcol4s(-q,0,0,0), mkcol4s(0,0,0,-1), mkcol4s(0,0,q,0)),
    mkmat4(mkcol4s(0,0,1,0), mkcol4s(0,0,0,1), mkcol4(gneg(p_pari),gen_0,gen_0,gen_0), mkcol4(gen_0,gneg(p_pari),gen_0,gen_0)),
    mkmat4(mkcol4s(0,0,0,1), mkcol4s(0,0,-q,0), mkcol4(gen_0,p_pari,gen_0,gen_0), mkcol4(gneg(gmul(stoi(q),p_pari)),gen_0,gen_0,gen_0))

    );

    GEN B = alg_csa_table(nf, multtable, var,0);

    GEN B1 = mkcol4(gen_1,gen_0,gen_0,gen_0);
    GEN B2 = mkcol4(ghalf,ghalf,gen_0,gen_0);
    GEN B3 = mkcol4(gen_0,gen_0,ghalf,gneg(ghalf));
    GEN B4 = mkcol4(gen_0,gdiv(gen_1,stoi(q)),gen_0,gneg(gdiv(stoi(c),stoi(q))));

    GEN B_O0 = alglathnf(B,mkmat4(B1,B2,B3,B4), gen_0); // the cannonical maximal order


    precomp_setup.p = p_pari;
    precomp_setup.B = B; // the quaternion algebra
    precomp_setup.qf = mkmat4(mkcol4s(1,0,0,0),
                             mkcol4s(0,q,0,0),
                             mkcol4(gen_0,gen_0,p_pari,gen_0),
                             mkcol4(gen_0,gen_0,gen_0,gmul(p_pari,stoi(q)))); // quadratic form defined by the reduced norm

    precomp_setup.torsion_fm = Z_factor_limit(strtoi(
        all_the_torsion_str
        ), 70000000);

    precomp_setup.O0 = B_O0; // the cannonical maximal order
    precomp_setup.one = B_1;
    precomp_setup.i = B_i;
    precomp_setup.j = B_j;
    precomp_setup.ji = B_ji;

    precomp_setup.O0_b1 = B1;
    precomp_setup.O0_b2 = B2;
    precomp_setup.O0_b3 = B3;
    precomp_setup.O0_b4 = B4;
    precomp_setup.O0_to_standard = mkmat4(B1, B2, B3, B4);
    precomp_setup.standard_to_O0 = RgM_inv(precomp_setup.O0_to_standard);







    proj E, E_twist;
    mont_to_ted(&E, &precomp_setup.E0, false);
    mont_to_ted(&E_twist, &precomp_setup.E0, true);








    GEN m_frob_two, m_dist_two;
    GEN gelle, inv2, invq;

    point basis_ted[on_curve_len][3];
    proj basis[on_curve_len][3];
















    // two-torsion
    
    proj basis_two[3];
    random_basis_2e(&basis_two[0], &basis_two[1], &precomp_setup.E0, two_tors_height, false);

    mont_add(&basis_two[2], &precomp_setup.E0, &basis_two[0], &basis_two[1]);
    assert(is_on_curve(&basis_two[2], &precomp_setup.E0));




    proj2 basisxy_two[3];

    xtoxy(&basisxy_two[0],&precomp_setup.E0,&basis_two[0]);
    xtoxy(&basisxy_two[1],&precomp_setup.E0,&basis_two[1]);

    xyADD(&basisxy_two[2], &precomp_setup.E0, &basisxy_two[0], &basisxy_two[1]);

    proj pt;
    xytox(&pt, &basisxy_two[2]);

    assert(mont_equal(&pt, &basis_two[2]));




    bool correct_sign;
    compute_action_2e(&m_dist_two, mont0_dist, &basis_two[0], &basis_two[1], &basis_two[2], &precomp_setup.E0, two_tors_height);

    do {
        // printf("check_action_2e call\n");
        correct_sign = check_action_2e(m_dist_two, mont0_dist, montxy0_dist, &basisxy_two[0], &basisxy_two[1], &precomp_setup.E0, two_tors_height);
        // printf("check_action_2e done\n");
        if (!correct_sign) { m_dist_two = gmod(gneg(m_dist_two), powuu(2,two_tors_height));}
    } while (!correct_sign);


    compute_action_2e(&m_frob_two, mont0_frob, &basis_two[0], &basis_two[1], &basis_two[2], &precomp_setup.E0, two_tors_height);
    do {
        correct_sign = check_action_2e(m_frob_two, mont0_frob, montxy0_frob, &basisxy_two[0], &basisxy_two[1], &precomp_setup.E0, two_tors_height);
        if (!correct_sign) { m_frob_two = gmod(gneg(m_frob_two), powuu(2,two_tors_height));}
    } while (!correct_sign);

    gelle = powuu(2,two_tors_height);

    // printf("\tglobal_setup.action_two_2 = %s;\n", pari_2x2_matrix_code(m_dist_two));
    // action_two_3_4(m_dist_two, m_frob_two, two_tors_height);


    // output(m1);
    // output(m_dist_two);
    // output(m_frob_two);
    // output(gmod(gmul(m_frob_two,m_dist_two),gelle));
    // output(gmod(gsqr(gadd(gmul(m_dist_two,gen_2),m1)),gelle));
    // output(gmod(gmul(m_frob_two,m_frob_two),gelle));

        // dist = -(1+i)/2, with i^2 = -3

    GEN action_two_2, action_two_3, action_two_4;
        action_two_2 = gmod(gneg(m_dist_two),gelle); //(1+i)/2
        action_two_3 = gmod(gmul(action_two_2,m_frob_two),gelle); //(j-ji)/2
        action_two_4 = gmod(gadd(action_two_2,gmul(action_two_3,stoi(c))),gelle); //(1+i+cj-cji)/2
        action_two_4 = gmod(gmul(action_two_4,gen_2),gelle); // (1+i+cj-cji)
        action_two_4 = gmod(gsub(gsub(action_two_4,gmul(m_frob_two,stoi(c))),m1),gelle); // (i-ji)
        action_two_4 = gmod(gmul(action_two_4,Fp_inv(stoi(q), gelle)),gelle); //(i-ji)/q


























    GEN m_frob, m_dist;


    GEN action_2[on_curve_len],action_3[on_curve_len],action_4[on_curve_len];


    // const long *on_curve_fact, *on_curve_mult;
    // long on_curve_len;

    // if (curve_order_is_p_plus_one) {
    //     on_curve_fact = p_plus_fact;
    //     on_curve_mult = p_plus_mult;
    //     on_curve_len = p_plus_len;
    // }
    // else {
    //     on_curve_fact = p_minus_fact;
    //     on_curve_mult = p_minus_mult;
    //     on_curve_len = p_minus_len;
    // }













    proj basis_twist[on_twist_len][3];
    point basis_twist_ted[on_twist_len][3];

    GEN action_twist_2[on_twist_len],action_twist_3[on_twist_len],action_twist_4[on_twist_len];

    for (int i = 0; i < on_twist_len; i++) {
        long ell = on_twist_fact[i], e = on_twist_mult[i];

        if (ell == q) e++; // a correction for the fact that constants.h only counts the q^(e-1)-torsion

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);
        random_basis(&basis_twist[i][0], &basis_twist[i][1], &basis_twist_ted[i][0], &basis_twist_ted[i][1], &precomp_setup.E0, ell, e, true);


        point test;
        proj test2;
        uintbig ell_big;
        gentobig(&ell_big, powuu(ell,e));
        ted_mul(&test, &basis_twist_ted[i][0], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_mul(&test, &basis_twist_ted[i][1], &E_twist, &ell_big);
        assert(ted_iszero(&test));

        ted_to_mont_point(&test2, &basis_twist_ted[i][0]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));

        ted_to_mont_point(&test2, &basis_twist_ted[i][1]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));


        ted_add(&basis_twist_ted[i][2], &E_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1]);
        ted_to_mont_point(&basis_twist[i][2], &basis_twist_ted[i][2]);

        ted_mul(&test, &basis_twist_ted[i][2], &E_twist, &ell_big);
        assert(ted_iszero(&test));
        ted_to_mont_point(&test2, &basis_twist_ted[i][2]);
        assert(!is_on_curve(&test2, &precomp_setup.E0));
        xMUL(&test2, &precomp_setup.E0, &test2, &ell_big);
        assert(mont_iszero(&test2));
        

        fprintf(stderr, "\\\\ computing m_dist %ld^%ld\n",ell,e);
        compute_action(&m_dist, ted0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        // printf("dist check\n");
        check_action(m_dist, ted0_dist, mont0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e);
        // printf("dist ok\n");

        fprintf(stderr, "\\\\ computing m_frob %ld^%ld\n",ell,e);
        compute_action(&m_frob, ted0_frob_twist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, ell, e);
        // printf("frob check\n");
        check_action(m_frob,       ted0_frob_twist, mont0_frob, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e);
        // printf("frob ok\n");



        action_twist_2[i] = gmod(gneg(m_dist),gelle); //(1+i)/2
        action_twist_3[i] = gmod(gmul(action_twist_2[i],m_frob),gelle); //(j-ji)/2
        action_twist_4[i] = gmod(gadd(action_twist_2[i],gmul(action_twist_3[i],stoi(c))),gelle); //(1+i+cj-cji)/2
        action_twist_4[i] = gmod(gmul(action_twist_4[i],gen_2),gelle); // (1+i+cj-cji)
        action_twist_4[i] = gmod(gsub(gsub(action_twist_4[i],gmul(m_frob,stoi(c))),m1),gelle); // (i-cji)
        if (ell != q) {
            action_twist_4[i] = gmod(gmul(action_twist_4[i], Fp_inv(stoi(q), gelle) ),gelle); //(i-cji)/q
        }
        else { // only compute the action on ell^(e-1)
            action_twist_2[i] = gmod(action_twist_2[i],powuu(ell,e-1));
            action_twist_3[i] = gmod(action_twist_3[i],powuu(ell,e-1));
            //output(gmod(action_twist_4[i],stoi(ell)));
            assert(gisexactzero(gmod(action_twist_4[i],stoi(ell))));
            action_twist_4[i] = gdiv(action_twist_4[i],stoi(ell));
            action_twist_4[i] = gmod(action_twist_4[i],powuu(ell,e-1));

            uintbig ell_big;
            uintbig_set(&ell_big,ell);
            ted_mul(&basis_twist_ted[i][0], &basis_twist_ted[i][0], &E_twist, &ell_big);
            ted_mul(&basis_twist_ted[i][1], &basis_twist_ted[i][1], &E_twist, &ell_big);
            ted_mul(&basis_twist_ted[i][2], &basis_twist_ted[i][2], &E_twist, &ell_big);
            xMUL(&basis_twist[i][0], &precomp_setup.E0, &basis_twist[i][0], &ell_big);
            xMUL(&basis_twist[i][1], &precomp_setup.E0, &basis_twist[i][1], &ell_big);
            xMUL(&basis_twist[i][2], &precomp_setup.E0, &basis_twist[i][2], &ell_big);

            check_action(gmod(m_frob,powuu(ell,e-1)), ted0_frob_twist, mont0_frob, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e-1);
            check_action(gmod(m_dist,powuu(ell,e-1)), ted0_dist, mont0_dist, &basis_twist_ted[i][0], &basis_twist_ted[i][1], &E_twist, &precomp_setup.E0, ell, e-1);

        }

    }









    for (int i = 0; i < on_curve_len; i++) {
        long ell = on_curve_fact[i], e = on_curve_mult[i];

        gelle = powuu(ell,e);
        inv2 = Fp_inv(gen_2, gelle);
        if (ell != q) invq = Fp_inv(stoi(q), gelle);









        random_basis(&basis[i][0], &basis[i][1], &basis_ted[i][0], &basis_ted[i][1], &precomp_setup.E0, ell, e, false);
        ted_add(&basis_ted[i][2], &E, &basis_ted[i][0], &basis_ted[i][1]);
        ted_to_mont_point(&basis[i][2], &basis_ted[i][2]);


        point test;
        proj test2, test3;

        // TEST that mont_add works properly
        mont_add(&test2, &precomp_setup.E0, &basis[i][0], &basis[i][1]);
        ted_neg(&test, &basis_ted[i][1]);
        ted_add(&test, &E, &basis_ted[i][0], &test);
        ted_to_mont_point(&test3, &test);
        assert(mont_equal(&test2,&basis[i][2]) || mont_equal(&test2,&test3));
        // END TEST



        fprintf(stderr, "\\\\ computing m_frob %ld^%ld\n",ell,e);

        compute_action(&m_frob, ted0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);



        check_action(m_frob, ted0_frob, mont0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e);



        fprintf(stderr, "\\\\ computing m_dist %ld^%ld\n",ell,e);
        compute_action(&m_dist, ted0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, ell, e);
        check_action(m_dist, ted0_dist, mont0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e);


        action_2[i] = gmod(gneg(m_dist),gelle); //(1+i)/2
        action_3[i] = gmod(gmul(action_2[i],m_frob),gelle); //(j-ji)/2
        action_4[i] = gmod(gadd(action_2[i],gmul(action_3[i],stoi(c))),gelle); //(1+i+cj-cji)/2
        action_4[i] = gmod(gmul(action_4[i],gen_2),gelle); // (1+i+cj-cji)
        action_4[i] = gmod(gsub(gsub(action_4[i],gmul(m_frob,stoi(c))),m1),gelle); // (i-ji)
        if (ell != q)
            action_4[i] = gmod(gmul(action_4[i],invq),gelle); //(i-ji)/q
        else { // only compute the action on ell^(e-1)
            action_2[i] = gmod(action_2[i],powuu(ell,e-1));
            action_3[i] = gmod(action_3[i],powuu(ell,e-1));
            assert(gisexactzero(gmod(action_4[i],stoi(ell))));
            action_4[i] = gdiv(action_4[i],stoi(ell));
            action_4[i] = gmod(action_4[i],powuu(ell,e-1));

            uintbig ell_big;
            uintbig_set(&ell_big,ell);
            ted_mul(&basis_ted[i][0], &basis_ted[i][0], &E, &ell_big);
            ted_mul(&basis_ted[i][1], &basis_ted[i][1], &E, &ell_big);
            ted_mul(&basis_ted[i][2], &basis_ted[i][2], &E, &ell_big);
            xMUL(&basis[i][0], &precomp_setup.E0, &basis[i][0], &ell_big);
            xMUL(&basis[i][1], &precomp_setup.E0, &basis[i][1], &ell_big);
            xMUL(&basis[i][2], &precomp_setup.E0, &basis[i][2], &ell_big);


            check_action(gmod(m_frob,powuu(ell,e-1)), ted0_frob, mont0_frob, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e-1);
            check_action(gmod(m_dist,powuu(ell,e-1)), ted0_dist, mont0_dist, &basis_ted[i][0], &basis_ted[i][1], &E, &precomp_setup.E0, ell, e-1);

        }

    }
























    point basis_ted_sum[3], basis_twist_ted_sum[3];
    proj basis_sum[3], basis_twist_sum[3];

    ted_add(&basis_ted_sum[0], &E, &basis_ted[0][0], &basis_ted[1][0]);
    ted_add(&basis_ted_sum[1], &E, &basis_ted[0][1], &basis_ted[1][1]);
    for (int i = 2; i < on_curve_len; i++) {
        ted_add(&basis_ted_sum[0], &E, &basis_ted_sum[0], &basis_ted[i][0]);
        ted_add(&basis_ted_sum[1], &E, &basis_ted_sum[1], &basis_ted[i][1]);
    }
    ted_add(&basis_ted_sum[2], &E, &basis_ted_sum[0], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[0], &basis_ted_sum[0]);
    ted_to_mont_point(&basis_sum[1], &basis_ted_sum[1]);
    ted_to_mont_point(&basis_sum[2], &basis_ted_sum[2]);

    ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted[0][0], &basis_twist_ted[1][0]);
    ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted[0][1], &basis_twist_ted[1][1]);
    for (int i = 2; i < on_twist_len; i++) {
        ted_add(&basis_twist_ted_sum[0], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted[i][0]);
        ted_add(&basis_twist_ted_sum[1], &E_twist, &basis_twist_ted_sum[1], &basis_twist_ted[i][1]);
    }
    ted_add(&basis_twist_ted_sum[2], &E_twist, &basis_twist_ted_sum[0], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[0], &basis_twist_ted_sum[0]);
    ted_to_mont_point(&basis_twist_sum[1], &basis_twist_ted_sum[1]);
    ted_to_mont_point(&basis_twist_sum[2], &basis_twist_ted_sum[2]);











    // printing stuff

    /*
   printf("#include <assert.h>\n");
   printf("#include \"precomputed.h\"\n");
   printf("#include \"curve.h\"\n");
   printf("#include \"uintbig.h\"\n\n");
    */

    printf("// each basis entry is a triple of the form P,Q,P+Q\n");
    printf("// this is initializing the point using the classical representation {0,...,p-1} for elements in GF(p).\n");
    printf("// We don't use this representation for actual computation but rather the montgomery representation (the conversion is made in init_precomputations using the fp_enc function)\n");
    printf("// hence the ***_uintbig[] defined below should not be used in any actual piece of code.\n\n");

    printf("const uintbig torsion_basis_sum_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_sum[0]), proj_code(&basis_sum[1]), proj_code(&basis_sum[2]));

    printf("const uintbig torsion_basis_twist_sum_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_twist_sum[0]), proj_code(&basis_twist_sum[1]), proj_code(&basis_twist_sum[2]));

    printf("\n");


    printf("const uintbig torsion_basis_ted_sum_uintbig[3][4][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", ted_code(&basis_ted_sum[0]), ted_code(&basis_ted_sum[1]), ted_code(&basis_ted_sum[2]));

    printf("const uintbig torsion_basis_twist_ted_sum_uintbig[3][4][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", ted_code(&basis_twist_ted_sum[0]), ted_code(&basis_twist_ted_sum[1]), ted_code(&basis_twist_ted_sum[2]));

    printf("\n");


    printf("const uintbig torsion_basis_uintbig[%ld][3][2][2] = {\n", on_curve_len);
    for (int i = 0; i < on_curve_len; i++) {
        printf("{ %s,\n %s,\n %s },\n", proj_code(&basis[i][0]), proj_code(&basis[i][1]), proj_code(&basis[i][2]));
    }
    printf("};\n\n");


    printf("const uintbig torsion_basis_twist_uintbig[%ld][3][2][2] = {\n", on_twist_len);
    for (int i = 0; i < on_twist_len; i++) {
        printf("{ %s,\n %s,\n %s },\n", proj_code(&basis_twist[i][0]), proj_code(&basis_twist[i][1]), proj_code(&basis_twist[i][2]));
    }
    printf("};\n\n");


    printf("const uintbig torsion_basis_two_uintbig[3][2][2] = \n");
    printf("{ %s,\n %s,\n %s };\n", proj_code(&basis_two[0]), proj_code(&basis_two[1]), proj_code(&basis_two[2]));




    printf("\n\n");
    printf("proj torsion_basis[%ld][3];\n", on_curve_len);
    printf("proj torsion_basis_sum[3];\n");
    printf("point torsion_basis_ted_sum[3];\n");
    printf("proj torsion_basis_twist[%ld][3];\n", on_twist_len);
    printf("proj torsion_basis_twist_sum[3];\n");
    printf("point torsion_basis_twist_ted_sum[3];\n");
    printf("proj torsion_basis_two[3];\n");



    printf("\n\n");
    printf("static void init_precomputations_generated() {\n");


    printf("\tglobal_setup.action_2 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_3 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_4 = malloc(%ld*sizeof(GEN));\n", on_curve_len);
    printf("\tglobal_setup.action_twist_2 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\tglobal_setup.action_twist_3 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\tglobal_setup.action_twist_4 = malloc(%ld*sizeof(GEN));\n", on_twist_len);
    printf("\n");


    for (int i = 0; i < on_curve_len; i++) {
        printf("\tglobal_setup.action_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_2[i]));
        printf("\tglobal_setup.action_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_3[i]));
        printf("\tglobal_setup.action_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_4[i]));
    }
    for (int i = 0; i < on_twist_len; i++) {
        printf("\tglobal_setup.action_twist_2[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_2[i]));
        printf("\tglobal_setup.action_twist_3[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_3[i]));
        printf("\tglobal_setup.action_twist_4[%d] = %s;\n", i, pari_2x2_matrix_code(action_twist_4[i]));
    }

    printf("\tglobal_setup.action_two_2 = %s;\n", pari_2x2_matrix_code(action_two_2));
    printf("\tglobal_setup.action_two_3 = %s;\n", pari_2x2_matrix_code(action_two_3));
    printf("\tglobal_setup.action_two_4 = %s;\n", pari_2x2_matrix_code(action_two_4));
    // action_two_3_4(m_dist_two, m_frob_two, two_tors_height);





    print_loop("\t","torsion_basis_two","torsion_basis_two_uintbig");
    print_loop("\t","torsion_basis_sum","torsion_basis_sum_uintbig");
    print_loop("\t","torsion_basis_twist_sum","torsion_basis_twist_sum_uintbig");
    printf("\tfor (int j=0;j<%ld;j++){\n", on_curve_len);
    print_loop("\t\t","torsion_basis[j]","torsion_basis_uintbig[j]");
    printf("\t}\n");
    printf("\tfor (int j=0;j<%ld;j++){\n", on_twist_len);
    print_loop("\t\t","torsion_basis_twist[j]","torsion_basis_twist_uintbig[j]");
    printf("\t}\n");

  printf("\
    \tfor (int i=0;i<3;i++){\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->x)->re, &(torsion_basis_ted_sum_uintbig[i][0][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->x)->im, &(torsion_basis_ted_sum_uintbig[i][0][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->y)->re, &(torsion_basis_ted_sum_uintbig[i][1][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->y)->im, &(torsion_basis_ted_sum_uintbig[i][1][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->z)->re, &(torsion_basis_ted_sum_uintbig[i][2][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->z)->im, &(torsion_basis_ted_sum_uintbig[i][2][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->t)->re, &(torsion_basis_ted_sum_uintbig[i][3][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_ted_sum[i])->t)->im, &(torsion_basis_ted_sum_uintbig[i][3][1]) );\n\
  \t}\n\
  \tfor (int i=0;i<3;i++){\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->re, &(torsion_basis_twist_ted_sum_uintbig[i][0][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->im, &(torsion_basis_twist_ted_sum_uintbig[i][0][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->re, &(torsion_basis_twist_ted_sum_uintbig[i][1][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->im, &(torsion_basis_twist_ted_sum_uintbig[i][1][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->re, &(torsion_basis_twist_ted_sum_uintbig[i][2][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->im, &(torsion_basis_twist_ted_sum_uintbig[i][2][1]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->re, &(torsion_basis_twist_ted_sum_uintbig[i][3][0]) );\n\
    \t\tfp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->im, &(torsion_basis_twist_ted_sum_uintbig[i][3][1]) );\n\
  \t}\n\n");




  //printf("/***** INSERT HERE INITIALISATION OF PARI VALUES *****/\n\n");





    printf("}\n\n");



    return 0;
}



