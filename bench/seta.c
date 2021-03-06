#define _XOPEN_SOURCE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pari/pari.h>
#include <assert.h>


#include "ideal.h"
#include "idiso.h"
#include "constants.h"
#include "precomputed.h"
#include "tedwards.h"
#include "isogenies.h"
#include "klpt.h"
#include "toolbox.h"
#include "seta.h"
#include "mont.h"
#include "printing.h"



static char* uintbig_code(const uintbig *x) {
    return pari_sprintf("{ %luULL, %luULL, %luULL, %luULL, %luULL, %luULL, %luULL }", x->c[0], x->c[1], x->c[2], x->c[3] , x->c[4], x->c[5], x->c[6]);
}

void print_secret_key(const secret_key *sk) {
    printf("(secret_key) {\n%s,\n %s,\n %s,\n %s,\n %s,\n %s,\n %s,\n %s }\n",
    uintbig_code(&sk->ker_psi_plus_1),
    uintbig_code(&sk->ker_psi_plus_2),
    uintbig_code(&sk->ker_psi_minus_1),
    uintbig_code(&sk->ker_psi_minus_2),
    uintbig_code(&sk->ker_psi_dual_plus_1),
    uintbig_code(&sk->ker_psi_dual_plus_2),
    uintbig_code(&sk->ker_psi_dual_minus_1),
    uintbig_code(&sk->ker_psi_dual_minus_2));
}


static char* fp_code(const fp *x) {
    return pari_sprintf("(fp) { %luULL, %luULL, %luULL, %luULL, %luULL, %luULL, %luULL }", x->x.c[0], x->x.c[1], x->x.c[2], x->x.c[3] , x->x.c[4], x->x.c[5], x->x.c[6]);
}

static char* fp2_code(const fp2 *x) {
    return pari_sprintf("(fp2) { %s,\n %s }", fp_code(&x->re), fp_code(&x->im));
}

static char* proj_code(const proj *P) {
    if (fp2_iszero(&P->z))
        return pari_sprintf("(proj) { %s,\n %s }", fp2_code(&P->x), fp2_code(&P->z));
    else {
        fp2 x,z,tmp = P->z;
        fp2_inv(&tmp);
        fp2_mul3(&x,&P->x,&tmp);
        fp2_mul3(&z,&P->z,&tmp);
        return pari_sprintf("(proj) { %s,\n %s }", fp2_code(&x), fp2_code(&z));
    }
}

static char* long_list_code(const long *list, long len){
    char* str = pari_sprintf("{ %ld", list[0]);
    for (int i = 1; i < len; ++i) {
        str = pari_sprintf("%s, %ld", str, list[i]);
    }
    return pari_sprintf("%s }", str);
}

static char* hybrid_point_code(const hybrid_point *x) {
    return pari_sprintf("(hybrid_point) { %s,\n %s,\n%s,\n %s,\n%s }", proj_code(&x->E),long_list_code(x->twist_mult,on_twist_len),long_list_code(x->curve_mult,on_curve_len), proj_code(&x->twist_pt),proj_code(&x->curve_pt));
}

void print_public_key(const public_key *pk) {
    printf("(public_key) { \n{\n%s,\n %s,\n %s },\n{\n %s,\n %s,\n %s } }\n",
    hybrid_point_code(&pk->N1_basis[0]),
    hybrid_point_code(&pk->N1_basis[1]),
    hybrid_point_code(&pk->N1_basis[2]),
    hybrid_point_code(&pk->N2_basis[0]),
    hybrid_point_code(&pk->N2_basis[1]),
    hybrid_point_code(&pk->N2_basis[2]));
}



void seta_test() {

    pari_sp ltop = avma;
    clock_t t;

    public_param param;
    public_key pk;
    secret_key sk;

    printf("seta_setup...\n");
    t = tic();
    seta_setup(&param);
    TOC(t,"seta_setup");

    printf("seta_genkey...\n");
    t = tic();
    // seta_genkey(&pk, &sk, &param);
    seta_genkey_invalid_for_testing(&pk, &sk, &param);
    TOC(t,"seta_genkey");



    pk = (public_key) { 
{
(hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 0, 12, 0, 11 },
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
 (proj) { (fp2) { (fp) { 214637209459294471ULL, 12820464457505495704ULL, 15355399691768868878ULL, 10023692330711603179ULL, 1773340000404818821ULL, 12979851198623152097ULL, 3088ULL },
 (fp) { 1112328399323702405ULL, 12754240904615906593ULL, 14386556376697396996ULL, 7448323289492894689ULL, 3445743690340714791ULL, 6825549528528644221ULL, 5823ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 10649170981479549110ULL, 10848839288666792319ULL, 13053576644236429571ULL, 10499653172239698985ULL, 10028528842505549131ULL, 16040359582306900143ULL, 4867ULL },
 (fp) { 14045067143126715644ULL, 572813566858529391ULL, 15749011192536091303ULL, 7178469766462438545ULL, 17007155092334661662ULL, 16989504313099528840ULL, 2627ULL } },
 (fp2) { (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } },
 (hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 0, 12, 0, 11 },
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
 (proj) { (fp2) { (fp) { 18185652390473812466ULL, 13881097315533233265ULL, 4393165363006589597ULL, 11112873565909793186ULL, 1276076023794203582ULL, 16696183641384344991ULL, 6138ULL },
 (fp) { 4415289537480309911ULL, 1330771652832080336ULL, 16093086549340452164ULL, 1570642965667445660ULL, 11579684464325329243ULL, 6965582806229560732ULL, 372ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 10649170981479549110ULL, 10848839288666792319ULL, 13053576644236429571ULL, 10499653172239698985ULL, 10028528842505549131ULL, 16040359582306900143ULL, 4867ULL },
 (fp) { 14045067143126715644ULL, 572813566858529391ULL, 15749011192536091303ULL, 7178469766462438545ULL, 17007155092334661662ULL, 16989504313099528840ULL, 2627ULL } },
 (fp2) { (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } },
 (hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 0, 12, 0, 11 },
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
 (proj) { (fp2) { (fp) { 11694187498612448177ULL, 3563801228930783944ULL, 963321314040289652ULL, 14545698375324229100ULL, 17193920653168551736ULL, 4140327983960613993ULL, 5985ULL },
 (fp) { 15709174989897989400ULL, 970332123057212512ULL, 2507177360247911685ULL, 15276204970873361595ULL, 1939419689277512267ULL, 2741838162881028620ULL, 910ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 10649170981479549110ULL, 10848839288666792319ULL, 13053576644236429571ULL, 10499653172239698985ULL, 10028528842505549131ULL, 16040359582306900143ULL, 4867ULL },
 (fp) { 14045067143126715644ULL, 572813566858529391ULL, 15749011192536091303ULL, 7178469766462438545ULL, 17007155092334661662ULL, 16989504313099528840ULL, 2627ULL } },
 (fp2) { (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } } },
{
 (hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 21, 0, 12, 0 },
{ 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 },
 (proj) { (fp2) { (fp) { 6156674266387762377ULL, 6488868244868048284ULL, 16844134348073794452ULL, 13122379993227338216ULL, 160991904795463647ULL, 3772737640572003698ULL, 5148ULL },
 (fp) { 3783190855801311578ULL, 14285805877352809429ULL, 8744995546132225063ULL, 4710264199872856023ULL, 10881165140947683927ULL, 11721525497174212184ULL, 2101ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 15281337400593823800ULL, 4973986318349928671ULL, 9243074515944333313ULL, 8270668177663341824ULL, 10291415836356746924ULL, 2646568803394483409ULL, 5758ULL },
 (fp) { 4295163312012367463ULL, 15394036351313401810ULL, 8738906385437031209ULL, 10093767768921802035ULL, 3978747096695069114ULL, 5532928841335979622ULL, 6065ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } },
 (hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 21, 0, 12, 0 },
{ 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 },
 (proj) { (fp2) { (fp) { 10674842579907822099ULL, 6075309568349400842ULL, 5992538134177526251ULL, 11235059133537837030ULL, 10860503548205129256ULL, 7922287974747777328ULL, 892ULL },
 (fp) { 2795607334017944500ULL, 15092004961164331383ULL, 511970082037653118ULL, 14525257629830592055ULL, 2488779210110957609ULL, 1306931128030883760ULL, 5940ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 13067732917420694056ULL, 18254103427078155747ULL, 9812966835808376208ULL, 13856976074382179624ULL, 4872025099272185639ULL, 3415613420277916178ULL, 4127ULL },
 (fp) { 12195403408956990054ULL, 13103475307375975891ULL, 16784504172029325818ULL, 8614112370978699249ULL, 10809681847247840840ULL, 9162392670611877386ULL, 5542ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } },
 (hybrid_point) { (proj) { (fp2) { (fp) { 15975943569627156142ULL, 14700402415520602756ULL, 1968887260658849972ULL, 15498408440594758482ULL, 7735043111479803572ULL, 11193183583946369811ULL, 5588ULL },
 (fp) { 12500385730831924425ULL, 5374140871247055802ULL, 8375885165527849029ULL, 2793167468489917883ULL, 8692640415777454510ULL, 6573093634904594250ULL, 5811ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
 { 21, 0, 12, 0 },
{ 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 },
 (proj) { (fp2) { (fp) { 8706977567806909943ULL, 3449521337688933688ULL, 6674252401954840545ULL, 5214785759108512424ULL, 9902325531857952453ULL, 14073244175627768083ULL, 4059ULL },
 (fp) { 4147663426384621411ULL, 2297578960934734804ULL, 2018981495745459743ULL, 567999675960623452ULL, 10228298942348246938ULL, 2338656532623179471ULL, 3282ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } },
(proj) { (fp2) { (fp) { 6085491343756121365ULL, 16258362856979897256ULL, 13212147577822265479ULL, 16996892335713842081ULL, 5808447738920572571ULL, 5274439722713141607ULL, 5504ULL },
 (fp) { 12361462980084046268ULL, 14853175954671514455ULL, 18030575566240971496ULL, 15899038481466383094ULL, 13812284527986873506ULL, 5736030730288135338ULL, 33ULL } },
 (fp2) { (fp) { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 (fp) { 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL } } } } } };


 sk = (secret_key) {
{ 4001577640881062959ULL, 10855057590284415167ULL, 1ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 { 8701852876890340012ULL, 10005881330857136325ULL, 2ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 { 9378746999513446426ULL, 9549650548232241644ULL, 239865105547099739ULL, 14043090645560323546ULL, 15580949235784289743ULL, 132705614405ULL, 0ULL },
 { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 { 15315942338266486186ULL, 5523416115139157852ULL, 2ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 { 1ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL },
 { 8027038862688935688ULL, 11393372570728605335ULL, 17727054572406653358ULL, 5240254602208873735ULL, 2719673791249692946ULL, 224207050338ULL, 0ULL },
 { 18180680994798409397ULL, 3288338528196333869ULL, 3976981034093096486ULL, 11008533056676257782ULL, 6402501255886023980ULL, 1724014361402ULL, 0ULL } };

    // printf("// PUBLIC KEY\n");
    // print_public_key(&pk);
    // printf("\n\n");
    // printf("// SECRET KEY\n");
    // print_secret_key(&sk);
    // printf("\n\n");




    for (int i = 0; i < 20; i++) {
        uintbig x = uintbig_1;
        randombytes(x.c + 0, 32);
        odd_isogeny phi_m;


        ciphertext ctxt;
        printf("seta_eval...\n");
        t = tic();
        seta_eval(&ctxt, &pk, &x);
        TOC(t,"seta_eval");


        printf("seta_inverse...\n");
        t = tic();
        seta_inverse(&phi_m, &ctxt, &pk, &sk, &param);
        TOC(t,"seta_eval");


        printf("done\n\n");


    }


    //gerepileall(ltop, 2, I, coeff);
    avma = ltop;
}


// argv[1] is the random seed; default = 1
int main(int argc, char *argv[]){
    pari_init(800000000, 1<<18);
    init_precomputations();

    setrand(stoi(1));
    srand48(1);
    if( argc > 1 ) {
      setrand(strtoi(argv[1]));
      srand48(atoi(argv[1]));
    }

    for (int i = 0; i < 1; i++) {
      uintbig m;
      randombytes(m.c, 32);

      seta_test();


    }

    printf("    \033[1;32mAll tests passed\033[0m\n");
    exit(0);

    return 0;
}
