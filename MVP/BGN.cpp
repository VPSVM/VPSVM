#include "BGN.h"

void BGN_Gen(BGN_PK &pk, BGN_SK &sk) {
    core_init();
    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(EP_MTYPE);
    bn_null(pk.q);
    bn_new(pk.q);

    g1_get_ord(pk.q);
    g1_get_gen(pk.g);
    g2_get_gen(pk.h);

    bn_new(sk.s);
    bn_rand_mod(sk.s, pk.q);

    g1_mul_gen(pk.gs, sk.s);
    g1_rand(pk.g1);
    g1_rand(pk.g2);

    g1_t term1, term2;
    g1_mul(term1, pk.g1, sk.s);
    g1_neg(term2, pk.g2);
    g1_add(sk.PI, term1, term2);
}

void BGN_Enc(BGN_CT &ct, BGN_PK pk, bn_t m) {
    bn_t r;
    bn_new(r);
    bn_rand_mod(r, pk.q);

    g1_t term1, term2;
    g1_mul(term1, pk.g1, m);
    g1_mul(term2, pk.g, r);
    g1_add(ct.c1, term1, term2);

    g1_mul(term1, pk.g2, m);
    g1_mul(term2, pk.gs, r);
    g1_add(ct.c2, term1, term2);
}

void BGN_Dec(bn_t &m, BGN_SK sk, BGN_CT ct) {
    g1_t PI_C, tmp, t;
    g1_mul(PI_C, ct.c1, sk.s);
    g1_neg(tmp, ct.c2);
    g1_add(PI_C, PI_C, tmp);

    g1_copy(t, sk.PI);
    for (auto i = 0; i < 2147483647; i++) {
        if (g1_cmp(t, PI_C) == RLC_EQ) {
            bn_new(m);
            bn_set_dig(m, i + 1);
            break;
        }
        g1_add(t, t, sk.PI);
        g1_norm(t, t);
    }
}

// a^x = b
// Pollard's rho algorithm
void BGN_Dlog(bn_t &res, BGN_PK pk, g1_t base, g1_t power){
    bn_t a, b, i, A, B;
    bn_new(a);
    bn_new(b);
    bn_new(i);
    bn_new(A);
    bn_new(B);
    g1_t x, X;
    g1_new(x);
    g1_new(X);
    g1_sub(x, x, x);
    g1_sub(X, X, X);

    while (bn_cmp(i, pk.q) == RLC_LT) {
        BGN_xab(pk, x, a, b, base, power);
        BGN_xab(pk, X, A, B, base, power);
        BGN_xab(pk, X, A, B, base, power);

        std::cout << "index: ";
        bn_print(i);
        bn_print(a);
        bn_print(b);
        bn_print(A);
        bn_print(B);
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        if (g1_cmp(x, X) == RLC_EQ) {
            break;
        }
        bn_add_dig(i, i, 1);
    }


}

void BGN_xab(BGN_PK pk, g1_t &x, bn_t &a, bn_t &b, g1_t base, g1_t power){
    switch (BGN_partition(x)) {
        case 0:
            g1_add(x, x, base);
            bn_add_dig(a, a, 1);
            bn_mod(a, a, pk.q);
            break;
        case 1:
            g1_add(x, x, power);
            bn_add_dig(b, b, 1);
            bn_mod(b, b, pk.q);
            break;
        case 2:
            g1_add(x, x, x);
            bn_mul_dig(a, a, 2);
            bn_mod(a, a, pk.q);
            bn_mul_dig(b, b, 2);
            bn_mod(b, b, pk.q);
            break;
    }
}

unsigned long BGN_partition(g1_t x){
    size_t size = g1_size_bin(x, 0);
    bn_t tmp;
    bn_new(tmp);
    auto *t = new uint8_t[size];
    g1_write_bin(t, size, x, 0);

    bn_read_bin(tmp, t, size);

    auto *s = new dig_t;
    bn_mod_dig(s, tmp, 3);

    return *s;
}


void Sim_Dec(bn_t &res, bn_t target){
    bn_t bound;
    bn_new(bound);
    bn_set_2b(bound, 256);
    if (bn_cmp(target, bound) == RLC_GT){
        bn_copy(res, target);
        return;
    }

    // binary search
    bn_t low, high, mid;
    bn_new(low);
    bn_new(high);
    bn_new(mid);
    bn_zero(low);
    bn_copy(high, bound);
    bn_zero(res);
    while (bn_cmp(low, high) == RLC_LT) {
        bn_add(mid, low, high);
        bn_hlv(mid, mid);
//        std::cout << "mid: ";
//        bn_print(mid);
        if (bn_cmp(mid, target) == RLC_EQ) {
            bn_copy(res, mid);
            break;
        } else if (bn_cmp(mid, target) == RLC_LT) {
            bn_copy(low, mid);
        } else {
            bn_copy(high, mid);
        }
    }
}

//void GenTable(){
//    bn_t tmp, bound;
//    bn_new(tmp);

//
//
//    // write to file
//    std::ofstream table("table.txt");
//    for (bn_zero(tmp); bn_cmp(tmp, bound) == RLC_LT; bn_add_dig(tmp, tmp, 1)){
//        auto bn_str = new char[bn_size_str(tmp, 16)];
//        bn_write_str(bn_str, bn_size_str(tmp, 16),tmp, 16);
//        table << bn_str << std::endl;
//    }
//}

