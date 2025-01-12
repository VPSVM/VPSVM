//
// Created by ly on 24-3-27.
//

#include "SSEDC.h"

void SSEDC_KeyGen(SSEDC_PK &pk, SSEDC_SK &sk, int features){
//    std::cout << "KeyGen" << std::endl;
    bn_t tmp;
    bn_new(tmp);
    BGN_Gen(pk.bgn_pk, sk.bgn_sk);
    pk.features = features;

    sk.t = new bn_t[features];
    pk.g1t = new g1_t[features];
    pk.gt = new g1_t[features];
    pk.g1t2 = new g1_t[features];
    pk.ht = new g2_t[features];

    for (auto i = 0; i < features; i++){
        bn_new(sk.t[i]);
        bn_rand_mod(sk.t[i], pk.bgn_pk.q);
        bn_mul(tmp, sk.t[i], sk.t[i]);
        g1_mul(pk.gt[i], pk.bgn_pk.g, sk.t[i]);
        g1_mul(pk.g1t[i], pk.bgn_pk.g1, sk.t[i]);
        g1_mul(pk.g1t2[i], pk.bgn_pk.g1, tmp);
        g2_mul(pk.ht[i], pk.bgn_pk.h, sk.t[i]);
    }


}

void SSEDC_Setup(SSEDC_FK &fk, SSEDC_PK pk, bn_t *x) {
//    std::cout << "Setup" << std::endl;
    g1_sub(fk, fk, fk);
    for (auto i = 0; i < pk.features; i++){
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, x[i], x[i]);

        g1_t g1x2;
        g1_new(g1x2);
        g1_mul(g1x2, pk.bgn_pk.g1, tmp);
        g1_add(fk, fk, g1x2);

        bn_mul_dig(tmp, x[i], 2);
        bn_neg(tmp, tmp);
        g1_t g1t_2x;
        g1_new(g1t_2x);
        g1_mul(g1t_2x, pk.g1t[i], tmp);
        g1_add(fk, fk, g1t_2x);

        g1_add(fk, fk, pk.g1t2[i]);
    }

}

void SSEDC_Enc(BGN_CT *ct1, BGN_CT *ct2,SSEDC_PK pk, bn_t *z) {
//    std::cout << "Enc" << std::endl;
    for (auto i = 0; i < pk.features; i++){
//        if (bn_is_zero(z[i])) {
//            ct1[i].is_zero = true;
//            continue;
//        }
//        ct1[i].is_zero = false;
//        ct2[i].is_zero = false;
        BGN_Enc(ct1[i], pk.bgn_pk, z[i]);
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, z[i], z[i]);
        BGN_Enc(ct2[i], pk.bgn_pk, tmp);
    }
}

void SSEDC_Compute(BGN_CT &V, g2_t &sigma0, g1_t *sigma, SSEDC_PK pk, BGN_CT *ct1, BGN_CT *ct2, bn_t *x){
//    std::cout << "Compute" << std::endl;
    g1_sub(V.c1, V.c1, V.c1);
    g1_sub(V.c2, V.c2, V.c2);
    bn_t d;
    bn_new(d);
    bn_rand_mod(d, pk.bgn_pk.q);
    g2_mul(sigma0, pk.bgn_pk.h, d);
    g2_norm(sigma0, sigma0);
    for (auto i = 0; i < pk.features; i++){
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, x[i], x[i]);
        g1_t g1x2, g2x2;
        g1_new(g1x2);
        g1_new(g2x2);
        g1_mul(g1x2, pk.bgn_pk.g1, tmp);
        g1_mul(g2x2, pk.bgn_pk.g2, tmp);
        g1_add(V.c1, V.c1, g1x2);
        g1_add(V.c2, V.c2, g2x2);

        bn_mul_dig(tmp, x[i], 2);
        bn_neg(tmp, tmp);
        g1_t g12xz, g22xz;
        g1_new(g12xz);
        g1_new(g22xz);

        g1_mul(g12xz, ct1[i].c1, tmp);
        g1_mul(g22xz, ct1[i].c2, tmp);

        g1_add(V.c1, V.c1, g12xz);
        g1_add(V.c2, V.c2, g22xz);

        g1_add(V.c1, V.c1, ct2[i].c1);
        g1_add(V.c2, V.c2, ct2[i].c2);

        bn_new(tmp);
        bn_mul_dig(tmp, x[i], 2);
        bn_neg(tmp, tmp);
        g1_new(g1x2);
        g1_mul(g1x2, pk.bgn_pk.g1, tmp);

        g1_add(sigma[i], pk.g1t[i], g1x2);
        g1_add(sigma[i], sigma[i], ct1[i].c1);
        g1_mul(sigma[i], sigma[i], d);
    }
//    // sigma
//    for (auto i = 0; i < pk.features; i++){
//        bn_t tmp;
//        bn_new(tmp);
//        bn_mul_dig(tmp, x[i], 2);
//        bn_neg(tmp, tmp);
//        g1_t g1x2;
//        g1_new(g1x2);
//        g1_mul(g1x2, pk.bgn_pk.g1, tmp);
//
//        g1_add(sigma[i], pk.g1t[i], g1x2);
//        g1_add(sigma[i], sigma[i], ct1[i].c1);
//        g1_mul(sigma[i], sigma[i], d);
//    }
}

void SSEDC_Dec(bn_t &v, BGN_CT V, SSEDC_SK sk) {
//    std::cout << "Dec" << std::endl;
    BGN_Dec(v, sk.bgn_sk, V);
}

bool SSEDC_Verify(SSEDC_PK pk, SSEDC_FK fk, bn_t v, bn_t *z, g2_t sigma0, g1_t *sigma) {
    gt_t left_side, right_side;
    gt_new(left_side);gt_new(right_side);

    g1_t Hf, g1_ls;
    g1_new(Hf);
    g1_sub(Hf, Hf, Hf);

    bn_t r;
    bn_new(r);
    bn_rand_mod(r, pk.bgn_pk.q);

    bn_t sum;
    bn_new(sum);
    for (auto i = 0; i < pk.features; i++){
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, r, z[i]);
        bn_add(sum, sum, tmp);
        g1_t gtr;
        g1_new(gtr);
        g1_mul(gtr, pk.gt[i], r);
        g1_add(Hf, Hf, gtr);
    }
    g1_t gsum;
    g1_new(gsum);
    g1_mul(gsum, pk.bgn_pk.g, sum);

    g1_sub(Hf, Hf, gsum);
    g1_norm(Hf, Hf);

    g1_mul(g1_ls, pk.bgn_pk.g1, v);
    g1_sub(g1_ls, Hf, g1_ls);
    g1_add(g1_ls, g1_ls, fk);
    g1_norm(g1_ls, g1_ls);

    pc_map(left_side, g1_ls, sigma0);

    gt_exp_dig(right_side, right_side, 0);
    for (auto i = 0; i < pk.features; i++){
        gt_t lsi;
        gt_new(lsi);
        g2_t g2_rs;
        g2_new(g2_rs);
        g2_mul(g2_rs, pk.bgn_pk.h, z[i]); // h^zi
        g2_sub(g2_rs, pk.ht[i], g2_rs); // h^(ti-zi)
        g2_norm(g2_rs, g2_rs);
//        std::cout<< i << std::endl;
//        g2_print(g2_rs);
        pc_map(lsi, sigma[i], g2_rs);
        gt_mul(right_side, right_side, lsi);
    }
    if (gt_cmp(left_side, right_side) == RLC_EQ){
        return true;
    }
    std::cout<< "Left side:" << std::endl;
    gt_print(left_side);
    std::cout<< "Right side:" << std::endl;
    gt_print(right_side);
    return false;
}