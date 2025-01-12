#include "SDPC.h"

void SDPC_KeyGen(SDPC_PK &pk, SDPC_SK &sk, int features) {
//    std::cout << "KeyGen" << std::endl;
    BGN_Gen(pk.bgn_pk, sk.bgn_sk);
    pk.features = features;

    sk.t = new bn_t[features];
    pk.g1t = new g1_t[features];
    pk.ht = new g2_t[features];

    bn_t sum;
    bn_new(sum);
    for (auto i = 0; i < features; i++){
        bn_rand_mod(sk.t[i], pk.bgn_pk.q);
        g1_mul(pk.g1t[i], pk.bgn_pk.g1, sk.t[i]);
        g2_mul(pk.ht[i], pk.bgn_pk.h, sk.t[i]);
        bn_t tmp;
        bn_new(tmp);
        bn_mul(tmp, sk.t[i], sk.t[i]);
        bn_add(sum, sum, tmp);
    }
    bn_mod(sum, sum, pk.bgn_pk.q);
    g1_mul(pk.gt2, pk.bgn_pk.g1, sum);
}

void SDPC_Setup(SDPC_FK &fk, SDPC_PK pk, bn_t *x) {
//    std::cout << "Setup" << std::endl;
    g1_sub(fk, fk, fk);
    for (auto i = 0; i < pk.features; i++){
        g1_t tmp;
        g1_new(tmp);
        g1_mul(tmp, pk.g1t[i], x[i]);
        g1_add(fk, fk, tmp);
    }
}

void SDPC_Enc(BGN_CT *ct, SDPC_PK pk, bn_t *z) {
//    std::cout << "Enc" << std::endl;
    for (auto i = 0; i < pk.features; i++){
//        if (bn_is_zero(z[i])) {
//            ct[i].is_zero = true;
//            continue;
//        }
//        ct[i].is_zero = false;
        BGN_Enc(ct[i], pk.bgn_pk, z[i]);
    }
}



void SDPC_Compute(BGN_CT &V, g2_t &sigma0, g1_t *sigma, SDPC_PK pk, BGN_CT *ct, bn_t *x) {
//    std::cout << "Compute" << std::endl;
    g1_sub(V.c1, V.c1, V.c1);
    g1_sub(V.c2, V.c2, V.c2);
    bn_t d;
    bn_new(d);
    bn_rand_mod(d, pk.bgn_pk.q);
    g2_mul(sigma0, pk.bgn_pk.h, d);
    g2_norm(sigma0, sigma0);
    for (auto i = 0; i < pk.features; i++) {
//        if (bn_is_zero(x[i])) {
//            continue;
//        }
        g1_t t1, t2, t3, t4;
        g1_new(t1);
        g1_new(t2);
        g1_mul(t1, ct[i].c1, x[i]);
        g1_add(V.c1, V.c1, t1);
        g1_mul(t2, ct[i].c2, x[i]);
        g1_add(V.c2, V.c2, t2);
        bn_t tmp;
        bn_mul(tmp, d, x[i]);
        g1_mul(t3, pk.bgn_pk.g1, tmp);
        g1_norm(t3, t3);
        g1_mul(t4, pk.g1t[i], d);
        g1_norm(t4, t4);
        g1_add(sigma[i], t3, t4);
        g1_norm(sigma[i], sigma[i]);
    }
}

void SDPC_Dec(bn_t &v, BGN_CT V, SDPC_SK sk) {
//    std::cout << "Dec" << std::endl;
    BGN_Dec(v, sk.bgn_sk, V);
}

bool SDPC_Verify(SDPC_PK pk, SDPC_FK fk, bn_t v, bn_t *z, g2_t sigma0, g1_t *sigma) {
    gt_t left_side, right_side;
    gt_new(left_side);gt_new(right_side);

    g1_t Hf, g1_ls;
    g1_new(Hf);
    g1_sub(Hf, Hf, Hf);
    for (auto i = 0; i < pk.features; i++){
        g1_t tmp;
        g1_new(tmp);
        g1_mul(tmp, pk.g1t[i], z[i]);
        g1_add(Hf, Hf, tmp);
    }
    g1_norm(Hf, Hf);
    g1_sub(Hf, pk.gt2, Hf);
    g1_norm(Hf, Hf);
    g1_mul(g1_ls, pk.bgn_pk.g1, v);
    g1_sub(g1_ls, Hf, g1_ls);
    g1_add(g1_ls, g1_ls, fk);
    g1_norm(g1_ls, g1_ls);
//    g1_print(g1_ls);

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
//    std::cout<< "Left side:" << std::endl;
//    gt_print(left_side);
//    std::cout<< "Right side:" << std::endl;
//    gt_print(right_side);
    return false;
}



