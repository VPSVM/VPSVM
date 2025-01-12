// Secure Dot Product Computation

#include "BGN.h"

typedef struct {
    BGN_PK bgn_pk;
    g1_t* g1t;
    g2_t* ht;
    g1_t gt2;
    int features;
}SDPC_PK;

typedef struct {
    BGN_SK bgn_sk;
    bn_t* t;
}SDPC_SK;

typedef g1_t SDPC_FK;

void SDPC_KeyGen(SDPC_PK &pk, SDPC_SK &sk, int features);

void SDPC_Setup(SDPC_FK &fk, SDPC_PK pk, bn_t* x);

void SDPC_Enc(BGN_CT* ct, SDPC_PK pk, bn_t* z);

void SDPC_Compute(BGN_CT &V, g2_t &sigma0, g1_t *sigma, SDPC_PK pk, BGN_CT *ct, bn_t *x);

void SDPC_Dec(bn_t &v, BGN_CT V, SDPC_SK sk);

bool SDPC_Verify(SDPC_PK pk, SDPC_FK fk, bn_t v, bn_t* z, g2_t sigma0, g1_t *sigma);