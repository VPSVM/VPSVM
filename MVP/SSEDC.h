//
// Created by ly on 24-3-27.
//
#include "BGN.h"

#ifndef MVP_SSEDC_H
#define MVP_SSEDC_H

typedef struct {
    BGN_PK bgn_pk;
    g1_t* gt;
    g1_t* g1t;
    g1_t* g1t2;
    g2_t* ht;
    int features;
}SSEDC_PK;

typedef struct {
    BGN_SK bgn_sk;
    bn_t* t;
}SSEDC_SK;

typedef g1_t SSEDC_FK;

void SSEDC_KeyGen(SSEDC_PK &pk, SSEDC_SK &sk, int features);

void SSEDC_Setup(SSEDC_FK &fk, SSEDC_PK pk, bn_t *x);

void SSEDC_Enc(BGN_CT *ct1, BGN_CT *ct2, SSEDC_PK pk, bn_t *z);

void SSEDC_Compute(BGN_CT &V, g2_t &sigma0, g1_t *sigma, SSEDC_PK pk, BGN_CT *ct1, BGN_CT *ct2, bn_t *x);

void SSEDC_Dec(bn_t &v, BGN_CT V, SSEDC_SK sk);

bool SSEDC_Verify(SSEDC_PK pk, SSEDC_FK fk, bn_t v, bn_t *z, g2_t sigma0, g1_t *sigma);
#endif //MVP_SSEDC_H
