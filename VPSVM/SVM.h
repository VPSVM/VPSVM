#pragma once

#include "BKS.h"

typedef struct {
    vec_ZZ_p alpha; // (y_j alpha_j) j \in SV
    Vec<vec_ZZ_p> X; // support vectors
    ZZ_p b; // bias

    // kernel
    ZZ_p gamma;
    ZZ_p c;
    ZZ p;

    vec_ZZ_pX hat_X; // encoded support vectors
    ZZ_pX hat_alpha; // encoded alpha
    ZZ_pX hat_b; // encoded bias
    ZZ_pX hat_gamma; // encoded gamma
    ZZ_pX hat_c; // encoded c

    int m; // num of support vectors
    int n; // num of features
    int n_;
    ZZ_p scale; // scale factor
}ModelPara;

typedef ZZ_p VK;

typedef struct {
    BKS_EK bksEk;
    Ciphertext C_delta;
}EK;

typedef struct {
    PKE_Para pkePara;
    int m; // num of support vectors
    int n; // num of features
}PubPara;

vec_ZZ_p HomMVMult(int b, const EK &ek, const PKE_Para &pkePara, const Vec<MemoryV> &t_M, const Ciphertext &C_v, int m, int n);

void SVM_Gen(PubPara &para, PKE_PK &pk, EK &ek1, EK &ek2, VK &pvk);

void SVM_ModelEnc(Vec<Ciphertext> &C, Ciphertext &C_d, Ciphertext &C_b, Ciphertext &C_g, Ciphertext &C_c,
                    PubPara &para, const PKE_PK &pk, const ModelPara &modelPara);

void Compute_poly(ZZ_p &y_1, ZZ_p &y_2, ZZ_p &Phi_1, ZZ_p &Phi_2, const EK &ek1, const EK &ek2, PubPara &para,
                    const Vec<Ciphertext> &C, const Ciphertext &C_d, const Ciphertext &C_b, const Ciphertext &C_g, const Ciphertext &C_c,
                    const ZZ &p, const PKE_PK &pk, VK &vk, const vec_ZZ_p &z,
                    std::vector<double> &Time);

void Compute_rbf(ZZ_p& y_1, ZZ_p& y_2, ZZ_p &Phi_1, ZZ_p &Phi_2, const EK &ek1, const EK &ek2, PubPara &para,
                    const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b, const Ciphertext & C_g, // server inputs
                    const PKE_PK& pk, VK& vk, const vec_ZZ_p& z,
                    std::vector<double> &Time);


bool Verify(ZZ_p &y, const ZZ_p& y_1, const ZZ_p& y_2, const ZZ_p &phi_1, const ZZ_p &phi_2, VK& vk, PubPara &para);