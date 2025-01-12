#include <numeric>
#include "SDPC.h"
#include "SSEDC.h"

void test_dp_input(int features);
void test_sed_input(int features);
void test_dp_verify(int size);
void test_sed_verify(int size);

int main(){


//    test_dp_input(128);
//    test_dp_input(128);
//    test_dp_input(256);
//    test_dp_input(512);
//    test_dp_input(1024);
//
//    test_sed_input(64);
//    test_sed_input(128);
//    test_sed_input(256);
//    test_sed_input(512);
//    test_sed_input(1024);

    test_dp_verify(100);
    test_dp_verify(1000);
    test_dp_verify(10000);
    test_dp_verify(100000);
    test_dp_verify(1000000);
//
//    test_sed_verify(100);
//    test_sed_verify(1000);
//    test_sed_verify(10000);
//    test_sed_verify(100000);
//    test_sed_verify(1000000);

//    BGN_PK pk;
//    BGN_SK sk;
//    BGN_Gen(pk, sk);
//    g1_t gg;
//    g1_new(gg);
//    g1_mul_dig(gg, pk.g1, 512);
//
//    bn_t res;
//    bn_new(res);
//    BGN_Dlog(res, pk, pk.g1, gg);

//    bn_print(pk.q);


//    bn_t m, bound, res;
//    bn_new(m);
//    bn_new(bound);
//    bn_new(res);
//    bn_set_2b(bound, 256);
//    bn_rand_mod(m, bound);
//
//    auto time = GetTime();
//    Sim_Dec(res, m);
//    time = GetTime() - time;
//    std::cout << "Sim_Dec time: " << time * 1000 << " ms" << std::endl;
//
//    bn_print(res);
//    bn_print(m);


//    std::cout << "Hello World!" << std::endl;
    return 0;
}

void test_sed_input(int features){
    SSEDC_PK pk;
    SSEDC_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t size;

    SSEDC_KeyGen(pk, sk, features);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/predict_students_dropout_and_academic_success_para_rbf.csv");


    get_user_inputs(X, y, "../data/predict_students_dropout_and_academic_success_test_data_rbf.csv", size);

    std::vector<double> user_enc_time;
    std::vector<double> server_compute_time;

    double time;
    for (auto in = 0; in < 1000; in++){
        auto *ct1 = new BGN_CT[modelPara.features];
        auto *ct2 = new BGN_CT[modelPara.features];
        time = GetTime();
        SSEDC_Enc(ct1, ct2, pk, X[in]);
        time = GetTime() - time;
        user_enc_time.push_back(time * 1000);


        time = GetTime();
        for (int i = 0; i < modelPara.SV; i++) {
            BGN_CT V;
            g2_t sigma0;
            auto *sigma = new g1_t[modelPara.features];
            g1_new(V.c1);
            g1_new(V.c2);

            SSEDC_Compute(V, sigma0, sigma, pk, ct1, ct2, modelPara.sv[i]);
        }
        time = GetTime() - time;
        server_compute_time.push_back(time * 1000);
    }

    // average time
    double mean, stdev;
    // mean of user encryption time
    mean = std::accumulate(user_enc_time.begin(), user_enc_time.end(), 0.0) / user_enc_time.size();
    // stdev of user encryption time
    stdev = std::sqrt(std::inner_product(user_enc_time.begin(), user_enc_time.end(), user_enc_time.begin(), 0.0) / user_enc_time.size() - mean * mean);
    std::cout << "Feature:" << features << std::endl;
    std::cout << "Mean of user encryption time: " << mean << "ms" << std::endl;
    std::cout << "Standard deviation of user encryption time: " << stdev << std::endl;

    // mean of server computing time
    mean = std::accumulate(server_compute_time.begin(), server_compute_time.end(), 0.0) / server_compute_time.size();
    // stdev of server computing time
    stdev = std::sqrt(std::inner_product(server_compute_time.begin(), server_compute_time.end(), server_compute_time.begin(), 0.0) / server_compute_time.size() - mean * mean);
    std::cout << "Mean of server computing time: " << mean << "ms" << std::endl;
    std::cout << "Standard deviation of server computing time: " << stdev << std::endl;
//    for (auto i = 0; i < size; i++){
//        auto *ct1 = new BGN_CT[modelPara.features];
//        auto *ct2 = new BGN_CT[modelPara.features];
//        auto time = GetTime();
//        SSEDC_Enc(ct1, ct2, pk, X[i]);
//        time = GetTime() - time;
//        enc_time << time * 1000 << std::endl;
//    }

//    for (int i = 0; i < modelPara.SV; i++){
//        SSEDC_FK fk;
//
//        auto time = GetTime();
//        SSEDC_Setup(fk, pk, modelPara.sv[i]);
//        time = GetTime() - time;
//        setup_time << time * 1000 << std::endl;
//
//        auto *ct1 = new BGN_CT[modelPara.features];
//        auto *ct2 = new BGN_CT[modelPara.features];
////        time = GetTime();
//        SSEDC_Enc(ct1, ct2, pk, X[0]);
////        time = GetTime() - time;
////        enc_time << time * 1000 << std::endl;
//
//        BGN_CT V;
//        g2_t sigma0;
//        auto *sigma = new g1_t[modelPara.features];
//        g1_new(V.c1);
//        g1_new(V.c2);
//        // Set to 0
//
//        time = GetTime();
//        SSEDC_Compute(V, sigma0, sigma, pk, ct1, ct2, modelPara.sv[i]);
//        time = GetTime() - time;
//        compute_time << time * 1000 << std::endl;
//        std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;
//
//
//
//
//    }

    /* std::ofstream verify_time("../benchmark/verify_edc_time_"+ std::to_string(features) +".txt");

    for (auto k = 0; k < size; k++){
        auto time = GetTime();
        g1_t FK;
        g1_new(FK);
        g1_sub(FK, FK, FK);

        bn_t sum_x, sum_xx, bn_SV, d;
        bn_new(sum_x);
        bn_new(sum_xx);
        bn_new(bn_SV);
        bn_new(d);
        bn_zero(sum_x);
        bn_zero(sum_xx);
        bn_set_dig(bn_SV, modelPara.SV);
        bn_rand_mod(d, pk.bgn_pk.q);
        g2_t sigma0;
        g2_new(sigma0);
        g2_mul(sigma0, pk.bgn_pk.h, d);

        auto sigma = new g1_t[modelPara.features];

        bn_t sum;
        bn_new(sum);
        bn_zero(sum);
        g1_t Hf;
        g1_new(Hf);

        for (auto i = 0; i < pk.features; i++){
            g1_new(sigma[i]);

            // Generate FK
            g1_t tmp;
            g1_new(tmp);
            g1_mul(tmp, pk.g1t2[i], bn_SV);
            g1_add(FK, FK, tmp);
            for (auto j = 0; j < modelPara.SV; j++){
                bn_t xx;
                bn_new(xx);
                bn_mul(xx, modelPara.sv[j][i], modelPara.sv[j][i]);
                bn_add(sum_xx, sum_xx, xx);
                bn_add(sum_x, sum_x, modelPara.sv[j][i]);
            }
            g1_mul(tmp, pk.g1t[i], sum_x);
            g1_mul_dig(tmp, tmp, 2);
            g1_sub(FK, FK, tmp);

            g1_mul(tmp, pk.bgn_pk.g1, sum_xx);

            g1_add(FK, FK, tmp);

            // sigma_i
            bn_t r;
            bn_new(r);
            bn_rand_mod(r, pk.bgn_pk.q);
            g1_mul(tmp, pk.bgn_pk.g, r);
            g1_add(sigma[i], pk.g1t[i], tmp);
            g1_mul(tmp, pk.bgn_pk.g1, X[k][i]);
            g1_add(sigma[i], sigma[i], tmp);
            g1_mul(sigma[i], sigma[i], d);
            g1_mul(sigma[i], sigma[i], bn_SV);
            g1_mul(tmp, pk.bgn_pk.g1, sum_x);
            g1_mul_dig(tmp, tmp, 2);
            g1_sub(sigma[i], sigma[i], tmp);

            // Hf
            bn_t xr;
            bn_new(xr);
            bn_mul(xr, r, X[k][i]);
            bn_add(sum, sum, xr);
            g1_t gtr;
            g1_new(gtr);
            g1_mul(gtr, pk.gt[i], r);
            g1_add(Hf, Hf, gtr);
        }
        g1_t gsum;
        g1_new(gsum);
        g1_mul(gsum, pk.bgn_pk.g, sum);

        g1_sub(Hf, Hf, gsum);
        g1_mul(Hf, Hf, bn_SV);

        // sum v
        auto v = new bn_t[modelPara.SV];
        bn_t sum_v;
        bn_new(sum_v);
        bn_zero(sum_v);
        for (int i = 0; i < modelPara.SV; i++){
            bn_new(v[i]);
            bn_zero(v[i]);
            for (int j = 0; j < modelPara.features; j++){
                bn_t tmp;
                bn_new(tmp);
                bn_mul(tmp, modelPara.sv[i][j], X[k][j]);
                bn_add(v[i], v[i], tmp);
            }
            bn_add(sum_v, sum_v, v[i]);
            Sim_Dec(v[i], v[i]);
        }

        gt_t left_side, right_side;
        gt_new(left_side);gt_new(right_side);

        g1_t g1_ls;
        g1_new(g1_ls);
        g1_mul(g1_ls, pk.bgn_pk.g1, sum_v);
        g1_neg(g1_ls, g1_ls);
        g1_add(g1_ls, g1_ls, FK);
        g1_add(g1_ls, g1_ls, Hf);

        pc_map(left_side, g1_ls, sigma0);

        gt_exp_dig(right_side, right_side, 0);
        for (auto i = 0; i < pk.features; i++){
            gt_t lsi;
            gt_new(lsi);
            g2_t g2_rs;
            g2_new(g2_rs);
            g2_mul(g2_rs, pk.bgn_pk.h, X[k][i]); // h^zi
            g2_sub(g2_rs, pk.ht[i], g2_rs); // h^(ti-zi)
            g2_norm(g2_rs, g2_rs);
            pc_map(lsi, sigma[i], g2_rs);
            gt_mul(right_side, right_side, lsi);
        }

        if (gt_cmp(left_side, right_side) == RLC_EQ){
            std::cout << "Batch Verification: True" << std::endl;
        } else {
            std::cout << "Batch Verification: False" << std::endl;
        }
        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }*/



}

void test_dp_input(int features){
    SDPC_PK pk;
    SDPC_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t size;

    SDPC_KeyGen(pk, sk, features);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/SMS_para_128.csv");


    get_user_inputs(X, y, "../data/SMS_test_data_128.csv", size);



    std::vector<double> user_enc_time;
    std::vector<double> server_compute_time;

    double time;
    for (auto in = 0; in < 1000; in++){
        auto *ct = new BGN_CT[modelPara.features];
//        auto *ct2 = new BGN_CT[modelPara.features];
        time = GetTime();
        SDPC_Enc(ct, pk, X[in]);
        time = GetTime() - time;
        user_enc_time.push_back(time * 1000);


        time = GetTime();
        for (int i = 0; i < modelPara.SV; i++) {
            BGN_CT V;
            g2_t sigma0;
            auto *sigma = new g1_t[modelPara.features];
            g1_new(V.c1);
            g1_new(V.c2);
            SDPC_Compute(V, sigma0, sigma, pk, ct, modelPara.sv[i]);
        }
        time = GetTime() - time;
        server_compute_time.push_back(time * 1000);
    }

    // average time
    double mean, stdev;
    // mean of user encryption time
    mean = std::accumulate(user_enc_time.begin(), user_enc_time.end(), 0.0) / user_enc_time.size();
    // stdev of user encryption time
    stdev = std::sqrt(std::inner_product(user_enc_time.begin(), user_enc_time.end(), user_enc_time.begin(), 0.0) / user_enc_time.size() - mean * mean);
//    std::cout << "Feature:" << features << std::endl;
    std::cout << "Mean of user encryption time: " << mean << "ms" << std::endl;
    std::cout << "Standard deviation of user encryption time: " << stdev << std::endl;

    // mean of server computing time
    mean = std::accumulate(server_compute_time.begin(), server_compute_time.end(), 0.0) / server_compute_time.size();
    // stdev of server computing time
    stdev = std::sqrt(std::inner_product(server_compute_time.begin(), server_compute_time.end(), server_compute_time.begin(), 0.0) / server_compute_time.size() - mean * mean);
    std::cout << "Mean of server computing time: " << mean << "ms" << std::endl;
    std::cout << "Standard deviation of server computing time: " << stdev << std::endl;

/*    for (int i = 0; i < modelPara.SV; i++){
        SDPC_FK fk;

        auto time = GetTime();
        SDPC_Setup(fk, pk, modelPara.sv[i]);
        time = GetTime() - time;
        setup_time << time * 1000 << std::endl;

        auto *ct = new BGN_CT[modelPara.features];
        time = GetTime();
        SDPC_Enc(ct, pk, X[0]);
        time = GetTime() - time;
        enc_time << time * 1000 << std::endl;

        BGN_CT V;
        g2_t sigma0;
        auto *sigma = new g1_t[modelPara.features];
        g1_new(V.c1);
        g1_new(V.c2);
        // Set to 0

        time = GetTime();
        SDPC_Compute(V, sigma0, sigma, pk, ct, modelPara.sv[i]);
        time = GetTime() - time;
        compute_time << time * 1000 << std::endl;
        std::cout << "Online computing time: " << time * 1000 << " ms" << std::endl;

        // true result
        bn_t res;
        bn_new(res);
        bn_zero(res);
        for (int j = 0; j < modelPara.features; j++){
            bn_t tmp;
            bn_new(tmp);
            bn_mul(tmp, modelPara.sv[i][j], X[0][j]);
            bn_add(res, res, tmp);
        }
//        time = GetTime();
//        SDPC_Verify(pk, fk, res, X[0], sigma0, sigma);
//        time = GetTime() - time;
//        std::cout << "Verify time: " << time * 1000 << " ms" << std::endl;
    }*/

/*   std::ofstream verify_time("../benchmark/verify_dp_time_"+ std::to_string(features) +".txt");
    for (auto k = 0; k < size; k++){
        auto time = GetTime();
        // Batch Verification
        g1_t FK;
        g1_new(FK);
        g1_sub(FK, FK, FK);

        bn_t d;
        bn_new(d);
        g2_t sigma0;
        g2_new(sigma0);
        bn_rand_mod(d, pk.bgn_pk.q);
        g2_mul(sigma0, pk.bgn_pk.h, d);
        g2_norm(sigma0, sigma0);

        auto sigma = new g1_t[modelPara.features];

        g1_t Hf;
        g1_new(Hf);
        g1_sub(Hf, Hf, Hf);
        bn_t bn_SV;
        bn_new(bn_SV);
        bn_set_dig(bn_SV, modelPara.SV);
        bn_t sum_x;
        bn_new(sum_x);
        bn_zero(sum_x);
        for (int i = 0; i < modelPara.features; i++){
            g1_new(sigma[i]);
            g1_sub(sigma[i], sigma[i], sigma[i]);
            // Generate FK

            for (int j = 0; j < modelPara.SV; j++){
                bn_add(sum_x, sum_x, modelPara.sv[j][i]);
            }

            g1_t tmp;
            g1_new(tmp);
            g1_mul(tmp, pk.g1t[i], sum_x);
            g1_add(FK, FK, tmp);

            // sigma_i
            g1_mul(tmp, pk.bgn_pk.g1, d);
            g1_mul(sigma[i], tmp, sum_x);

            bn_t dSV;
            bn_new(dSV);
            bn_mul_dig(dSV, d, modelPara.SV);
            g1_mul(tmp, pk.g1t[i], dSV);
            g1_add(sigma[i], sigma[i], tmp);

            // Hf
            g1_mul(tmp, pk.g1t[i], X[k][i]);
            g1_add(Hf, Hf, tmp);
        }
//
        g1_mul(Hf, Hf, bn_SV);

        auto v = new bn_t[modelPara.SV];
        bn_t sum_v;
        bn_new(sum_v);
        bn_zero(sum_v);
        for (int i = 0; i < modelPara.SV; i++){
            bn_new(v[i]);
            bn_zero(v[i]);
            for (int j = 0; j < modelPara.features; j++){
                bn_t tmp;
                bn_new(tmp);
                bn_mul(tmp, modelPara.sv[i][j], X[k][j]);
                bn_add(v[i], v[i], tmp);
            }
            bn_add(sum_v, sum_v, v[i]);
            Sim_Dec(v[i], v[i]);
        }

        gt_t left_side, right_side;
        gt_new(left_side);gt_new(right_side);
        g1_t g1_ls;
        g1_new(g1_ls);
        g1_mul(g1_ls, pk.bgn_pk.g1, sum_v);
        g1_neg(g1_ls, g1_ls);
        g1_add(g1_ls, g1_ls, FK);
        g1_add(g1_ls, g1_ls, Hf);

        pc_map(left_side, g1_ls, sigma0);

        gt_exp_dig(right_side, right_side, 0);
        for (auto i = 0; i < pk.features; i++){
            gt_t lsi;
            gt_new(lsi);
            g2_t g2_rs;
            g2_new(g2_rs);
            g2_mul(g2_rs, pk.bgn_pk.h, X[0][i]); // h^zi
            g2_sub(g2_rs, pk.ht[i], g2_rs); // h^(ti-zi)
            g2_norm(g2_rs, g2_rs);
//        std::cout<< i << std::endl;
//        g2_print(g2_rs);
            pc_map(lsi, sigma[i], g2_rs);
            gt_mul(right_side, right_side, lsi);
        }

        if (gt_cmp(left_side, right_side) == RLC_EQ){
            std::cout << "Batch Verification: True" << std::endl;
        } else {
            std::cout << "Batch Verification: False" << std::endl;
        }
        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }*/

}

void test_sed_verify(int size){
    SSEDC_PK pk;
    SSEDC_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t s;

    SSEDC_KeyGen(pk, sk, 200);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/predict_students_dropout_and_academic_success_para_rbf.csv");


    get_user_inputs(X, y, "../data/predict_students_dropout_and_academic_success_test_data_rbf.csv", s);

    std::ofstream verify_time("../benchmark/verify_sed_time_size_"+ std::to_string(size) +".txt");

    for (auto l = 0; l< 1000; l++){
        auto time = GetTime();
        g1_t FK;
        g1_new(FK);
        g1_sub(FK, FK, FK);

        bn_t sum_x, sum_xx, bn_SV, d, bn_size;
        bn_new(sum_x);
        bn_new(sum_xx);
        bn_new(bn_SV);
        bn_new(d);
        bn_new(bn_size);
        bn_zero(sum_x);
        bn_zero(sum_xx);
        bn_set_dig(bn_SV, modelPara.SV);
        bn_rand_mod(d, pk.bgn_pk.q);
        bn_set_dig(bn_size, size);
        g2_t sigma0;
        g2_new(sigma0);
        g2_mul(sigma0, pk.bgn_pk.h, d);

        auto sigma = new g1_t[modelPara.features];

        bn_t sum;
        bn_new(sum);
        bn_zero(sum);
        g1_t Hf;
        g1_new(Hf);

        for (auto i = 0; i < pk.features; i++){
            g1_new(sigma[i]);

            // Generate FK
            g1_t tmp;
            g1_new(tmp);
            g1_mul(tmp, pk.g1t2[i], bn_SV);
            g1_add(FK, FK, tmp);
            for (auto j = 0; j < modelPara.SV; j++){
                bn_t xx;
                bn_new(xx);
                bn_mul(xx, modelPara.sv[j][i], modelPara.sv[j][i]);
                bn_add(sum_xx, sum_xx, xx);
                bn_add(sum_x, sum_x, modelPara.sv[j][i]);
            }
            g1_mul(tmp, pk.g1t[i], sum_x);
            g1_mul_dig(tmp, tmp, 2);
            g1_sub(FK, FK, tmp);

            g1_mul(tmp, pk.bgn_pk.g1, sum_xx);

            g1_add(FK, FK, tmp);

            g1_mul(FK, FK, bn_size);

            // sigma_i
            bn_t r;
            bn_new(r);
            bn_rand_mod(r, pk.bgn_pk.q);
            g1_mul(tmp, pk.bgn_pk.g, r);
            g1_add(sigma[i], pk.g1t[i], tmp);
            bn_t x;
            bn_new(x);
            bn_rand_mod(x, pk.bgn_pk.q);
            g1_mul(tmp, pk.bgn_pk.g1, x);
            g1_add(sigma[i], sigma[i], tmp);
            g1_mul(sigma[i], sigma[i], d);
            g1_mul(sigma[i], sigma[i], bn_SV);
            g1_mul(tmp, pk.bgn_pk.g1, sum_x);
            g1_mul_dig(tmp, tmp, 2);
            g1_sub(sigma[i], sigma[i], tmp);

            // Hf
            bn_t xr;
            bn_new(xr);
            bn_mul(xr, r, x);
            bn_add(sum, sum, xr);
            g1_t gtr;
            g1_new(gtr);
            g1_mul(gtr, pk.gt[i], r);
            g1_add(Hf, Hf, gtr);
        }
        g1_t gsum;
        g1_new(gsum);
        g1_mul(gsum, pk.bgn_pk.g, sum);

        g1_sub(Hf, Hf, gsum);
        g1_mul(Hf, Hf, bn_SV);

        // sum v
        auto v = new bn_t[modelPara.SV];
        bn_t sum_v;
        bn_new(sum_v);
        bn_zero(sum_v);
        for (int k = 0; k < size; k++){
            for (int i = 0; i < modelPara.SV; i++){
                bn_new(v[i]);
                bn_zero(v[i]);
                for (int j = 0; j < modelPara.features; j++){
                    bn_t tmp;
                    bn_new(tmp);
                    bn_t x;
                    bn_new(x);
                    bn_rand_mod(x, pk.bgn_pk.q);
                    bn_mul(tmp, modelPara.sv[i][j], x);
                    bn_add(v[i], v[i], tmp);
                }
                bn_add(sum_v, sum_v, v[i]);
                Sim_Dec(v[i], v[i]);
            }
        }

        gt_t left_side, right_side;
        gt_new(left_side);gt_new(right_side);

        g1_t g1_ls;
        g1_new(g1_ls);
        g1_mul(g1_ls, pk.bgn_pk.g1, sum_v);
        g1_neg(g1_ls, g1_ls);
        g1_add(g1_ls, g1_ls, FK);
        g1_add(g1_ls, g1_ls, Hf);

        pc_map(left_side, g1_ls, sigma0);

        gt_exp_dig(right_side, right_side, 0);
        for (auto k = 0; k < size; k++){
            for (auto i = 0; i < pk.features; i++){
                gt_t lsi;
                gt_new(lsi);
                g2_t g2_rs;
                g2_new(g2_rs);
                bn_t x;
                bn_new(x);
                bn_rand_mod(x, pk.bgn_pk.q);
                g2_mul(g2_rs, pk.bgn_pk.h, x); // h^zi
                g2_sub(g2_rs, pk.ht[i], g2_rs); // h^(ti-zi)
                g2_norm(g2_rs, g2_rs);
                pc_map(lsi, sigma[i], g2_rs);
                gt_mul(right_side, right_side, lsi);
            }
            gt_mul(right_side, right_side, right_side);
        }


        if (gt_cmp(left_side, right_side) == RLC_EQ){
            std::cout << "Batch Verification: True" << std::endl;
        } else {
            std::cout << "Batch Verification: False" << std::endl;
        }
        time = GetTime() - time;
        verify_time << time * 1000 << std::endl;
    }
}


void test_dp_verify(int size){
    SDPC_PK pk;
    SDPC_SK sk;
    bn_t **X = nullptr;
    bn_t *y = nullptr;
    size_t s;

    SDPC_KeyGen(pk, sk, 200);

    ModelPara modelPara;
    bn_new(modelPara.scale);
    bn_new(modelPara.gamma);
    bn_new(modelPara.b);
    bn_new(modelPara.c);
    set_model_paras(modelPara, 0.005, -1.4690208, 2, "../data/SMS_para_128.csv");


    get_user_inputs(X, y, "../data/SMS_test_data_128.csv", s);

//    std::ofstream verify_time("../benchmark/verify_dp_time_size_"+ std::to_string(size) +".txt");

    std::vector<double> verify_time_vec;
    for (auto l = 0; l < 1000; l++){
        // Batch Verification
        g1_t FK;
        g1_new(FK);
        g1_sub(FK, FK, FK);

        bn_t d;
        bn_new(d);
        g2_t sigma0;
        g2_new(sigma0);
        bn_rand_mod(d, pk.bgn_pk.q);
        g2_mul(sigma0, pk.bgn_pk.h, d);
        g2_norm(sigma0, sigma0);

        auto sigma = new g1_t[modelPara.features];

        g1_t Hf;
        g1_new(Hf);
        g1_sub(Hf, Hf, Hf);
        bn_t bn_SV, bn_size;
        bn_new(bn_SV);
        bn_new(bn_size);
        bn_set_dig(bn_SV, modelPara.SV);
        bn_set_dig(bn_size, size);

        bn_t sum_x;
        bn_new(sum_x);
        bn_zero(sum_x);
        for (int i = 0; i < modelPara.features; i++){
            g1_new(sigma[i]);
            g1_sub(sigma[i], sigma[i], sigma[i]);
            // Generate FK

            for (int j = 0; j < modelPara.SV; j++){
                bn_add(sum_x, sum_x, modelPara.sv[j][i]);
            }

            g1_t tmp;
            g1_new(tmp);
            g1_mul(tmp, pk.g1t[i], sum_x);
            g1_add(FK, FK, tmp);

            // sigma_i
            g1_mul(tmp, pk.bgn_pk.g1, d);
            g1_mul(sigma[i], tmp, sum_x);

            bn_t dSV;
            bn_new(dSV);
            bn_mul_dig(dSV, d, modelPara.SV);
            g1_mul(tmp, pk.g1t[i], dSV);
            g1_add(sigma[i], sigma[i], tmp);

            // Hf
            bn_t x, x_tmp;
            bn_new(x);
            bn_zero(x);
            bn_new(x_tmp);
            for (auto k = 0; k < size; k++){
                bn_rand_mod(x_tmp, pk.bgn_pk.q);
                bn_add(x, x, x_tmp);
            }
            g1_mul(tmp, pk.g1t[i], x);
            g1_add(Hf, Hf, tmp);
        }
//
        g1_mul(Hf, Hf, bn_SV);
        g1_mul(Hf, Hf, bn_size);

        auto v = new bn_t[modelPara.SV];
        bn_t sum_v;
        bn_new(sum_v);
        bn_set_2b(sum_v, 12);
//        for (auto k = 0; k< size; k++){
//            for (int i = 0; i < modelPara.SV; i++){
//                bn_new(v[i]);
//                bn_zero(v[i]);
//                for (int j = 0; j < modelPara.features; j++){
//                    bn_t tmp;
//                    bn_new(tmp);
//                    bn_t x;
//                    bn_new(x);
//                    bn_rand_mod(x, pk.bgn_pk.q);
//                    bn_mul(tmp, modelPara.sv[i][j], x);
//                    bn_add(v[i], v[i], tmp);
//                }
//                bn_add(sum_v, sum_v, v[i]);
//            }
////            Sim_Dec(sum_v, sum_v);
//        }
        bn_t mod;
        bn_new(mod);
//        bn_set_dig(mod, 524287);
        bn_rand_mod(sum_v, pk.bgn_pk.q);
        auto time = GetTime();
        for (auto k = 0; k< size; k++){
            for (auto l = 0; l < 20; l++) {
                Sim_Dec(sum_v, sum_v);
            }
        }
        g1_mul(FK, FK, bn_size);

        gt_t left_side, right_side;
        gt_new(left_side);gt_new(right_side);
        g1_t g1_ls;
        g1_new(g1_ls);
        g1_mul(g1_ls, pk.bgn_pk.g1, sum_v);
        g1_neg(g1_ls, g1_ls);
        g1_add(g1_ls, g1_ls, FK);
        g1_add(g1_ls, g1_ls, Hf);
//
        pc_map(left_side, g1_ls, sigma0);

        gt_exp_dig(right_side, right_side, 0);
       /* for (auto i = 0; i < 16; i++){
            gt_t lsi;
            gt_new(lsi);
            g2_t g2_rs;
            g2_new(g2_rs);
            bn_t x, x_tmp;
            bn_new(x);
            bn_zero(x);
            bn_new(x_tmp);
            for (auto k = 0; k < size; k++){
                bn_rand_mod(x_tmp, pk.bgn_pk.q);
                bn_add(x, x, x_tmp);
            }
            g2_mul(g2_rs, pk.bgn_pk.h, x); // h^zi
            g2_t htphi;
            g2_new(htphi);
            g2_mul(htphi, pk.ht[i], bn_size);
            g2_sub(g2_rs, htphi, g2_rs); // h^(ti-zi)
            g2_norm(g2_rs, g2_rs);
//        std::cout<< i << std::endl;
//        g2_print(g2_rs);
            pc_map(lsi, sigma[i], g2_rs);
            gt_mul(right_side, right_side, lsi);
        }*/
//
        if (gt_cmp(left_side, right_side) == RLC_EQ){
//            std::cout << "Batch Verification: True" << std::endl;
        } else {
//            std::cout << "Batch Verification: False" << std::endl;
        }
        time = GetTime() - time;
        verify_time_vec.push_back(time * 1000);
//        verify_time << time * 1000 << std::endl;
    }

    std::cout << "Size: " << size << std::endl;
    double mean, stdev;
    mean = std::accumulate(verify_time_vec.begin(), verify_time_vec.end(), 0.0) / verify_time_vec.size();
    stdev = std::sqrt(std::inner_product(verify_time_vec.begin(), verify_time_vec.end(), verify_time_vec.begin(), 0.0) / verify_time_vec.size() - mean * mean);
    std::cout << "Mean of verification time: " << mean << "ms" << std::endl;
    std::cout << "Standard deviation of verification time: " << stdev << std::endl;
}