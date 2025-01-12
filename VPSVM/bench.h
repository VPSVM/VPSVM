#ifndef PVPSVMC_BENCH_H
#define PVPSVMC_BENCH_H

#include "SVM.h"

void test_SVM_Gen();

void bench_ModelEnc(const string& para_file, const string& datafile);

void load_data_poly(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, const string& para_file, const string& datafile, int slots);

void load_data_rbf(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, const string& para_file, const string& datafile, int slots);

void bench_compute_poly(const string& para_file, const string& datafile);

void bench_compute_rbf(const string& para_file, const string& datafile);

void bench_verify(int size);

void bench_communication();



#endif //PVPSVMC_BENCH_H
