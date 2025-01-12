#pragma once
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <gmp.h>

extern "C"{
#include <relic/relic.h>
}

using namespace std;
using namespace NTL;


typedef struct {
    bn_t *alpha; // (y_j alpha_j) j \in SV
    bn_t ** sv; // support vectors
    bn_t b; // bias

    // poly kernel
    bn_t gamma;
    bn_t c;

    long SV; // num of support vectors
    long features; // num of features
    bn_t scale; // scale factor
} ModelPara;


// Read data from csv file.
void read_csv(std::vector<std::vector<std::string>> &data, std::string filename);
void data_preprocess(std::vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data);
void f_data2bn(bn_t **data, std::vector<std::vector<double>> f_data, double scale);

void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file);
void get_user_inputs(bn_t **&X, bn_t *&y, std::string in_file, size_t &size);