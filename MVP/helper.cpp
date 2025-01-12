#include "helper.h"


void read_csv(std::vector<std::vector<std::string>> &data, std::string filename) {
    std::ifstream file(filename);

    // check that the file open correctly
    if (!file.is_open()) {
        std::cerr << "Open file failed" << std::endl;
        return ;
    }

    // read csv file by rows
    std::string line;
    while (std::getline(file, line)) {
        std::vector<std::string> rowData;
        std::istringstream iss(line);
        std::string cell;

        while (std::getline(iss, cell, ',')) {
            rowData.push_back(cell);
        }

        data.push_back(rowData);
    }

    file.close();
}

void data_preprocess(vector<std::vector<double>> &f_data, std::vector<std::vector<std::string>> data) {
    for (const auto& r: data){
        std::vector<double> rows;
        rows.reserve(r.size());
        for (const auto& cell: r){
            rows.push_back(std::stod(cell));
        }
        f_data.push_back(rows);
    }
}


void f_data2bn(bn_t **data, vector<std::vector<double>> f_data, double scale){
    bool neg;
    for (auto i = 0; i < f_data.size(); i++){
        for (auto j = 0; j < f_data[0].size(); j++){
            auto tmp = f_data[i][j] * scale;
            auto int_value = static_cast<long long>(std::round(tmp));
            neg = false;
            if (int_value < 0) {
                int_value = - int_value;
                neg = true;
            }
            bn_read_str(data[i][j], std::to_string(int_value).c_str(), floor(log10(int_value)) + 1, 10);
            if (neg) {
                bn_neg(data[i][j], data[i][j]);
            }
        }
    }
}

void set_model_paras(ModelPara &modelPara, double gamma, double b, double c, std::string para_file) {
    std::vector<std::vector<std::string>> para_data;

    read_csv(para_data, para_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, para_data);

    modelPara.SV = f_data.size();
    modelPara.features = pow(2, ceil(log2(f_data[0].size() - 1)));

    auto scale = pow(2, 16);

    auto int_scale = static_cast<long long>(std::round(scale));

    bn_read_str(modelPara.scale, std::to_string(int_scale).c_str(), floor(log10(int_scale))+1, 10);

//    bn_print(modelPara.scale);

    modelPara.alpha = new bn_t[modelPara.SV];
    for (int i = 0; i < modelPara.SV; i++){
        bn_new(modelPara.alpha[i]);
    }
    modelPara.sv = new bn_t*[modelPara.SV];
    for (int i = 0; i < modelPara.SV; i++){
        modelPara.sv[i] = new bn_t[modelPara.features];
        for (auto j = 0; j < modelPara.features; j++){
            bn_new(modelPara.sv[i][j]);
        }
    }

    auto bn_data = new bn_t*[f_data.size()];
    for (int i = 0; i < f_data.size(); i++){
        bn_data[i] = new bn_t[f_data[0].size()];
        for (auto j = 0; j < f_data[0].size(); j++){
            bn_new(bn_data[i][j]);
        }
    }

    f_data2bn(bn_data, f_data, scale);
    for(auto i = 0; i < modelPara.SV; i++){
        bn_copy(modelPara.alpha[i], bn_data[i][0]);
        for (auto j = 1; j < modelPara.features + 1; j++){
            if (j < f_data[0].size()){
                bn_copy(modelPara.sv[i][j-1], bn_data[i][j]);
            } else {
                bn_zero(modelPara.sv[i][j-1]);
            }
        }
    }
    auto int_gamma = static_cast<long long>(std::round(scale * gamma));
    auto int_b = static_cast<long long>(std::round(scale * b));
    auto int_c = static_cast<long long>(std::round(scale * c));

    bn_read_str(modelPara.gamma, std::to_string(int_gamma).c_str(), floor(log10(int_gamma))+1, 10);
    bn_read_str(modelPara.b, std::to_string(-int_b).c_str(), floor(log10(-int_b))+1, 10);
    bn_neg(modelPara.b, modelPara.b);
    bn_read_str(modelPara.c, std::to_string(int_c).c_str(), floor(log10(int_c))+1, 10);

    for (auto i = 0; i < modelPara.SV; i++){
        delete []bn_data[i];
    }
    delete[] bn_data;
}

void get_user_inputs(bn_t ** &X, bn_t* &y, std::string in_file, size_t &size) {
    std::vector<std::vector<std::string>> input_data;

    read_csv(input_data, in_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, input_data);

    auto scale = pow(10, 10);
    size = f_data.size();
    long features = pow(2, ceil(log2(f_data[0].size() - 1)));


    y = new bn_t[size];
    for (auto i = 0; i < size; i++){
        bn_new(y[i]);
    }
    X = new bn_t*[size];
    for (auto i = 0; i < size; i++){
        X[i] = new bn_t[features];
        for (auto j = 0; j < features; j++){
            bn_new(X[i][j]);
        }
    }

    auto bn_data = new bn_t*[f_data.size()];
    for (int i = 0; i < f_data.size(); i++){
        bn_data[i] = new bn_t[f_data[0].size()];
        for (auto j = 0; j < f_data[0].size(); j++){
            bn_new(bn_data[i][j]);
        }
    }


    f_data2bn(bn_data, f_data, scale);
    for(auto i = 0; i < size; i++){
        bn_copy(y[i], bn_data[i][0]);
        for (auto j = 1; j < features + 1 ; j++){
            if (j < f_data[0].size()){
                bn_copy(X[i][j-1], bn_data[i][j]);
            } else {
                bn_zero(X[i][j-1]);
            }
        }
    }

    for (auto i = 0; i < size; i++){
        delete []bn_data[i];
    }
    delete[] bn_data;
}

