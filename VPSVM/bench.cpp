#include "bench.h"


void test_SVM_Gen(){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 128;

    // compute average time for running SVM_Gen 100 times
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 1000; i++){
        SVM_Gen(para, pk, ek1, ek2, vk);
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Average time for running SVM_Gen 1000 times: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 100 << "ms" << std::endl;
}


void load_data_rbf(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z,  const string& para_file, const string& datafile, int slots) {
    std::vector<std::vector<std::string>> data;
    read_csv(data, "../data/"+para_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, data);
    long features = pow(2, ceil(log2(f_data[0].size() - 1)));

    auto scale = pow(2, 13);
    Vec<vec_ZZ_p> ZZ_p_data;
    f_data2ZZ_p(ZZ_p_data, f_data, scale);
    modelPara.m = ZZ_p_data.length();
    modelPara.n = features;
    modelPara.n_ = ceil(modelPara.m / floor(slots / modelPara.n));
    vec_ZZ_p pad;
    pad.SetLength(features + 1, ZZ_p::zero());
    if (ZZ_p_data.length() < modelPara.n_ * int(floor(slots / modelPara.n))){
        // padding zero
        ZZ_p_data.SetLength(modelPara.n_ * int(floor(slots / modelPara.n)), pad);
    }

    for (auto row: ZZ_p_data){
        modelPara.alpha.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        row_t.SetLength(features, ZZ_p::zero());
        modelPara.X.append(row_t);
    }

    modelPara.gamma = ZZ_p(64);
    modelPara.b = ZZ_p(-5723);

    read_csv(data, "../data/" + datafile);
    data_preprocess(f_data, data);
    f_data2ZZ_p(ZZ_p_data, f_data, scale);

    for (auto row: ZZ_p_data){
        y.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        row_t.SetLength(features, ZZ_p::zero());
        Z.append(row_t);
    }
}

void load_data_poly(ModelPara &modelPara, vec_ZZ_p &y, Vec<vec_ZZ_p> &Z, const string& para_file, const string& datafile, int slots) {
    std::vector<std::vector<std::string>> data;
    read_csv(data, "../data/"+para_file);
    std::vector<std::vector<double>> f_data;
    data_preprocess(f_data, data);

    int features = pow(2, ceil(log2(f_data[0].size() - 1)));

    auto scale = pow(2, 13);
    Vec<vec_ZZ_p> ZZ_p_data;
    f_data2ZZ_p(ZZ_p_data, f_data, scale);
    modelPara.m = ZZ_p_data.length();
    modelPara.n = features;
    modelPara.n_ = ceil(modelPara.m / floor(slots / modelPara.n));
    vec_ZZ_p pad;
    pad.SetLength(features + 1, ZZ_p::zero());
    if (ZZ_p_data.length() < modelPara.n_ * int(floor(slots / modelPara.n))){
        ZZ_p_data.SetLength(modelPara.n_ * int(floor(slots / modelPara.n)), pad);
    }

    for (auto row: ZZ_p_data){
        modelPara.alpha.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        row_t.SetLength(features, ZZ_p::zero());
        modelPara.X.append(row_t);
    }

    modelPara.gamma = ZZ_p(64);
    modelPara.c = ZZ_p(134217728);
    modelPara.p = ZZ(3);
    modelPara.b = ZZ_p(-12529);

    read_csv(data, "../data/" + datafile);
    data_preprocess(f_data, data);
    f_data2ZZ_p(ZZ_p_data, f_data, scale);

    for (auto row: ZZ_p_data){
        y.append(row[0]);
        vec_ZZ_p row_t;
        for (int i = 1; i < row.length(); i++){
            row_t.append(row[i]);
        }
        row_t.SetLength(features, ZZ_p::zero());
        Z.append(row_t);
    }
}

void bench_ModelEnc(const string& para_file, const string& datafile){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 64;

    SVM_Gen(para, pk, ek1, ek2, vk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_poly(modelPara, y, Z, para_file, datafile, para.pkePara.slots);

        std::cout << "n_:" << modelPara.n_ << std::endl;

        Vec<vec_ZZ_p> X_;
        vec_ZZ_p b_, g_, c_;
        b_.SetLength(para.pkePara.slots, ZZ_p::zero());
        b_[0] = modelPara.b;
        g_.SetLength(para.pkePara.slots, modelPara.gamma);
        c_.SetLength(para.pkePara.slots, modelPara.c);

        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);

        for (auto i = 0; i < modelPara.n_; i++) {
            for (auto j = 0; j < ceil(para.pkePara.slots / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.slots / modelPara.n + j]);
            }
        }

        for (auto i = 0; i < modelPara.n_; i++) {
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        modelPara.hat_b = Encode(para.pkePara, b_);
        modelPara.hat_gamma = Encode(para.pkePara, g_);
        modelPara.hat_c = Encode(para.pkePara, c_);
    }
    Vec<Ciphertext> C;
    Ciphertext C_d, C_b, C_g, C_c;

    std::vector<double> Time;
    double time;
    for (auto i = 0; i < 1000; i++){
        time = GetTime();
        SVM_ModelEnc(C, C_d, C_b, C_g, C_c, para, pk, modelPara);
        time = GetTime() - time;
        Time.push_back(time);
    }

    double mean, stdev;
    // mean of time
    mean = std::accumulate(Time.begin(), Time.end(), 0.0) / Time.size();
    // stdev of time
    stdev = std::sqrt(std::inner_product(Time.begin(), Time.end(), Time.begin(), 0.0) / Time.size() - mean * mean);

    std::cout << "Mean of time: " << mean*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time: " << stdev << std::endl;
}

void bench_compute_poly(const string& para_file, const string& datafile){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 64;

    SVM_Gen(para, pk, ek1, ek2, vk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_poly(modelPara, y, Z, para_file, datafile, para.pkePara.slots);

        Vec<vec_ZZ_p> X_;
        vec_ZZ_p b_, g_, c_;
        b_.SetLength(para.pkePara.slots, ZZ_p::zero());
        b_[0] = modelPara.b;
        g_.SetLength(para.pkePara.slots, modelPara.gamma);
        c_.SetLength(para.pkePara.slots, modelPara.c);

        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);

        for (auto i = 0; i < modelPara.n_; i++) {
            for (auto j = 0; j < ceil(para.pkePara.slots / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.slots / modelPara.n + j]);
            }
        }

        for (auto i = 0; i < modelPara.n_; i++) {
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        modelPara.hat_b = Encode(para.pkePara, b_);
        modelPara.hat_gamma = Encode(para.pkePara, g_);
        modelPara.hat_c = Encode(para.pkePara, c_);
    }
    Vec<Ciphertext> C;
    Ciphertext C_d, C_b, C_g, C_c;

    SVM_ModelEnc(C, C_d, C_b, C_g, C_c, para, pk, modelPara);

    std::vector<double> s1_offline;
    std::vector<double> s2_offline;
    std::vector<double> s1_online;
    std::vector<double> s2_online;
    std::vector<double> user_enc;
    std::vector<double> user_compute;

    ZZ_p y_1, y_2, Phi_1, Phi_2;
    for (auto i = 0; i < 1000; i ++){
        // 0 time offline for server 1
        // 1 time offline for server 2
        // 2 time online for server 1
        // 3 time online for server 2
        // 4 time user enc
        // 5 time user compute
        std::vector<double> Time;
        Time.resize(6);
        Compute_poly(y_1, y_2, Phi_1, Phi_2, ek1, ek2, para, C, C_d, C_b, C_g, C_c,
                            modelPara.p, pk, vk, Z[i], Time);
        s1_offline.push_back(Time[0]);
        s2_offline.push_back(Time[1]);
        s1_online.push_back(Time[2]);
        s2_online.push_back(Time[3]);
        user_enc.push_back(Time[4]);
        user_compute.push_back(Time[5]);
    }
//    std::cout << "Feature: " << features <<std::endl;
    // mean s1 offline
    double mean_s1_offline = std::accumulate(s1_offline.begin(), s1_offline.end(), 0.0) / s1_offline.size();
    // stdev s1 offline
    double stdev_s1_offline = std::sqrt(std::inner_product(s1_offline.begin(), s1_offline.end(), s1_offline.begin(), 0.0) / s1_offline.size() - mean_s1_offline * mean_s1_offline);
    std::cout << "Mean of time for s1 offline: " << mean_s1_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 offline: " << stdev_s1_offline << std::endl;

    // mean s2 offline
    double mean_s2_offline = std::accumulate(s2_offline.begin(), s2_offline.end(), 0.0) / s2_offline.size();
    // stdev s2 offline
    double stdev_s2_offline = std::sqrt(std::inner_product(s2_offline.begin(), s2_offline.end(), s2_offline.begin(), 0.0) / s2_offline.size() - mean_s2_offline * mean_s2_offline);
    std::cout << "Mean of time for s2 offline: " << mean_s2_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 offline: " << stdev_s2_offline << std::endl;

    // mean s1 online
    double mean_s1_online = std::accumulate(s1_online.begin(), s1_online.end(), 0.0) / s1_online.size();
    // stdev s1 online
    double stdev_s1_online = std::sqrt(std::inner_product(s1_online.begin(), s1_online.end(), s1_online.begin(), 0.0) / s1_online.size() - mean_s1_online * mean_s1_online);
    std::cout << "Mean of time for s1 online: " << mean_s1_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 online: " << stdev_s1_online << std::endl;

    // mean s2 online
    double mean_s2_online = std::accumulate(s2_online.begin(), s2_online.end(), 0.0) / s2_online.size();
    // stdev s2 online
    double stdev_s2_online = std::sqrt(std::inner_product(s2_online.begin(), s2_online.end(), s2_online.begin(), 0.0) / s2_online.size() - mean_s2_online * mean_s2_online);
    std::cout << "Mean of time for s2 online: " << mean_s2_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 online: " << stdev_s2_online << std::endl;

    // mean user enc
    double mean_user_enc = std::accumulate(user_enc.begin(), user_enc.end(), 0.0) / user_enc.size();
    // stdev user enc
    double stdev_user_enc = std::sqrt(std::inner_product(user_enc.begin(), user_enc.end(), user_enc.begin(), 0.0) / user_enc.size() - mean_user_enc * mean_user_enc);
    std::cout << "Mean of time for user enc: " << mean_user_enc*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user enc: " << stdev_user_enc << std::endl;

    // mean user compute
    double mean_user_compute = std::accumulate(user_compute.begin(), user_compute.end(), 0.0) / user_compute.size();
    // stdev user compute
    double stdev_user_compute = std::sqrt(std::inner_product(user_compute.begin(), user_compute.end(), user_compute.begin(), 0.0) / user_compute.size() - mean_user_compute * mean_user_compute);
    std::cout << "Mean of time for user compute: " << mean_user_compute*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user compute: " << stdev_user_compute << std::endl;


}

void bench_compute_rbf(const string& para_file, const string& datafile){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 64;

    SVM_Gen(para, pk, ek1, ek2, vk);

    ModelPara modelPara;
    vec_ZZ_p y;
    Vec<vec_ZZ_p> Z;

    {
        ZZ_pPush push(para.pkePara.r_context);

        load_data_rbf(modelPara, y, Z, para_file, datafile, para.pkePara.N);

        Vec<vec_ZZ_p> X_;
//        modelPara.n_ = ceil((modelPara.m * modelPara.n) / para.pkePara.N);
        X_.SetLength(modelPara.n_);
        modelPara.hat_X.SetLength(modelPara.n_);
        for (auto i = 0; i < modelPara.n_; i++) {
//        X_[i].SetLength(para.pkePara.N);
            for (auto j = 0; j < ceil(para.pkePara.N / modelPara.n); j++) {
                X_[i].append(modelPara.X[i * para.pkePara.N / modelPara.n + j]);
            }
            modelPara.hat_X[i] = Encode(para.pkePara, X_[i]);
        }
        modelPara.hat_alpha = Encode(para.pkePara, modelPara.alpha);
        vec_ZZ_p b_;
        b_.SetLength(para.pkePara.N, ZZ_p::zero());
        b_[0] = modelPara.b;
        modelPara.hat_b = Encode(para.pkePara, b_);
    }

    Vec<Ciphertext> C;
    Ciphertext C_d, C_b, C_g, C_c;

    SVM_ModelEnc(C, C_d, C_b, C_g, C_c, para, pk, modelPara);

    std::vector<double> s1_offline;
    std::vector<double> s2_offline;
    std::vector<double> s1_online;
    std::vector<double> s2_online;
    std::vector<double> user_enc;
    std::vector<double> user_compute;

    ZZ_p y_1, y_2;
    ZZ_p Phi_1, Phi_2;
    for (auto i = 0; i < 1000; i ++){
        // 0 time offline for server 1
        // 1 time offline for server 2
        // 2 time online for server 1
        // 3 time online for server 2
        // 4 time user enc
        // 5 time user compute
        std::vector<double> Time;
        Time.resize(6);
        Compute_rbf(y_1, y_2, Phi_1, Phi_2, ek1, ek2, para, C, C_d, C_b,
                            C_g, pk, vk, Z[i], Time);
        s1_offline.push_back(Time[0]);
        s2_offline.push_back(Time[1]);
        s1_online.push_back(Time[2]);
        s2_online.push_back(Time[3]);
        user_enc.push_back(Time[4]);
        user_compute.push_back(Time[5]);
    }
//    std::cout << "Feature: " << features <<std::endl;
    // mean s1 offline
    double mean_s1_offline = std::accumulate(s1_offline.begin(), s1_offline.end(), 0.0) / s1_offline.size();
    // stdev s1 offline
    double stdev_s1_offline = std::sqrt(std::inner_product(s1_offline.begin(), s1_offline.end(), s1_offline.begin(), 0.0) / s1_offline.size() - mean_s1_offline * mean_s1_offline);
    std::cout << "Mean of time for s1 offline: " << mean_s1_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 offline: " << stdev_s1_offline << std::endl;

    // mean s2 offline
    double mean_s2_offline = std::accumulate(s2_offline.begin(), s2_offline.end(), 0.0) / s2_offline.size();
    // stdev s2 offline
    double stdev_s2_offline = std::sqrt(std::inner_product(s2_offline.begin(), s2_offline.end(), s2_offline.begin(), 0.0) / s2_offline.size() - mean_s2_offline * mean_s2_offline);
    std::cout << "Mean of time for s2 offline: " << mean_s2_offline*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 offline: " << stdev_s2_offline << std::endl;

    // mean s1 online
    double mean_s1_online = std::accumulate(s1_online.begin(), s1_online.end(), 0.0) / s1_online.size();
    // stdev s1 online
    double stdev_s1_online = std::sqrt(std::inner_product(s1_online.begin(), s1_online.end(), s1_online.begin(), 0.0) / s1_online.size() - mean_s1_online * mean_s1_online);
    std::cout << "Mean of time for s1 online: " << mean_s1_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s1 online: " << stdev_s1_online << std::endl;

    // mean s2 online
    double mean_s2_online = std::accumulate(s2_online.begin(), s2_online.end(), 0.0) / s2_online.size();
    // stdev s2 online
    double stdev_s2_online = std::sqrt(std::inner_product(s2_online.begin(), s2_online.end(), s2_online.begin(), 0.0) / s2_online.size() - mean_s2_online * mean_s2_online);
    std::cout << "Mean of time for s2 online: " << mean_s2_online*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for s2 online: " << stdev_s2_online << std::endl;

    // mean user enc
    double mean_user_enc = std::accumulate(user_enc.begin(), user_enc.end(), 0.0) / user_enc.size();
    // stdev user enc
    double stdev_user_enc = std::sqrt(std::inner_product(user_enc.begin(), user_enc.end(), user_enc.begin(), 0.0) / user_enc.size() - mean_user_enc * mean_user_enc);
    std::cout << "Mean of time for user enc: " << mean_user_enc*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user enc: " << stdev_user_enc << std::endl;

    // mean user compute
    double mean_user_compute = std::accumulate(user_compute.begin(), user_compute.end(), 0.0) / user_compute.size();
    // stdev user compute
    double stdev_user_compute = std::sqrt(std::inner_product(user_compute.begin(), user_compute.end(), user_compute.begin(), 0.0) / user_compute.size() - mean_user_compute * mean_user_compute);
    std::cout << "Mean of time for user compute: " << mean_user_compute*1000 << "ms" << std::endl;
    std::cout << "Standard deviation of time for user compute: " << stdev_user_compute << std::endl;


}


void bench_communication(){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 32;

    SVM_Gen(para, pk, ek1, ek2, vk);

    ZZ_p f;
    random(f);
    std::cout << NumBits(rep(f)) <<std::endl;

}

void bench_verify(int size){
    PubPara para;
    PKE_PK pk;
    EK ek1, ek2;
    VK vk;

    para.pkePara.msg_bit = 256;

    SVM_Gen(para, pk, ek1, ek2, vk);
    ZZ_p y, y_1, y_2, Phi_1, Phi_2;
    {
        ZZ_pPush push(para.pkePara.r_context);
        auto time = GetTime();
        for (auto i = 0; i < size; i++) {
            random(y_1);
            random(y_2);
            random(Phi_1);
            random(Phi_2);
            Verify(y, y_1, y_2, Phi_1, Phi_2, vk, para);
        }
        time = GetTime() - time;
        std::cout << "Verify " << size << " times: " << time*1000 << "ms" << std::endl;
    }
}
