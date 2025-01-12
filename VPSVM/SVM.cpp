#include "SVM.h"

vec_ZZ_p HomMVMult(int b, const EK &ek, const PKE_Para &pkePara, const Vec<MemoryV> &t_M, const Ciphertext &C_v, int m, int n) {
    vec_ZZ_p u;
    MemoryV t_P;
    Vec<vec_ZZ_p> p_b;
    vec_ZZ_p p;
    u.SetLength(m);
    p_b.SetLength(t_M.length());
    double time;
    for (auto i = 0; i < t_M.length(); i++){
        BKS_Mult(t_P, b, ek.bksEk, pkePara, t_M[i], C_v);
        {
            // switch to mod r
            ZZ_pPush push(pkePara.r_context);
            p_b[i] = Decode(pkePara, BKS_Output(b, ek.bksEk, t_P, pkePara));
            p.append(p_b[i]);
        }
    }

    ZZ_pPush push(pkePara.r_context);
    for (auto i = 0; i < m; i++){
        for (auto j = 0; j < n; j++){
            u[i] += p[i * n + j];
        }
    }

    return u;
}

void SVM_Gen(PubPara &para, PKE_PK &pk, EK &ek1, EK &ek2, VK &vk) {
    BKS_Gen(ek1.bksEk, ek2.bksEk, pk, para.pkePara);

    ZZ_p::init(para.pkePara.r);
    ZZ_p alpha;
    random(alpha);

    vk = alpha;

//    std::cout << "start to compute vk" << std::endl;

    vec_ZZ_p vec_alpha;
    vec_alpha.SetLength(para.pkePara.slots, alpha);

    ZZ_p::init(para.pkePara.q);

    BKS_Enc(ek1.C_delta, para.pkePara, pk, Encode(para.pkePara, vec_alpha));
    ek2.C_delta = ek1.C_delta;
}


void SVM_ModelEnc(Vec<Ciphertext> &C, Ciphertext &C_d, Ciphertext &C_b, Ciphertext &C_g, Ciphertext &C_c,
                            PubPara &para, const PKE_PK &pk, const ModelPara &modelPara) {
    C.SetLength(modelPara.n_);
    BKS_Enc(C_d, para.pkePara, pk, modelPara.hat_alpha);
    BKS_Enc(C_b, para.pkePara, pk, modelPara.hat_b);
    BKS_Enc(C_g, para.pkePara, pk, modelPara.hat_gamma);
    BKS_Enc(C_c, para.pkePara, pk, modelPara.hat_c);

    for (auto i = 0; i < modelPara.n_; i++) {
        BKS_Enc(C[i], para.pkePara, pk, modelPara.hat_X[i]);
    }
    para.m = modelPara.m;
    para.n = modelPara.n;
}

// compute protocol of basic scheme
void Compute_poly(ZZ_p &y_1, ZZ_p &y_2, ZZ_p &Phi_1, ZZ_p &Phi_2, const EK &ek1, const EK &ek2, PubPara &para,
                    const Vec<Ciphertext> &C, const Ciphertext &C_d, const Ciphertext &C_b, const Ciphertext &C_g, const Ciphertext &C_c,
                    const ZZ &p, const PKE_PK &pk, VK &vk, const vec_ZZ_p &z,
                    std::vector<double> &Time) { // user inputs
    // Offline
    Vec<MemoryV> t_X_1, t_X_2;
    Vec<MemoryV> t_X_gamma_1, t_X_gamma_2;
    MemoryV t_b_1, t_b_2, t_d_1, t_d_2, t_g_1, t_g_2, t_c_1, t_c_2, t_delta_b_1, t_delta_b_2, t_delta_c_1, t_delta_c_2;
    t_X_1.SetLength(C.length());
    t_X_gamma_1.SetLength(C.length());
    t_X_2.SetLength(C.length());
    t_X_gamma_2.SetLength(C.length());
    // Server 1
    BKS_Load(t_b_1, 1, ek1.bksEk, para.pkePara, C_b);
    BKS_Load(t_d_1, 1, ek1.bksEk, para.pkePara, C_d);
    BKS_Load(t_g_1, 1, ek1.bksEk, para.pkePara, C_g);
    BKS_Load(t_c_1, 1, ek1.bksEk, para.pkePara, C_c);
    BKS_Mult(t_delta_b_1, 1, ek1.bksEk, para.pkePara, t_b_1, ek1.C_delta);
    BKS_Mult(t_delta_c_1, 1, ek1.bksEk, para.pkePara, t_c_1, ek1.C_delta);

    vec_ZZ_p vec_c_1, vec_tau_1;
    {
        ZZ_pPush push(para.pkePara.r_context);
        vec_c_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_c_1, para.pkePara));
        vec_tau_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_delta_c_1, para.pkePara));
    }

    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_X_gamma_1[i], 1, ek1.bksEk, para.pkePara, t_g_1, C[i]);
    }

//    std::cout << "Server 1 offline time: " << Time[0] << std::endl;

    // Server 2
    BKS_Load(t_b_2, 2, ek2.bksEk, para.pkePara, C_b);
    BKS_Load(t_d_2, 2, ek2.bksEk, para.pkePara, C_d);
    BKS_Load(t_g_2, 2, ek2.bksEk, para.pkePara, C_g);
    BKS_Load(t_c_2, 2, ek2.bksEk, para.pkePara, C_c);
    BKS_Mult(t_delta_b_2, 2, ek2.bksEk, para.pkePara, t_b_2, ek2.C_delta);
    BKS_Mult(t_delta_c_2, 2, ek2.bksEk, para.pkePara, t_c_2, ek2.C_delta);

    vec_ZZ_p vec_c_2, vec_tau_2;
    {
        ZZ_pPush push(para.pkePara.r_context);
        vec_c_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_c_2, para.pkePara));
        vec_tau_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_delta_c_2, para.pkePara));
    }

    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_X_gamma_2[i], 2, ek1.bksEk, para.pkePara, t_g_2, C[i]);
    }

//    std::cout << "Server 2 offline time: " << Time[1] << std::endl;

    // Online
    // User input
    vec_ZZ_p z_1, z_2, theta, theta_1, theta_2;
    z_1.SetLength(para.n);
    z_2.SetLength(para.n);
    theta.SetLength(para.n);
    theta_1.SetLength(para.n);
    theta_2.SetLength(para.n);
    auto time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.n; i++) {
            theta[i] = vk * z[i];
            random(z_1[i]);
            random(theta_1[i]);
            z_2[i] = z[i] - z_1[i];
            theta_2[i] = theta[i] - theta_1[i];
        }
    }
    // user share time
    Time[4] = GetTime() - time;
//    for (auto i = 0 ; i < para.pkePara.N / z.length(); i++){
//        z_.append(z);
//    }
//    auto hat_z = Encode(para.pkePara, z_);
//    time = GetTime();
//    Ciphertext C_z;
//    BKS_Enc(C_z, para.pkePara, pk, hat_z);
//    Time[4] = GetTime() - time;
//    std::cout << "User enc time: " << Time[4] << std::endl;

    // Sending z_1 theta_1 z_2 theta_2 to server

    // Server 1 re encrypt
    vec_ZZ_p z_1_, theta_1_, z_2_, theta_2_;
    for (auto i = 0 ; i < para.pkePara.slots / z.length(); i++){
        z_1_.append(z_1);
        theta_1_.append(theta_1);
    }
    auto hat_z_1 = Encode(para.pkePara, z_1_);
    auto hat_theta_1 = Encode(para.pkePara, theta_1_);
    time = GetTime();
    Ciphertext C_z_1, C_theta_1;
    BKS_Enc(C_z_1, para.pkePara, pk, hat_z_1);
    BKS_Enc(C_theta_1, para.pkePara, pk, hat_theta_1);
    // server 1 enc time
    Time[2] = GetTime() - time;

    // Server 2 re encrypt
    for (auto i = 0 ; i < para.pkePara.slots / z.length(); i++){
        z_2_.append(z_2);
        theta_2_.append(theta_2);
    }
    auto hat_z_2 = Encode(para.pkePara, z_2_);
    auto hat_theta_2 = Encode(para.pkePara, theta_2_);
    time = GetTime();
    Ciphertext C_z_2, C_theta_2;
    BKS_Enc(C_z_2, para.pkePara, pk, hat_z_2);
    BKS_Enc(C_theta_2, para.pkePara, pk, hat_theta_2);
    // server 2 enc time
    Time[3] = GetTime() - time;

    // Server 1 reconstruct
    time = GetTime();
    Ciphertext C_z, C_theta;
    BKS_ADD2(C_z, 1, ek1.bksEk, C_z_1, C_z_2);
    BKS_ADD2(C_theta, 1, ek1.bksEk, C_theta_1, C_theta_2);
    vec_ZZ_p u_1 = HomMVMult(1, ek1, para.pkePara, t_X_gamma_1, C_z, para.m, para.n);
    vec_ZZ_p mu_1 = HomMVMult(1, ek1, para.pkePara, t_X_gamma_1, C_theta, para.m, para.n);
    vec_ZZ_p w_1, omega_1;
    {
        ZZ_pPush push(para.pkePara.r_context);
        w_1.SetLength(para.m);
        omega_1.SetLength(para.m);
        for (auto i = 0; i < para.m; i++) {
            w_1[i] = u_1[i] + vec_c_1[i];
            omega_1[i] = mu_1[i] + vec_tau_1[i];
        }
    }
    Time[2] += GetTime() - time;

    // Server 2 reconstruct
    time = GetTime();
    BKS_ADD2(C_z, 2, ek2.bksEk, C_z_1, C_z_2);
    BKS_ADD2(C_theta, 2, ek2.bksEk, C_theta_1, C_theta_2);
    vec_ZZ_p u_2 = HomMVMult(2, ek2, para.pkePara, t_X_gamma_2, C_z, para.m, para.n);
    vec_ZZ_p mu_2 = HomMVMult(2, ek2, para.pkePara, t_X_gamma_2, C_theta, para.m, para.n);
    vec_ZZ_p w_2, omega_2;
    {
        ZZ_pPush push(para.pkePara.r_context);
        w_2.SetLength(para.m);
        omega_2.SetLength(para.m);
        for (auto i = 0; i < para.m; i++) {
            w_2[i] = u_2[i] + vec_c_2[i];
            omega_2[i] = mu_2[i] + vec_tau_2[i];
        }
    }
    Time[3] += GetTime() - time;

    // Sending w and omega to the user.
    // User check
    vec_ZZ_p k, k_1, k_2, phi, phi_1, phi_2;
    {
        ZZ_pPush push(para.pkePara.r_context);

        // User compute kernel function
        time = GetTime();
        ZZ_p w, omega;
        for (auto i = 0; i < para.m; i++) {
            w += w_1[i] + w_2[i];
            omega += omega_1[i] + omega_2[i];
        }
        if (w * vk == omega){
            std::cout << "true" << std::endl;
        }else{
            std::cout << "false" << std::endl;
        }



        k.SetLength(para.m);
        k_1.SetLength(para.m);
        k_2.SetLength(para.m);
        phi.SetLength(para.m);
        phi_1.SetLength(para.m);
        phi_2.SetLength(para.m);
        for (auto i = 0; i < para.m; i++) {
            k[i] = (w_1[i] + w_2[i]);
            power(k[i], k[i], p);
            phi[i] = vk * k[i];
            random(k_1[i]);
            random(phi_1[i]);
            k_2[i] = k[i] - k_1[i];
            phi_2[i] = phi[i] - phi_1[i];
        }

        Time[5] = GetTime() - time;
    }

    // server 1 re encrypt 2
    Ciphertext C_k_1, C_phi_1;
    auto hat_k_1 = Encode(para.pkePara, k_1);
    auto hat_phi_1 = Encode(para.pkePara, phi_1);
    time = GetTime();
    BKS_Enc(C_k_1, para.pkePara, pk, hat_k_1);
    BKS_Enc(C_phi_1, para.pkePara, pk, hat_phi_1);
    Time[2] += GetTime() - time;

    // server 2 re encrypt 2
    Ciphertext C_k_2, C_phi_2;
    auto hat_k_2 = Encode(para.pkePara, k_2);
    auto hat_phi_2 = Encode(para.pkePara, phi_2);
    time = GetTime();
    BKS_Enc(C_k_2, para.pkePara, pk, hat_k_2);
    BKS_Enc(C_phi_2, para.pkePara, pk, hat_phi_2);
    Time[3] += GetTime() - time;


    // Sending C_k to servers
    Ciphertext C_k, C_phi;
    MemoryV t_kd_1, t_kd_2, t_phi_d_1, t_phi_d_2, t_y_1, t_y_2, t_varphi_1, t_varphi_2;
    // Server 1 compute
    time = GetTime();
    BKS_ADD2(C_k, 1, ek1.bksEk, C_k_1, C_k_2);
    BKS_ADD2(C_phi, 1, ek1.bksEk, C_phi_1, C_phi_2);
    BKS_Mult(t_kd_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_k);
    BKS_Mult(t_phi_d_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_phi);
    BKS_ADD1(t_y_1, 1, ek1.bksEk, t_kd_1, t_b_1);
    BKS_ADD1(t_varphi_1, 1, ek1.bksEk, t_phi_d_1, t_delta_b_1);
    Time[2] += GetTime() - time;
    vec_ZZ_p y_vec_1, varphi_vec_1;
    {
        ZZ_pPush push(para.pkePara.r_context);
        y_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_y_1, para.pkePara));
        varphi_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_varphi_1, para.pkePara));
    }
    time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        Phi_1 = ZZ_p::zero();
        for (auto i = 0; i < para.m; i++) {
            y_1 += y_vec_1[i];
            Phi_1 += varphi_vec_1[i];
        }
    }
    Time[2] += GetTime() - time;
//    std::cout << "Server 1 compute time: " << Time[2] << std::endl;

    // Server 2 compute
    time = GetTime();
    BKS_ADD2(C_k, 2, ek2.bksEk, C_k_1, C_k_2);
    BKS_ADD2(C_phi, 2, ek2.bksEk, C_phi_1, C_phi_2);
    BKS_Mult(t_kd_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_k);
    BKS_Mult(t_phi_d_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_phi);
    BKS_ADD1(t_y_2, 2, ek2.bksEk, t_kd_2, t_b_2);
    BKS_ADD1(t_varphi_2, 2, ek2.bksEk, t_phi_d_2, t_delta_b_2);
    Time[3] += GetTime() - time;
    vec_ZZ_p y_vec_2, varphi_vec_2;

    {
        ZZ_pPush push(para.pkePara.r_context);
        y_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_y_2, para.pkePara));
        varphi_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_varphi_2, para.pkePara));
    }
    time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        Phi_2 = ZZ_p::zero();
        for (auto i = 0; i < para.m; i++) {
            y_2 += y_vec_2[i];
            Phi_2 += varphi_vec_2[i];
        }
    }
    Time[3] += GetTime() - time;
}

void Compute_rbf(ZZ_p& y_1, ZZ_p& y_2, ZZ_p &Phi_1, ZZ_p &Phi_2, const EK &ek1, const EK &ek2, PubPara &para,
                        const Vec<Ciphertext> &C, const Ciphertext& C_d, const Ciphertext &C_b, const Ciphertext & C_g, // server inputs
                        const PKE_PK& pk, VK& vk, const vec_ZZ_p& z,
                        std::vector<double> &Time) {
    // Offline
    Vec<MemoryV> t_X_1, t_X_2;
    Vec<MemoryV> t_X_gamma_1, t_X_gamma_2, t_XX_gamma_1, t_XX_gamma_2, t_delta_XX_gamma_1, t_delta_XX_gamma_2;
    MemoryV t_b_1, t_b_2, t_d_1, t_d_2, t_g_1, t_g_2,  t_delta_b_1, t_delta_b_2;
//    Vec<vec_ZZ_p> vec_h_1, vec_h_2, vec_eta_1, vec_eta_2;
    vec_ZZ_p h_1, h_2, eta_1, eta_2;
    vec_ZZ_p H_1, H_2, Eta_1, Eta_2;
    vec_ZZ_p h_i_1, h_i_2, eta_i_1, eta_i_2;
    t_X_1.SetLength(C.length());
    t_X_gamma_1.SetLength(C.length());
    t_XX_gamma_1.SetLength(C.length());
    t_delta_XX_gamma_1.SetLength(C.length());
    t_X_2.SetLength(C.length());
    t_X_gamma_2.SetLength(C.length());
    t_XX_gamma_2.SetLength(C.length());
    t_delta_XX_gamma_2.SetLength(C.length());
    H_1.SetLength(para.m);
    H_2.SetLength(para.m);
    Eta_1.SetLength(para.m);
    Eta_2.SetLength(para.m);
//    vec_h_1.SetLength(C.length());
//    vec_h_2.SetLength(C.length());
//    vec_eta_1.SetLength(C.length());
//    vec_eta_2.SetLength(C.length());
    // Server 1
    BKS_Load(t_b_1, 1, ek1.bksEk, para.pkePara, C_b);
    BKS_Load(t_d_1, 1, ek1.bksEk, para.pkePara, C_d);
    BKS_Load(t_g_1, 1, ek1.bksEk, para.pkePara, C_g);
    BKS_Mult(t_delta_b_1, 1, ek1.bksEk, para.pkePara, t_b_1, ek1.C_delta);

    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_X_gamma_1[i], 1, ek1.bksEk, para.pkePara, t_g_1, C[i]);
        BKS_Mult(t_XX_gamma_1[i], 1, ek1.bksEk, para.pkePara, t_X_gamma_1[i], C[i]);
        BKS_Mult(t_delta_XX_gamma_1[i], 1, ek1.bksEk, para.pkePara, t_XX_gamma_1[i], ek1.C_delta);

        {
            ZZ_pPush push(para.pkePara.r_context);
            h_i_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_XX_gamma_1[i], para.pkePara));
            eta_i_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_delta_XX_gamma_1[i], para.pkePara));
            h_1.append(h_i_1);
            eta_1.append(eta_i_1);
        }
    }
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            H_1[i] = 0;
            for (auto j = 0; j < para.n; j++) {
                H_1[i] += h_1[i * para.n + j];
                Eta_1[i] += eta_1[i * para.n + j];
            }
        }
    }

//    std::cout << "Server 1 offline time: " << Time[0] << std::endl;

    // Server 2
    BKS_Load(t_b_2, 2, ek2.bksEk, para.pkePara, C_b);
    BKS_Load(t_d_2, 2, ek2.bksEk, para.pkePara, C_d);
    BKS_Load(t_g_2, 2, ek2.bksEk, para.pkePara, C_g);
    BKS_Mult(t_delta_b_2, 2, ek2.bksEk, para.pkePara, t_b_2, ek2.C_delta);

    for (auto i = 0; i < C.length(); i++){
        BKS_Mult(t_X_gamma_2[i], 2, ek2.bksEk, para.pkePara, t_g_2, C[i]);
        BKS_Mult(t_XX_gamma_2[i], 2, ek2.bksEk, para.pkePara, t_X_gamma_2[i], C[i]);
        BKS_Mult(t_delta_XX_gamma_2[i], 2, ek2.bksEk, para.pkePara, t_XX_gamma_2[i], ek2.C_delta);

        {
            ZZ_pPush push(para.pkePara.r_context);
            h_i_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_XX_gamma_2[i], para.pkePara));
            eta_i_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_delta_XX_gamma_2[i], para.pkePara));
            h_2.append(h_i_2);
            eta_2.append(eta_i_2);
        }
    }
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.m; i++) {
            H_2[i] = 0;
            for (auto j = 0; j < para.n; j++) {
                H_2[i] += h_2[i * para.n + j];
                Eta_2[i] += eta_2[i * para.n + j];
            }
        }
    }

    // Online
    // Online
    // User input
    vec_ZZ_p z_1, z_2, theta, theta_1, theta_2;
    z_1.SetLength(para.n);
    z_2.SetLength(para.n);
    theta.SetLength(para.n);
    theta_1.SetLength(para.n);
    theta_2.SetLength(para.n);
    auto time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        for (auto i = 0; i < para.n; i++) {
            theta[i] = vk * z[i];
            random(z_1[i]);
            random(theta_1[i]);
            z_2[i] = z[i] - z_1[i];
            theta_2[i] = theta[i] - theta_1[i];
        }
    }
    Time[4] = GetTime() - time;


    // Server 1 re encrypt
    vec_ZZ_p z_1_, theta_1_, z_2_, theta_2_;
    for (auto i = 0 ; i < para.pkePara.slots / z.length(); i++){
        z_1_.append(z_1);
        theta_1_.append(theta_1);
    }
    auto hat_z_1 = Encode(para.pkePara, z_1_);
    auto hat_theta_1 = Encode(para.pkePara, theta_1_);
    time = GetTime();
    Ciphertext C_z_1, C_theta_1;
    BKS_Enc(C_z_1, para.pkePara, pk, hat_z_1);
    BKS_Enc(C_theta_1, para.pkePara, pk, hat_theta_1);
    Time[2] = GetTime() - time;

    // Server 2 re encrypt
    for (auto i = 0 ; i < para.pkePara.slots / z.length(); i++){
        z_2_.append(z_2);
        theta_2_.append(theta_2);
    }
    auto hat_z_2 = Encode(para.pkePara, z_2_);
    auto hat_theta_2 = Encode(para.pkePara, theta_2_);
    time = GetTime();
    Ciphertext C_z_2, C_theta_2;
    BKS_Enc(C_z_2, para.pkePara, pk, hat_z_2);
    BKS_Enc(C_theta_2, para.pkePara, pk, hat_theta_2);
    Time[3] = GetTime() - time;

    // Server 1 reconstruct
    Ciphertext C_z, C_theta;
    MemoryV t_z_1, t_z_2, t_zz_1, t_zz_2, t_theta_z_1, t_theta_z_2, t_gzz_1, t_gzz_2, t_gtz_1, t_gtz_2;
    vec_ZZ_p v_1, v_2, nu_1, nu_2, w_1, w_2, omega_1, omega_2;
    w_1.SetLength(para.m);
    w_2.SetLength(para.m);
    omega_1.SetLength(para.m);
    omega_2.SetLength(para.m);

    time = GetTime();
    BKS_ADD2(C_z, 1, ek1.bksEk, C_z_1, C_z_2);
    BKS_ADD2(C_theta, 1, ek1.bksEk, C_theta_1, C_theta_2);
    BKS_Load(t_z_1, 1, ek1.bksEk, para.pkePara, C_z);
    BKS_Mult(t_zz_1, 1, ek1.bksEk, para.pkePara, t_z_1, C_z);
    BKS_Mult(t_theta_z_1, 1, ek1.bksEk, para.pkePara, t_z_1, ek1.C_delta);
    BKS_Mult(t_gzz_1, 1, ek1.bksEk, para.pkePara, t_zz_1, C_g);
    BKS_Mult(t_gtz_1, 1, ek1.bksEk, para.pkePara, t_theta_z_1, C_g);
    Time[2] = GetTime() - time;
    {
        ZZ_pPush push(para.pkePara.r_context);
        v_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_gzz_1, para.pkePara));
        nu_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_gtz_1, para.pkePara));
    }
    time = GetTime();
    vec_ZZ_p u_1 = HomMVMult(1, ek1, para.pkePara, t_X_gamma_1, C_z, para.m, para.n);
    vec_ZZ_p mu_1 = HomMVMult(1, ek1, para.pkePara, t_X_gamma_1, C_theta, para.m, para.n);
    {
        ZZ_pPush push(para.pkePara.r_context);
        ZZ_p sum_v, sum_nu;
        sum_v = 0;
        sum_nu = 0;
        for (auto j = 0; j < para.n; j++){
            sum_v += v_1[j];
            sum_nu += nu_1[j];
        }
        for (auto i = 0; i < para.m; i++) {
            w_1[i] = H_1[i] - 2 * u_1[i] + sum_v;
            omega_1[i] = Eta_1[1] - 2*mu_1[i] + sum_nu;
        }
    }
    Time[2] += GetTime() - time;

    // Server 2 reconstruct
    time = GetTime();
    BKS_ADD2(C_z, 2, ek2.bksEk, C_z_1, C_z_2);
    BKS_ADD2(C_theta, 2, ek2.bksEk, C_theta_1, C_theta_2);
    BKS_Load(t_z_2, 2, ek2.bksEk, para.pkePara, C_z);
    BKS_Mult(t_zz_2, 2, ek2.bksEk, para.pkePara, t_z_2, C_z);
    BKS_Mult(t_theta_z_2, 2, ek2.bksEk, para.pkePara, t_z_2, ek1.C_delta);
    BKS_Mult(t_gzz_2, 2, ek2.bksEk, para.pkePara, t_zz_2, C_g);
    BKS_Mult(t_gtz_2, 2, ek2.bksEk, para.pkePara, t_theta_z_2, C_g);
    Time[3] = GetTime() - time;
    {
        ZZ_pPush push(para.pkePara.r_context);
        v_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_gzz_2, para.pkePara));
        nu_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_gtz_2, para.pkePara));
    }
    time = GetTime();
    vec_ZZ_p u_2 = HomMVMult(2, ek2, para.pkePara, t_X_gamma_2, C_z, para.m, para.n);
    vec_ZZ_p mu_2 = HomMVMult(2, ek2, para.pkePara, t_X_gamma_2, C_theta, para.m, para.n);
    {
        ZZ_pPush push(para.pkePara.r_context);
        ZZ_p sum_v, sum_nu;
        sum_v = 0;
        sum_nu = 0;
        for (auto j = 0; j < para.n; j++){
            sum_v += v_2[j];
            sum_nu += nu_2[j];
        }
        for (auto i = 0; i < para.m; i++) {
            w_2[i] = H_2[i] - 2 * u_2[i] + sum_v;
            omega_2[i] = Eta_2[1] - 2*mu_2[i] + sum_nu;
        }
    }
    Time[3] += GetTime() - time;



    // Sending w and omega to the user.
    // User check
    vec_ZZ_p k, k_1, k_2, phi, phi_1, phi_2;
    {
        ZZ_pPush push(para.pkePara.r_context);

        // User compute kernel function
        time = GetTime();
        ZZ_p w, omega;
        for (auto i = 0; i < para.m; i++) {
            w += w_1[i] + w_2[i];
            omega += omega_1[i] + omega_2[i];
        }
        if (w * vk == omega){
            std::cout << "true" << std::endl;
        }else{
            std::cout << "false" << std::endl;
        }



        k.SetLength(para.m);
        k_1.SetLength(para.m);
        k_2.SetLength(para.m);
        phi.SetLength(para.m);
        phi_1.SetLength(para.m);
        phi_2.SetLength(para.m);
        for (auto i = 0; i < para.m; i++) {
            k[i] = (w_1[i] + w_2[i]);
            // TODO: exp
            phi[i] = vk * k[i];
            random(k_1[i]);
            random(phi_1[i]);
            k_2[i] = k[i] - k_1[i];
            phi_2[i] = phi[i] - phi_1[i];
        }

        Time[5] = GetTime() - time;
    }
    // server 1 re encrypt 2
    Ciphertext C_k_1, C_phi_1;
    auto hat_k_1 = Encode(para.pkePara, k_1);
    auto hat_phi_1 = Encode(para.pkePara, phi_1);
    time = GetTime();
    BKS_Enc(C_k_1, para.pkePara, pk, hat_k_1);
    BKS_Enc(C_phi_1, para.pkePara, pk, hat_phi_1);
    Time[2] += GetTime() - time;

    // server 2 re encrypt 2
    Ciphertext C_k_2, C_phi_2;
    auto hat_k_2 = Encode(para.pkePara, k_2);
    auto hat_phi_2 = Encode(para.pkePara, phi_2);
    time = GetTime();
    BKS_Enc(C_k_2, para.pkePara, pk, hat_k_2);
    BKS_Enc(C_phi_2, para.pkePara, pk, hat_phi_2);
    Time[3] += GetTime() - time;


    // Sending C_k to servers
    Ciphertext C_k, C_phi;
    MemoryV t_kd_1, t_kd_2, t_phi_d_1, t_phi_d_2, t_y_1, t_y_2, t_varphi_1, t_varphi_2;
    // Server 1 compute
    time = GetTime();
    BKS_ADD2(C_k, 1, ek1.bksEk, C_k_1, C_k_2);
    BKS_ADD2(C_phi, 1, ek1.bksEk, C_phi_1, C_phi_2);
    BKS_Mult(t_kd_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_k);
    BKS_Mult(t_phi_d_1, 1, ek1.bksEk, para.pkePara, t_d_1, C_phi);
    BKS_ADD1(t_y_1, 1, ek1.bksEk, t_kd_1, t_b_1);
    BKS_ADD1(t_varphi_1, 1, ek1.bksEk, t_phi_d_1, t_delta_b_1);
    Time[2] += GetTime() - time;
    vec_ZZ_p y_vec_1, varphi_vec_1;
    {
        ZZ_pPush push(para.pkePara.r_context);
        y_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_y_1, para.pkePara));
        varphi_vec_1 = Decode(para.pkePara, BKS_Output(1, ek1.bksEk, t_varphi_1, para.pkePara));
    }
    time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        Phi_1 = ZZ_p::zero();
        for (auto i = 0; i < para.m; i++) {
            y_1 += y_vec_1[i];
            Phi_1 += varphi_vec_1[i];
        }
    }
    Time[2] += GetTime() - time;
//    std::cout << "Server 1 compute time: " << Time[2] << std::endl;

    // Server 2 compute
    time = GetTime();
    BKS_ADD2(C_k, 2, ek2.bksEk, C_k_1, C_k_2);
    BKS_ADD2(C_phi, 2, ek2.bksEk, C_phi_1, C_phi_2);
    BKS_Mult(t_kd_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_k);
    BKS_Mult(t_phi_d_2, 2, ek2.bksEk, para.pkePara, t_d_2, C_phi);
    BKS_ADD1(t_y_2, 2, ek2.bksEk, t_kd_2, t_b_2);
    BKS_ADD1(t_varphi_2, 2, ek2.bksEk, t_phi_d_2, t_delta_b_2);
    Time[3] += GetTime() - time;
    vec_ZZ_p y_vec_2, varphi_vec_2;

    {
        ZZ_pPush push(para.pkePara.r_context);
        y_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_y_2, para.pkePara));
        varphi_vec_2 = Decode(para.pkePara, BKS_Output(2, ek2.bksEk, t_varphi_2, para.pkePara));
    }
    time = GetTime();
    {
        ZZ_pPush push(para.pkePara.r_context);
        Phi_2 = ZZ_p::zero();
        for (auto i = 0; i < para.m; i++) {
            y_2 += y_vec_2[i];
            Phi_2 += varphi_vec_2[i];
        }
    }
    Time[3] += GetTime() - time;
}

bool Verify(ZZ_p &y, const ZZ_p& y_1, const ZZ_p& y_2, const ZZ_p &phi_1, const ZZ_p &phi_2, VK& vk, PubPara &para){
    ZZ_pPush push(para.pkePara.r_context);
    y = y_1 + y_2;

    auto Phi = phi_1 + phi_2;

    if (y * vk != Phi){
        return false;
    }
    return true;
}
