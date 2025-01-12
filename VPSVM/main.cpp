#include "bench.h"

int main() {
//    test_SVM_Gen();

//    bench_ModelEnc("breast_cancer_wisconsin_diagnostic_para_.csv", "breast_cancer_wisconsin_diagnostic_test_data_.csv");
//    bench_ModelEnc("breast_cancer_wisconsin_diagnostic_para_rbf.csv", "breast_cancer_wisconsin_diagnostic_test_data_rbf.csv");
//
//    bench_ModelEnc("car_evaluation_para.csv", "car_evaluation_test_data.csv");
//    bench_ModelEnc("car_evaluation_para_rbf.csv", "car_evaluation_test_data_rbf.csv");
//
//    bench_ModelEnc("mushroom_para.csv", "mushroom_test_data.csv");
//    bench_ModelEnc("mushroom_para_rbf.csv", "mushroom_test_data_rbf.csv");

    std::cout << "Model Enc" << std::endl;
    std::cout << "WQ:" << std::endl;
    bench_ModelEnc("wine_quality_para.csv", "wine_quality_test_data.csv");
    std::cout << "SSC:" << std::endl;
    bench_ModelEnc("SMS_para_128.csv", "SMS_test_data_128.csv");
    std::cout << "PSDAS:" << std::endl;
    bench_ModelEnc("predict_students_dropout_and_academic_success_para_rbf.csv", "predict_students_dropout_and_academic_success_test_data_rbf.csv");
    std::cout << "MUS:" << std::endl;
    bench_ModelEnc("mushroom_para_rbf.csv", "mushroom_test_data_rbf.csv");



//    bench_compute_poly("breast_cancer_wisconsin_diagnostic_para_.csv", "breast_cancer_wisconsin_diagnostic_test_data_.csv");
    std::cout << "Compute:" << std::endl;
    std::cout << "WQ:" << std::endl;
    bench_compute_poly("wine_quality_para.csv", "wine_quality_test_data.csv");
//    bench_compute_poly("mushroom_para.csv", "mushroom_test_data.csv");
    std::cout << "SSC:" << std::endl;
    bench_compute_poly("SMS_para_128.csv", "SMS_test_data_128.csv");



//    bench_compute_poly(128);
//    bench_compute_poly(256);
//    bench_compute_poly(512);
//    bench_compute_poly(1024);
    std::cout << "PSDAS:" << std::endl;
    bench_compute_rbf("predict_students_dropout_and_academic_success_para_rbf.csv", "predict_students_dropout_and_academic_success_test_data_rbf.csv");
//    bench_compute_rbf("car_evaluation_para_rbf.csv", "car_evaluation_test_data_rbf.csv");
    std::cout << "MUS:" << std::endl;
    bench_compute_rbf("mushroom_para_rbf.csv", "mushroom_test_data_rbf.csv");
//    bench_compute_rbf("SMS_para_rbf_128.csv", "SMS_test_data_128.csv");


//    bench_compute_rbf(64);
//    bench_compute_rbf(128);
//    bench_compute_rbf(256);
//    bench_compute_rbf(512);
//    bench_compute_rbf(1024);



//    bench_verify(100);
//    bench_verify(1000);
//    bench_verify(10000);
//    bench_verify(100000);
//    bench_verify(1000000);

    return 0;
}









