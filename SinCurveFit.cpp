#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense> // 使用Eigen库进行矩阵计算
#define PI 3.1415926


using namespace std;
using namespace Eigen;


// 定义sin函数模型，A是振幅，omega是频率，phi是相位，b是偏移量
double sinFunc(double A, double omega, double phi, double b, double x) {
    return A * sin(omega * x + phi) + b;
}

// 最小二乘拟合sin函数
void fitSinCurve(const vector<double>& xData, const vector<double>& yData, double omega, double& A, double& phi, double& b) {
    int n = xData.size();
    MatrixXd AMatrix(n, 3);
    VectorXd bVector(n);

    // 构造系数矩阵和右侧向量
    for (int i = 0; i < n; i++) {
        AMatrix(i, 0) = sin(omega * xData[i]);
        AMatrix(i, 1) = cos(omega * xData[i]);
        AMatrix(i, 2) = 1.0;
        bVector(i) = yData[i];
    }

    // 使用最小二乘法计算拟合参数
    Vector3d params = AMatrix.colPivHouseholderQr().solve(bVector);

    A = sqrt(params(0) * params(0) + params(1) * params(1));
    phi = atan2(params(1), params(0));
    b = params(2); // 偏移量b
}

int main() {
    // 指定测试参数
    double true_A = 1.5;
    double true_omega = 2.0;
    double true_phi = PI / 4.0;
    double true_b = 0.5;

    // 生成一组测试数据
    vector<double> xData;
    vector<double> yData;
    for (double x = 0.0; x <= 10.0; x += 0.1) {
        double y = sinFunc(true_A, true_omega, true_phi, true_b, x);
        xData.push_back(x);
        yData.push_back(y);
    }

    // 执行sin函数拟合，拟合结果存储在以下变量中
    double fitted_A, fitted_phi, fitted_b;

    // 调用拟合函数
    fitSinCurve(xData, yData, true_omega, fitted_A, fitted_phi, fitted_b);

    // 输出拟合结果和指定的参数结果
    cout << "指定参数：" << endl;
    cout << "振幅 A: " << true_A << endl;
    cout << "角频率 omega: " << true_omega << endl;
    cout << "相位 phi: " << true_phi << endl;
    cout << "偏移量 b: " << true_b << endl;
    cout << "---------------------------" << endl;
    cout << "拟合结果：" << endl;
    cout << "振幅 A: " << fitted_A << endl;
    cout << "相位 phi: " << fitted_phi << endl;
    cout << "偏移量 b: " << fitted_b << endl;

    return 0;
}
