// OptMeth_lab2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
using namespace std;

const double eps1 = 0.01;
const double eps2 = 0.01;

double f(double x1, double x2, double x3, double x4)
{
    return x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4;
}

double df1(double x1)
{
    return 2 * x1;
}

double df2(double x2)
{
    return 2 * x2;
}

double df3(double x3)
{
    return 2 * x3;
}

double df4(double x4)
{
    return 2 * x4;
}

double h1(double x1, double x2, double x3, double x4)
{
    return x1 + 2 * x2 + 3 * x3 + 5 * x4 - 10;
}

double h2(double x1, double x2, double x3, double x4)
{
    return x1 + 2 * x2 + 5 * x3 + 6 * x4 - 15;
}

double P(double x1, double x2, double x3, double x4, double rk)
{
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return f(x1, x2, x3, x4) + rk * (h_12 * h_12 + h_22 * h_22) / 2.;
}

double Phi(double x1, double x2, double x3, double x4, double rk)
{
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return rk * (h_12 * h_12 + h_22 * h_22) / 2.;
}

double dPdx1(double x1, double x2, double x3, double x4, double rk) {
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return 2. * x1 + rk * (h_12 + h_22);
}

double dPdx2(double x1, double x2, double x3, double x4, double rk) {
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return 2. * x2 + 2. * rk * (h_12 + h_22);
}

double dPdx3(double x1, double x2, double x3, double x4, double rk) {
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return 2. * x3 + rk * (3. * h_12 + 5. * h_22);
}

double dPdx4(double x1, double x2, double x3, double x4, double rk) {
    double h_12 = h1(x1, x2, x3, x4);
    double h_22 = h2(x1, x2, x3, x4);
    return 2. * x4 + rk * (5. * h_12 + 6. * h_22);
}
#pragma region d^2P/dx_i*x_j
double dPdx1x1(double rk) {
    return 2. + 2. * rk;
}
double dPdx1x2(double rk) {
    return 4. * rk;
}
double dPdx1x3(double rk) {
    return 8. * rk;
}
double dPdx1x4(double rk) {
    return 11. * rk;
}
double dPdx2x2(double rk) {
    return 2. + 8. * rk;
}
double dPdx2x3(double rk) {
    return 16. * rk;
}
double dPdx2x4(double rk) {
    return 22. * rk;
}
double dPdx3x3(double rk) {
    return 2. + 34. * rk;
}
double dPdx3x4(double rk) {
    return 45. * rk;
}
double dPdx4x4(double rk) {
    return 2. + 61 * rk;
}
#pragma endregion



vector<double> Markvardt(double x1, double x2, double x3, double x4, double rk, int step)
{
    double h_inv[4][4];
    cout.precision(9);
    int k = 0;

    double l = 10000.;

    #pragma region Start value
        double pr_x1 = x1;
        double pr_x2 = x2;
        double pr_x3 = x3;
        double pr_x4 = x4;
    #pragma endregion 
    #pragma region End value
        double fut_x1 = 0.;
        double fut_x2 = 0.;
        double fut_x3 = 0.;
        double fut_x4 = 0.;
    #pragma endregion

    double d1, d2, d3, d4;

    while (true)
    {
        double P1 = dPdx1(pr_x1, pr_x2, pr_x3, pr_x4, rk);
        double P2 = dPdx2(pr_x1, pr_x2, pr_x3, pr_x4, rk);
        double P3 = dPdx3(pr_x1, pr_x2, pr_x3, pr_x4, rk);
        double P4 = dPdx4(pr_x1, pr_x2, pr_x3, pr_x4, rk);

        if (pow(P1 * P1 + P2 * P2 + P3 * P3 + P4 * P4, 0.5) < eps2) { //if( ||grad(P)|| < eps2 )
            
            //cout << "Markvardt :: step = " << step << endl;
            //cout << "(" << pr_x1 << " , " << pr_x2 << " , " << pr_x3 << " , " << pr_x4 << "),\tP_min = " << P(pr_x1, pr_x2, pr_x4, pr_x4, rk) << endl;
            //cout << "k = " << k << endl;
            vector<double> res = { pr_x1, pr_x2, pr_x3, pr_x4 };
            return res;
        }
        else {
            // 
            double det = (l * l * l + 105. * l * l * rk + 74. * l * rk * rk + 6. * l * l + 420. * l * rk + 148. * rk * rk + 12. * l + 420. * rk + 8.);
            #pragma region inversed Hesse
            h_inv[0][0] = (l * l + 103 * l * rk + 69 * rk * rk + 4 * l + 206 * rk + 4) / det;
            h_inv[0][1] = h_inv[1][0] = -(2 * (2 * l + 5 * rk + 4)) * rk / det;
            h_inv[0][2] = h_inv[2][0] = -(8 * l - 7 * rk + 16) * rk / det;
            h_inv[0][3] = h_inv[3][0] = -(11 * l + 14 * rk + 22) * rk / det;
            h_inv[1][1] = (l * l + 97 * l * rk + 54 * rk * rk + 4 * l + 194 * rk + 4) / det;
            h_inv[1][2] = h_inv[2][1] = -(2 * (8 * l - 7 * rk + 16)) * rk / det;
            h_inv[1][3] = h_inv[3][1] = -(2 * (11 * l + 14 * rk + 22)) * rk / det;
            h_inv[2][2] = (l * l + 71 * l * rk + 5 * rk * rk + 4 * l + 142 * rk + 4) / det;
            h_inv[2][3] = h_inv[3][2] = -(5 * (9 * l + 2 * rk + 18)) * rk / det;
            h_inv[3][3] = (l * l + 44 * l * rk + 20 * rk * rk + 4 * l + 88 * rk + 4) / det;
            #pragma endregion           
            #pragma region vector d_k
            d1 = -(h_inv[0][0] * P1 + h_inv[0][1] * P2 + h_inv[0][2] * P3 + h_inv[0][3] * P4);
            d2 = -(h_inv[1][0] * P1 + h_inv[1][1] * P2 + h_inv[1][2] * P3 + h_inv[1][3] * P4);
            d3 = -(h_inv[2][0] * P1 + h_inv[2][1] * P2 + h_inv[2][2] * P3 + h_inv[2][3] * P4);
            d4 = -(h_inv[3][0] * P1 + h_inv[3][1] * P2 + h_inv[3][2] * P3 + h_inv[3][3] * P4);
            #pragma endregion

            fut_x1 = pr_x1 + d1;
            fut_x2 = pr_x2 + d2;
            fut_x3 = pr_x3 + d3;
            fut_x4 = pr_x4 + d4;

            if (P(fut_x1, fut_x2, fut_x3, fut_x4, rk) < P(pr_x1, pr_x2, pr_x3, pr_x4, rk)) {
                l /= 2;
                k++;
                pr_x1 = fut_x1;
                pr_x2 = fut_x2;
                pr_x3 = fut_x3;
                pr_x4 = fut_x4;
                continue;
            }
            else
            {
                l *= 2;
                continue;
            }
            
        }
    }
}

void Shtraf() {
    cout.precision(9);
    double x1, x2, x3,x4;
    double rk = 0.1;
    double C = 5.;
    int k = 0;
#pragma region start value
    x1 = x2 = x3 = x4 = 5.;
#pragma endregion

    vector<double> x;
M:  x = Markvardt(x1, x2, x3, x4, rk, k);
    double phi = Phi(x[0], x[1], x[2], x[3], rk);
    
    if (abs(phi) <= eps1)
    {
        cout << "step = " << k << "\tPhi = " << phi << endl;
        x1 = x[0];
        x2 = x[1];
        x3 = x[2];
        x4 = x[3];
        cout << "(" << x1 << " , " << x2 << " , " << x3 << " , " << x4 << ")\nf_min = " << f(x1, x2, x3, x4) << endl;;
    }
    else
    {
        rk /= C;
        x1 = x[0];
        x2 = x[1];
        x3 = x[2];
        x4 = x[3];
        k++;
        goto M;
    }

}

int main()
{
    Shtraf();
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
