#include <complex>

int main()
{
    double N = 3.0;
    double U = 2.2;
    double E0 = 0.5;
// ну типа задаем матрицы
    std::complex<double> matrixHd1[3] = {{3, -2}, {3, -2}, {3, -2}};
    std::complex<double> matrixHd2[3] = {{0, 0}, {0, 0}, {0, 0}};
    std::complex<double> matrixHd3[3] = {{3, -2}, {3, -2}, {3, -2}};

    std::complex<double> matrixH[9] = {{0, 0}, {3, -2}, {0, 0}, {3, -2}, {0, 0}, {3, -2}, {0, 0}, {3, -2}, {0, 0}};

    for (int i = 0; i < N*N; ++i) {
        if (i == 0 || i == 4 || i == 8) {
            matrixHd2[i] = (matrixH[1] * E0 + matrixHd1[1]*matrixHd1[1] * U / N);
        }
    }

    std::complex<double> h[8] = {};
    for (int i = 1; i < 9; ++i) {

    }


// ну типа генерим матрицы базиса


}




// ну типа считаем h и l



// ну типа считаем Fmns и Zmns

