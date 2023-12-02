#include <complex>

int main()
{
    double N = 3.0;
    double U = 2.2;
    double E0 = 0.5;

    std::complex<double> trace = {};

// ну типа задаем матрицы
    //std::complex<double> matrixHd1[3] = {{3, -2}, {3, -2}, {3, -2}};
    //std::complex<double> matrixHd2[3] = {{0, 0}, {0, 0}, {0, 0}};
    //std::complex<double> matrixHd3[3] = {{3, -2}, {3, -2}, {3, -2}};

    std::complex<double> matrixH[9] = {{0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0}};

    std::complex<double> matrixL[9] = {{0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0}};

    for (int i = 0; i < N*N; i+=4) {
        if (i == 0 || i == 4 || i == 8) {
            matrixH[i] = (matrixH[1] * E0 + matrixH[1] * matrixH[1] * U / N);
            trace += matrixH[i];
        }
    }

    matrixH[1] = {2, -4};
    matrixH[5] = {3, 1};
    matrixH[3] = {2,4};
    matrixH[7] = {3, -1};

    matrixL[4] = {1, 0};
    matrixL[8] = {2, 0};





    trace /= N;

// ну типа генерим матрицы базиса

    std::complex<double> fBasis1[9] = {{}, {1, 0}, {},
                                       {1, 0}, {}, {},
                                       {},{}, {}};
    std::complex<double> fBasis2[9] = {{}, {}, {1, 0},
                                       {}, {}, {},
                                       {1, 0},{}, {}};
    std::complex<double> fBasis3[9] = {{}, {}, {},
                                       {}, {}, {1, 0},
                                       {},{1, 0}, {}};

    std::complex<double> fBasis4[9] = {{}, {0, -1}, {},
                                       {0, 1}, {}, {},
                                       {},{}, {}};
    std::complex<double> fBasis5[9] = {{}, {}, {0, -1},
                                       {}, {}, {},
                                       {0, 1},{}, {}};
    std::complex<double> fBasis6[9] = {{}, {}, {},
                                       {}, {}, {0, -1},
                                       {},{0, 1}, {}};

    std::complex<double> fBasis7[9] = {{1, 0}, {}, {},
                                       {}, {-1, 0}, {},
                                       {},{}, {}};
    std::complex<double> fBasis8[9] = {{1, 0}, {}, {},
                                       {}, {1, 0}, {},
                                       {},{}, {-2, 0}};

// ну типа считаем h и l

    std::complex<double> h[9] = {};
    for (int i = 1; i < 9; ++i) {
        if (i == 2 || i == 4 || i == 5 || i == 6) {
            h[i] = {0, 0};
        } else if (i == 1) {
            for (int j = 0; j < N*N; ++j) {
                h[i] += matrixH[j]*fBasis1[j];
            }
        } else if (i == 3) {
            for (int j = 0; j < N*N; ++j) {
                h[i] += matrixH[j]*fBasis3[j];
            }
        } else if (i == 7) {
            for (int j = 0; j < N*N; ++j) {
                h[i] += matrixH[j]*fBasis7[j];
            }
        } else if (i == 8) {
            for (int j = 0; j < N * N; ++j) {
                h[i] += matrixH[j] * fBasis8[j];
            }
        }
    }

    std::complex<double> l[9] = {};
    for (int i = 1; i < 9; ++i) {
        if (i == 1 || i == 2 || i == 3 || i == 5) {
            l[i] = {0, 0};
        } else if (i == 4) {
            for (int j = 0; j < N*N; ++j) {
                l[i] += matrixL[j]*fBasis4[j];
            }
        } else if (i == 6) {
            for (int j = 0; j < N*N; ++j) {
                l[i] += matrixL[j]*fBasis6[j];
            }
        } else if (i == 7   ) {
            for (int j = 0; j < N*N; ++j) {
                l[i] += matrixL[j]*fBasis7[j];
            }
        } else if (i == 8) {
            for (int j = 0; j < N * N; ++j) {
                l[i] += matrixL[j] * fBasis8[j];
            }
        }
    }

// ну типа считаем Fmns и Zmns

    std::complex<double> Fmns = {0};
}
