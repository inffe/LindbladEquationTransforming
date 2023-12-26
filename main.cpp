#include <complex>
#include <vector>
#include <tuple>
#include <algorithm>

std::vector<std::complex<double>> matrixMul(int size, std::complex<double> matrixA[], std::complex<double> matrixB[]) {
    std::vector<std::complex<double>> matrixC;
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            std::complex<double> filler = 0;
            for (int k = 0; k < size; ++k) {
                filler += matrixA[row * size + k] * matrixB[k * size + col];
            }
            matrixC.push_back(filler);
        }
    }
    return matrixC;
}

std::complex<double> trace(std::complex<double> matrix[]) {
    std::complex<double> result = {};
    result += matrix[0];
    result += matrix[5];
    result += matrix[8];
    return result;
}

int main()
{
    int size = 3;
    double N = 3.0;
    double U = 2.2;
    double E0 = 0.5;


// ну типа задаем матрицы
    std::complex<double> matrixH[9] = {{0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0}};

    std::complex<double> matrixL[9] = {{0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0}};


    matrixH[1] = {2, 4};
    matrixH[4] = {1, 0};
    matrixH[5] = {3, 1};
    matrixH[3] = {2,4};
    matrixH[7] = {3, 1};
    matrixH[8] = {2, 0};


    matrixL[1] = {4, 1};
    matrixL[4] = {2, 0};
    matrixL[5] = {2, 3};
    matrixL[3] = {4,-1};
    matrixL[7] = {2, -3};
    matrixL[8] = {4, 0};

// ну типа генерим матрицы базиса

    std::complex<double> fBasis[81] = {{1,0}, {1,0}, {1,0},
                                       {1,0}, {1,0}, {1,0},
                                       {1,0},{1,0}, {1,0},

                                       {}, {1, 0}, {},
                                       {1, 0}, {}, {},
                                       {},{}, {},

                                       {}, {}, {1, 0},
                                       {}, {}, {},
                                       {1, 0},{}, {},

                                       {}, {}, {},
                                       {}, {}, {1, 0},
                                       {},{1, 0}, {},

                                       {}, {0, -1}, {},
                                       {0, 1}, {}, {},
                                       {},{}, {},

                                       {}, {}, {0, -1},
                                       {}, {}, {},
                                       {0, 1},{}, {},

                                       {}, {}, {},
                                       {}, {}, {0, -1},
                                       {},{0, 1}, {},

                                       {1, 0}, {}, {},
                                       {}, {-1, 0}, {},
                                       {},{}, {},

                                       {1, 0}, {}, {},
                                       {}, {1, 0}, {},
                                       {},{}, {-2, 0}};

// ну типа считаем h и l

    std::complex<double> h[9] = {};

    for (int i = 1; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            h[i] += matrixH[j]*fBasis[i * 9 + j];
        }
    }

    // сократить, зачем тут считать только нужные
    std::complex<double> l[9] = {};
    for (int i = 1; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            l[i] += matrixL[j]*fBasis[i * 9 + j];
        }
    }

// ну типа считаем Fmns и Zmns

    std::complex<double> Fmns[8][8][8] = {0};
    std::complex<double> Zmns[8][8][8] = {0};
    std::complex<double> filler = {0 , -1};
    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            for (int k = 1; k < 9; ++k) {
                std::vector<std::complex<double>> vector = {};
                std::vector<std::complex<double>> vector1 = {};
                std::vector<std::complex<double>> vector2 = {};

                std::complex<double> Fm[9] = {};
                std::complex<double> Fn[9] = {};
                std::complex<double> Fs[9] = {};
                for (int o = 0; o < 9; ++o) {
                    Fm[o] = fBasis[i * 9 + o];
                    Fn[o] = fBasis[k * 9 + o];
                    Fs[o] = fBasis[j * 9 + o];
                }

                vector1 = matrixMul(3, Fm, Fn);

                vector2 = matrixMul(3, Fn, Fm);

                std::transform(vector1.begin(), vector1.end(), vector2.begin(), std::back_inserter(vector), [&](std::complex<double> one, std::complex<double> two){
                    return one - two;
                });

                vector1 = matrixMul(3, Fs, vector);

                Fmns[k][j][i] = filler * trace(vector);
            }
        }
    }

    // part C
    std::complex<double> matrixQ[8][8] = {0};
    for (int i = 1; i < 9; ++i){
        if (arg(h[i]) != 0 || imag(h[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 8; ++n) {
                    matrixQ[m][n] += h[i]*Fmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }

    // part D
    std::complex<double> matrixK[8] = {};
    std::complex<double> fillerMatrixK[8][8] = {};
    for (int i = 1; i < 9; ++i) {
        if (arg(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 8; ++n) {
                    fillerMatrixK[m][n] += l[i]*Fmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }
    for (int i = 1; i < 9; ++i) {
        for (int m = 1; m < 9; ++m) {
            for (int n = 1; n < 9; ++n) {
                fillerMatrixK[m][n] = l[i]*fillerMatrixK[m][n];
            }
        }
    }
    for (int i = 1; i < 9; ++i) {
        for (int m = 1; m < 9; ++m) {
            matrixK[i] += fillerMatrixK[i][m];
        }
    }

    //part E
    std::complex<double> fillerMatrixEone[8][8] = {0};
    for (int i = 1; i < 9; ++i){
        if (arg(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 8; ++n) {
                    fillerMatrixEone[m][n] += l[i]*Zmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }
}
