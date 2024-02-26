#include <complex>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>

<<<<<<< HEAD
//matrix mul

//check

=======
>>>>>>> 4495cea560a32078e2e4011add6f682ec1f69447
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

<<<<<<< HEAD
=======

>>>>>>> 4495cea560a32078e2e4011add6f682ec1f69447
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
    matrixH[3] = {2, 4};
    matrixH[7] = {3, 1};
    matrixH[8] = {2, 0};


    matrixL[1] = {4, 1};
    matrixL[4] = {2, 0};
    matrixL[5] = {2, 3};
    matrixL[3] = {4, -1};
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
    std::complex<double> lPodvox[9] = {};
    for (int i = 1; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            l[i] += matrixL[j]*fBasis[i * 9 + j];
        }
    }

    for (int i = 1; i < 9; ++i) {
        lPodvox[i] = std::conj(l[i]);
    }

// ну типа считаем Fmns и Zmns

    int countFmns = 0;
    int countZmns = 0;

    std::complex<double> Fmns[9][9][9] = {0};
    std::complex<double> Zmns[9][9][9] = {0};
    std::complex<double> filler = {0 , 1};
    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            for (int k = 1; k < 9; ++k) {
                std::vector<std::complex<double>> vectorComm = {};
                std::vector<std::complex<double>> vector1 = {};
                std::vector<std::complex<double>> vector2 = {};
                std::vector<std::complex<double>> vectorResult = {};

                std::complex<double> Fm[9] = {};
                std::complex<double> Fn[9] = {};
                std::complex<double> Fs[9] = {};

                std::complex<double> comm[9] = {};
                std::complex<double> result[9] = {};

                for (int o = 0; o < 9; ++o) {
                    Fm[o] = fBasis[i * 9 + o];
                    Fn[o] = fBasis[k * 9 + o];
                    Fs[o] = fBasis[j * 9 + o];
                }

                // расчет коммутатора
                vector1 = matrixMul(3, Fm, Fn);
                vector2 = matrixMul(3, Fn, Fm);
                std::transform(vector1.begin(), vector1.end(), vector2.begin(), std::back_inserter(vectorComm), [&](std::complex<double> one, std::complex<double> two){
                    return one - two;
                });

                for (int p = 0; p < 9; ++p) {
                    comm[p] = vectorComm[p];
                }

                vectorResult = matrixMul(3, Fs, comm);

                for (int p = 0; p < 9; ++p) {
                    result[p] = vectorResult[p];
                }

                Fmns[i][j][k] = -filler * trace(result);
                if (arg(Fmns[i][j][k]) != 0 || imag(Fmns[i][j][k]) != 0) {
                    ++countFmns;
                }
            }
        }
    }

    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            for (int k = 1; k < 9; ++k) {
                std::vector<std::complex<double>> vectorComm = {};
                std::vector<std::complex<double>> vector1 = {};
                std::vector<std::complex<double>> vector2 = {};
                std::vector<std::complex<double>> vectorResult = {};

                std::complex<double> Fm[9] = {};
                std::complex<double> Fn[9] = {};
                std::complex<double> Fs[9] = {};

                std::complex<double> comm[9] = {};
                std::complex<double> result[9] = {};

                for (int o = 0; o < 9; ++o) {
                    Fm[o] = fBasis[i * 9 + o];
                    Fn[o] = fBasis[k * 9 + o];
                    Fs[o] = fBasis[j * 9 + o];
                }

                // расчет коммутатора
                vector1 = matrixMul(3, Fm, Fn);
                vector2 = matrixMul(3, Fn, Fm);
                std::transform(vector1.begin(), vector1.end(), vector2.begin(), std::back_inserter(vectorComm), [&](std::complex<double> one, std::complex<double> two){
                    return one + two;
                });

                for (int p = 0; p < 9; ++p) {
                    comm[p] = vectorComm[p];
                }

                vectorResult = matrixMul(3, Fs, comm);

                for (int p = 0; p < 9; ++p) {
                    result[p] = vectorResult[p];
                }

                Zmns[i][j][k] = Fmns[i][j][k] + filler * trace(result);
                if (arg(Zmns[i][j][k]) != 0 || imag(Zmns[i][j][k]) != 0) {
                    ++countZmns;
                }
            }
        }
    }

    // part C
    std::complex<double> matrixQ[9][9] = {};
    for (int i = 1; i < 9; ++i){
        if (arg(h[i]) != 0 || imag(h[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    matrixQ[m][n] += h[i]*Fmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }

    // part D
    std::complex<double> matrixK[9] = {};
    std::complex<double> fillerMatrixK[9][9] = {};
    for (int i = 1; i < 9; ++i) {
        if (arg(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixK[m][n] += l[i]*Fmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }
    for (int i = 1; i < 9; ++i) {
        for (int m = 1; m < 9; ++m) {
            if (arg(l[m]) != 0 || imag(l[m]) != 0) {
                fillerMatrixK[i][m] = lPodvox[m]*fillerMatrixK[m][i];
            } else {
                fillerMatrixK[i][m] = {0};
            }
        }
    }
    for (int i = 1; i < 9; ++i) {
        for (int m = 1; m < 9; ++m) {
            matrixK[i] += fillerMatrixK[i][m];
        }
    }

    //part E
    std::complex<double> fillerMatrixEone[9][9] = {0};
    for (int i = 1; i < 9; ++i){
        if (arg(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixEone[m][n] += l[i]*Zmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }

    std::complex<double> fillerMatrixEtwo[9][9] = {0};
    for (int i = 1; i < 9; ++i){
        if (arg(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixEtwo[m][n] += lPodvox[i]*Fmns[i][m][n]; //проверить порядок индексов
                }
            }
        }
    }


    std::complex<double> matrixR[9][9] = {};
    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            std::complex<double> qqq = {};
            for (int u = 1; u < 9; ++u) {
                qqq += fillerMatrixEone[i][u] + fillerMatrixEtwo[j][u];
            }
            matrixR[i][j] = qqq;
        }
    }

    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            std::cout << matrixR[i][j].real() << " " << i << " " << j << std::endl;
        }
    }

    int qweqweqwe = 1;

<<<<<<< HEAD
}
=======
}
>>>>>>> 4495cea560a32078e2e4011add6f682ec1f69447
