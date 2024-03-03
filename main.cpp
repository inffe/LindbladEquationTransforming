#include <complex>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>

#define OST (1/sqrt(2))
#define OSS (1/sqrt(6))

std::complex<double> *matrixMul(int size, std::complex<double> matrixA[], std::complex<double> matrixB[], std::complex<double> resultMatrix[]) {
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            std::complex<double> filler = 0;
            for (int k = 0; k < size; ++k) {
                filler += matrixA[row * size + k] * matrixB[k * size + col];
            }
            resultMatrix[row * size + col] = filler;
        }
    }
    return resultMatrix;
}

std::complex<double> trace(std::complex<double> matrix[], int size) {
    std::complex<double> result;
    for (int i = 0; i < size; ++i) {
        result += matrix[i * size + i];
    }
    return result;
}

int main()
{

    const int size = 3;

    double N = 3.0;
    double U = 2.2;
    double E0 = 0.5;


    // matrix init
    std::complex<double> matrixH[9] = { {0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0} };

    std::complex<double> matrixL[9] = { {0}, {0}, {0},
                                       {0}, {0}, {0},
                                       {0}, {0}, {0} };

    //hamiltonian 

    matrixH[0] = { 1, 0 };
    matrixH[1] = { -3, 0 };
    matrixH[3] = { -3, 0 };
    matrixH[4] = { 2, 0 };
    matrixH[5] = { 5, 0 };
    matrixH[7] = { 5, 0 };
    matrixH[8] = { 3, 0 };

    //lindbladian
    
    matrixL[0] = { 1, 0 };
    matrixL[1] = { 4, 0 };
    matrixL[3] = { -4, 0 };
    matrixL[4] = { 2, 0 };
    matrixL[5] = { -2, 0 };
    matrixL[7] = { 2, 0 };
    matrixL[8] = { 4, 0 };

    // init basis matrix

    std::complex<double> fBasis[81] = { {1,0}, {1,0}, {1,0},
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
                                       {},{}, {-2, 0} };

    // h and l vectors

    std::complex<double> h[size * size] = {};
    for (int i = 1; i < size * size; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i != size * size - 1) {
                h[i] += matrixH[j] * (OST * fBasis[i * size * size + j]);
            }
            else
            {
                h[i] += matrixH[j] * (OSS * fBasis[i * size * size + j]);
            }   
        }
    }

    std::complex<double> l[size*size] = {};
    std::complex<double> lConj[size*size] = {}; 

    for (int i = 1; i < size*size; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i != size * size - 1) {
                l[i] += matrixL[j] * (OST * fBasis[i * size*size + j]);
            }
            else
            {
                l[i] += matrixL[j] * (OSS * fBasis[i * size*size + j]);
            }
        }
    }

    for (int i = 1; i < size*size; ++i) {
        lConj[i] = std::conj(l[i]);
    }

    //Fmns and Zmns

    int countFmns = 0;
    int countZmns = 0;

    std::complex<double> Fmns[size * size][size * size][size * size] = { 0 };
    std::complex<double> Zmns[size * size][size * size][size * size] = { 0 };
    std::complex<double> filler = { 0 , 1 };

    // Fmns
    for (int i = 1; i < size*size; ++i) {
        for (int j = 1; j < size*size; ++j) {
            for (int k = 1; k < size*size; ++k) {
                std::complex<double> ab[size * size];
                std::complex<double> ba[size * size];

                std::complex<double> Fm[size * size];
                std::complex<double> Fn[size * size];
                std::complex<double> Fs[size * size];

                std::complex<double> commutator[size * size];
                std::complex<double> result[size * size];

                for (int l = 0; l < size*size; ++l) { // filling basis matrix (need to fix OST and OSS coeffs)
                    if (l != 8) {
                        Fm[l] = OST * fBasis[i * 9 + l];
                        Fn[l] = OST * fBasis[k * 9 + l];
                        Fs[l] = OST * fBasis[j * 9 + l];
                    }
                    else {
                        Fm[l] = OSS * fBasis[i * 9 + l];
                        Fn[l] = OSS * fBasis[k * 9 + l];
                        Fs[l] = OSS * fBasis[j * 9 + l];
                    }
                    
                }

                matrixMul(3, Fm, Fn, ab); // AB
                matrixMul(3, Fn, Fm, ba); // BA

                for (int z = 0; z < size * size; ++z) { // AB - BA or [Fm, Fn]
                    commutator[z] = ab[z] - ba[z];
                }

                matrixMul(3, Fs, commutator, result); // Fs[Fm, Fn]

                Fmns[i][j][k] = -filler * trace(result, size); // -iTr(Fs[Fm, Fn])

                if (real(Fmns[i][j][k]) != 0 || imag(Fmns[i][j][k]) != 0) { 
                    ++countFmns;
                }
            }
        }
    }

    // Zmns
    for (int i = 1; i < size * size; ++i) {
        for (int j = 1; j < size * size; ++j) { 
            for (int k = 1; k < size * size; ++k) {
                std::complex<double> ab[size * size];
                std::complex<double> ba[size * size];

                std::complex<double> Fm[size * size];
                std::complex<double> Fn[size * size];
                std::complex<double> Fs[size * size];

                std::complex<double> antiCommutator[size * size];
                std::complex<double> result[size * size];

                for (int l = 0; l < size * size; ++l) { // filling basis matrix
                    if (l != 8) {
                        Fm[l] = OST * fBasis[i * 9 + l];
                        Fn[l] = OST * fBasis[k * 9 + l];
                        Fs[l] = OST * fBasis[j * 9 + l];
                    }
                    else {
                        Fm[l] = OSS * fBasis[i * 9 + l];
                        Fn[l] = OSS * fBasis[k * 9 + l];
                        Fs[l] = OSS * fBasis[j * 9 + l];
                    }
                }

                matrixMul(3, Fm, Fn, ab); // AB
                matrixMul(3, Fn, Fm, ba); // BA

                for (int z = 0; z < size * size; ++z) { // AB + BA or {Fm, Fn}
                    antiCommutator[z] = ab[z] + ba[z];
                }

                matrixMul(3, Fs, antiCommutator, result); // Fs{Fm, Fn}

                Zmns[i][j][k] = Fmns[i][j][k] + filler * trace(result, size); // +iTr(Fs{Fm, Fn})

                if (real(Zmns[i][j][k]) != 0 || imag(Zmns[i][j][k]) != 0) { 
                    ++countZmns;
                }
            }
        }
    }

    // part C
    std::complex<double> matrixQ[9][9] = { 0 }; // checked
    for (int i = 1; i < 9; ++i) {
        if (real(h[i]) != 0 || imag(h[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    matrixQ[n][m] += h[i] * Fmns[i][m][n]; 
                }
            }
        }
    }

    // part D
    std::complex<double> matrixK[9] = { 0 };
    std::complex<double> fillerMatrixK[9][9] = { 0 };
    for (int i = 1; i < 9; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixK[m][n] += l[i] * Fmns[i][m][n]; 
                }
            }
        }
    }
    for (int i = 1; i < 9; ++i) {
        for (int m = 1; m < 9; ++m) {
            fillerMatrixK[i][m] = lConj[m] * fillerMatrixK[i][m]; // need to proof the type 
        }
    }
    for (int i = 1; i < 9; ++i) { // merge with previous loop
        for (int m = 1; m < 9; ++m) {
            matrixK[i] += fillerMatrixK[i][m]; 
        }
    }

    //part E
    std::complex<double> fillerMatrixEone[9][9] = { 0 };
    for (int i = 1; i < 9; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixEone[m][n] += l[i] * Zmns[i][m][n];
                }
            }
        }
    }

    std::complex<double> fillerMatrixEtwo[9][9] = { 0 };
    for (int i = 1; i < 9; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < 9; ++m) {
                for (int n = 1; n < 9; ++n) {
                    fillerMatrixEtwo[m][n] += lConj[i] * Fmns[i][m][n]; 
                }
            }
        }
    }

    std::complex<double> matrixR[9][9] = {};
    for (int i = 1; i < 9; ++i) {
        for (int j = 1; j < 9; ++j) {
            std::complex<double> filler = {};
            for (int u = 1; u < 9; ++u) {
                filler += fillerMatrixEone[i][u] + fillerMatrixEtwo[j][u]; 
            }
            matrixR[i][j] = filler;
        }
    }

    // rk4 for initial system


    // rk4 for result system

}
