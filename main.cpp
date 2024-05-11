#include <complex>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <array>

#define OST (1/sqrt(2))
#define OSS (1/sqrt(6))
#define ODT (1/3.0)
#define M_PI 3.14159265358979323846

void matrixMul(int size, double matrixA[], double matrixB[], double resultMatrix[]) {
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            double filler = 0;
            for (int k = 0; k < size; ++k) {
                filler += matrixA[row * size + k] * matrixB[k * size + col];
            }
            resultMatrix[row * size + col] = filler;
        }
    }
}

void cmplmatrixMul(int size, std::complex<double> matrixA[], std::complex<double> matrixB[], std::complex<double> resultMatrix[]) {
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            std::complex<double> filler = 0;
            for (int k = 0; k < size; ++k) {
                filler += matrixA[row * size + k] * matrixB[k * size + col];
            }
            resultMatrix[row * size + col] = filler;
        }
    }
}

void matvecMul(int size, double matrix[][8], std::vector<double> vectV, std::vector<double> tmp) {
    for (int i = 0; i < size * size - 1; ++i) {
        for (int j = 0; j < size * size - 1; ++j) {
            tmp[i] += matrix[i][j] * vectV[j];
        }
    }
}

void matrixAdd(int size, double matrixOne[][9], double matrixTwo[][9], double resultMatrix[][8]) {
    for (int i = 0; i < size * size - 1; ++i) {
        for (int j = 0; j < size * size - 1; ++j) {
            resultMatrix[i][j] = matrixOne[i + 1][j + 1] + matrixTwo[i + 1][j + 1];
        }
    }
}

void matrixTripleMul(int size, double matrixOne[], double matrixTwo[], double matrixThree[], double resultMatrix[]) {

    double tmpArray[9] = { 0 };

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double tmp = 0;
            for (int k = 0; k < size; ++k) {
                tmp += matrixOne[i * size + k] * matrixTwo[k * size + j];
            }
            tmpArray[i * size + j] = tmp;
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double tmp = 0;
            for (int k = 0; k < size; ++k) {
                tmp += tmpArray[i * size + k] * matrixThree[k * size + j];
            }
            resultMatrix[i * size + j] = tmp;
        }
    }
}

std::vector<double> vectorAdd(int size, std::vector<double> vectorV, double vectorK[]) {
    std::vector<double> res(8);
    for (int i = 0; i < size * size - 1; ++i) {
        res[i] = vectorV[i] + vectorK[i + 1];
    }
    return res;
}

std::complex<double> trace(std::complex<double> matrix[], int size) {
    std::complex<double> result = { 0, 0 };
    for (int i = 0; i < size; ++i) {
        result += matrix[i * size + i];
    }
    return result;
}

std::vector<double> func(const int size, double t, double matrixQOne[][9], double matrixQTwo[][9], double matrixR[][9], std::vector<double> vectorV, double matrixK[]) { //разобраться с форматами хранения данных (формула 10)
    double filler[8][8] = { 0 };
    std::vector<double> tmp(8);
    std::vector<double> vvv(8);
    if (t < M_PI) {

        matrixAdd(size, matrixQOne, matrixR, filler); // 8x8 + 8x8 = 8x8

        matvecMul(size, filler, vectorV, tmp); // 8x8 * 8x1 = 8x1

        vvv = vectorAdd(size, tmp, matrixK); // 8x1 + 8x1

    }
    else {

        matrixAdd(size, matrixQTwo, matrixR, filler); // 8x8 + 8x8 = 8x8

        matvecMul(size, filler, vectorV, tmp);

        vvv = vectorAdd(size, tmp, matrixK); // 8x1 + 8x1
    }

    return vvv;
}

double *funcInit(const int size, double t, double hamiltonian[], double hamiltonian2[],  double lindbladian[], double lindbladianT[], double ro[]) {

    double tmp1[9] = { 0 };
    double tmp2[9] = { 0 };
    double tmp3[9] = { 0 };

    double comtmp[9] = { 0 };
    double comtmp2[9] = { 0 };

    if (t < M_PI) {

        //anti-commutator [H(t), ro] + dissipation part of linbladian

        //commutator AB - BA

        matrixMul(size, hamiltonian, ro, comtmp);
        matrixMul(size, ro, hamiltonian, comtmp2);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                comtmp[i * size + j] -= comtmp2[i * size + j];
            }
        }

        // dissipation part 

        matrixTripleMul(size, lindbladian, ro, lindbladianT, tmp1); // LpL+
        matrixTripleMul(size, lindbladianT, lindbladian, ro,  tmp2); // L+Lp
        matrixTripleMul(size, ro, lindbladianT, lindbladian, tmp3); // pL+L

        // 

        for (int i = 0; i < size; ++i) { 
            for (int j = 0; j < size; ++j) {
                tmp1[i * size + j] -= 0.5 * (tmp2[i * size + j] + tmp3[i * size + j]);
                tmp1[i * size + j] = tmp1[i * size + j] * 0.05;
            }
        }
    
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                tmp1[i * size + j] += comtmp[i * size + j];
            }
        }
    }
    else {

        matrixMul(size, hamiltonian2, ro, comtmp);
        matrixMul(size, ro, hamiltonian2, comtmp2);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                comtmp[i * size + j] -= comtmp2[i * size + j];
            }
        }

        // dissipation part 

        matrixTripleMul(size, lindbladian, ro, lindbladianT, tmp1); // LpL+
        matrixTripleMul(size, lindbladianT, lindbladian, ro, tmp2); // L+Lp
        matrixTripleMul(size, ro, lindbladianT, lindbladian, tmp3); // pL+L

        // 

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                tmp1[i * size + j] -= 0.5 * (tmp2[i * size + j] + tmp3[i * size + j]);
                tmp1[i * size + j] = tmp1[i * size + j] * 0.05;
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                tmp1[i * size + j] += comtmp[i * size + j];
            }
        }
    }
    
    return tmp1;
}

std::vector<double> rk4(int size, double t, double matrixQOne[][9], double matrixQTwo[][9], double matrixR[][9], std::vector<double> vectorV, double matrixK[]) {

    std::vector<double> result(size * size - 1);

    double h = 0.05;

    std::vector <double> k1(size * size - 1);

    std::vector <double> k1shift(size * size - 1);
    std::vector <double> k2shift(size * size - 1);
    std::vector <double> k3shift(size * size - 1);

    k1 = func(size, t, matrixQOne, matrixQTwo, matrixR, vectorV, matrixK);

    for (int i = 0; i < size * size - 1; ++i) {
        k1shift[i] = vectorV[i] + (h / 2) * k1[i];
    }

    auto k2 = func(size, t + h / 2, matrixQOne, matrixQTwo, matrixR, k1shift, matrixK);

    for (int i = 0; i < size * size - 1; ++i) {
        k2shift[i] = vectorV[i] + (h / 2) * k2[i];
    }

    auto k3 = func(size, t + h / 2, matrixQOne, matrixQTwo, matrixR, k2shift, matrixK);

    for (int i = 0; i < size * size - 1; ++i) {
        k3shift[i] = vectorV[i] + h * k3[i];
    }

    auto k4 = func(size, t + h, matrixQOne, matrixQTwo, matrixR, k3shift, matrixK);

    for (int i = 0; i < size * size - 1; ++i) {
        result[i] = vectorV[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return result;
}

double* rk4init(const int size, double t, double hamiltonian1[], double hamiltonian2[], double lindbladian[], double lindbladianT[], double ro[]) {

    double result[9];

    double *k1;
    double *k2;
    double* k3;
    double* k4;
    double kShift[9] = { 0 };

    double h = 0.05;

    k1 = funcInit(size, t, hamiltonian1, hamiltonian2, lindbladian, lindbladianT, ro);

    for (int i = 0; i < size * size - 1; ++i) {
        kShift[i] = ro[i] + (h / 2) * k1[i];
    }

    k2 = funcInit(size, t + h / 2, hamiltonian1, hamiltonian2, lindbladian, lindbladianT, kShift);

    for (int i = 0; i < size * size - 1; ++i) {
        kShift[i] = ro[i] + (h / 2) * k2[i];
    }

    k3 = funcInit(size, t + h / 2, hamiltonian1, hamiltonian2, lindbladian, lindbladianT, kShift);

    for (int i = 0; i < size * size - 1; ++i) {
        kShift[i] = ro[i] + h * k3[i];
    }

    k4 = funcInit(size, t + h, hamiltonian1, hamiltonian2, lindbladian, lindbladianT, kShift);

    for (int i = 0; i < size * size - 1; ++i) {
        result[i] = ro[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return result;

}

int main()
{

    const int size = 3;
    const int m = size * size - 1;

    double J = -1.0;
    double N = 3.0;
    double U = 2.2;
    double E0 = 0.5;
    double gamma = 0.1;

    // rewrite matrix init

    //hamiltonian

    double hamiltonian1[size * size] = { 0 };

    hamiltonian1[0] = -9;
    hamiltonian1[1] = sqrt(3);
    hamiltonian1[3] = sqrt(3);
    hamiltonian1[4] = 9.2;
    hamiltonian1[5] = 2;
    hamiltonian1[7] = 2;
    hamiltonian1[8] = -0.2;

    double hamiltonian2[size * size] = { 0 };

    hamiltonian2[0] = -7;
    hamiltonian2[1] = sqrt(3);
    hamiltonian2[3] = sqrt(3);
    hamiltonian2[4] = 5.2;
    hamiltonian2[5] = 2;
    hamiltonian2[7] = 2;
    hamiltonian2[8] = 1.8;

    //добавить нормировку главной диагонали

    //lindbladian

    double linbladian[size * size] = { 0 };

    linbladian[0] = -1.6667;
    linbladian[1] = sqrt(3);
    linbladian[3] = -sqrt(3);
    linbladian[4] = 1.3333;
    linbladian[5] = 2;
    linbladian[7] = -2;
    linbladian[8] = 0.3333;

    double lindbladianT[size * size] = { 0 };

    lindbladianT[0] = -1.6667;
    lindbladianT[1] = -sqrt(3);
    lindbladianT[3] = sqrt(3);
    lindbladianT[4] = 1.3333;
    lindbladianT[5] = -2;
    lindbladianT[7] = 2;
    lindbladianT[8] = 0.3333;


    // init basis matrix

    std::complex<double> fBasis[81] = { {ODT,0}, {ODT,0}, {ODT,0},
                                       {ODT,0}, {ODT,0}, {ODT,0},
                                       {ODT,0},{ODT,0}, {ODT,0},

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
    std::complex<double> h2[size * size] = {};

    for (int i = 1; i < size * size; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i != size * size - 1) {
                h[i] += hamiltonian1[j] * (OST * fBasis[i * size * size + j]);
            }
            else
            {
                h[i] += hamiltonian1[j] * (OSS * fBasis[i * size * size + j]);
            }
        }
    }

    for (int i = 1; i < size * size; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i != size * size - 1) {
                h2[i] += hamiltonian2[j] * (OST * fBasis[i * size * size + j]);
            }
            else
            {
                h2[i] += hamiltonian2[j] * (OSS * fBasis[i * size * size + j]);
            }
        }
    }

    std::complex<double> l[size * size] = {};
    std::complex<double> lConj[size * size] = {};

    for (int i = 1; i < size * size; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i != size * size - 1) {
                l[i] += linbladian[j] * (OST * fBasis[i * size * size + j]);
            }
            else
            {
                l[i] += linbladian[j] * (OSS * fBasis[i * size * size + j]);
            }
        }
    }

    for (int i = 1; i < size * size; ++i) {
        lConj[i] = std::conj(l[i]);
    }

    //Fmns and Zmns

    int countFmns = 0;
    int countZmns = 0;

    std::complex<double> Fmns[size * size][size * size][size * size] = { 0 };
    std::complex<double> Zmns[size * size][size * size][size * size] = { 0 };
    std::complex<double> filler = { 0 , 1 };

    // Fmns
    for (int i = 1; i < size * size; ++i) {
        for (int j = 1; j < size * size; ++j) {
            for (int k = 1; k < size * size; ++k) {
                std::complex<double> ab[size * size];
                std::complex<double> ba[size * size];

                std::complex<double> Fm[size * size];
                std::complex<double> Fn[size * size];
                std::complex<double> Fs[size * size];

                std::complex<double> commutator[size * size];
                std::complex<double> result[size * size];

                for (int l = 0; l < size * size; ++l) { // filling basis matrix (need to fix OST and OSS coeffs)
                    if (l != 8) {
                        Fm[l] = OST * fBasis[i * size * size + l];
                        Fn[l] = OST * fBasis[k * size * size + l];
                        Fs[l] = OST * fBasis[j * size * size + l];
                    }
                    else {
                        Fm[l] = OSS * fBasis[i * size * size + l];
                        Fn[l] = OSS * fBasis[k * size * size + l];
                        Fs[l] = OSS * fBasis[j * size * size + l];
                    }

                }

                cmplmatrixMul(3, Fm, Fn, ab); // AB
                cmplmatrixMul(3, Fn, Fm, ba); // BA

                for (int z = 0; z < size * size; ++z) { // AB - BA or [Fm, Fn]
                    commutator[z] = ab[z] - ba[z];
                }

                cmplmatrixMul(3, Fs, commutator, result); // Fs[Fm, Fn]

                Fmns[i][j][k] = -filler * trace(result, size); // -iTr(Fs[Fm, Fn])

                if (imag(Fmns[i][j][k]) != 0) {
                    ++countFmns;
                }
            }
        }
    }

    // Zmns
    for (int i = 1; i < size * size; ++i) {
        for (int j = 1; j < size * size; ++j) {
            for (int k = 1; k < size * size; ++k) {
                std::complex<double> ab[size * size] = { 0 };
                std::complex<double> ba[size * size] = { 0 };

                std::complex<double> Fm[size * size] = { 0 };
                std::complex<double> Fn[size * size] = { 0 };
                std::complex<double> Fs[size * size] = { 0 };

                std::complex<double> antiCommutator[size * size] = { 0 };
                std::complex<double> result[size * size] = { 0 };

                for (int l = 0; l < size * size; ++l) { // filling basis matrix
                    if (l != 8) {
                        Fm[l] = OST * fBasis[i * size * size + l];
                        Fn[l] = OST * fBasis[k * size * size + l];
                        Fs[l] = OST * fBasis[j * size * size + l];
                    }
                    else {
                        Fm[l] = OSS * fBasis[i * size * size + l];
                        Fn[l] = OSS * fBasis[k * size * size + l];
                        Fs[l] = OSS * fBasis[j * size * size + l];
                    }
                }

                cmplmatrixMul(3, Fm, Fn, ab); // AB
                cmplmatrixMul(3, Fn, Fm, ba); // BA

                for (int z = 0; z < size * size; ++z) { // AB + BA or {Fm, Fn}
                    antiCommutator[z] = ab[z] + ba[z];
                }

                cmplmatrixMul(3, Fs, antiCommutator, result); // Fs{Fm, Fn}

                Zmns[i][j][k] = Fmns[i][j][k] + filler * trace(result, size); // +iTr(Fs{Fm, Fn})

                if (imag(Zmns[i][j][k]) != 0) {
                    ++countZmns;
                }
            }
        }
    }

    // part C
    double matrixQ[size * size][size * size] = { 0 }; // checked
    double matrixQ2[size * size][size * size] = { 0 };

    for (int i = 1; i < size * size; ++i) {
        if (real(h[i]) != 0 || imag(h[i]) != 0) {
            for (int m = 1; m < size * size; ++m) {
                for (int n = 1; n < size * size; ++n) {
                    matrixQ[n][m] += real(h[i] * Fmns[i][m][n]);
                }
            }
        }
    }

    for (int i = 1; i < size * size; ++i) {
        if (real(h2[i]) != 0 || imag(h2[i]) != 0) {
            for (int m = 1; m < size * size; ++m) {
                for (int n = 1; n < size * size; ++n) {
                    matrixQ2[n][m] += real(h2[i] * Fmns[i][m][n]);
                }
            }
        }
    }

    // part D
    double matrixK[size * size] = { 0 };
    std::complex<double> fillerMatrixK[size * size][size * size] = { 0 };
    for (int i = 1; i < size * size; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < size * size; ++m) {
                for (int n = 1; n < size * size; ++n) {
                    fillerMatrixK[m][n] += l[i] * Fmns[i][m][n]; // checked
                }
            }
        }
    }

    for (int i = 1; i < size * size; ++i) {
        for (int m = 1; m < size * size; ++m) { // написать проверку чтобы не считать лишнее
            fillerMatrixK[i][m] = lConj[m] * fillerMatrixK[i][m]; // need to proof the type 
        }
    }

    for (int i = 1; i < size * size; ++i) { // merge with previous loop
        for (int m = 1; m < size * size; ++m) {
            matrixK[i] += real(fillerMatrixK[i][m] * filler) / 3.0; // check matrixK[5] and add i/N
        }
    }

    //part E
    std::complex<double> fillerMatrixEone[size * size][size * size] = { 0 };
    for (int i = 1; i < size * size; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < size * size; ++m) {
                for (int n = 1; n < size * size; ++n) {
                    fillerMatrixEone[m][n] += l[i] * Zmns[i][m][n]; // checked
                }
            }
        }
    }

    std::complex<double> fillerMatrixEtwo[size * size][size * size] = { 0 };
    for (int i = 1; i < size * size; ++i) {
        if (real(l[i]) != 0 || imag(l[i]) != 0) {
            for (int m = 1; m < size * size; ++m) {
                for (int n = 1; n < size * size; ++n) {
                    fillerMatrixEtwo[m][n] += lConj[i] * Fmns[i][m][n]; // checked
                }
            }
        }
    }

    double matrixR[size * size][size * size] = { 0 };
    for (int i = 1; i < size * size; ++i) {
        for (int j = 1; j < size * size; ++j) {
            std::complex<double> fillerR = { 0 };
            for (int u = 1; u < size * size; ++u) {
                matrixR[j][i] += -gamma * real(fillerMatrixEone[i][u] * fillerMatrixEtwo[j][u]);
            }
        }
    }

    int rty = 15;

    // rk4 for initial system

    double ro[size * size] = { 0 };

    for (double t = 0; t < 2 * M_PI; t += 0.05) {

        double* inittmp;
        
        inittmp = rk4init(size, t, hamiltonian1, hamiltonian2, linbladian, lindbladianT, ro);

        for (int i = 0; i < size * size; ++i) {
            ro[i] = inittmp[i];
        }
    }
    

    // rk4 for result system

    std::vector<double> vectorV(8);


    for (double t = 0; t < 2 * M_PI; t += 0.05) {   
        std::vector<double> tmpp(8);

        tmpp = rk4(size, t, matrixQ, matrixQ2, matrixR, vectorV, matrixK);

        for (int i = 0; i < 8; ++i) {
            vectorV[i] = tmpp[i];
        }
    }

    double final[9] = { 0 };

    for (int i = 0; i < size * size - 1; ++i) {
        for (int j = 0; j < size * size; ++j) {
            if (i == 0) {
                final[j] += real(fBasis[i * size * size + j]);
            }
            else {
                final[j] += vectorV[i - 1] * real(fBasis[i * size * size + j]);
            }
            
        }
    }    
}
