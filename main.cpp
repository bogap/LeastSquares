//bogapova aliia

#include <iomanip>
#include <valarray>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cstring>
#include <random>

using namespace std;

class ColumnVector {
private:
    int n;
    vector<double> cvector;
public:
    ColumnVector(int _n) {
        n = _n;
        cvector = vector<double>(_n, 0);
    }

    ColumnVector(ColumnVector &other) {
        n = other.n;
        cvector = other.cvector;
    }

    int get_n() const { return n; }

    double &operator()(int i) { return cvector[i]; }

    const double &operator()(int i) const { return cvector[i]; }

    ColumnVector &operator=(const ColumnVector &other) {
        if (this != &other) {
            n = other.n;
            cvector = other.cvector;
        }
        return *this;
    }

    ColumnVector operator+(const ColumnVector &other) const {
        ColumnVector result(n);
        if (n != other.n) {
            return result;
        }
        for (int i = 0; i < n; ++i) {

            result(i) = (*this)(i) + other(i);
        }
        return result;
    }

    ColumnVector operator-(const ColumnVector &other) const {
        ColumnVector result(n);
        if (n != other.n) {
            return result;
        }
        for (int i = 0; i < n; ++i) {
            result(i) = (*this)(i) - other(i);
        }
        return result;
    }

    ColumnVector operator*(const double other) {
        ColumnVector result(n);
        for (int i = 0; i < n; ++i) {
            result(i) = (*this)(i) * other;
        }
        return result;
    }

    double operator*(const ColumnVector &other) {
        double result = 0;
        for (int i = 0; i < n; ++i) {
            result += (*this)(i) * other(i);
        }
        return result;
    }

    friend istream &operator>>(istream &inp, ColumnVector &cvec) {
        for (int i = 0; i < cvec.n; ++i) {
            inp >> cvec(i);
        }
        return inp;
    }

    friend ostream &operator<<(ostream &out, ColumnVector &cvec) {
        for (int i = 0; i < cvec.n; ++i) {
            if (to_string(cvec(i)).substr(0, 7) == "-0.0000") {
                out << fixed << setprecision(4) << 0.0000;
            } else {
                out << fixed << setprecision(4) << cvec(i);
            }
            out << "\n";
        }
        return out;
    }

    double norm() {
        double result = 0;
        for (int i = 0; i < n; ++i) {
            result += (*this)(i) * (*this)(i);
        }
        return sqrt(result);
    }
};


class Matrix {
private:
    int n, m;
    vector<vector<double>> matrix;
public:
    Matrix(int _n, int _m) {
        n = _n;
        m = _m;
        matrix = vector<vector<double> >(_n, vector<double>(_m, 0));
    }

    Matrix(Matrix &other) {
        n = other.n;
        m = other.m;
        matrix = other.matrix;
    }

    string s;

    int get_n() const { return n; }

    int get_m() const { return m; }

    double &operator()(int i, int j) { return matrix[i][j]; }

    const double &operator()(int i, int j) const { return matrix[i][j]; }

    Matrix &operator=(const Matrix &other) {
        if (this != &other) {
            n = other.n;
            m = other.m;
            matrix = other.matrix;
        }
        return *this;
    }

    Matrix operator+(const Matrix &other) {
        Matrix result(n, m);
        if (n != other.n || m != other.m) {
            return result;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix &other) {
        Matrix result(n, m);
        if (n != other.n || m != other.m) {
            return result;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &other) {
        Matrix result(n, other.m);
        if (m != other.n) {
            return result;
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < other.m; ++j) {
                for (int k = 0; k < m; ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }

    Matrix operator*(const ColumnVector &other) {
        Matrix result(n, 1);
        if (m != other.get_n()) {
            return result;
        }

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int k = 0; k < m; ++k) {
                sum += (*this)(i, k) * other(k);
            }
            result(i, 0) = sum;
        }
        return result;
    }

    Matrix transpose() {
        Matrix result(m, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    friend istream &operator>>(istream &inp, Matrix &matrixa) {
        for (int i = 0; i < matrixa.n; ++i) {
            for (int j = 0; j < matrixa.m; ++j) {
                inp >> matrixa(i, j);
            }
        }
        return inp;
    }

    friend ostream &operator<<(ostream &out, Matrix &matrixa) {
        if (matrixa.n == 0) {
            out << ("Error: the dimensional problem occurred\n");
            return out;
        }
        for (int i = 0; i < matrixa.n; ++i) {
            for (int j = 0; j < matrixa.m; ++j) {
                if (to_string(matrixa(i, j)).substr(0, 7) == "-0.0000") {
                    out << fixed << setprecision(4) << 0.0000 << " ";
                } else out << fixed << setprecision(4) << matrixa(i, j) << " ";
            }
            out << "\n";
        }
        return out;
    }
};

class SquareMatrix : public Matrix {
public:

    explicit SquareMatrix(int n) : Matrix(n, n) {
    }

    SquareMatrix(const SquareMatrix &other) : Matrix((Matrix &) other) {}


    SquareMatrix &operator=(const SquareMatrix &other) {
        if (this != &other) {
            Matrix::operator=(other);
        }
        return *this;
    }

    SquareMatrix operator+(const SquareMatrix &other) {
        SquareMatrix result(get_n());
        if (get_n() != other.get_n() || get_m() != other.get_m()) {
            return result;
        }
        for (int i = 0; i < get_n(); ++i) {
            for (int j = 0; j < get_m(); ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }

    SquareMatrix operator-(const SquareMatrix &other) {
        SquareMatrix result(get_n());
        if (get_n() != other.get_n() || get_m() != other.get_m()) {
            return result;
        }
        for (int i = 0; i < get_n(); ++i) {
            for (int j = 0; j < get_m(); ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }

    SquareMatrix operator*(const SquareMatrix &other) {
        SquareMatrix result(get_n());
        if (get_m() != other.get_n()) {
            return result;
        }
        for (int i = 0; i < get_n(); ++i) {
            for (int j = 0; j < other.get_m(); ++j) {
                for (int k = 0; k < get_m(); ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return result;
    }

    SquareMatrix transpose() {
        SquareMatrix result(get_n());
        for (int i = 0; i < get_n(); i++) {
            for (int j = 0; j < get_m(); j++) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    double det() {
        int step = 1;
        double det_val = 1;
        for (int k = 0; k < (*this).get_n(); k++) {
            int max_row = k;
            double max_val = abs((*this)(k, k));
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double abs_val = abs((*this)(i, k));
                if (abs_val > max_val) {
                    max_val = abs_val;
                    max_row = i;
                }
            }
            if (max_row != k) {
                swap(max_row, k);
                cout << "step #" << step << ": permutation\n";
                step++;
                //cout << (*this);
                det_val *= -1;
            }
            det_val *= (*this)(k, k);
            if (abs(det_val) == 0) {
                return 0;
            }
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double factor = (*this)(i, k) / (*this)(k, k);
                for (int j = k + 1; j < (*this).get_n(); j++) {
                    (*this)(i, j) -= factor * (*this)(k, j);
                }
                (*this)(i, k) = 0;
                cout << "step #" << step << ": elimination\n";

                cout << (*this);
                step++;
            }
        }
        cout << "result:\n";
        return det_val;
    }

    SquareMatrix inverse() {
        SquareMatrix aug((*this).get_n());
        for (int i = 0; i < (*this).get_n(); ++i) {
            for (int j = 0; j < (*this).get_n(); ++j) {
                aug(i, j) = 0.00;
                if (i == j) {
                    aug(i, j) = 1.00;
                }
            }
        }

        int step = 1;
        double det_val = 1;
        for (int k = 0; k < (*this).get_n(); k++) {
            int max_row = k;
            double max_val = abs((*this)(k, k));
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double abs_val = abs((*this)(i, k));
                if (abs_val > max_val) {
                    max_val = abs_val;
                    max_row = i;
                }
            }
            if (max_row != k) {
                double tmp;
                for (int m = 0; m < (*this).get_n(); m++) {
                    tmp = (*this)(max_row, m);
                    (*this)(max_row, m) = (*this)(k, m);
                    (*this)(k, m) = tmp;

                    tmp = (aug)(max_row, m);
                    (aug)(max_row, m) = (aug)(k, m);
                    (aug)(k, m) = tmp;
                }
                step++;
                det_val *= -1;
            }

            det_val *= (*this)(k, k);
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double factor = (*this)(i, k) / (*this)(k, k);
                if (abs(factor) == 0 || isnan(factor)) {
                    continue;
                }
                for (int j = 0; j < (*this).get_n(); j++) {
                    (*this)(i, j) -= factor * (*this)(k, j);
                }
                for (int j = 0; j < (aug).get_n(); j++) {
                    (aug)(i, j) -= factor * (aug)(k, j);
                }
                step++;
            }
        }
        int n = (*this).get_n() - 1;
        for (int k = 0; k < (*this).get_n(); k++) {
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double factor = (*this)(n - i, n - k) / (*this)(n - k, n - k);
                if (abs((*this)(n - k, n - k)) == 0) {
                    factor = 0;
                }
                if (abs(factor) == 0 || isnan(factor)) {
                    continue;
                }
                for (int j = 0; j < (*this).get_n(); j++) {
                    (*this)(n - i, j) -= factor * (*this)(n - k, j);
                }
                for (int j = 0; j < (aug).get_n(); j++) {
                    (aug)(n - i, j) -= factor * (aug)(n - k, j);
                }
                step++;
            }
        }
        for (int k = 0; k < (*this).get_n(); k++) {
            double pivot = (*this)(k, k);
            for (int j = 0; j < (*this).get_n(); j++) {
                if (abs(pivot) == 0 || isnan(pivot)) {
                    continue;
                }
                (*this)(k, j) /= pivot;
            }
            for (int j = 0; j < (aug).get_n(); j++) {
                if (abs(pivot) == 0 || isnan(pivot)) {
                    continue;
                }
                (aug)(k, j) /= pivot;
            }
            if (abs(pivot) == 0 || isnan(pivot)) {
                continue;
            }
        }
        return aug;
    }

    SquareMatrix solve(ColumnVector cvec) {
        ColumnVector vec = cvec;
        SquareMatrix aug((*this).get_n());
        for (int i = 0; i < (*this).get_n(); ++i) {
            for (int j = 0; j < (*this).get_n(); ++j) {
                aug(i, j) = 0.00;
                if (i == j) {
                    aug(i, j) = 1.00;
                }
            }
        }

        cout << "step #0:\n";
        cout << (*this);
        cout << vec;
        int step = 1;
        double det_val = 1;
        for (int k = 0; k < (*this).get_n(); k++) {
            int max_row = k;
            double max_val = abs((*this)(k, k));
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double abs_val = abs((*this)(i, k));
                if (abs_val > max_val) {
                    max_val = abs_val;
                    max_row = i;
                }
            }
            if (max_row != k) {
                double tmp;
                for (int m = 0; m < (*this).get_n(); m++) {
                    tmp = (*this)(max_row, m);
                    (*this)(max_row, m) = (*this)(k, m);
                    (*this)(k, m) = tmp;

                    tmp = (aug)(max_row, m);
                    (aug)(max_row, m) = (aug)(k, m);
                    (aug)(k, m) = tmp;
                }
                tmp = vec(max_row);
                vec(max_row) = vec(k);
                vec(k) = tmp;
                cout << "step #" << step << ": permutation\n";
                step++;
                cout << (*this);
                cout << vec;
                det_val *= -1;
            }

            det_val *= (*this)(k, k);
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double factor = (*this)(i, k) / (*this)(k, k);
                if (abs(factor) == 0 || isnan(factor)) {
                    continue;
                }
                for (int j = 0; j < (*this).get_n(); j++) {
                    (*this)(i, j) -= factor * (*this)(k, j);
                }
                for (int j = 0; j < (aug).get_n(); j++) {
                    (aug)(i, j) -= factor * (aug)(k, j);
                }
                vec(i) -= factor * vec(k);
                cout << "step #" << step << ": elimination\n";
                cout << (*this);
                cout << vec;
                step++;
            }
        }

        int n = (*this).get_n() - 1;
        for (int k = 0; k < (*this).get_n(); k++) {
            for (int i = k + 1; i < (*this).get_n(); i++) {
                double factor = (*this)(n - i, n - k) / (*this)(n - k, n - k);
                if (abs((*this)(n - k, n - k)) == 0) {
                    factor = 0;
                }
                if (abs(factor) == 0 || isnan(factor)) {
                    continue;
                }
                for (int j = 0; j < (*this).get_n(); j++) {
                    (*this)(n - i, j) -= factor * (*this)(n - k, j);
                }
                for (int j = 0; j < (aug).get_n(); j++) {
                    (aug)(n - i, j) -= factor * (aug)(n - k, j);
                }
                vec(n - i) -= factor * vec(n - k);
                cout << "step #" << step << ": elimination\n";
                cout << (*this);
                cout << vec;
                step++;
            }
        }
        cout << "Diagonal normalization:\n";
        for (int k = 0; k < (*this).get_n(); k++) {
            double pivot = (*this)(k, k);
            for (int j = 0; j < (*this).get_n(); j++) {
                if (abs(pivot) == 0 || isnan(pivot)) {
                    continue;
                }
                (*this)(k, j) /= pivot;
            }
            for (int j = 0; j < (aug).get_n(); j++) {
                if (abs(pivot) == 0 || isnan(pivot)) {
                    continue;
                }
                (aug)(k, j) /= pivot;
            }
            if (abs(pivot) == 0 || isnan(pivot)) {
                continue;
            }
            vec(k) /= pivot;
        }
        cout << (*this);
        cout << vec;
        cout << "result:\n";

        for (int i = 0; i < vec.get_n(); ++i) {
            if (abs((*this)(i, i)) == 0) {
                if (abs(vec(i)) == 0) {
                    cout << "INF\n";
                    return aug;
                }
                cout << "NO\n";
                return aug;
            }
        }
        for (int i = 0; i < vec.get_n(); ++i) {
            if (vec(i) != cvec(i)) {
                cout << vec;
                break;
            }
        }
        return aug;
    }

    void print_aug(SquareMatrix aug) const {
        for (int i = 0; i < (*this).get_n(); ++i) {
            for (int j = 0; j < 2 * (*this).get_n(); ++j) {
                if (j < (*this).get_n()) {
                    cout << std::fixed << setprecision(2) << (*this)(i, j) << " ";
                } else {
                    cout << std::fixed << setprecision(2) << aug(i, j - (*this).get_n()) << " ";
                }
            }
            cout << endl;
        }
    }

};

class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            (*this)(i, i) = 1;
        }
    }
};

class EliminationMatrix : public SquareMatrix {
public:

    EliminationMatrix(int n, int i, int j, int c) : SquareMatrix(n) {

        for (int k = 0; k < n; k++) {
            (*this)(k, k) = 1;
        }

        (*this)(i, j) = -c;
    }
};

class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix(int n, int i, int j) : SquareMatrix(n) {
        for (int k = 0; k < n; k++) {
            if (k == i) {
                (*this)(j, k) = 1;
            } else if (k == j) {
                (*this)(i, k) = 1;
            } else {
                (*this)(k, k) = 1;
            }
        }
    }
};

Matrix leastSquareApproximation(Matrix A, ColumnVector b) {
    cout << "A:\n" << A;
    cout << "A_T*A:\n";
    Matrix B(A.get_m(), A.get_n());
    B = A.transpose();
    Matrix D(A.get_n(), A.get_m());
    D = B * A;
    cout << D;
    cout << "(A_T*A)^-1:\n";
    Matrix C(D.get_n(), D.get_n());
    C = (*(SquareMatrix *) (&D)).inverse();
    cout << C;
    cout << "A_T*b:\n";
    Matrix Atb(B.get_n(), 1);
    Atb = B * b;
    cout << Atb;
    cout << "x~:\n";
    Matrix E(C.get_n(), Atb.get_m());
    E = C * Atb;
    cout << E;
    return E;
}

#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main() {
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    default_random_engine _random{std::random_device{}()};
    uniform_real_distribution<double> interval(0.0, 10.0);
    int n;
    n = interval(_random);
    vector<double> vectorX;
    ColumnVector vectorY(n);
    double x, y;

    for (int i = 0; i < n; ++i) {
        x = interval(_random);
        y = interval(_random);
        vectorX.push_back(x);
        vectorY(i) = y;
    }

    int degree;
    degree = interval(_random);
    Matrix A(n, degree + 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < degree + 1; ++j) {
            A(i, j) = pow(vectorX[i], j);
        }
    }
    Matrix X(n, 1);
    X = leastSquareApproximation(A, vectorY);
    string tmp;
    tmp += "plot [0:10] [0:10] ";
    for (int i = 0; i <= degree; i++) {
        tmp += "%lf*x**";
        tmp += to_string(i);
        if (i == degree)break;
        tmp += " + ";
    }
    tmp += ", '-' with lines\n";
    double vectr[X.get_n()];
    for (int i = 0; i < X.get_n(); ++i) {
        vectr[i] = X(i, 0);
    }

    fprintf(pipe, tmp.c_str(), vectr[0], vectr[1],
            vectr[2], vectr[3], vectr[4], vectr[5], vectr[6], vectr[7], vectr[8], vectr[9], vectr[10]);

    fprintf(pipe, "e\n");
    fflush(pipe);
#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return 0;
}
