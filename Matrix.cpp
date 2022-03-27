//
// Created by 79009 on 12.12.2021.
//

#include "Matrix.h"

char comma_punct::do_decimal_point() const { return ','; }


// Перегрузка опреатора сложения
Matrix Matrix::operator+(const Matrix &plus) const {

    if (this->row != plus.row || this->col != plus.col) {
        throw std::logic_error("Error plus");
    }

    Matrix out_tmp(this->row, this->col);

    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            out_tmp.val[i][j] = plus.val[i][j] + this->val[i][j];
        }
    }

    return out_tmp;
}


// Перегрузка оператора вычитания
Matrix Matrix::operator-(const Matrix &minus) const {

    if (this->row != minus.row || this->col != minus.col) {
        throw std::logic_error("Error minus");
    }

    Matrix out_tmp(this->row, this->col);

    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            out_tmp.val[i][j] = minus.val[i][j] - this->val[i][j];
        }
    }

    return out_tmp;
}


// перегрузка оператора для произведения матрицы на число
Matrix Matrix::operator*(float coef) const {

    Matrix out_tmp(this->row, this->col);

    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            out_tmp.val[i][j] = this->val[i][j] * coef;
        }
    }

    return out_tmp;
}


// Перегрузка оператора умножения для матриц
Matrix Matrix::operator*(const Matrix &mult) const {

    if (this->col != mult.row) {
        throw std::logic_error("Error mult");
    }

    Matrix out_tmp(this->row, mult.col);

    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < mult.col; ++j) {
            float sum_here = 0;

            for (int p = 0; p < this->col; ++p) {
                sum_here += this->val[i][p] * mult.val[p][j];
            }

            out_tmp.val[i][j] = sum_here;
        }
    }

    return out_tmp;
}


//перегрузка оператора вывода
std::ostream &operator<<(std::ostream &s, const Matrix &out_tmp) {
    s.imbue(std::locale(s.getloc(), new comma_punct()));;
    for (int i = 0; i < out_tmp.row; ++i) {
        for (int j = 0; j < out_tmp.col; ++j) {
            s << round(out_tmp.val[i][j]*1000) / 1000 << "\t";
        }
        s << std::endl;
    }
    return s;
}

std::istream &operator>>(std::istream &s, Matrix &in_tmp) {
    s.imbue(std::locale(s.getloc(), new comma_punct()));
    for (int i = 0; i < in_tmp.row; ++i)
        for (int j = 0; j < in_tmp.col; ++j) {
            float v;
            s >> v;
            in_tmp.val[i][j] = v;
        }
    return s;
}

void Matrix::outBin(std::ostream &s) const {
    s << *this;
}

void Matrix::inBin(std::istream &s) {
    s >> *this;
}


// Произведение Адамара
Matrix Matrix::adamar(const Matrix &tmp) const {
    if (this->row != tmp.row || this->col != tmp.col) {
        throw std::logic_error("Error adamar");
    }

    Matrix out_tmp(this->row, this->col);

    for (int i = 0; i < this->row; ++i) {
        for (int j = 0; j < this->col; ++j) {
            out_tmp.val[i][j] = tmp.val[i][j] * this->val[i][j];
        }
    }


    return out_tmp;
}


// След матрицы
float Matrix::trace() const {
    if (col != row) {
        throw std::logic_error("Error. Matrix isn't square");
    }

    float sum = 0;
    for (int i = 0; i < row; ++i) {
        sum += val[i][i];
    }
    return sum;
}


// Приведение матрицы методом Гаусса
void Matrix::gauss() {
    for (int d = 0; d < col && d < row; ++d) {
        int k = d;
        for (; k < row && val[k][d] == 0; ++k) {}
        if (k > d && k < row)
            for (int xx = 0; xx < col; ++xx)
                val[d][xx] = val[d][xx] + val[k][xx] / val[k][d];
        for (int y = d + 1; y < row; ++y)
            if (val[y][d] != 0) {
                float c = val[y][d] / val[d][d];
                for (int xx = 0; xx < col; ++xx)
                    val[y][xx] = val[y][xx] - val[d][xx] * c;
            }
    }
}


// Связная матрица
Matrix Matrix::adj() const {
    Matrix a(col, row);
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j) {
            Matrix dm(row - 1, col - 1);
            for (int ii = 0; ii < row - 1; ++ii)
                for (int jj = 0; jj < col - 1; ++jj)
                    dm.val[ii][jj] = val[ii + (ii >= i)][jj + (jj >= j)];
            a.val[j][i] = dm.determ() * static_cast<float>(1 - 2 * ((i + j) % 2));
        }
    return a;
}

// Детерминант методом Гаусса
float Matrix::determ() const {
    if (col != row) {
        throw std::logic_error("Error. Matrix isn't square");
    }

    Matrix t(row, col);
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            t.val[i][j] = val[i][j];

    t.gauss();
    float r = 1;
    for (int d = 0; d < row && d < col; ++d)
        r *= t.val[d][d];
    return r;
}


// Скалярное произведение векторов
float Vector::scalar(const Vector &tmp) const {
    if (this->col != tmp.col) {
        throw std::logic_error("Error scalar");
    }

    float result = 0;
    for (int i = 0; i < col; ++i) {
        result += this->val[0][i] * tmp.val[0][i];
    }
    return result;
}


// Евклидова норма вектора
float Vector::norm_evklid_vector() const {
    float sum = 0;
    for (int i = 0; i < col; ++i) {
        sum += val[0][i] * val[0][i];
    }

    return sqrtf(sum);
}


// максимальная норма вектора
float Vector::norm_max_vector() const {
    float max = 0;
    for (int i = 0; i < col; ++i) {
        if (i == 0 || (fabsf(val[0][i]) > max))
            max = fabsf(val[0][i]);
    }

    return max;
}


// норма матрицы
float Matrix::norm_matrix() const {
    if (col != row) {
        throw std::logic_error("Error. Matrix isn't square");
    }

    float sum = 0;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            sum += val[i][j] * val[i][j];
        }
    }

    return sqrtf(sum);
}


// угол между векторами
float Vector::corner(const Vector &tmp) const {
    double cos = this->scalar(tmp) / (this->norm_evklid_vector() * tmp.norm_evklid_vector());
    auto result = static_cast<float>(acos(cos) * 180.0 / M_PI);
    return result;
}


// транспонирование матрицы
Matrix Matrix::transpos() const {
    Matrix tmp(this->col, this->row);

    for (int i = 0; i < tmp.row; ++i)
        for (int j = 0; j < tmp.col; ++j)
            tmp.val[i][j] = this->val[j][i];

    return tmp;
}


// Ранг матрицы методом Гауса
int Matrix::rank() const {
    Matrix t(row, col);
    for (int i = 0; i < row; ++i)
        for (int j = 0; j < col; ++j)
            t.val[i][j] = val[i][j];
    t.gauss();

    int r = 0;
    for (int d = 0; d < row && d < col; ++d)
        r += (t.val[d][d] != 0);
    return r;
}


// Обратная матрица
Matrix Matrix::reverse() const {
    float d = determ();
    if (d == 0)
        throw std::logic_error("Determinant is zero.");

    Matrix a = adj();
    for (int i = 0; i < a.row; ++i)
        for (int j = 0; j < a.col; ++j)
            a.val[i][j] /= d;
    return a;
}


// Чтобы переопределить какие-то значения matrix[i][j] = a;
std::vector<float> &Matrix::operator[](int i) {
    return val[i];
}

