#ifndef LABS_1_4_SECOND_MODULE_MATRIX_H
#define LABS_1_4_SECOND_MODULE_MATRIX_H


#include <vector>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <iostream>
#include <stdexcept>  // for errors
#include <cmath>

class comma_punct : public std::numpunct<char> {
protected:
    char do_decimal_point() const override;
};

class Matrix {

protected:
    std::vector<std::vector<float>> val;
    int row;
    int col;

public:

    Matrix() {        // для случайных нажатий
        row = col = 0;
    }

    Matrix(int m, int n) {
        row = m;
        col = n;

        for (int i = 0; i < m; ++i) {
            std::vector<float> tmp;
            val.push_back(tmp);
            for (int j = 0; j < n; ++j) {
                val[i].push_back(0);
            }
        }
    }

    std::vector<float> &operator[](int i);

    Matrix operator+(const Matrix &plus) const;

    Matrix operator-(const Matrix &minus) const;

    Matrix operator*(float coef) const;

    Matrix operator*(const Matrix &mult) const;

    Matrix adamar(const Matrix &second) const;

    friend std::ostream &operator<<(std::ostream &s, const Matrix &out_tmp);

    friend std::istream &operator>>(std::istream &s, Matrix &in_tmp);

    void outBin(std::ostream &s) const;

    void inBin(std::istream &s);

    float trace() const;

    float determ() const;

    float norm_matrix() const;

    int rank() const;

    Matrix transpos() const;

    Matrix reverse() const;

    void gauss();

    Matrix adj() const;
};


class Vector : public Matrix {
public:

    Vector() {
        row = col = 0;
    }

    Vector(int a) : Matrix(1, a) {}

    float scalar(const Vector &tmp) const;

    float norm_evklid_vector() const;

    float norm_max_vector() const;

    float corner(const Vector &tmp) const;


};

// подкласс квадратная матрица
class Square : public Matrix {        // все подклассы для квадратных матриц,
public:                            // поэтому создадим такой подкласс

    Square() {
        row = col = 0;
    }

    explicit Square(int a) : Matrix(a, a) {}

};


// подкласс единичная матрица
class E : public Square {
public:

    E() {
        row = col = 0;
    }

    explicit E(int a) : Square(a) {

        for (int i = 0; i < row; ++i) {    // нулями остальное заполнили
            val[i][i] = 1;                // еще раньше, тк класс наследуется
        }
    }

};


// подкласс диагональная матрица
class Diagonal : public Square {
public:

    Diagonal() {
        col = row = 0;
    }

    explicit Diagonal(int a) : Square(a) {
        for (int i = 0; i < row; ++i) {
            val[i][i] = static_cast<float>(i + 1);
        }
    }


};


// подкласс верхняя треугольная матрица
class Top_tri : public Square {
public:

    Top_tri() {
        col = row = 0;
    }

    explicit Top_tri(int a) : Square(a) {
        for (int i = 0; i < row; ++i) {
            for (int j = i; j < row; ++j) {
                val[i][j] = 1;
            }
        }
    }

};


// подкласс нижняя треугольная матрица
class Low_tri : public Square {
public:

    Low_tri() {
        col = row = 0;
    }

    explicit Low_tri(int a) : Square(a) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j <= i; ++j) {
                val[i][j] = 1;
            }
        }
    }
};


// подкласс симметричная матрица
class Symmetr : public Square {
public:

    Symmetr() {
        col = row = 0;
    }

    explicit Symmetr(int a) : Square(a) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                val[i][j] = static_cast<float>(i + j);
            }
        }
    }

};


#endif //LABS_1_4_SECOND_MODULE_MATRIX_H
