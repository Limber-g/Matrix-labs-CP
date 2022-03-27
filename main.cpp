#include <iostream>
#include <fstream>
#include <locale>

#include "Matrix.h"


int main() {
    std::cin.imbue(std::locale(std::cin.getloc(), new comma_punct()));
    std::cout.imbue(std::locale(std::cout.getloc(), new comma_punct()));

    Symmetr s(5);
    Matrix p = s;

    Matrix nekv(2, 3);
    nekv[0][1] = 1;
    nekv[1][2] = 2;
    nekv[1][1] = 3;



    p[1][3] = 9;
    p[2][4] = -1;
    p[0][5] = -12;

    s[0][3] = 6;
    s[2][4] = 10;
    s[3][1] = -3;

    std::cout << "Transpos nekv\n" << nekv.transpos() << std::endl;
    std::cout << "Our Martix:\n"<< s << std::endl;
    std::cout << "Norm matrix:\n" << s.norm_matrix() << std::endl;
    std::cout << "Trace:\n" << s.trace() << std::endl;
    std::cout << "Rank:\n" << s.rank() << std::endl;
    std::cout << "Determinant:\n" << s.determ() << std::endl;
    std::cout << "Transpose:\n" << s.transpos() << std::endl;
    std::cout << "Adjacent matrix:\n" << s.adj() << std::endl;
    std::cout << "Inverse matrix:\n" << s.reverse() << std::endl;
    std::cout << "Matrix * Inverse:\n" << s * s.reverse() << std::endl;
    std::cout << "S+P:\n" << s+p << std::endl;
    std::cout << "S-P:\n" << s-p << std::endl;
    std::cout << "S*P:\n" << s*p << std::endl;
    std::cout << "Adamar S*P:\n" << s.adamar(p) << std::endl;

    std::ofstream of;
    of.open("out1.txt");
    of << s;
    of.close();
    of.open("out2.txt", std::ios::binary);
    s.outBin(of);
    of.close();

    Matrix tmp(5, 5);
    std::ifstream f;
    f.open("out1.txt");
    f >> tmp;
    std::cout << "Read matrix 1:\n" << tmp << std::endl;
    f.close();
    f.open("out2.txt", std::ios::binary);
    tmp.inBin(f);
    std::cout << "Read matrix 2:\n" << tmp << std::endl;
    f.close();

    Matrix data(32, 12);
    f.open("data.txt");
    f >> data;
    f.close();
    std::cout << "data.txt matrix:\n" << data << std::endl;
    Matrix loading(12, 12);
    f.open("loading.txt");
    f >> loading;
    f.close();
    std::cout << "loading.txt matrix:\n" << loading << std::endl;
    Matrix scores(32, 12);
    f.open("scores.txt");
    f >> scores;
    f.close();
    std::cout << "scores.txt matrix:\n" << scores << std::endl;

    return 0;
}





