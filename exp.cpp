/*
 *    This file is part of etLAP - the Expression Templated Linear Algebra Package
 *    Copyright (C) 2002 by Norbert Nemec
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <etLAP/Matrix.h>
#include <etLAP/Operator.h>

#include <etLAP/Output.h>

#include <iostream>
#include <string>
#include <ctime>

using namespace etLAP;
using namespace std;

const int D = 3;
const int MAX = 10000000;

inline void fixed() {
    clock_t cpptime;

    Matrix<double,D,D> a = 3;
    for(int i=0;i<D-1;i++)
        a(i,i+1)=1;
    for(int i=0;i<D-1;i++)
        a(i+1,i)=1;

    Matrix<double,D,D> exp = 1;
    Matrix<double,D,D> tmp = 1;

    cpptime = clock();

    for(int n=1;n<MAX;n++) {
#if ASSIGN_POLICY == AP_manual
        tmp = buf(tmp * a / n);
#else
        tmp = tmp * a / n;
#endif
        exp += tmp;
    };

    cpptime = clock() - cpptime;

    cout << "FIX"
         << " D: " << setw(3) << D
         << " N: " << setw(3) << MAX
         << " tr: " << setw(10) << setprecision(4) << setiosflags(ios::fixed) << trace(exp)
         << " time: " << setw(10) << setprecision(4) << setiosflags(ios::fixed) << double(cpptime)/double(CLOCKS_PER_SEC) << "\n";
};

inline void var() {
    clock_t cpptime;

    Matrix<double> a(D,D); a = 3;
    for(int i=0;i<D-1;i++)
        a(i,i+1)=1;
    for(int i=0;i<D-1;i++)
        a(i+1,i)=1;

    Matrix<double> exp(D,D); exp = 1;
    Matrix<double> tmp(D,D); tmp = 1;

    cpptime = clock();

    for(int n=1;n<MAX;n++) {
#if ASSIGN_POLICY == AP_manual
        tmp = buf(tmp * a / n);
#else
        tmp = tmp * a / n;
#endif
        exp += tmp;
    };

    cpptime = clock() - cpptime;

    cout << "VAR"
         << " D: " << setw(3) << D
         << " N: " << setw(3) << MAX
         << " tr: " << setw(10) << setprecision(4) << setiosflags(ios::fixed) << trace(exp)
         << " time: " << setw(10) << setprecision(4) << setiosflags(ios::fixed) << double(cpptime)/double(CLOCKS_PER_SEC) << "\n";
};

int main() {
    fixed();
    var();
};
