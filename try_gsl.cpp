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
#include <etLAP/gsl.h>
#include <etLAP/Output.h>

#include <iostream>
#include <string>
#include <complex>
#include <typeinfo>

using namespace etLAP;
using namespace std;

int main() {
    Matrix<double,3,3> A = 2;
/*
    A(0,0) = 4;
    A(1,1) = 16;
    A(2,2) = 25;
    A(1,0) = A(0,1) = 2;
    A(2,0) = A(0,2) = 3;
    A(2,1) = A(1,2) = 4;
*/
//    Matrix<double,3,3> B;

    cout << A << "\n";

    for(int i=0; i<10000000; i++)
        A = gsl_inv_sqrt(A);

    cout << A << "\n";
//    cout << B << "\n";
//    cout << B * A * B << "\n";
}
