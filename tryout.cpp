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

#define etLAP_DEFAULT_ENGINE Packed

#include <etLAP/Vector.h>
#include <etLAP/Matrix.h>
#include <etLAP/Operator.h>
#include <etLAP/Output.h>
#include <etLAP/Tupel.h>

#include <iostream>
#include <string>
#include <complex>
#include <typeinfo>

using namespace etLAP;
using namespace std;

int main() {
    Vector<double,2> L = tupel(1,1);
    Vector<double,2> G = tupel(1,1);
    Matrix<double,2,2> Q = 1;
    Vector<double,2> res = Q*buf(G+Q*L);
    cout << res << "\n";

/*
    cout << sizeof(Matrix<char,1,1,Packed>) << "\n";
    Vector<Matrix<complex<double>,3,3>,4> v;
    cout << sizeof(v) << " " << typeid(v).name() << "\n";
    cout << &v << "\n";
    cout << sizeof(v(0)) << "\n";
    cout << &(v(0)) << "\n";
    cout << &(v(1)) << "\n";
//    cout << &(v(0)(0,0)) << "\n";
//    Vector<int,3> v = tupel(3,4,5);
*/
};
