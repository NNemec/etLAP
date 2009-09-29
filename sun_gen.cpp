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

#include <etLAP/SU_N.h>
#include <etLAP/Operator.h>
#include <etLAP/Output.h>

#define Pi 3.1415926535897932384626433832795029

using namespace etLAP;

typedef double REAL;
typedef std::complex<REAL> COMPLEX;

const int NCOL = 2;

const uint MIN_A = (NCOL==1)?0:1;
typedef Vector<REAL,NCOL*NCOL-MIN_A> suN_VECTOR;
typedef Matrix<COMPLEX,NCOL,NCOL> SUN_MATRIX;

int main() {
    for(int i=0;i<NCOL*NCOL;i++) {
        std::cout << SUNgenerator<NCOL,REAL>(i) << "\n";
        std::cout << trace(SUNgenerator<NCOL,REAL>(i) * SUNgenerator<NCOL,REAL>(i)) << "\n";
    };
    
    suN_VECTOR v;
    v.clear(); v(0) = 1.0;
    for(int i=0;i<=16;i++)
	std::cout << SUN_from_suN(2*Pi*i/16.0 * v,*((SUN_MATRIX*)0)) << "\n";
}
