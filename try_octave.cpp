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
#include <etLAP/octave.h>
#include <etLAP/Output.h>
#include <etLAP/SU_N.h>

#include <iostream>
#include <string>
#include <complex>
#include <typeinfo>

using namespace etLAP;
using namespace std;

const int NDIM = 5;

int main(int argc,char **argv) {
    complex<double> I = complex<double>(0,1);
    double zero = argc==-17?1.0:0.0;
    complex<double> res;
    
    etLAP::Matrix<complex<double>,NDIM,NDIM> B = 0.5 * SUNgenerator<NDIM,double>(2) + 0.5 * SUNgenerator<NDIM,double>(7);
    
    std::cerr << B << "\n";
    std::cerr << octave_exp(-I*B) << "\n";
    std::cerr << exp(-I*B) << "\n";
    
    std::cerr << "start\n";

    res = zero;    
    for(int i=1;i<10000;i++) {
        I += zero;
        res += octave_exp(-I*B)(0,0);
    };
    
    std::cerr << "middle: " << res << "\n";
    
    res = zero;    
    for(int i=1;i<10000;i++) {
        I += zero;
        res += exp(-I*B)(0,0);
    };
    
    std::cerr << "end: " << res << "\n";
}
