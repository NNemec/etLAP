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

#ifndef _ETLAP_SIGMA_H_
#define _ETLAP_SIGMA_H_

#include "Common.h"
#include "Costs.h"
#include <complex>
#include <cassert>

namespace etLAP {

class TagSigma;

template <> struct Costs<TagSigma> { enum { c = 0 }; };

template <int N,typename REAL>
class Matrix<std::complex<REAL>,N,N,TagSigma>
: public Common<Matrix<std::complex<REAL>,N,N,TagSigma> > {
    typedef TagSigma E;
    typedef std::complex<REAL> COMPLEX;
    typedef COMPLEX T;
    typedef Matrix<T,N,N,E> SAME;

    enum { unit, diag, real, imag } type;
    int a;
    int b;

  public:
    typedef SAME Ref_t;

    Matrix<T,N,N,E>(int idx) {
        assert(idx < N*N);
        int x;
        if(idx == 0)
            type = unit;
        else {
            for(a = 2;a*a <= idx;a++);
            x = a*a - idx - 1; // x >= 0
            a--;
            if(x == 0)
                type = diag;
            else if(x % 2) {
                type = imag;
                b = a-((x+1)/2); // b < a
            } else {
                type = real;
                b = a-(x/2); // b < a
            };
        };
//std::cout << "idx: " << idx << " type: " << type << " x: " << x << " a: " << a << " b: " << b << "\n";
    };

    const T operator() (int r,int c) const {
        switch(type) {
        case unit:
            return r==c?COMPLEX(sqrt(2.0/N)):COMPLEX(0.0);
        case diag:
            if(r!=c || r>a)
                return COMPLEX(0.0);
            if(r==a)
                return COMPLEX(-sqrt(2.0*a/(a+1)));
            else
                return COMPLEX(sqrt(2.0/((a+1)*a)));
        case real:
            if(r==a && c==b)
                return COMPLEX(1.0);
            if(r==b && c==a)
                return COMPLEX(1.0);
            return COMPLEX(0.0);
        case imag:
            if(r==a && c==b)
                return COMPLEX(0.0,   1.0);
            if(r==b && c==a)
                return COMPLEX(0.0, - 1.0);
            return COMPLEX(0.0);
        default:
            assert(false);
            throw 0;
        };
    };

    int rows() const { return N; };
    int cols() const { return N; };
};

template <int N,typename REAL>
inline Matrix<std::complex<REAL>,N,N,TagSigma> Sigma(int idx) {
    return Matrix<std::complex<REAL>,N,N,TagSigma>(idx);
};

}; // namespace etLAP

#endif // _ETLAP_MATRIX_H_
