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

#ifndef _ETLAP_SU_N_H_
#define _ETLAP_SU_N_H_

#include "Common.h"
#include "Costs.h"
#include <complex>
#include <cassert>

namespace etLAP {

class TagSUNgenerator;

template <> struct Costs<TagSUNgenerator> { enum { c = 0 }; };

template <int N,typename REAL>
class Matrix<std::complex<REAL>,N,N,TagSUNgenerator>
: public Common<Matrix<std::complex<REAL>,N,N,TagSUNgenerator> > {
    typedef TagSUNgenerator E;
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
            return r==c?COMPLEX(sqrt(0.5/N)):COMPLEX(0.0);
        case diag:
            if(r!=c || r>a)
                return COMPLEX(0.0);
            if(r==a)
                return COMPLEX(-sqrt(0.5*a/(a+1)));
            else
                return COMPLEX(sqrt(0.5/((a+1)*a)));
        case real:
            if(r==a && c==b)
                return COMPLEX(0.5);
            if(r==b && c==a)
                return COMPLEX(0.5);
            return COMPLEX(0.0);
        case imag:
            if(r==a && c==b)
                return COMPLEX(0.0,   0.5);
            if(r==b && c==a)
                return COMPLEX(0.0, - 0.5);
            return COMPLEX(0.0);
        default:
            assert(false);
            throw 0;
        };
    };

    int rows() const { return N; };
    int cols() const { return N; };

    bool is_locked() const { return false; };
};

template <int N,typename REAL>
inline Matrix<std::complex<REAL>,N,N,TagSUNgenerator> SUNgenerator(int idx) {
    return Matrix<std::complex<REAL>,N,N,TagSUNgenerator>(idx);
};


template <int N,typename REAL>
class SUNstructure_values {
    bool _nonzero[N*N*N*N*N*N];
    REAL _data[N*N*N*N*N*N];
    REAL &data(uint a,uint b,uint c) { return _data[a*(N*N*N*N) + b*(N*N) + c]; }
    REAL data(uint a,uint b,uint c) const { return _data[a*(N*N*N*N) + b*(N*N) + c]; }
    bool &nonzero(uint a,uint b,uint c) { return _nonzero[a*(N*N*N*N) + b*(N*N) + c]; }
    bool nonzero(uint a,uint b,uint c) const { return _nonzero[a*(N*N*N*N) + b*(N*N) + c]; }

  public:
    SUNstructure_values() {
        for(uint i=0;i<N*N*N*N*N*N;i++) {
            _nonzero[i] = false;
            _data[i] = 0.0;
        };

        bool is_diag[N*N];
        for(uint i=0;i<N*N;i++)
            is_diag[i] = false;
        for(uint i=1;i<=N;i++)
            is_diag[i*i-1] = true;

        for(uint a=1;a<N*N-2;a++)
        for(uint b=a+1;b<N*N-1;b++)
        for(uint c=b+1;c<N*N;c++)
        if(!((is_diag[a] && is_diag[b]) || (is_diag[b] && is_diag[c]) || (is_diag[c] && is_diag[a]))) {
            REAL val = 2.0 * trace(imag(commute(SUNgenerator<N,REAL>(a),SUNgenerator<N,REAL>(b)) * SUNgenerator<N,REAL>(c)));
            if(val != 0.0) {
                data(a,b,c) = val;
                data(b,c,a) = val;
                data(c,a,b) = val;
                data(c,b,a) = -val;
                data(b,a,c) = -val;
                data(a,c,b) = -val;
                nonzero(a,b,c) = true;
                nonzero(b,c,a) = true;
                nonzero(c,a,b) = true;
                nonzero(c,b,a) = true;
                nonzero(b,a,c) = true;
                nonzero(a,c,b) = true;
            };
        };
    };

    REAL value(uint a,uint b,uint c) const {
        assert(a > 0 && a < N*N);
        assert(b > 0 && b < N*N);
        assert(c > 0 && c < N*N);
        return data(a,b,c);
    };

    bool is_nonzero(uint a,uint b,uint c) {
        assert(a > 0 && a < N*N);
        assert(b > 0 && b < N*N);
        assert(c > 0 && c < N*N);
        return nonzero(a,b,c);
    };
};

template <typename REAL>
struct SUNstructure_values<1,REAL> {
    REAL value(uint a,uint b,uint c) {
        assert(a == 0);
        assert(b == 0);
        assert(c == 0);
        return 0.0;
    };

    bool is_nonzero(uint a,uint b,uint c) {
        assert(a == 0);
        assert(b == 0);
        assert(c == 0);
        return false;
    };
};

template <typename REAL>
struct SUNstructure_values<2,REAL> {
    REAL value(uint a,uint b,uint c) {
        assert(a > 0 && a < 4);
        assert(b > 0 && b < 4);
        assert(c > 0 && c < 4);
        if(a==b || b==c || c==a)
            return 0.0;
        return REAL((int)((4+b-a)%3)-1);
    };

    bool is_nonzero(uint a,uint b,uint c) {
        assert(a > 0 && a < 4);
        assert(b > 0 && b < 4);
        assert(c > 0 && c < 4);
        return !(a==b || b==c || c==a);
    };
};

template <int N,typename REAL>
inline SUNstructure_values<N,REAL> &SUNstructure_val() {
    static SUNstructure_values<N,REAL> value;
    return value;
};

template <int N,typename REAL>
inline REAL SUNstructure(uint a,uint b,uint c) {
    return SUNstructure_val<N,REAL>().value(a,b,c);
};

template <int N,typename REAL>
inline bool SUNstructure_nonzero(uint a,uint b,uint c) {
    return SUNstructure_val<N,REAL>().is_nonzero(a,b,c);
};

}; // namespace etLAP

#endif // _ETLAP_SU_N_H_
