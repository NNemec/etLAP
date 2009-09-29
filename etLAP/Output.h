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

#ifndef _ETLAP_OUTPUT_H_
#define _ETLAP_OUTPUT_H_

#include "Vector.h"
#include "Matrix.h"

#include <iostream>
#include <iomanip>

namespace etLAP {

using std::ostream;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::ios;

template <typename T,int N,class E>
ostream &operator<<(ostream &o,const Vector<T,N,E> &v) {
    o << "[" << setw(10) << setprecision(8) << setiosflags(ios::fixed) << v(0);
    for(int n=1;n<v.size();n++)
        o << "," << setw(10) << setprecision(8) << setiosflags(ios::fixed) << v(n);
    o << "]";
    return o;
};

template <typename T,int R,int C,class E>
ostream &operator<<(ostream &o,const Matrix<T,R,C,E> &m) {
    for(int r=0;r<m.rows();r++) {
        o << ((r==0)?"[":"],\n ");
        for(int c=0;c<m.cols();c++)
            o << ((c==0)?"[":",") << setw(10) << setprecision(8) << setiosflags(ios::fixed) << m(r,c);
    };
    o << "]]\n";
    return o;
};

template <typename T,int N>
ostream &operator<<(ostream &o,const Tuple<T,N> &t) {
    o << "(" << setw(10) << setprecision(8) << setiosflags(ios::fixed) << t[0];
    for(int n=1;n<N;n++)
        o << "," << setw(10) << setprecision(8) << setiosflags(ios::fixed) << t[n];
    o << ")";
    return o;
};

}; // namespace etLAP

#endif // _ETLAP_OUTPUT_H_
