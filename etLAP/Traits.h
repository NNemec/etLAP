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

#pragma once

namespace etLAP {

template <class C> struct Traits;

} // namespace etLAP

#include "Tags.h"

namespace etLAP {

template <> struct Traits<int> { static const int costs = 0; static const bool flataccess = true; };

template <> struct Traits<Smart> { static const int costs = 0; static const bool flataccess = true; };
template <> struct Traits<Packed> { static const int costs = 0; static const bool flataccess = true; };
template <> struct Traits<ExtPtr> { static const int costs = 0; static const bool flataccess = true; };
template <> struct Traits<Diagonal> { static const int costs = 0; static const bool flataccess = false; };
template <> struct Traits<Scalar> { static const int costs = 0; static const bool flataccess = false; };
template <> struct Traits<Zero> { static const int costs = 0; static const bool flataccess = true; };
template <> struct Traits<One> { static const int costs = 0; static const bool flataccess = false; };

template <class A,class Op>
struct Traits<ElemUnOp<A,Op> > {
    static const int costs = Traits<A>::costs + Traits<Op>::costs;
    static const bool flataccess = Traits<A>::flataccess;
};

template <class A,class B,class Op>
struct Traits<ElemBinOp<A,B,Op> > {
    static const int costs = Traits<A>::costs + Traits<B>::costs + Traits<Op>::costs;
    static const bool flataccess = Traits<A>::flataccess && Traits<B>::flataccess;
};

template <class A,typename T,class Op>
struct Traits<ScalarOp<A,T,Op> > {
    static const int costs = Traits<A>::costs + Traits<Op>::costs;
    static const bool flataccess = Traits<A>::flataccess;
};

template <class A>
struct Traits<Transposed<A> > { static const int costs = Traits<A>::costs; static const bool flataccess = false; };

template <class T>
struct Traits<Row<T> > { static const int costs = Traits<T>::costs; static const bool flataccess = true; };
template <class T>
struct Traits<Column<T> > { static const int costs = Traits<T>::costs; static const bool flataccess = true; };
template <> struct Traits<Buffer> { static const int costs = 0; static const bool flataccess = true; };
template <class A>
struct Traits<NoBuffer<A> > { static const int costs = 0; static const bool flataccess = Traits<A>::flataccess; };

template <> struct Traits<OpIdent> { static const int costs = 0; };
template <> struct Traits<OpNeg> { static const int costs = 1; };
template <> struct Traits<OpConj> { static const int costs = 1; };
template <> struct Traits<OpAbs> { static const int costs = 1; };
template <> struct Traits<OpReal> { static const int costs = 0; };
template <> struct Traits<OpImag> { static const int costs = 0; };

template <> struct Traits<OpAdd> { static const int costs = 1; };
template <> struct Traits<OpSub> { static const int costs = 1; };
template <> struct Traits<OpMul> { static const int costs = 1; };
template <> struct Traits<OpDiv> { static const int costs = 1; };
template <> struct Traits<OpRevMul> { static const int costs = 1; };

} // namespace etLAP

//#include "Vector.h"
//#include "Matrix.h"

namespace etLAP {

template <typename T,int N,class E>
struct Traits<Vector<T,N,E> > { static const int costs = Traits<E>::costs; static const bool flataccess = false; };
template <int R,int C,typename T,class E>
struct Traits<Matrix<T,R,C,E> > { static const int costs = Traits<E>::costs; static const bool flataccess = Traits<E>::flataccess; };

template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct Traits<Multiplied<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > > {
    static const int M = SizeCombine<C1,R2>::n;
    static const int factor = M?M:5;
    static const int costs = factor * (Traits<E1>::costs + Traits<E2>::costs + Traits<OpMul>::costs);
    static const bool flataccess = false;
};

template <typename T1,int R1,int C1,class E1,typename T2,int N2,class E2>
struct Traits<Multiplied<Matrix<T1,R1,C1,E1>,Vector<T2,N2,E2> > > {
    static const int M = SizeCombine<C1,N2>::n;
    static const int factor = M?M:5;
    static const int costs = factor * (Traits<E1>::costs + Traits<E2>::costs + Traits<OpMul>::costs);
};

template <typename T1,int N1,class E1,typename T2,int R2,int C2,class E2>
struct Traits<Multiplied<Vector<T1,N1,E1>,Matrix<T2,R2,C2,E2> > > {
    static const int M = SizeCombine<N1,R2>::n;
    static const int factor = M?M:5;
    static const int costs = factor * (Traits<E1>::costs + Traits<E2>::costs + Traits<OpMul>::costs);
};

} // namespace etLAP
