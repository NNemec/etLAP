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

#ifndef _ETLAP_COSTS_H_
#define _ETLAP_COSTS_H_

namespace etLAP {

template <class C> struct Costs;// { enum { c = 0 }; };

}; // namespace etLAP

#include "Tags.h"

namespace etLAP {

template <> struct Costs<int> { enum { c = 0 }; };

template <> struct Costs<Smart> { enum { c = 0 }; };
template <> struct Costs<Packed> { enum { c = 0 }; };
template <> struct Costs<ExtPtr> { enum { c = 0 }; };
template <> struct Costs<Diagonal> { enum { c = 0 }; };
template <> struct Costs<Scalar> { enum { c = 0 }; };
template <> struct Costs<Zero> { enum { c = 0 }; };
template <> struct Costs<One> { enum { c = 0 }; };

template <class A,class Op>
struct Costs<ElemUnOp<A,Op> > { enum { c = Costs<A>::c + Costs<Op>::c }; };

template <class A,class B,class Op>
struct Costs<ElemBinOp<A,B,Op> > { enum { c = Costs<A>::c + Costs<B>::c + Costs<Op>::c }; };

template <class A,typename T,class Op>
struct Costs<ScalarOp<A,T,Op> > { enum { c = Costs<A>::c + Costs<Op>::c }; };

template <class A>
struct Costs<Transposed<A> > { enum { c = Costs<A>::c }; };

template <class T>
struct Costs<Row<T> > { enum { c = Costs<T>::c }; };
template <class T>
struct Costs<Column<T> > { enum { c = Costs<T>::c }; };
template <> struct Costs<Buffer> { enum { c = 0 }; };
template <class A>
struct Costs<NoBuffer<A> > { enum { c = 0 }; };

template <> struct Costs<OpIdent> { enum { c = 0 }; };
template <> struct Costs<OpNeg> { enum { c = 1 }; };
template <> struct Costs<OpConj> { enum { c = 1 }; };
template <> struct Costs<OpAbs> { enum { c = 1 }; };
template <> struct Costs<OpReal> { enum { c = 0 }; };
template <> struct Costs<OpImag> { enum { c = 0 }; };

template <> struct Costs<OpAdd> { enum { c = 1 }; };
template <> struct Costs<OpSub> { enum { c = 1 }; };
template <> struct Costs<OpMul> { enum { c = 1 }; };
template <> struct Costs<OpDiv> { enum { c = 1 }; };
template <> struct Costs<OpRevMul> { enum { c = 1 }; };

}; // namespace etLAP

//#include "Vector.h"
//#include "Matrix.h"

namespace etLAP {

template <typename T,int N,class E>
struct Costs<Vector<T,N,E> > { enum { c = Costs<E>::c }; };
template <int R,int C,typename T,class E>
struct Costs<Matrix<T,R,C,E> > { enum { c = Costs<E>::c }; };

template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct Costs<Multiplied<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > > {
    enum { M = SizeCombine<C1,R2>::n };
    enum { factor = M?M:5 };
    enum { c = factor * (Costs<E1>::c + Costs<E2>::c + Costs<OpMul>::c) };
};

template <typename T1,int R1,int C1,class E1,typename T2,int N2,class E2>
struct Costs<Multiplied<Matrix<T1,R1,C1,E1>,Vector<T2,N2,E2> > > {
    enum { M = SizeCombine<C1,N2>::n };
    enum { factor = M?M:5 };
    enum { c = factor * (Costs<E1>::c + Costs<E2>::c + Costs<OpMul>::c) };
};

template <typename T1,int N1,class E1,typename T2,int R2,int C2,class E2>
struct Costs<Multiplied<Vector<T1,N1,E1>,Matrix<T2,R2,C2,E2> > > {
    enum { M = SizeCombine<N1,R2>::n };
    enum { factor = M?M:5 };
    enum { c = factor * (Costs<E1>::c + Costs<E2>::c + Costs<OpMul>::c) };
};

}; // namespace etLAP

#endif // _ETLAP_COSTS_H_
