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

#ifndef _ETLAP_OPERATOR_H_
#define _ETLAP_OPERATOR_H_

#include "Tags.h"
#include "Types.h"
#include "Vector.h"
#include "Matrix.h"

namespace etLAP {

/******************************************************************************
 * Unary operators
 */

#define DefineVectorElemUnOp(OPTAG)                                                                 \
template <typename T,int N,class E>                                                                 \
struct UnaryOp<OPTAG,Vector<T,N,E> > {                                                              \
    typedef Vector<typename UnaryOp<OPTAG,T>::Result_t,N,ElemUnOp<Vector<T,N,E>,OPTAG> > Result_t;  \
    static Result_t apply(const Vector<T,N,E> &v) {                                                 \
        return Result_t(v);                                                                         \
    };                                                                                              \
};

DefineVectorElemUnOp(OpIdent)
DefineVectorElemUnOp(OpNeg)
DefineVectorElemUnOp(OpConj)
DefineVectorElemUnOp(OpAbs)
DefineVectorElemUnOp(OpReal)
DefineVectorElemUnOp(OpImag)

#define DefineMatrixElemUnOp(OPTAG)                                                                     \
template <typename T,int R,int C,class E>                                                               \
struct UnaryOp<OPTAG,Matrix<T,R,C,E> > {                                                                \
    typedef Matrix<typename UnaryOp<OPTAG,T>::Result_t,R,C,ElemUnOp<Matrix<T,R,C,E>,OPTAG> > Result_t;  \
    static Result_t apply(const Matrix<T,R,C,E> &m) {                                                   \
        return Result_t(m);                                                                             \
    };                                                                                                  \
};

DefineMatrixElemUnOp(OpIdent)
DefineMatrixElemUnOp(OpNeg)
DefineMatrixElemUnOp(OpConj)
DefineMatrixElemUnOp(OpAbs)
DefineMatrixElemUnOp(OpReal)
DefineMatrixElemUnOp(OpImag)

// sqr(Vector)
template <typename T,int N,class E>
struct UnaryOp<OpSqr,Vector<T,N,E> > {
    typedef typename UnaryOp<OpSqr,T>::Result_t Result_t;
    static Result_t apply(const Vector<T,N,E> &a) {
        Result_t res = 0;
        for (int n = a.size(); n-- > 0;) {
            res += sqr(a(n));
        };
        return res;
    };
};


// transpose(Matrix)
template <int R,int C,typename T,class E>
inline const Matrix<T,R,C,Transposed<Matrix<T,C,R,E> > >
transpose(const Matrix<T,C,R,E> &m) {
    return Matrix<T,R,C,Transposed<Matrix<T,C,R,E> > >(m);
};

// adj(Matrix)
template <typename T,int R,int C,class E>
inline const Matrix<T,C,R,Transposed<Matrix<T,R,C,ElemUnOp<Matrix<T,R,C,E>,OpConj> > > >
//inline typeof(transpose(conj(Matrix<T,R,C,E>)))
adj(const Matrix<T,R,C,E> &m) {
    return transpose(conj(m));
};

// trace(Matrix)
template <int N,typename T,class E>
inline T trace(const Matrix<T,N,N,E> &m) {
    assert(m.rows() == m.cols());
    T res = (T)0;
    for(int i=m.rows();i-->0;) res += m(i,i);
    return res;
};

// det(Matrix)
template <typename T,class E>
inline T det(const Matrix<T,2,2,E> &m) {
    assert(m.rows() == m.cols());
    return m(0,0)*m(1,1)-m(0,1)*m(1,0);
};

// inv(Matrix)
template <int N,typename T,class E>
inline const Matrix<T,N,N> inv(Matrix<T,N,N,E> a) {
    assert(a.cols() == a.rows());
    int size = a.cols();
    Matrix<T,N,N> res(size,size);
    res = 1;

    for (int c = size; c-- > 0;) {
        int rmax = c;
        T pivot = a(c, c);
        for (int r = c; r-- > 0;)
            if (abs(a(r, c)) > abs(pivot)) {
                rmax = r;
                pivot = a(r, c);
            };
        for (int i = size; i-- > 0;) {
            T x = a(rmax, i);
            a(rmax, i) = a(c, i);
            a(c, i) = x / pivot;

            x = res(rmax, i);
            res(rmax, i) = res(c, i);
            res(c, i) = x / pivot;
        };
        for (int r = size; r-- > 0;)
            if (r != c) {
                pivot = a(r, c);
                for (int i = size; i-- > 0;) {
                    a(r, i) -= pivot * a(c, i);
                    res(r, i) -= pivot * res(c, i);
                };
            };
    };
    return res;
};


// max(Matrix)
template <typename T,int R,int C,class E>
inline const typename UnaryOp<OpAbs,T>::Result_t max(Matrix<T,R,C,E> a) {
    typename UnaryOp<OpAbs,T>::Result_t res = 0;
//inline const typeof(abs(T(0)) max(Matrix<T,R,C,E> a) {
//    typeof(abs(T(0)) res = 0;
    for (int r = a.rows(); r-- > 0;)
    for (int c = a.cols(); c-- > 0;)
        res = res >? abs(a(r,c));
    return res;
};

// simplemax(Matrix)
template <typename T,int R,int C,class E>
inline const typename UnaryOp<OpAbs,T>::Result_t simplemax(Matrix<T,R,C,E> a) {
    typename UnaryOp<OpAbs,T>::Result_t res = 0;
//inline const typeof(abs(T(0)) max(Matrix<T,R,C,E> a) {
//    typeof(abs(T(0)) res = 0;
    for (int r = a.rows(); r-- > 0;)
    for (int c = a.cols(); c-- > 0;)
        res = res >? fabs(a(r,c));
    return res;
};

// simplemax(Matrix)
template <typename T,int R,int C,class E>
inline const typename UnaryOp<OpAbs,T>::Result_t simplemax(Matrix<std__complex<T>,R,C,E> a) {
    typename UnaryOp<OpAbs,T>::Result_t res = 0;
//inline const typeof(abs(T(0)) max(Matrix<T,R,C,E> a) {
//    typeof(abs(T(0)) res = 0;
    for (int r = a.rows(); r-- > 0;)
    for (int c = a.cols(); c-- > 0;)
        res = res >? fabs(real(a(r,c))) >? fabs(imag(a(r,c)));
    return res;
};


// exp(Matrix)
template <int N,typename T,class E>
inline const Matrix<T,N,N> exp(Matrix<T,N,N,E> a) {
    assert(a.cols() == a.rows());

    typedef typename UnaryOp<OpAbs,T>::Result_t REAL;
    REAL max = simplemax(a);
    if(max > 10.0) {
        // use product expansion exp(a) = lim_[n->inf] (1+a/n)^n
        int lb_n = ilogb(max)+((TypeEqual<REAL,double>::res)?26:13);
        Matrix<T,N,N> res = 1 + a/pow(2,lb_n);
        while(lb_n-- > 0)
            res = buf(res*res);
        return res;
    } else {
        // use product expansion exp(a) = sum[n] a^n/n!
        int size = a.cols();
        Matrix<T,N,N> res1(size,size);
        Matrix<T,N,N> res2(size,size);
        Matrix<T,N,N> tmp1(size,size);
        Matrix<T,N,N> tmp2(size,size);
        tmp1 = a*a/2;
        res1 = 1 + a + tmp1;
        for(int n=3;n<1000;n++) {
            tmp2 = tmp1 * a / n;
            res2 = res1+tmp2;
            if(res1 == res2)
                return res1;
            n++;
            tmp1 = tmp2 * a / n;
            res1 = res2+tmp1;
            if(res1 == res2)
                return res1;
        };
        std::cerr << "exp(Matrix) failed:\n"
                     "a=" << a << "\n"
                     "res1=" << res1 << "\n"
                     "res2=" << res2 << "\n";
/*
        std::cerr << a;
        std::cerr << res1;
        std::cerr << res2;
        std::cerr << tmp1;
        std::cerr << tmp2;
*/
        abort();
        throw "exp(Matrix) failed!!";
    };
};

// buf(Vector)
template <typename T,int N,class E>
inline const Vector<T,N,Buffer>
buf(const Vector<T,N,E> &v) {
    return Vector<T,N,Buffer>(v);
};

// buf(Vector<Buffer>)
template <typename T,int N>
inline const Vector<T,N,Buffer>
buf(const Vector<T,N,Buffer> &v) {
    return v;
};

// buf(Matrix)
template <typename T,int R,int C,class E>
inline const Matrix<T,R,C,Buffer>
buf(const Matrix<T,R,C,E> &m) {
    return Matrix<T,R,C,Buffer>(m);
};

// buf(Matrix<Buffer>)
template <typename T,int R,int C,class E>
inline const Matrix<T,R,C,Buffer>
buf(const Matrix<T,R,C,Buffer> &m) {
    return m;
};

// nobuf(Vector)
template <typename T,int N,class E>
inline const Vector<T,N,NoBuffer<Vector<T,N,E> > >
nobuf(const Vector<T,N,E> &m) {
    return Vector<T,N,NoBuffer<Vector<T,N,E> > >(m);
};

// nobuf(Matrix)
template <typename T,int R,int C,class E>
inline const Matrix<T,R,C,NoBuffer<Matrix<T,R,C,E> > >
nobuf(const Matrix<T,R,C,E> &m) {
    return Matrix<T,R,C,NoBuffer<Matrix<T,R,C,E> > >(m);
};

/******************************************************************************
 * Binary operator structures
 */

// Vector + Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
struct BinaryOp<OpAdd,Vector<T1,N1,E1>,Vector<T2,N2,E2> > {
    typedef Vector<T1,N1,E1> X1;
    typedef Vector<T2,N2,E2> X2;
    typedef typename BinaryOp<OpAdd,T1,T2>::Result_t T;
    static const int N = SizeCombine<N1,N2>::n;
    typedef ElemBinOp<X1,X2,OpAdd> E;
    typedef Vector<T,N,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) {
        assert(x1.size() == x2.size());
        return Result_t(x1,x2);
    };
};

// Vector - Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
struct BinaryOp<OpSub,Vector<T1,N1,E1>,Vector<T2,N2,E2> > {
    typedef Vector<T1,N1,E1> X1;
    typedef Vector<T2,N2,E2> X2;
    typedef typename BinaryOp<OpSub,T1,T2>::Result_t T;
    static const int N = SizeCombine<N1,N2>::n;
    typedef ElemBinOp<X1,X2,OpSub> E;
    typedef Vector<T,N,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) {
        assert(x1.size() == x2.size());
        return Result_t(x1,x2);
    };
};

// Vector * Scalar
template <int N,typename T1,class E1,typename T2>
struct BinaryOp<OpMul,Vector<T1,N,E1>,T2> {
    typedef Vector<T1,N,E1> V;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef ScalarOp<V,T2,OpMul> E;
    typedef Vector<T,N,E> Result_t;
    static Result_t apply(const V &v,const T2 &s) { return Result_t(v,s); };
};

// Vector / Scalar
template <int N,typename T1,class E1,typename T2>
struct BinaryOp<OpDiv,Vector<T1,N,E1>,T2> {
    typedef Vector<T1,N,E1> V;
    typedef typename BinaryOp<OpDiv,T1,T2>::Result_t T;
    typedef ScalarOp<V,T2,OpDiv> E;
    typedef Vector<T,N,E> Result_t;
    static Result_t apply(const V &v,const T2 &s) { return Result_t(v,s); };
};

// Scalar * Vector
template <typename T1,int N,typename T2,class E2>
struct BinaryOp<OpMul,T1,Vector<T2,N,E2> > {
    typedef Vector<T2,N,E2> V;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef ScalarOp<V,T1,OpRevMul> E;
    typedef Vector<T,N,E> Result_t;
    static Result_t apply(const T1 &s,const V &v) { return Result_t(v,s); };
};

// Vector * Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
struct BinaryOp<OpMul,Vector<T1,N1,E1>,Vector<T2,N2,E2> > {
    typedef Vector<T1,N1,E1> X1;
    typedef Vector<T2,N2,E2> X2;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t Result_t;
    static Result_t apply(const X1 &v1,const X2 &v2) {
        CTAssert(N1 == 0 || N2 == 0 || N1 == N2);
        assert(v1.size() == v2.size());
        Result_t res = (Result_t)0;
        for(int i=v1.size();i-->0;) res += BinaryOp<OpMul,T1,T2>::apply(v1(i),v2(i));
        return res;
    };
};

// Vector == Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
struct BinaryOp<OpIsEq,Vector<T1,N1,E1>,Vector<T2,N2,E2> > {
    typedef Vector<T1,N1,E1> X1;
    typedef Vector<T2,N2,E2> X2;
    typedef bool Result_t;
    static Result_t apply(const X1 &v1,const X2 &v2) {
        CTAssert(N1 == 0 || N2 == 0 || N1 == N2);
        assert(v1.size() == v2.size());
        for(int i=v1.size();i-->0;)
            if(! v1(i) == v2(i))
                return false;
        return true;
    };
};

// Matrix + Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct BinaryOp<OpAdd,Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > {
    typedef Matrix<T1,R1,C1,E1> X1;
    typedef Matrix<T2,R2,C2,E2> X2;
    typedef typename BinaryOp<OpAdd,T1,T2>::Result_t T;
    static const int R = SizeCombine<R1,R2>::n;
    static const int C = SizeCombine<C1,C2>::n;
    typedef ElemBinOp<X1,X2,OpAdd> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) {
        assert(x1.rows() == x2.rows() && x1.cols() == x2.cols());
        return Result_t(x1,x2);
    };
};

// Matrix - Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct BinaryOp<OpSub,Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > {
    typedef Matrix<T1,R1,C1,E1> X1;
    typedef Matrix<T2,R2,C2,E2> X2;
    typedef typename BinaryOp<OpSub,T1,T2>::Result_t T;
    static const int R = SizeCombine<R1,R2>::n;
    static const int C = SizeCombine<C1,C2>::n;
    typedef ElemBinOp<X1,X2,OpSub> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) {
        assert(x1.rows() == x2.rows() && x1.cols() == x2.cols());
        return Result_t(x1,x2);
    };
};

// Matrix * Scalar
template <typename T1,int R,int C,class E1,typename T2>
struct BinaryOp<OpMul,Matrix<T1,R,C,E1>,T2> {
    typedef Matrix<T1,R,C,E1> M;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef ScalarOp<M,T2,OpMul> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const M &m,const T2 &s) { return Result_t(m,s); };
};

// Matrix / Scalar
template <typename T1,int R,int C,class E1,typename T2>
struct BinaryOp<OpDiv,Matrix<T1,R,C,E1>,T2> {
    typedef Matrix<T1,R,C,E1> M;
    typedef typename BinaryOp<OpDiv,T1,T2>::Result_t T;
    typedef ScalarOp<M,T2,OpDiv> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const M &m,const T2 &s) { return Result_t(m,s); };
};

// Scalar * Matrix
template <typename T1,typename T2,int R,int C,class E2>
struct BinaryOp<OpMul,T1,Matrix<T2,R,C,E2> > {
    typedef Matrix<T2,R,C,E2> M;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef ScalarOp<M,T1,OpRevMul> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const T1 &s,const M &m) { return Result_t(m,s); };
};

// Matrix + Scalar
template <typename T1,int R,int C,class E1,typename T2>
struct BinaryOp<OpAdd,Matrix<T1,R,C,E1>,T2> {
    typedef Matrix<T1,R,C,E1> X1;
    typedef Matrix<T2,R,C,Scalar> X2;
    typedef typename BinaryOp<OpAdd,T1,T2>::Result_t T;
    typedef ElemBinOp<X1,X2,OpAdd> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const X1 &m,const T2 &s) {
        CTAssert(R == C);
        assert(m.rows() == m.cols());
        return Result_t(m,X2(m.rows(),s));
    };
};

// Matrix - Scalar
template <typename T1,int R,int C,class E1,typename T2>
struct BinaryOp<OpSub,Matrix<T1,R,C,E1>,T2> {
    typedef Matrix<T1,R,C,E1> X1;
    typedef Matrix<T2,R,C,Scalar> X2;
    typedef typename BinaryOp<OpSub,T1,T2>::Result_t T;
    typedef ElemBinOp<X1,X2,OpSub> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const X1 &m,const T2 &s) {
        CTAssert(R == C);
        assert(m.rows() == m.cols());
        return Result_t(m,X2(m.rows(),s));
    };
};

// Scalar + Matrix
template <typename T1,typename T2,int R,int C,class E2>
struct BinaryOp<OpAdd,T1,Matrix<T2,R,C,E2> > {
    typedef Matrix<T1,R,C,Scalar> X1;
    typedef Matrix<T2,R,C,E2> X2;
    typedef typename BinaryOp<OpAdd,T1,T2>::Result_t T;
    typedef ElemBinOp<X1,X2,OpAdd> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const T1 &s,const X2 &m) {
        CTAssert(R == C);
        assert(m.rows() == m.cols());
        return Result_t(X1(m.rows(),s),m);
    };
};

// Scalar - Matrix
template <typename T1,typename T2,int R,int C,class E2>
struct BinaryOp<OpSub,T1,Matrix<T2,R,C,E2> > {
    typedef Matrix<T1,R,C,Scalar> X1;
    typedef Matrix<T2,R,C,E2> X2;
    typedef typename BinaryOp<OpSub,T1,T2>::Result_t T;
    typedef ElemBinOp<X1,X2,OpSub> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const T1 &s,const X1 &m) {
        CTAssert(R == C);
        assert(m.rows() == m.cols());
        return Result_t(X1(m.rows(),s),m);
    };
};

// Matrix * Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct BinaryOp<OpMul,Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > {
    static const int M=SizeCombine<C1,R2>::n;
    static const int factor = (M==0?1:M-1);
    typedef Matrix<T1,R1,C1,E1> X1;
    typedef typename IF<(factor*Traits<E1>::costs > 2),Matrix<T1,R1,C1,Buffer>,X1>::T X1_buf;
    typedef Matrix<T2,R2,C2,E2> X2;
    typedef typename IF<(factor*Traits<E2>::costs > 2),Matrix<T2,R2,C2,Buffer>,X2>::T X2_buf;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    static const int R = C2?R1:0;
    static const int C = R1?C2:0;
    typedef Multiplied<X1_buf,X2_buf> E;
    typedef Matrix<T,R,C,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) { assert(x1.cols() == x2.rows()); return Result_t(x1,x2); };
};

// Vector * Matrix
template <typename T1,int N1,class E1,typename T2,int R2,int C2,class E2>
struct BinaryOp<OpMul,Vector<T1,N1,E1>,Matrix<T2,R2,C2,E2> > {
    static const int factor = (N1==0?1:N1-1);
    typedef Vector<T1,N1,E1> X1;
    typedef typename IF<(factor*Traits<E1>::costs > 2),Vector<T1,N1,Buffer>,X1>::T X1_buf;
    typedef Matrix<T2,R2,C2,E2> X2;
    typedef typename IF<(factor*Traits<E2>::costs > 2),Matrix<T2,R2,C2,Buffer>,X2>::T X2_buf;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef Multiplied<X1_buf,X2_buf> E;
    typedef Vector<T,C2,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) { assert(x1.size() == x2.rows()); return Result_t(x1,x2); };
};

// Matrix * Vector
template <typename T1,int R1,int C1,class E1,typename T2,int N2,class E2>
struct BinaryOp<OpMul,Matrix<T1,R1,C1,E1>,Vector<T2,N2,E2> > {
    static const int factor = (N2==0?1:N2-1);
    typedef Matrix<T1,R1,C1,E1> X1;
    typedef typename IF<(factor*Traits<E1>::costs > 2),Matrix<T1,R1,C1,Buffer>,X1>::T X1_buf;
    typedef Vector<T2,N2,E2> X2;
    typedef typename IF<(factor*Traits<E2>::costs > 2),Vector<T2,N2,Buffer>,X2>::T X2_buf;
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    typedef Multiplied<X1_buf,X2_buf> E;
    typedef Vector<T,R1,E> Result_t;
    static Result_t apply(const X1 &x1,const X2 &x2) { assert(x1.cols() == x2.size()); return Result_t(x1,x2); };
};

// commute(Matrix,Matrix)
template <int N,typename T1,class E1,typename T2,class E2>
inline const Matrix<typename BinaryOp<OpMul,T1,T2>::Result_t,N,N> commute(const Matrix<T1,N,N,E1> &m1,const Matrix<T2,N,N,E2> &m2) {
    typedef typename BinaryOp<OpMul,T1,T2>::Result_t T;
    if(Traits<E1>::costs > 0)
        if(Traits<E2>::costs > 0) {
            Matrix<T,N,N> tmp1 = m1;
            Matrix<T,N,N> tmp2 = m2;
            return tmp1*tmp2-tmp2*tmp1;
        } else {
            Matrix<T,N,N> tmp1 = m1;
            return tmp1*m2-m2*tmp1;
        }
    else
        if(Traits<E2>::costs > 0) {
            Matrix<T,N,N> tmp1 = m1;
            return tmp1*m2-m2*tmp1;
        } else {
            return m1*m2-m2*m1;
        }
};

// Matrix == Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
struct BinaryOp<OpIsEq,Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > {
    typedef Matrix<T1,R1,C1,E1> X1;
    typedef Matrix<T2,R2,C2,E2> X2;
    typedef bool Result_t;
    static Result_t apply(const X1 &m1,const X2 &m2) {
        CTAssert((R1 == 0 && C1 == 0) || (R2 == 0 && C2 == 0) || (R1 == R2 && C1 == C2));
        assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
        for(int r=m1.rows();r-->0;)
        for(int c=m1.cols();c-->0;)
            if(! (m1(r,c) == m2(r,c)))
                return false;
        return true;
    };
};



/*******************************************************************************
 * Operator routines
 */

#define Xheader1 typename X1
#define Xheader2 typename X2
#define Xtype1 X1
#define Xtype2 X2
#define Vheader1 typename T1,int N1,class E1
#define Vheader2 typename T2,int N2,class E2
#define Vtype1 Vector<T1,N1,E1>
#define Vtype2 Vector<T2,N2,E2>
#define Mheader1 typename T1,int R1,int C1,class E1
#define Mheader2 typename T2,int R2,int C2,class E2
#define Mtype1 Matrix<T1,R1,C1,E1>
#define Mtype2 Matrix<T2,R2,C2,E2>

#define DefineUnOp(NAME,FNAME,ARG)                         \
template <ARG##header1>                                     \
inline const typename UnaryOp<NAME,ARG##type1>::Result_t    \
FNAME(const ARG##type1 &x) {                                \
    return UnaryOp<NAME,ARG##type1>::apply(x);             \
}

#define DefineUnOps(NAME,FNAME) \
DefineUnOp(NAME,FNAME,V)      \
DefineUnOp(NAME,FNAME,M)

DefineUnOps(OpIdent,ident)
DefineUnOps(OpNeg,operator-)
DefineUnOps(OpConj,conj)
DefineUnOps(OpAbs,abs)
DefineUnOps(OpReal,real)
DefineUnOps(OpImag,imag)
DefineUnOp(OpSqr,sqr,V)

#undef DefineUnOps
#undef DefineUnOp


#define DefineBinOp(NAME,FNAME,ARG1,ARG2)                              \
template <ARG1##header1,ARG2##header2>                                  \
inline const typename BinaryOp<NAME,ARG1##type1,ARG2##type2>::Result_t  \
FNAME(const ARG1##type1 &x1,const ARG2##type2 &x2) {                    \
    return BinaryOp<NAME,ARG1##type1,ARG2##type2>::apply(x1,x2);        \
}

#define DefineBinOps(NAME,FNAME)     \
DefineBinOp(NAME,FNAME,V,X)        \
DefineBinOp(NAME,FNAME,X,V)        \
DefineBinOp(NAME,FNAME,V,V)        \
DefineBinOp(NAME,FNAME,M,X)        \
DefineBinOp(NAME,FNAME,X,M)        \
DefineBinOp(NAME,FNAME,M,M)        \
DefineBinOp(NAME,FNAME,V,M)        \
DefineBinOp(NAME,FNAME,M,V)

DefineBinOps(OpAdd,operator+)
DefineBinOps(OpSub,operator-)
DefineBinOps(OpMul,operator*)
DefineBinOps(OpDiv,operator/)
DefineBinOps(OpIsEq,operator==)
DefineBinOps(OpIsNeq,operator!=)

#undef DefineBinOps
#undef DefineBinOp

/******************************************************************************
 * Assignment operators
 */

// Vector += Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
inline const Vector<T1,N1,E1> &
operator+=(Vector<T1,N1,E1> &v1,const Vector<T2,N2,E2> &v2) {
    assign_add(v1,v2);
    return v1;
};

// Vector -= Vector
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
inline const Vector<T1,N1,E1> &
operator-=(Vector<T1,N1,E1> &v1,const Vector<T2,N2,E2> &v2) {
    assign_sub(v1,v2);
    return v1;
};

// Vector *= Scalar
template <typename T1,int N1,class E1,typename T2>
inline const Vector<T1,N1,E1> &
operator*=(Vector<T1,N1,E1> &v,const T2 &s) {
    assign_mul(v,s);
    return v;
};

// Vector /= Scalar
template <typename T1,int N1,class E1,typename T2>
inline const Vector<T1,N1,E1> &
operator/=(Vector<T1,N1,E1> &v,const T2 &s) {
    assign_div(v,s);
    return v;
};

// Vector *= Matrix
template <typename T1,int N1,class E1,typename T2,int N2,class E2>
inline const Vector<T1,N1,E1> &
operator*=(Vector<T1,N1,E1> &v,const Matrix<T2,N2,N2,E2> &m) {
    CTAssert(N1==0 || N2==0 || N1 == N2);
    assert(m.cols() == m.rows() && v.size() == m.rows());
    return v = buf(v * m);
};

// Matrix += Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
inline const Matrix<T1,R1,C1,E1> &
operator+=(Matrix<T1,R1,C1,E1> &m1,const Matrix<T2,R2,C2,E2> &m2) {
    assign_add(m1,m2);
    return m1;
};

// Matrix -= Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
inline const Matrix<T1,R1,C1,E1> &
operator-=(Matrix<T1,R1,C1,E1> &m1,const Matrix<T2,R2,C2,E2> &m2) {
    assign_sub(m1,m2);
    return m1;
};

// Matrix *= Scalar
template <typename T1,int R1,int C1,class E1,typename T2>
inline const Matrix<T1,R1,C1,E1> &
operator*=(Matrix<T1,R1,C1,E1> &m,const T2 &s) {
    assign_mul(m,s);
    return m;
};

// Matrix /= Scalar
template <typename T1,int R1,int C1,class E1,typename T2>
inline const Matrix<T1,R1,C1,E1> &
operator/=(Matrix<T1,R1,C1,E1> &m,const T2 &s) {
    assign_div(m,s);
    return m;
};

// Matrix *= Matrix
template <typename T1,int R1,int C1,class E1,typename T2,int N2,class E2>
inline const Matrix<T1,R1,C1,E1> &
operator*=(Matrix<T1,R1,C1,E1> &m1,const Matrix<T2,N2,N2,E2> &m2) {
    CTAssert(C1==0 || N2==0 || C1 == N2);
    assert(m2.cols() == m2.rows() && m1.cols() == m2.rows());
    return m1 = buf(m1 * m2);
};

/******************************************************************************
 * Special "clear" operator
 */

template <typename T>
void clear(T &t) {
    t = (T)0;
};

template <typename T,int N,class E>
void clear(Vector<T,N,E> &v) {
    v.clear();
};

template <typename T,int R,int C,class E>
void clear(Matrix <T,R,C,E> &m) {
    m.clear();
};

}; // namespace etLAP

#endif // _ETLAP_OPERATOR_H_
