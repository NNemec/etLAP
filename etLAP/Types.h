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

#ifndef _ETLAP_TYPES_H_
#define _ETLAP_TYPES_H_

#if __GNUC__==2
template <typename T>
class complex;
#define std__complex complex
#else
namespace std {
template <typename T>
class complex;
}
#define std__complex std::complex
#endif

#include <complex>
#include "Util.h"
#include "Tags.h"

namespace etLAP {


struct AnyType;


/*******************************************************************************
 * TypeCast
 */
class NoCast {};
typedef class NoClass {} *NoType;
template <typename Dest,typename Src> struct TypeCast {
    static const bool valid = false;
    typedef NoCast Cast_t;
    typedef NoType Dest_t;
};

template <typename T> struct TypeCast<T,T> {
    static const bool valid = true;
    typedef TypeCast<T,T> Cast_t;
    typedef T      Dest_t;
    static T cast(T x) { return x; }
};

#define DEFINE_NUMERIC_CAST(TO,FROM)     \
template <> struct TypeCast<TO,FROM> {   \
    static const bool valid = true;      \
    typedef TypeCast<TO,FROM> Cast_t;    \
    typedef TO Dest_t;                   \
    static TO cast(FROM x) { return x; } \
};

DEFINE_NUMERIC_CAST(float,int)
DEFINE_NUMERIC_CAST(double,int)
DEFINE_NUMERIC_CAST(double,float)
DEFINE_NUMERIC_CAST(long double,int)
DEFINE_NUMERIC_CAST(long double,float)
DEFINE_NUMERIC_CAST(long double,double)

#undef DEFINE_NUMERIC_CAST

template <typename Dest,typename Src>
struct TypeCast<std__complex<Dest>,std__complex<Src> > {
    static const bool valid = TypeCast<Dest,Src>::valid;
    typedef typename IF<valid,TypeCast<std__complex<Dest>,std__complex<Src> >,NoCast>::T Cast_t;
    typedef typename IF<valid,std__complex<typename TypeCast<Dest,Src>::Dest_t>,NoType>::T Dest_t;
    static std__complex<Dest> cast(std__complex<Src> x) { return std__complex<Dest>(TypeCast<Dest,Src>::cast(real(x)),TypeCast<Dest,Src>::cast(imag(x))); }
};

template <typename Dest,typename Src>
struct TypeCast<std__complex<Dest>,Src > {
    static const bool valid = TypeCast<Dest,Src>::valid;
    typedef typename IF<valid,TypeCast<std__complex<Dest>,Src >,NoCast>::T Cast_t;
    typedef typename IF<valid,std__complex<typename TypeCast<Dest,Src>::Dest_t>,NoType>::T Dest_t;
    static std__complex<Dest> cast(Src x) { return std__complex<Dest>(TypeCast<Dest,Src>::cast(x)); }
};

template <typename T>
struct TypeCast<std__complex<T>,std__complex<T> > {
    static const bool valid = true;
    typedef TypeCast<std__complex<T>,std__complex<T> > Cast_t;
    typedef std__complex<T> Dest_t;
    static std__complex<T> cast(std__complex<T> x) { return x; }
};

/*******************************************************************************
 * Additional routines for numeric types
 */

int sqr(int x) { return x*x; }
uint sqr(uint x) { return x*x; }
float sqr(float x) { return x*x; }
double sqr(double x) { return x*x; }
long double sqr(long double x) { return x*x; }
template <typename T>
T sqr(std__complex<T> x) { T r=real(x),i=imag(x); return r*r+i*i; }

int max(int a,int b) { return a>b?a:b; }
uint max(uint a,uint b) { return a>b?a:b; }
float max(float a,float b) { return a>b?a:b; }
double max(double a,double b) { return a>b?a:b; }
long double max(long double a,long double b) { return a>b?a:b; }

/*******************************************************************************
 * UnaryOp
 */


template <class OpTag,typename T=AnyType> struct UnaryOp;

template <class OpTag> struct UnaryOp<OpTag,AnyType> {
    template <typename T>
    static typename UnaryOp<OpTag,T>::Result_t apply(T a) { return UnaryOp<OpTag,T>::apply(a); }
};

template <class OpTag,typename RES,typename T> struct _RawUnOp;
template <typename RES,typename T> struct _RawUnOp<OpIdent,RES,T> { static RES apply(T a) { return a; } };
template <typename RES,typename T> struct _RawUnOp<OpNeg  ,RES,T> { static RES apply(T a) { return -a; } };
template <typename RES,typename T> struct _RawUnOp<OpConj ,RES,T> { static RES apply(T a) { return conj(a); } };
template <typename RES,typename T> struct _RawUnOp<OpExp  ,RES,T> { static RES apply(T a) { return exp(a); } };
template <typename RES,typename T> struct _RawUnOp<OpAbs  ,RES,T> { static RES apply(T a) { return abs(a); } };
template <typename RES,typename T> struct _RawUnOp<OpReal ,RES,T> { static RES apply(T a) { return real(a); } };
template <typename RES,typename T> struct _RawUnOp<OpImag ,RES,T> { static RES apply(T a) { return imag(a); } };
template <typename RES,typename T> struct _RawUnOp<OpSqr  ,RES,T> { static RES apply(T a) { return sqr(a); } };


#define DEFINE_UNARYOP(REAL) \
template <class OpTag> \
struct UnaryOp<OpTag,REAL> { \
    typedef REAL Result_t; \
    static Result_t apply(REAL a) { return _RawUnOp<OpTag,Result_t,REAL>::apply(a); } \
};

DEFINE_UNARYOP(int)
DEFINE_UNARYOP(float)
DEFINE_UNARYOP(double)
DEFINE_UNARYOP(long double)
#undef DEFINE_UNARYOP


template <typename T> struct UnaryOp<OpIdent,std__complex<T> > { typedef std__complex<T> Result_t; static Result_t apply(std__complex<T> a) { return a; } };
template <typename T> struct UnaryOp<OpNeg,std__complex<T> >   { typedef std__complex<T> Result_t; static Result_t apply(std__complex<T> a) { return -a; } };
template <typename T> struct UnaryOp<OpConj,std__complex<T> >  { typedef std__complex<T> Result_t; static Result_t apply(std__complex<T> a) { return conj(a); } };
template <typename T> struct UnaryOp<OpExp,std__complex<T> >   { typedef std__complex<T> Result_t; static Result_t apply(std__complex<T> a) { return exp(a); } };

template <typename T> struct UnaryOp<OpAbs,std__complex<T> >   { typedef T Result_t; static Result_t apply(std__complex<T> a) { return abs(a); } };
template <typename T> struct UnaryOp<OpReal,std__complex<T> >  { typedef T Result_t; static Result_t apply(std__complex<T> a) { return real(a); } };
template <typename T> struct UnaryOp<OpImag,std__complex<T> >  { typedef T Result_t; static Result_t apply(std__complex<T> a) { return imag(a); } };
template <typename T> struct UnaryOp<OpSqr,std__complex<T> >   { typedef T Result_t; static Result_t apply(std__complex<T> a) { return sqr(a); } };


/*******************************************************************************
 * BinaryOp
 */

template <class OpTag,typename T1=AnyType,typename T2=AnyType> struct BinaryOp;

template <class OpTag> struct BinaryOp<OpTag,AnyType,AnyType> {
    template <typename T1,typename T2>
    static typename BinaryOp<OpTag,T1,T2>::Result_t apply(T1 a,T2 b) { return BinaryOp<OpTag,T1,T2>::apply(a,b); }
};

template <> struct BinaryOp<OpRevMul,AnyType,AnyType> {
    template <typename T1,typename T2>
    static typename BinaryOp<OpMul,T2,T1>::Result_t apply(T1 a,T2 b) { return BinaryOp<OpMul,T2,T1>::apply(b,a); }
};

template <typename T1,typename T2>
struct BinaryOp<OpRevMul,T1,T2> {
    typedef typename BinaryOp<OpMul,T2,T1>::Result_t Result_t;
    static Result_t apply(T1 a,T2 b) { return BinaryOp<OpMul,T2,T1>::apply(cast1(b),cast2(a)); }
};

template <typename T>
struct BinaryOp<OpRevMul,T,T> {
    typedef typename BinaryOp<OpMul,T,T>::Result_t Result_t;
    static Result_t apply(T a,T b) { return BinaryOp<OpMul,T,T>::apply(cast1(b),cast2(a)); }
};



template <class OpTag,typename RES,typename T1,typename T2> struct _RawBinOp;
template <typename RES,typename T1,typename T2> struct _RawBinOp<OpAdd,RES,T1,T2> { static RES apply(T1 a,T2 b) { return a + b; } };
template <typename RES,typename T1,typename T2> struct _RawBinOp<OpSub,RES,T1,T2> { static RES apply(T1 a,T2 b) { return a - b; } };
template <typename RES,typename T1,typename T2> struct _RawBinOp<OpMul,RES,T1,T2> { static RES apply(T1 a,T2 b) { return a * b; } };
template <typename RES,typename T1,typename T2> struct _RawBinOp<OpDiv,RES,T1,T2> { static RES apply(T1 a,T2 b) { return a / b; } };
template <typename RES,typename T1,typename T2> struct _RawBinOp<OpRevMul,RES,T1,T2> { static RES apply(T1 a,T2 b) { return b * a; } };

#define DEFINE_IDTYPE_BINARYOP(TYPE)                                                             \
template <class OpTag>                                                                           \
struct BinaryOp<OpTag,TYPE,TYPE> {                                                               \
    typedef TYPE Result_t;                                                                       \
    static Result_t apply(TYPE a,TYPE b) { return _RawBinOp<OpTag,TYPE,TYPE,TYPE>::apply(a,b); } \
};

#define DEFINE_BINARYOP(RES,T1,T2)                                                              \
template <typename OpTag>                                                                       \
struct BinaryOp<OpTag,T1,T2> {                                                                  \
    typedef RES Result_t;                                                                       \
    static Result_t apply(T1 a,T2 b) { return BinaryOp<OpTag,RES,RES>::apply((RES)a,(RES)b); }  \
};

#define DEFINE_IDTYPE_BINARYOPS(TYPE)                       \
DEFINE_IDTYPE_BINARYOP(TYPE)                                \
DEFINE_BINARYOP(std__complex<TYPE>,std__complex<TYPE>,TYPE) \
DEFINE_BINARYOP(std__complex<TYPE>,TYPE,std__complex<TYPE>)

#define DEFINE_BINARYOPS(HIGH,LOW)                          \
DEFINE_BINARYOP(HIGH,HIGH,LOW)                              \
DEFINE_BINARYOP(HIGH,LOW,HIGH)                              \
DEFINE_BINARYOP(std__complex<HIGH>,std__complex<HIGH>,LOW)  \
DEFINE_BINARYOP(std__complex<HIGH>,std__complex<LOW>,HIGH)  \
DEFINE_BINARYOP(std__complex<HIGH>,HIGH,std__complex<LOW>)  \
DEFINE_BINARYOP(std__complex<HIGH>,LOW,std__complex<HIGH>)

DEFINE_IDTYPE_BINARYOPS(int)
DEFINE_IDTYPE_BINARYOPS(float)
DEFINE_IDTYPE_BINARYOPS(double)
DEFINE_IDTYPE_BINARYOPS(long double)

DEFINE_BINARYOPS(float,int)
DEFINE_BINARYOPS(double,int)
DEFINE_BINARYOPS(double,float)
DEFINE_BINARYOPS(long double,int)
DEFINE_BINARYOPS(long double,float)
DEFINE_BINARYOPS(long double,double)

#undef DEFINE_BINARYOPS
#undef DEFINE_IDTYPE_BINARYOPS
#undef DEFINE_BINARYOP
#undef DEFINE_IDTYPE_BINARYOP

/*
template <class OpTag,typename T1,typename T2>
struct BinaryOp<OpTag,std__complex<T1>,T2> {
    typedef typename BinaryOp<OpTag,T1,T2>::Result_t Real_t;
    typedef std__complex<Real_t> Result_t;
    static Result_t cast1(std__complex<T1> a){return TypeCast<Result_t,std__complex<T1> >::cast(a);}
    static Real_t cast2(T2 a){return TypeCast<Real_t,T2>::cast(a);}
    static Result_t apply(std__complex<T1> a,T2 b) { return _RawBinOp<OpTag,Result_t,Result_t,Real_t>::apply(cast1(a),cast2(b)); }
};
template <class OpTag,typename T1,typename T2>
struct BinaryOp<OpTag,T1,std__complex<T2> > {
    typedef typename BinaryOp<OpTag,T1,T2>::Result_t Real_t;
    typedef std__complex<Real_t> Result_t;
    static Real_t cast1(T1 a){return TypeCast<Real_t,T1>::cast(a);}
    static Result_t cast2(std__complex<T2> a){return TypeCast<Result_t,std__complex<T2> >::cast(a);}
    static Result_t apply(T1 a,std__complex<T2> b) { return _RawBinOp<OpTag,Result_t,Real_t,Result_t>::apply(cast1(a),cast2(b)); }
};
*/
template <class OpTag,typename T1,typename T2>
struct BinaryOp<OpTag,std__complex<T1>,std__complex<T2> > {
    typedef std__complex<typename BinaryOp<OpTag,T1,T2>::Result_t> Result_t;
    static Result_t cast1(std__complex<T1> a){return TypeCast<Result_t,std__complex<T1> >::cast(a);}
    static Result_t cast2(std__complex<T2> a){return TypeCast<Result_t,std__complex<T2> >::cast(a);}
    static Result_t apply(std__complex<T1> a,std__complex<T2> b) { return _RawBinOp<OpTag,Result_t,Result_t,Result_t>::apply(cast1(a),cast2(b)); }
};
template <class OpTag,typename T>
struct BinaryOp<OpTag,std__complex<T>,std__complex<T> > {
    typedef std__complex<T> Result_t;
    static Result_t apply(std__complex<T> a,std__complex<T> b) { return _RawBinOp<OpTag,Result_t,Result_t,Result_t>::apply(a,b); }
};

/*******************************************************************************
 * SizeCombine
 */
template <int N1,int N2>
struct SizeCombine;

template <int N>
struct SizeCombine<N,N> { static const int n = N; };

template <int N>
struct SizeCombine<0,N> { static const int n = N; };

template <int N>
struct SizeCombine<N,0> { static const int n = N; };

template <>
struct SizeCombine<0,0> { static const int n = 0; };




#define DEFINE_ASSIGNMENT(TD,OP,OPNAME)        \
template <typename TS>                         \
inline void OPNAME(TD &restrict dest,TS src) { \
    dest OP TypeCast<TD,TS>::cast(src);        \
}

#define DEFINE_ASSIGNMENTS(TD)      \
DEFINE_ASSIGNMENT(TD,=,assign)      \
DEFINE_ASSIGNMENT(TD,+=,assign_add) \
DEFINE_ASSIGNMENT(TD,-=,assign_sub) \
DEFINE_ASSIGNMENT(TD,*=,assign_mul) \
DEFINE_ASSIGNMENT(TD,/=,assign_div)

DEFINE_ASSIGNMENTS(int)
DEFINE_ASSIGNMENTS(float)
DEFINE_ASSIGNMENTS(double)
DEFINE_ASSIGNMENTS(long double)

#undef DEFINE_ASSIGNMENTS
#undef DEFINE_ASSIGNMENT

#define DEFINE_CPX_ASSIGNMENT(OP,OPNAME)        \
template <typename TD,typename TS>                         \
inline void OPNAME(std__complex<TD> &restrict dest,TS src) { \
    dest OP TypeCast<std__complex<TD>,TS>::cast(src);        \
}

DEFINE_CPX_ASSIGNMENT(=,assign)      \
DEFINE_CPX_ASSIGNMENT(+=,assign_add) \
DEFINE_CPX_ASSIGNMENT(-=,assign_sub) \
DEFINE_CPX_ASSIGNMENT(*=,assign_mul) \
DEFINE_CPX_ASSIGNMENT(/=,assign_div)

#undef DEFINE_CPX_ASSIGNMENT

} // namespace etLAP

#endif // _ETLAP_TYPES_H_
