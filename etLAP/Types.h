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
};
#define std__complex std::complex
#endif

#include "Util.h"
#include "Tags.h"

namespace etLAP {

#if 0
template <typename T> SafeCast {};
template <> struct SafeCast<int> {
    static int cast(int a) { return a; }
};
template <> struct SafeCast<float> {
    static float cast(int a) { return a; }
    static float cast(double a) { return a; }
};
template <> struct SafeCast<double> {
    static double cast(int a) { return a; }
    static double cast(float a) { return a; }
    static double cast(double a) { return a; }
};
template <T> struct SafeCast<std__complex<T> > {
    static std__complex<T> cast(std_complex<T> a) { return a; }
    template <T1>
    static std__complex<T> cast(std_complex<T1> a) { return std__complex<T>(SafeCast<T>::cast(real(a)),SafeCast<T>::cast(real(a)));
    template <T1>
    static std__complex<T> cast(T1 a) { return std__complex<T>(SafeCast<T>::cast(a)); };
};
#endif

/*******************************************************************************
 * TypeCast
 */
class NoCast {};
typedef class NoClass {} *NoType;
template <typename Dest,typename Src> struct TypeCast {
    enum { valid = false };
    typedef NoCast Cast_t;
    typedef NoType Dest_t;
};

template <typename T> struct TypeCast<T,T> {
    enum { valid = true };
    typedef TypeCast<T,T> Cast_t;
    typedef T      Dest_t;
    static T cast(T x) { return x; }
};

#define DEFINE_NUMERIC_CAST(TO,FROM)     \
template <> struct TypeCast<TO,FROM>     { \
    enum { valid = true };                 \
    typedef TypeCast<TO,FROM> Cast_t;      \
    typedef TO Dest_t;                     \
    static TO cast(FROM x) { return x; }   \
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
    enum { valid = TypeCast<Dest,Src>::valid };
    typedef typename IF<valid,TypeCast<std__complex<Dest>,std__complex<Src> >,NoCast>::T Cast_t;
    typedef typename IF<valid,std__complex<typename TypeCast<Dest,Src>::Dest_t>,NoType>::T Dest_t;
    static std__complex<Dest> cast(std__complex<Src> x) { return std__complex<Dest>(TypeCast<Dest,Src>::cast(real(x)),TypeCast<Dest,Src>::cast(imag(x))); };
};

template <typename Dest,typename Src>
struct TypeCast<std__complex<Dest>,Src > {
    enum { valid = TypeCast<Dest,Src>::valid };
    typedef typename IF<valid,TypeCast<std__complex<Dest>,Src >,NoCast>::T Cast_t;
    typedef typename IF<valid,std__complex<typename TypeCast<Dest,Src>::Dest_t>,NoType>::T Dest_t;
    static std__complex<Dest> cast(Src x) { return std__complex<Dest>(TypeCast<Dest,Src>::cast(x)); };
};

template <typename T>
struct TypeCast<std__complex<T>,std__complex<T> > {
    enum { valid = true };
    typedef TypeCast<std__complex<T>,std__complex<T> > Cast_t;
    typedef std__complex<T> Dest_t;
    static std__complex<T> cast(std__complex<T> x) { return x; };
};

/*******************************************************************************
 * UnaryResult
 */

template <typename T,class OpTag> struct UnaryResult { typedef T Result_t; };

template <typename REAL> struct UnaryResult<std__complex<REAL>,OpAbs> { typedef REAL Result_t; };

/*******************************************************************************
 * TypeCombine
 */
template <typename T1,typename T2,class OpTag> struct TypeCombine;

template <typename T,class OpTag> struct TypeCombine<T,T,OpTag> { typedef T Result_t; static T cast(T a) {return a;};};

#define DEFINE_NUMERIC_TYPE_COMBINE(RES,T1,T2) \
template <class OpTag>                         \
struct TypeCombine<T1,T2,OpTag> {              \
    typedef RES Result_t;                    \
    static RES cast(T1 a) { return a; };       \
    static RES cast(T2 a) { return a; };       \
};

#define DEFINE_NUMERIC_TYPE_COMBINES(HIGH,LOW) \
DEFINE_NUMERIC_TYPE_COMBINE(HIGH,HIGH,LOW)     \
DEFINE_NUMERIC_TYPE_COMBINE(HIGH,LOW,HIGH)

DEFINE_NUMERIC_TYPE_COMBINES(float,int)
DEFINE_NUMERIC_TYPE_COMBINES(double,int)
DEFINE_NUMERIC_TYPE_COMBINES(double,float)
DEFINE_NUMERIC_TYPE_COMBINES(long double,int)
DEFINE_NUMERIC_TYPE_COMBINES(long double,float)
DEFINE_NUMERIC_TYPE_COMBINES(long double,double)

#undef DEFINE_NUMERIC_TYPE_COMBINES
#undef DEFINE_NUMERIC_TYPE_COMBINE

template <typename T1,typename T2,class OpTag>
struct TypeCombine<std__complex<T1>,T2,OpTag> {
    typedef typename TypeCombine<T1,T2,OpTag>::Result_t Real_t;
    typedef std__complex<Real_t> Result_t;
    static Result_t cast(std__complex<T1> a){return TypeCast<Result_t,std__complex<T1> >::cast(a);};
    static Real_t cast(T2 a){return TypeCast<Real_t,T2>::cast(a);};
};
template <typename T1,typename T2,class OpTag>
struct TypeCombine<T1,std__complex<T2>,OpTag> {
    typedef typename TypeCombine<T1,T2,OpTag>::Result_t Real_t;
    typedef std__complex<Real_t> Result_t;
    static Real_t cast(T1 a){return TypeCast<Real_t,T1>::cast(a);};
    static Result_t cast(std__complex<T2> a){return TypeCast<Result_t,std__complex<T2> >::cast(a);};
};
template <typename T1,typename T2,class OpTag>
struct TypeCombine<std__complex<T1>,std__complex<T2>,OpTag> {
    typedef std__complex<typename TypeCombine<T1,T2,OpTag>::Result_t> Result_t;
    static Result_t cast(std__complex<T1> a){return TypeCast<Result_t,std__complex<T1> >::cast(a);};
    static Result_t cast(std__complex<T2> a){return TypeCast<Result_t,std__complex<T2> >::cast(a);};
};
template <typename T,class OpTag>
struct TypeCombine<std__complex<T>,std__complex<T>,OpTag> {
    typedef std__complex<T> Result_t;
    static std__complex<T> cast(std__complex<T> a){return TypeCast<Result_t,std__complex<T> >::cast(a);};
};

/*******************************************************************************
 * SizeCombine
 */
template <int N1,int N2>
struct SizeCombine;

template <int N>
struct SizeCombine<N,N> { enum { n = N }; };

template <int N>
struct SizeCombine<0,N> { enum { n = N }; };

template <int N>
struct SizeCombine<N,0> { enum { n = N }; };

template <>
struct SizeCombine<0,0> { enum { n = 0 }; };


}; // namespace etLAP

#endif // _ETLAP_TYPES_H_
