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

#ifndef _ETLAP_MATRIX_H_
#define _ETLAP_MATRIX_H_

#include "Common.h"

namespace etLAP {

/*****************************************************************************
 *  Smart Arrays
 */

template <int R,int C,typename T>
class Matrix<T,R,C,Smart>
: public Common_Smart<Matrix<T,R,C,Smart>,T> {
    typedef Smart E;
    typedef Matrix<T,R,C,E> SAME;

    T data[R*C];

  public:
#if ASSIGN_POLICY == AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Matrix<T,R,C,E>() {};
    Matrix<T,R,C,E>(int rw_,int cl_) { assert(rw_ == R && cl_ == C); };

    template<typename X>
    Matrix<T,R,C,E>(const X &restrict src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { assert(R == rw_ && C == cl_); };
    void clear() { for(int n=R*C;n-->0;) data[n] = (T)0; };

    const T operator() (int r,int c) const { return data[r*C+c]; };
    T &operator() (int r,int c,bool unsafe=false) { return data[r*C+c]; };

//    const T operator() (int n) const { return data[n]; };
//    T &operator() (int n,bool unsafe=false) { return data[n]; };

    T *rawdata() { return &data[0]; };
    // be careful with this one -- only use it if you know exactly what you are doing!

    int rows() const { return R; };
    int cols() const { return C; };

    void prepare_write_clone() {}
    void prepare_write_noclone() {}

    void dimension(int rw_,int cl_) { resize(rw_,cl_); };  // introduced for MDP
    int rowmax() { return rows(); };  // introduced for MDP
    T *address() { return data; };  // introduced for MDP
    int size() { return rows()*cols(); };  // introduced for MDP
};

/*****************************************************************************
 *  Smart Arrays - Variable sized
 */

template <typename T>
class Matrix<T,0,0,Smart>
: public Common_Smart<Matrix<T,0,0,Smart>,T> {
    typedef Smart E;
    typedef Matrix<T,0,0,E> SAME;

    int rw,cl;
    Storage<T> data;

  public:
#if ASSIGN_POLICY==AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Matrix<T,0,0,E>(): rw(0),cl(0) {};
    Matrix<T,0,0,E>(int rw_,int cl_): rw(rw_),cl(cl_),data(rw*cl) {};

    template<typename X>
    Matrix<T,0,0,E>(const X &restrict src) { rw=-1; cl=-1; assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { rw=rw_;cl=cl_;data.resize(rw*cl); };
    void clear() { data.clear(); };

    const T operator() (int r,int c) const { return data[r*cl + c]; };
    T &operator() (int r,int c,bool unsafe=false) { if(!unsafe)prepare_write_clone(); return data[r*cl + c]; };

//    const T operator() (int n) const { return data[n]; };
//    T &operator() (int n,bool unsafe) { if(!unsafe)prepare_write_clone(); return data[n]; };

    T *rawdata() { return data.rawdata(); };
    // be careful with this one -- only use it if you know exactly what you are doing!
    
    int rows() const { return rw; }
    int cols() const { return cl; }

    void prepare_write_clone() { data.prepare_write_clone(); }
    void prepare_write_noclone() { data.prepare_write_noclone(); }

    void dimension(int rw_,int cl_) { resize(rw_,cl_); };  // introduced for MDP
    int rowmax() { return rows(); };  // introduced for MDP
    T *address() { data.prepare_write_clone(); return &data[0]; };  // introduced for MDP
    int size() { return rows()*cols(); };  // introduced for MDP
};

/*****************************************************************************
 *  Packed Arrays
 */

template <int R,int C,typename T>
class Matrix<T,R,C,Packed>
: public Common_Packed<Matrix<T,R,C,Packed>,T> {
    typedef Packed E;
    typedef Matrix<T,R,C,E> SAME;

    T data[R*C];

  public:
#if ASSIGN_POLICY == AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Matrix<T,R,C,E>() {};
    Matrix<T,R,C,E>(int rw_,int cl_) { assert(rw_ == R && cl_ == C); };

    template<typename X>
    Matrix<T,R,C,E>(const X &restrict src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { assert(R == rw_ && C == cl_); };
    void clear() { for(int n=R*C;n-->0;) data[n] = (T)0; };

    const T operator() (int r,int c) const { return data[r*C+c]; };
    T &operator() (int r,int c,bool unsafe=false) { return data[r*C+c]; };

//    const T operator() (int n) const { return data[n]; };
//    T &operator() (int n,bool unsafe=false) { return data[n]; };

    T *rawdata() { return &data[0]; };
    // be careful with this one -- only use it if you know exactly what you are doing!

    int rows() const { return R; };
    int cols() const { return C; };

    void prepare_write_clone() {}
    void prepare_write_noclone() {}
};

/*****************************************************************************
 *  Packed Arrays - Variable sized - don't exist!!
 */

template <typename T>
class Matrix<T,0,0,Packed> {
    Matrix<T,0,0,Packed>();
};

/*****************************************************************************
 *  ExtPtr - fixed size
 */

template <int R,int C,typename T>
class Matrix<T,R,C,ExtPtr>
: public Common_Smart<Matrix<T,R,C,ExtPtr>,T> {
    typedef ExtPtr E;
    typedef Matrix<T,R,C,E> SAME;

    T *ptr;

  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(T *ptr_): ptr(ptr_) { };
    Matrix<T,R,C,E>(T *ptr_,int rw_,int cl_): ptr(ptr_) { assert(rw_==R && cl_ == C); }; // for usage in mdp

    Matrix<T,R,C,E>(const SAME &restrict other): ptr(other.ptr) {};
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { assert(R == rw_ && C == cl_); };
    void clear() { for(int i=R*C;i-->0;) ptr[i] = (T)0; };

    const T operator() (int r,int c) const { return ptr[r*C + c]; };
    T &operator() (int r,int c,bool unsafe=false) { return ptr[r*C + c]; };

    const T operator() (int n) const { return ptr[n]; };
    T &operator() (int n,bool unsafe) { return ptr[n]; };

    void prepare_write_clone() {}
    void prepare_write_noclone() {}

    int rows() const { return R; };
    int cols() const { return C; };
};


/*****************************************************************************
 *  ExtPtr - Variable size
 */

template <typename T>
class Matrix<T,0,0,ExtPtr>
: public Common_Smart<Matrix<T,0,0,ExtPtr>,T> {
    typedef ExtPtr E;
    typedef Matrix<T,0,0,E> SAME;

    int rw,cl;
    T *ptr;

  public:
    typedef SAME Ref_t;

    Matrix<T,0,0,E>(T *ptr_,int rw_,int cl_): rw(rw_),cl(cl_),ptr(ptr_) {};
    Matrix<T,0,0,E>(const SAME &restrict other): rw(other.rw),cl(other.cl),ptr(other.ptr) {};
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { rw = rw_; cl = cl_; };
    void clear() { for(int i=rw*cl;i-->0;) ptr[i] = (T)0; };

    const T operator() (int r,int c) const { return ptr[r*cl + c]; };
    T &operator() (int r,int c,bool unsafe=false) { return ptr[r*cl + c]; };

    const T operator() (int n) const { return ptr[n]; };
    T &operator() (int n,bool unsafe) { return ptr[n]; };

    void prepare_write_clone() {}
    void prepare_write_noclone() {}

    int rows() const { return rw; }
    int cols() const { return cl; }
};



// Matrix = Matrix

template <typename TD,int RD,int CD,class ED,typename TS,int RS,int CS,class ES>
inline void assign(Matrix<TD,RD,CD,ED> &restrict dest,const Matrix<TS,RS,CS,ES> &restrict src,NoCast) {
    CTAssert((RD == RS && CD == CS) || (RD==0 && CD==0) || (RS==0 && CS==0));
    dest.resize(src.rows(),src.cols());
    dest.prepare_write_noclone();
#ifndef FLATLOOP
    for(int r=dest.rows();r-->0;)
    for(int c=dest.cols();c-->0;)
        dest(r,c,true) = TypeCast<TD,TS>::cast(src(r,c));
#else
    for(int n=dest.rows()*dest.cols();n-->0;)
        dest(n,true) = TypeCast<TD,TS>::cast(src(n));
#endif
};

// Matrix += Matrix

template <typename TD,int RD,int CD,class ED,typename TS,int RS,int CS,class ES>
inline void assign_add(Matrix<TD,RD,CD,ED> &restrict dest,const Matrix<TS,RS,CS,ES> &restrict src) {
    CTAssert((RD == RS && CD == CS) || (RD==0 && CD==0));
    assert(dest.rows() == src.rows() && dest.cols() == src.cols());
    dest.prepare_write_clone();
#ifndef FLATLOOP
    for(int r=dest.rows();r-->0;)
    for(int c=dest.cols();c-->0;)
        dest(r,c,true) += TypeCast<TD,TS>::cast(src(r,c));
#else
    for(int n=dest.rows()*dest.cols();n-->0;)
        dest(n,true) += TypeCast<TD,TS>::cast(src(n);
#endif
};

// Matrix -= Matrix

template <typename TD,int RD,int CD,class ED,typename TS,int RS,int CS,class ES>
inline void assign_sub(Matrix<TD,RD,CD,ED> &restrict dest,const Matrix<TS,RS,CS,ES> &restrict src) {
    CTAssert((RD == RS && CD == CS) || (RD==0 && CD==0));
    assert(dest.rows() == src.rows() && dest.cols() == src.cols());
    dest.prepare_write_clone();
#ifndef FLATLOOP
    for(int r=dest.rows();r-->0;)
    for(int c=dest.cols();c-->0;)
        dest(r,c,true) -= TypeCast<TD,TS>::cast(src(r,c));
#else
    for(int n=dest.rows()*dest.cols();n-->0;)
        dest(n,true) -= TypeCast<TD,TS>::cast(src(n);
#endif
};

// Matrix = Scalar

template <int N,typename TD,class ED,typename TS>
inline void assign(Matrix<TD,N,N,ED> &restrict dest,const TS &restrict s,TypeCast<TD,TS>) {
    assert(dest.rows() == dest.cols());
    dest.clear();
    for(int n=dest.rows();n-->0;)
        dest(n,n,true) = TypeCast<TD,TS>::cast(s);
};

// Matrix *= Scalar

template <typename TD,int R,int C,class E,typename TS>
inline void assign_mul(Matrix<TD,R,C,E> &restrict dest,TS s) {
    dest.prepare_write_clone();
#ifndef FLATLOOP
    for(int r=dest.rows();r-->0;)
    for(int c=dest.cols();c-->0;)
        dest(r,c,true) *= TypeCast<TD,TS>::cast(s);
#else
    for(int n=dest.rows()*dest.cols();n-->0;)
        dest(n,true) *= TypeCast<TD,TS>::cast(s);
#endif
};

// Matrix /= Scalar

template <typename TD,int R,int C,class E,typename TS>
inline void assign_div(Matrix<TD,R,C,E> &restrict dest,TS s) {
    dest.prepare_write_clone();
#ifndef FLATLOOP
    for(int r=dest.rows();r-->0;)
    for(int c=dest.cols();c-->0;)
        dest(r,c,true) /= TypeCast<TD,TS>::cast(s);
#else
    for(int n=dest.rows()*dest.cols()-1;n>=0;n--)
        dest(n,true) /= TypeCast<TD,TS>::cast(s);
#endif
};

// Matrix<N,N> = Vector<N>

template <typename TD,int ND,class ED,typename TS,int NS,class ES>
inline void assign(Matrix<TD,ND,ND,ED> &restrict dest,const Vector<TS,NS,ES> &restrict src,NoCast) {
    CTAssert(ND == NS || ND==0);
    dest.resize(src.size(),src.size());
    dest.clear();
    for(int n=src.size();n-->0;)
        dest(n,n,true) = TypeCast<TD,TS>::cast(src(n));
};

// Matrix<N,1> = Vector<N>

template <typename TD,int RD,class ED,typename TS,int RS,class ES>
inline void assign(Matrix<TD,RD,1,ED> &restrict dest,const Vector<TS,RS,ES> &restrict src,NoCast) {
    CTAssert(RS == 0 || RD == RS);
    assert(dest.rows() == src.size());
    for(int r=dest.rows();r-->0;)
#ifndef FLATLOOP
        dest(r,0,true) = TypeCast<TD,TS>::cast(src(r));
#else
        dest(r,true) = TypeCast<TD,TS>::cast(src(r));
#endif
};

// Matrix<1,N> = Vector<N>

template <typename TD,int CD,class ED,typename TS,int CS,class ES>
inline void assign(Matrix<TD,1,CD,ED> &restrict dest,const Vector<TS,CS,ES> &restrict src,NoCast) {
    CTAssert(CS == 0 || CD == CS);
    assert(dest.cols() == src.size());
    for(int c=dest.cols();c-->0;)
#ifndef FLATLOOP
        dest(0,c,true) = TypeCast<TD,TS>::cast(src(c));
#else
        dest(c,true) = TypeCast<TD,TS>::cast(src(c));
#endif
};



#if 0
/*****************************************************************************
 *  Diagonal
 */
template <int N,typename T>
class Matrix<N,N,T,Diagonal>
: public Common_Smart<Matrix<N,N,T,Diagonal>,T> {
    typedef Diagonal E;

    T data[N];
  public:
//    typedef const SAME &Ref_t;
    typedef SAME Ref_t;

    Matrix<N,N,T,E>() {};
//    Matrix<N,N,T,E>(const SAME &restrict other): data(other.data) {};
    template<typename X>
    Matrix<N,N,T,E>(const X &restrict src) { assign_from(src); };
    template<typename X>
    SAME &operator=(const X &src) { assign_from(src); return *(SAME *)this; };

    T operator() (int r,int c) const { assert(!locked); return r==c?data[r]:(T)0; }
    T &operator() (int n) { return data[n]; }
};

template <int N,typename TD,typename TS,class E>
inline void assign(Matrix<N,N,TD,Diagonal> &dest,const Vector<N,TS,E> &src,NoCast) {
    for(int n=0;n<N;n++)
        dest(n) = TypeCast<TD,TS>::cast(src(n));
};

template <int N,typename TD,typename TS,class E>
inline void assign(Matrix<N,N,TD,Diagonal> &dest,const Matrix<N,N,TS,Diagonal> &src,NoCast) {
    for(int n=0;n<N;n++)
        dest(n) = TypeCast<TD,TS>::cast(src(n,n));
};

template <int N,typename TD,typename TS>
inline void assign(Matrix<N,N,TD,Diagonal> &dest,TS s,TypeCast<TD,TS>) {
    for(int n=0;n<N;n++)
        dest(n) = TypeCast<TD,TS>::cast(s);
};
#endif

/*****************************************************************************
 *  Scalar
 */
template <typename T,int N>
class Matrix<T,N,N,Scalar>
: public Common_Smart<Matrix<T,N,N,Scalar>,T> {
    typedef Scalar E;
    typedef Matrix<T,N,N,E> SAME;

    T s;
  public:
    typedef SAME Ref_t;

//    Matrix<T,N,N,E>() {};
//    Matrix<T,N,N,E>(int sz_) { assert(sz_ == N); };
    Matrix<T,N,N,E>(T s_): s(s_) {};
    Matrix<T,N,N,E>(int sz_,T s_): s(s_) { assert(sz_ == N); };

    template<typename X>
    Matrix<T,N,N,E>(const X &restrict src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    const T operator() (int r,int c) const { return r==c?s:(T)0; }
    T &operator() () { return s; }
    const T operator() (int n) const { return (n%(N+1)==1)?s:(T)0; }

    int rows() const { return N; }
    int cols() const { return N; }
};

template <typename T>
class Matrix<T,0,0,Scalar>
: public Common_Smart<Matrix<T,0,0,Scalar>,T> {
    typedef Scalar E;
    typedef Matrix<T,0,0,E> SAME;

    T s;
    int sz;
  public:
    typedef SAME Ref_t;

    Matrix<T,0,0,E>(): sz(0) {};
//    Matrix<T,0,0,E>(int sz_): sz(sz_) {};
//    Matrix<T,0,0,E>(T s_): sz(sz_), s(s_) {};
    Matrix<T,0,0,E>(int sz_,T s_): sz(sz_),s(s_) {};

    template<typename X>
    Matrix<T,0,0,E>(const X &restrict src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    const T operator() (int r,int c) const { return r==c?s:(T)0; }
    T &operator() () { return s; }
    const T operator() (int n) const { return (n%(N+1)==1)?s:(T)0; }

    int rows() const { return sz; }
    int cols() const { return sz; }
};

template <int N,typename TD,typename TS>
inline void assign(Matrix<TD,N,N,Scalar> &restrict dest,TS s,TypeCast<TD,TS>*) {
    dest() = TypeCast<TD,TS>::cast(s);
};

/*****************************************************************************
 *  Trivials
 */

template <int R,int C,typename T>
class Matrix<T,R,C,Zero>
: public Common<Matrix<T,R,C,Zero> > {
    typedef Matrix<T,R,C,Zero> SAME;
  public:
    typedef SAME Ref_t;
    Matrix<T,R,C,Zero>() {};
    const T operator() (int r,int c) const { return (T)0; }
    const T operator() (int n) const { return (T)0; }
    int rows() const { return R; };
    int cols() const { return C; };
    bool is_locked() const { return false; };
};

template <typename T>
class Matrix<T,0,0,Zero>
: public Common<Matrix<T,0,0,Zero> > {
    typedef Matrix<T,0,0,Zero> SAME;

    int rw,cl;
  public:
    typedef SAME Ref_t;
    Matrix<T,0,0,Zero>(int rw_,int cl_): rw(rw_),cl(cl_) {};
    const T operator() (int r,int c) const { return (T)0; }
    const T operator() (int n) const { return (T)0; }
    int rows() const { return rw; };
    int cols() const { return cl; };
    bool is_locked() const { return false; };
};

template <int N,typename T>
class Matrix<T,N,N,One>
: public Common<Matrix<T,N,N,One> > {
    typedef Matrix<T,N,N,One> SAME;
  public:
    typedef SAME Ref_t;
    Matrix<T,N,N,One>() {};
    const T operator() (int r,int c) const { return r==c?(T)1:(T)0; }
    const T operator() (int n) const { return (n%(N+1)==1)?(T)1:(T)0; }
    int rows() const { return N; };
    int cols() const { return N; };
    bool is_locked() const { return false; };
};

template <typename T>
class Matrix<T,0,0,One>
: public Common<Matrix<T,0,0,One> > {
    typedef Matrix<T,0,0,One> SAME;

    int sz;
  public:
    typedef SAME Ref_t;
    Matrix<T,0,0,One>(int sz_): sz(sz_) {};
    const T operator() (int r,int c) const { return r==c?(T)1:(T)0; }
    const T operator() (int n) const { return (n%(sz+1)==1)?(T)1:(T)0; }
    int rows() const { return sz; };
    int cols() const { return sz; };
    bool is_locked() const { return false; };
};


/*****************************************************************************
 *  Expressions
 */

template <int R,int C,typename T,typename T1,class E1,class Op>
class Matrix<T,R,C,ElemUnOp<Matrix<T1,R,C,E1>,Op> >
: public Common<Matrix<T,R,C,ElemUnOp<Matrix<T1,R,C,E1>,Op> > > {
    typedef Matrix<T1,R,C,E1> M;
    typedef ElemUnOp<M,Op> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M::Ref_t m;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M &restrict m_) : m(m_) {}
    const T operator() (int r,int c) const { return apply((Op*)0,m(r,c)); }
    const T operator() (int n) const { return apply((Op*)0,m(n)); }
    int rows() const { return m.rows(); };
    int cols() const { return m.cols(); };
    bool is_locked() const { return m.is_locked(); };
};

template <typename T,int R,int C,typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2,class Op>
class Matrix<T,R,C,ElemBinOp<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2>,Op> >
: public Common<Matrix<T,R,C,ElemBinOp<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2>,Op> > > {
    typedef Matrix<T1,R1,C1,E1> M1;
    typedef Matrix<T2,R2,C2,E2> M2;
    typedef ElemBinOp<M1,M2,Op> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M1::Ref_t m1;
    const typename M2::Ref_t m2;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M1 &restrict m1_, const M2 &restrict m2_) : m1(m1_), m2(m2_) {}
    const T operator() (int r,int c) const { return apply((Op*)0,m1(r,c), m2(r,c)); }
    const T operator() (int n) const { return apply((Op*)0,m1(n), m2(n)); }
    int rows() const { return m1.rows(); };
    int cols() const { return m1.cols(); };
    bool is_locked() const { return m1.is_locked() || m2.is_locked(); };
};

template <int R,int C,typename T,class E1>
class Matrix<T,R,C,Transposed<Matrix<T,C,R,E1> > >
: public Common<Matrix<T,R,C,Transposed<Matrix<T,C,R,E1> > > > {
    typedef Matrix<T,C,R,E1> M;
    typedef Transposed<M> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M::Ref_t m;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M &restrict m_) : m(m_) {}
    const T operator() (int r,int c) const { return m(c,r); }
    const T operator() (int n) const { return m(n); }
    int rows() const { return m.cols(); };
    int cols() const { return m.rows(); };
    bool is_locked() const { return m.is_locked(); };
};


template <int R,int C,typename T,typename T1,class E1,typename T2,class Op>
class Matrix<T,R,C,ScalarOp<Matrix<T1,R,C,E1>,T2,Op> >
: public Common<Matrix<T,R,C,ScalarOp<Matrix<T1,R,C,E1>,T2,Op> > > {
    typedef Matrix<T1,R,C,E1> M;
    typedef ScalarOp<Matrix<T1,R,C,E1>,T2,Op> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M::Ref_t m;
    T2 s;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M &restrict m_, T2 s_) : m(m_), s(s_) {}
    const T operator() (int r,int c) const { return apply((Op*)0, m(r,c), s); }
    const T operator() (int n) const { return apply((Op*)0, m(n), s); }
    int rows() const { return m.rows(); };
    int cols() const { return m.cols(); };
    bool is_locked() const { return m.is_locked(); };
};

/*****************************************************************************
 *  Row/Column
 */

template <int R,typename T,typename T1,class E1>
class Matrix<T,R,1,Column<Vector<T1,R,E1> > >
: public Common<Matrix<T,R,1,Column<Vector<T1,R,E1> > > > {
    typedef Vector<T1,R,E1> V;
    typedef Column<V> E;
    typedef Matrix<T,R,1,E> SAME;

    const typename V::Ref_t v;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,1,E>(const V &restrict v_) : v(v_) {}
    const T operator() (int r,int c) const { return v(r); }
    const T operator() (int n) const { return v(r); }
    int rows() const { return v.size(); };
    int cols() const { return 1; };
    bool is_locked() const { return v.is_locked(); };
};

template <int C,typename T,typename T1,class E1>
class Matrix<T,1,C,Row<Vector<T1,C,E1> > >
: public Common<Matrix<T,1,C,Row<Vector<T1,C,E1> > > > {
    typedef Vector<T1,C,E1> V;
    typedef Row<V> E;
    typedef Matrix<T,1,C,E> SAME;

    const typename V::Ref_t v;
  public:
    typedef SAME Ref_t;

    Matrix<T,1,C,E>(const V &restrict v_) : v(v_) {};
    const T operator() (int r,int c) const { return v(c); };
    const T operator() (int n) const { return v(c); };
    int rows() const { return 1; };
    int cols() const { return v.size(); };
    bool is_locked() const { return v.is_locked(); };
};

/*****************************************************************************
 *  Buffer
 */

template <int R,int C,typename T>
class Matrix<T,R,C,Buffer>
: public Common<Matrix<T,R,C,Buffer> > {
    typedef Buffer E;
    typedef Matrix<T,R,C,E> SAME;

    const Matrix<T,R,C,Packed> buf;
  public:
    typedef SAME Ref_t;

    template <class E1>
    Matrix<T,R,C,E>(const Matrix<T,R,C,E1> &restrict m_): buf(m_) {};
    const T operator() (int r,int c) const { return buf(r,c); };
    const T operator() (int n) const { return buf(n); };
    int rows() const { return buf.rows(); };
    int cols() const { return buf.cols(); };
    bool is_locked() const { return false; };
};

/*****************************************************************************
 *  Buffer
 */

template <typename T>
class Matrix<T,0,0,Buffer>
: public Common<Matrix<T,0,0,Buffer> > {
    typedef Buffer E;
    typedef Matrix<T,0,0,E> SAME;

    const Matrix<T,0,0,Smart> buf;
  public:
    typedef SAME Ref_t;

    template <class E1>
    Matrix<T,0,0,E>(const Matrix<T,0,0,E1> &restrict m_): buf(m_) {};
    const T operator() (int r,int c) const { return buf(r,c); };
    const T operator() (int n) const { return buf(n); };
    int rows() const { return buf.rows(); };
    int cols() const { return buf.cols(); };
    bool is_locked() const { return false; };
};

/*****************************************************************************
 *  NoBuffer
 */

template <int R,int C,typename T,class E1>
class Matrix<T,R,C,NoBuffer<Matrix<T,R,C,E1> > >
: public Common<Matrix<T,R,C,NoBuffer<Matrix<T,R,C,E1> > > > {
    typedef Matrix<T,R,C,E1> M;
    typedef NoBuffer<M> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M::Ref_t m;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M &restrict m_) : m(m_) {}
    const T operator() (int r,int c) const { return m(r,c); }
    const T operator() (int n) const { return m(n); }
    int rows() const { return m.rows(); };
    int cols() const { return m.cols(); };
    bool is_locked() const { return m.is_locked(); };
};

/*****************************************************************************
 *  Multiplied
 */

template <typename T,int R,int C,typename T1,int R1,int C1,class E1,typename T2,int R2,int C2,class E2>
class Matrix<T,R,C,Multiplied<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > >
: public Common<Matrix<T,R,C,Multiplied<Matrix<T1,R1,C1,E1>,Matrix<T2,R2,C2,E2> > > > {
    typedef Matrix<T1,R1,C1,E1> M1;
    typedef Matrix<T2,R2,C2,E2> M2;
    typedef Multiplied<M1,M2> E;
    typedef Matrix<T,R,C,E> SAME;

    const typename M1::Ref_t m1;
    const typename M2::Ref_t m2;
  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(const M1 &restrict m1_, const M2 &restrict m2_) : m1(m1_), m2(m2_) {}

    const T operator() (int r,int c) const {
        T res = 0;
        for(int x=m1.cols();x-->0;)
            res += m1(r,x)*m2(x,c);
        return res;
    }
    const T operator() (int n) const {
        T res = 0;
        int cs=cols();
        int c=n;
        while((c-=cs)>=0);
        int n1=c+cs*rows();
        int n2=n-c;
        do {
            --n2;
            res += m1(n1)*m2(n2);
            n1-=cs;
        } while(n1>=0);
        return res;
    }

    int rows() const { return m1.rows(); };
    int cols() const { return m2.cols(); };
    bool is_locked() const { return m1.is_locked() || m2.is_locked(); };
};

}; // namespace etLAP

#endif // _ETLAP_MATRIX_H_
