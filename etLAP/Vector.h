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

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "Common.h"
#include "Tupel.h"

namespace etLAP {

/*
template <int N,typename T,class E>
class Vector_Common<Vector<N,T,E> > {
  public:
    template <typename TD>
    operator Matrix<N,1,TD,Column<Vector<N,T,E> > >() {
        return Matrix<N,1,typename TypeCast<TD,T>::Dest_t,Column<Vector<N,T,E> > >(*this);
    };
    template <typename TD>
    operator Matrix<1,N,TD,Row<Vector<N,T,E> > >() {
        return Matrix<1,N,typename TypeCast<TD,T>::Dest_t,Row<Vector<N,T,E> > >(*this);
    };
};
*/

/*****************************************************************************
 *  Smart Array - Fixed size
 */

template <int N,typename T>
class Vector<T,N,Smart>
: public Common_Smart<Vector<T,N,Smart>,T> {
    typedef Smart E;
    typedef Vector<T,N,E> SAME;

    T data[N];
  public:
#if ASSIGN_POLICY == AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Vector<T,N,E>() {};
    template<typename X>
    Vector<T,N,E>(const X &src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &src) { assign_from(src); return *(SAME *)this; };

    void resize(int sz_) { assert(sz_ == N); };
    void clear() { for(int n=N;n-->0;) data[n] = (T)0; };

    const T operator() (int n) const { return data[n]; };
    T &operator() (int n,bool unsafe=false) { return data[n]; };

    T *rawdata() { return &data[0]; };
    // be careful with this one -- only use it if you know exactly what you are doing!

    int size() const { return N; };

    void prepare_write_clone() {};
    void prepare_write_noclone() {};
};

/*****************************************************************************
 *  Smart Array - Variable size
 */

template <typename T>
class Vector<T,0,Smart>
: public Common_Smart<Vector<T,0,Smart>,T> {
    typedef Smart E;
    typedef Vector<T,0,E> SAME;

    int sz;
    Storage<T> data;

  public:
#if ASSIGN_POLICY==AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Vector<T,0,E>(): sz(0) {};
    Vector<T,0,E>(int sz_): sz(sz_),data(sz) {};

    template<typename X>
    Vector<T,0,E>(const X &src) { sz=-1; assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &src) { assign_from(src); return *(SAME *)this; };

    void resize(int sz_) { sz=sz_;data.resize(sz); };
    void clear() { data.clear(); };

    const T operator() (int n) const { return data[n]; };
    T &operator() (int n,bool unsafe=false) { if(!unsafe)prepare_write_clone(); return data[n]; };

    T *rawdata() { return data.rawdata(); };
    // be careful with this one -- only use it if you know exactly what you are doing!

    int size() const { return sz; }

    void prepare_write_clone() { data.prepare_write_clone(); }
    void prepare_write_noclone() { data.prepare_write_noclone(); }
};

/*****************************************************************************
 *  Packed Array - Fixed size
 */

template <int N,typename T>
class Vector<T,N,Packed>
: public Common_Packed<Vector<T,N,Packed>,T> {
    typedef Packed E;
    typedef Vector<T,N,E> SAME;

    T data[N];
  public:
#if ASSIGN_POLICY == AP_norefs
    typedef SAME Ref_t;
#else
    typedef const SAME &Ref_t;
#endif

    Vector<T,N,E>() {};
    template<typename X>
    Vector<T,N,E>(const X &src) { assign_from(src); };
    template<typename X>
    const SAME &operator=(const X &src) { assign_from(src); return *(SAME *)this; };

    void resize(int sz_) { assert(sz_ == N); };
    void clear() { for(int n=N;n-->0;) data[n] = (T)0; };

    const T operator() (int n) const { return data[n]; };
    T &operator() (int n,bool unsafe=false) { return data[n]; };

    T *rawdata() { return &data[0]; };
    // be careful with this one -- only use it if you know exactly what you are doing!

    int size() const { return N; };

    void prepare_write_clone() {};
    void prepare_write_noclone() {};
};

/*****************************************************************************
 *  Packed Array - Variable size - don't exist!!
 */

template <typename T>
class Vector<T,0,Packed> {
    Vector<T,0,Packed>();
  public:
    void dummy();
};

// Vector = Vector

template <typename TD,int ND,class ED,typename TS,int NS,class ES,class CAST_TAG>
inline void assign(Vector<TD,ND,ED> &dest,const Vector<TS,NS,ES> &src,CAST_TAG) {
    CTAssert(ND == NS || ND == 0 || NS == 0);
    dest.resize(src.size());
    dest.prepare_write_noclone();
    for(int i=dest.size();i-->0;)
        dest(i,true) = TypeCast<TD,TS>::cast(src(i));
};

// Vector = Tupel

template <typename TD,int ND,class E,typename TS,int NS,class CAST_TAG>
inline void assign(Vector<TD,ND,E> &dest,const Tupel<TS,NS> &src,CAST_TAG) {
    CTAssert(ND == NS || ND == 0);
    dest.resize(NS);
    dest.prepare_write_noclone();
    for(int i=dest.size();i-->0;)
        dest(i,true) = TypeCast<TD,TS>::cast(src[i]);
};

/*****************************************************************************
 *  Trivials
 */

template <int N,typename T>
class Vector<T,N,Zero>
: public Common<Vector<T,N,Zero> > {
    typedef Vector<T,N,Zero> SAME;
  public:
    typedef SAME Ref_t;

    Vector<T,N,Zero>() {};
    T operator() (int n) const { return (T)0; }
    int size() const { return N; }
    bool is_locked() const { return false; };
};

template <typename T>
class Vector<T,0,Zero>
: public Common<Vector<T,N,Zero> > {
    typedef Vector<T,0,Zero> SAME;

    int sz;
  public:
    typedef SAME Ref_t;

    Vector<T,0,Zero>(int sz_): sz(sz_) {};
    T operator() (int n) const { return (T)0; }
    int size() const { return sz; }
    bool is_locked() const { return false; };
};

/*****************************************************************************
 *  Expressions
 */

template <int N,typename T,typename T1,class E1,class Op>
class Vector<T,N,ElemUnOp<Vector<T1,N,E1>,Op> >
: public Common<Vector<T,N,ElemUnOp<Vector<T1,N,E1>,Op> > > {
    typedef Vector<T1,N,E1> V;
    typedef ElemUnOp<V,Op> E;
    typedef Vector<T,N,E> SAME;

    const typename V::Ref_t v;
  public:
    typedef SAME Ref_t;

    Vector<T,N,E>(const V &v_) : v(v_) {}
    T operator() (int n) const { return apply( (Op*)0,v(n) ); }
    int size() const { return v.size(); }
    bool is_locked() const { return v.is_locked(); };
};

template <int N,typename T,typename T1,class E1,typename T2,class E2,class Op>
class Vector<T,N,ElemBinOp<Vector<T1,N,E1>,Vector<T2,N,E2>,Op> >
: public Common<Vector<T,N,ElemBinOp<Vector<T1,N,E1>,Vector<T2,N,E2>,Op> > > {
    typedef Vector<T1,N,E1> V1;
    typedef Vector<T2,N,E2> V2;
    typedef ElemBinOp<V1,V2,Op> E;
    typedef Vector<T,N,E> SAME;

    const typename V1::Ref_t v1;
    const typename V2::Ref_t v2;
  public:
    typedef SAME Ref_t;

    Vector<T,N,E>(const V1 &v1_, const V2 &v2_) : v1(v1_), v2(v2_) {}
    T operator() (int n) const { return apply( (Op*)0,v1(n),v2(n) ); }
    int size() const { return v1.size(); }
    bool is_locked() const { return v1.is_locked() || v2.is_locked(); };
};

template <int N,typename T,typename T1,class E1,typename T2,class Op>
class Vector<T,N,ScalarOp<Vector<T1,N,E1>,T2,Op> >
: public Common<Vector<T,N,ScalarOp<Vector<T1,N,E1>,T2,Op> > > {
    typedef Vector<T1,N,E1> V;
    typedef ScalarOp<V,T2,Op> E;
    typedef Vector<T,N,E> SAME;

    const typename V::Ref_t v;
    T2 s;
  public:
    typedef SAME Ref_t;

    Vector<T,N,E>(const V &v_, T2 s_) : v(v_), s(s_) {}
    T operator() (int n) const { return apply((Op*)0, v(n), s); }
    int size() const { return v.size(); }
    bool is_locked() const { return v.is_locked(); };
};

/*****************************************************************************
 *  Buffer
 */

template <int N,typename T>
class Vector<T,N,Buffer>
: public Common<Vector<T,N,Buffer> > {
    typedef Buffer E;
    typedef Vector<T,N,E> SAME;

    Vector<T,N,Packed> buf;
  public:
    typedef SAME Ref_t;

    template <class E1>
    Vector<T,N,E>(const Vector<T,N,E1> &v_) : buf(v_) {};
    T operator() (int n) const { return buf(n); };
    int size() const { return buf.size(); };
    bool is_locked() const { return false; };
};

template <typename T>
class Vector<T,0,Buffer>
: public Common<Vector<T,0,Buffer> > {
    typedef Buffer E;
    typedef Vector<T,0,E> SAME;

    Vector<T,0,Smart> buf;
  public:
    typedef SAME Ref_t;

    template <class E1>
    Vector<T,0,E>(const Vector<T,0,E1> &v_) : buf(v_) {};
    T operator() (int n) const { return buf(n); };
    int size() const { return buf.size(); };
    bool is_locked() const { return false; };
};

template <int N,typename T,class E1>
class Vector<T,N,NoBuffer<Vector<T,N,E1> > >
: public Common<Vector<T,N,NoBuffer<Vector<T,N,E1> > > > {
    typedef Vector<T,N,E1> V;
    typedef NoBuffer<V> E;
    typedef Vector<T,N,E> SAME;

    const typename V::Ref_t v;
  public:
    typedef SAME Ref_t;

    Vector<T,N,E>(const V &v_) : v(v_) {}
    T operator() (int n) const { return v(n); }
    int size() const { return v.size(); };
    bool is_locked() const { return v.is_locked(); };
};

/*****************************************************************************
 *  Multiplied
 */

template <int R,int C,typename T,typename T1,class E1,typename T2,class E2>
class Vector<T,R,Multiplied<Matrix<T1,R,C,E1>,Vector<T2,C,E2> > >
: public Common<Vector<T,R,Multiplied<Matrix<T1,R,C,E1>,Vector<T2,C,E2> > > > {
    typedef Matrix<T1,R,C,E1> M;
    typedef Vector<T2,C,E2> V;
    typedef Multiplied<M,V> E;
    typedef Vector<T,R,E> SAME;

    const typename M::Ref_t m;
    const typename V::Ref_t v;
  public:
    typedef SAME Ref_t;

    Vector<T,R,E>(const M &m_, const V &v_) : m(m_), v(v_) {}
    T operator() (int r) const {
        T res = 0;
        for(int c=m.cols();c-->0;)
            res += m(r,c)*v(c);
        return res;
    }
    int size() const { return m.rows(); }
    bool is_locked() const { return m.is_locked() || v.is_locked(); };
};

template <int R,int C,typename T,typename T1,class E1,typename T2,class E2>
class Vector<T,C,Multiplied<Vector<T1,R,E1>, Matrix<T2,R,C,E2> > >
: public Common<Vector<T,C,Multiplied<Vector<T1,R,E2>, Matrix<T2,R,C,E2> > > > {
    typedef Vector<T1,R,E1> V;
    typedef Matrix<T2,R,C,E2> M;
    typedef Multiplied<V,M> E;
    typedef Vector<T,N,E> SAME;

    const typename V::Ref_t v;
    const typename M::Ref_t m;
  public:
    typedef SAME Ref_t;

    Vector<T,C,E>(const V &v_, const M &m_) : v(v_), m(m_) {}
    T operator() (int c) const {
        T res = 0;
        for(int r=m.rows();r-->0;)
            res += v(r)*m(r,c);
        return res;
    }
    int size() const { return m.cols(); }
    bool is_locked() const { return v.is_locked() || m.is_locked(); };
};

}; // namespace etLAP

#endif // _VECTOR_H_
