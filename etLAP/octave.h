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

#include "Matrix.h"
#include "Vector.h"
#include "Operator.h"

#define GCC_ATTR_NORETURN

extern "C" {
void kpse_clear_dir_cache() {};
};

#include <octave/octave/CMatrix.h>

namespace etLAP {

class Octave;

/*****************************************************************************
 *  Octave - fixed size
 */

template <> struct Traits<Octave> { static const int costs = 0; static const bool flataccess = false; };


template <int R,int C>
class Matrix<std__complex<double>,R,C,Octave>
: public Common_Smart<Matrix<std__complex<double>,R,C,Octave>,std__complex<double> > {
    typedef std__complex<double> T;
    typedef Octave E;
    typedef Matrix<T,R,C,E> SAME;

    ::ComplexMatrix mat;

  public:
    typedef SAME Ref_t;

    Matrix<T,R,C,E>(): mat(R,C) { };
    Matrix<T,R,C,E>(::ComplexMatrix mat_): mat(mat_) { assert(mat_.rows() == R && mat_.cols() == C); };
    Matrix<T,R,C,E>(int rw_,int cl_): mat(R,C) { assert(rw_==R && cl_ == C); }; // for usage in mdp

    Matrix<T,R,C,E>(const SAME &restrict other): mat(other.mat) {};
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { assert(R == rw_ && C == cl_); };
    void clear() { mat.fill(0); };

    const T operator() (int r,int c) const { return mat.elem(r,c); };
    T &operator() (int r,int c,bool unsafe=false) { return mat.elem(r,c); };

//    const T operator[] (int n) const { return mat.elem(n); };
//    T &operator[] (int n) { return mat.elem(n); };
    const T operator[] (int n) const { const Array<T> *m = &mat; return m->elem(n); };
    T &operator[] (int n) { Array<T> *m = &mat; return m->elem(n); };
//    T &operator[] (int n) { return mat.data()[n]; };
//    T operator[] (int n) const { return mat.data()[n]; };

    void prepare_write_clone() {}
    void prepare_write_noclone() {}

    int rows() const { return R; };
    int cols() const { return C; };

    ::ComplexMatrix &octave_matrix() { return mat; }
};


/*****************************************************************************
 *  ExtPtr - Variable size
 */

template <>
class Matrix<std__complex<double>,0,0,Octave>
: public Common_Smart<Matrix<std__complex<double>,0,0,Octave>,std__complex<double> > {
    typedef std__complex<double> T;
    typedef Octave E;
    typedef Matrix<T,0,0,E> SAME;

    ::ComplexMatrix mat;

  public:
    typedef SAME Ref_t;

    Matrix<T,0,0,E>(int rw_,int cl_): mat(rw_,cl_) {};
    Matrix<T,0,0,E>(const SAME &restrict other): mat(other.mat) {};
    template<typename X>
    const SAME &operator=(const X &restrict src) { assign_from(src); return *(SAME *)this; };

    void resize(int rw_,int cl_) { mat.resize(rw_,cl_); };
    void clear() { mat.fill(0); };

    const T operator() (int r,int c) const { return mat.elem(r,c); };
    T &operator() (int r,int c,bool unsafe=false) { return mat.elem(r,c); };
/*
    const T operator[] (int n) const { return mat.elem(n); };
    T &operator[] (int n) { return mat.elem(n); };
*/
    void prepare_write_clone() {}
    void prepare_write_noclone() {}

    int rows() const { return mat.rows(); }
    int cols() const { return mat.cols(); }
    
    ::ComplexMatrix &octave_matrix() { return mat; }
};

// exp(Matrix)
template <int N,class E>
inline const Matrix<std__complex<double>,N,N,Octave> octave_exp(Matrix<std__complex<double>,N,N,E> a) {
    assert(a.cols() == a.rows());

    Matrix<std__complex<double>,N,N,Octave> tmp(a.rows(),a.cols());
    tmp = a;
    return Matrix<std__complex<double>,N,N,Octave>(tmp.octave_matrix().expm());
};



}; // namespace etLAP
