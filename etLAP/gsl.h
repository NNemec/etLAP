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

#ifndef _GSL_H_
#define _GSL_H_

#include "Matrix.h"
#include "Vector.h"
#include "Operator.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>


namespace etLAP {

template <int N>
struct _gsl_eigen_symmv_workspace {
    gsl_eigen_symmv_workspace *gsl_workspace;

    gsl_matrix gsl_res;

    Vector<double,N,Packed> D;
    gsl_vector gsl_D;

    Matrix<double,N,N,Packed> U;
    gsl_matrix gsl_U;

    _gsl_eigen_symmv_workspace() {
        gsl_workspace = gsl_eigen_symmv_alloc(N);

        gsl_res.size1 = N;
        gsl_res.size2 = N;
        gsl_res.tda = N;
        gsl_res.block = 0;
        gsl_res.owner = 0;

        gsl_D.size = N;
        gsl_D.stride = 1;
        gsl_D.data = D.rawdata();
        gsl_D.block = 0;
        gsl_D.owner = 0;

        gsl_U.size1 = N;
        gsl_U.size2 = N;
        gsl_U.tda = N;
        gsl_U.data = U.rawdata();
        gsl_U.block = 0;
        gsl_U.owner = 0;
    };

    ~_gsl_eigen_symmv_workspace() {
        gsl_eigen_symmv_free(gsl_workspace);
    };

    void resize(int _size) {
	assert(_size == N);
    };
};

template <>
struct _gsl_eigen_symmv_workspace<0> {
    int size;

    gsl_eigen_symmv_workspace *gsl_workspace;

    gsl_matrix gsl_res;

    Vector<double,0> D;
    gsl_vector gsl_D;

    Matrix<double,0,0> U;
    gsl_matrix gsl_U;

    _gsl_eigen_symmv_workspace() {
        size = 0;
    };

    ~_gsl_eigen_symmv_workspace() {
        if(size)
            gsl_eigen_symmv_free(gsl_workspace);
    };

    void resize(int _size) {
        if(size == _size)
            return;
    
        if(size > 0)
            gsl_eigen_symmv_free(gsl_workspace);

        size = _size;
        D.resize(size);
        U.resize(size,size);

        gsl_workspace = gsl_eigen_symmv_alloc(size);

        gsl_res.size1 = size;
        gsl_res.size2 = size;
        gsl_res.tda = size;
        gsl_res.block = 0;
        gsl_res.owner = 0;

        gsl_D.size = size;
        gsl_D.stride = 1;
        gsl_D.data = D.rawdata();
        gsl_D.block = 0;
        gsl_D.owner = 0;

        gsl_U.size1 = size;
        gsl_U.size2 = size;
        gsl_U.tda = size;
        gsl_U.data = U.rawdata();
        gsl_U.block = 0;
        gsl_U.owner = 0;
    };
};

// inv(Matrix)
template <int N,class E>
inline const Matrix<double,N,N> gsl_inv(Matrix<double,N,N,E> a) {
    assert(a.cols() == a.rows());
    assert(N > 0);

    static _gsl_eigen_symmv_workspace<N> ws;
    
    ws.resize(a.cols());

    Matrix<double,N,N> res = a;
    ws.gsl_res.data = res.rawdata();

    gsl_eigen_symmv(&ws.gsl_res,&ws.gsl_D,&ws.gsl_U,ws.gsl_workspace);
    // res = U * diag(D) * transpose(U)

    res.clear();

/*
    for(uint x=N;x-->0;)
    for(uint r=N;r-->0;)
    for(uint c=N;c-->0;)
        res(r,c) += ws.U(r,x) * 1/ws.D(x) * transpose(ws.U)(x,c);
*/

    for(uint x=N;x-->0;) {
        assert(ws.D(x) > 0);
        double inv = 1/ws.D(x);
        for(uint r=N;r-->0;) {
            double p1 = ws.U(r,x) * inv;
            for(uint c=N;c-->0;)
                res(r,c) += p1 * transpose(ws.U)(x,c);
        };
    };

    return res;
};

struct ENonPositiveEigenvalue {};

// inv_sqrt(Matrix)
template <int N,class E>
inline const Matrix<double,N,N> gsl_inv_sqrt(Matrix<double,N,N,E> a) {
    assert(a.cols() == a.rows());
    assert(N > 0);

    static _gsl_eigen_symmv_workspace<N> ws;
    
    ws.resize(a.cols());

    Matrix<double,N,N> res = a;
    ws.gsl_res.data = res.rawdata();

    gsl_eigen_symmv(&ws.gsl_res,&ws.gsl_D,&ws.gsl_U,ws.gsl_workspace);
    // res = U * diag(D) * transpose(U)

    res.clear();

/*
    for(uint x=N;x-->0;)
    for(uint r=N;r-->0;)
    for(uint c=N;c-->0;)
        res(r,c) += ws.U(r,x) * 1/sqrt(ws.D(x)) * transpose(ws.U)(x,c);
*/

    for(uint x=N;x-->0;) {
        double D = ws.D(x);
        if(D<=0)
            throw ENonPositiveEigenvalue();
        double sqrt_inv = 1/sqrt(D);
        for(uint r=N;r-->0;) {
            double p1 = ws.U(r,x) * sqrt_inv;
            for(uint c=N;c-->0;)
                res(r,c) += p1 * transpose(ws.U)(x,c);
        };
    };

    return res;
};

}; // namespace etLAP

#endif // _GSL_H_
