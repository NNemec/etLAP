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

#ifndef _ETLAP_TAGS_H_
#define _ETLAP_TAGS_H_

namespace etLAP {

/*******************************************************************************
 * Engine tags
 */
class Smart;
class Packed;
class ExtPtr;
class Diagonal;
class Scalar;
class Zero;
class One;

template <class A,class Op>
class ElemUnOp;

template <class A,class B,class Op>
class ElemBinOp;

template <class A,typename T,class Op>
class ScalarOp;

template <class A>
class Transposed;

template <class A,class B>
class Multiplied;

template <class T>
class Row;

template <class T>
class Column;

class Buffer;

template <class A>
class NoBuffer;



/*******************************************************************************
 * Operator tags
 */

class OpIdent;
class OpNeg;
class OpConj;
class OpAbs;
class OpReal;
class OpImag;
class OpExp;
class OpTranspose;
class OpAdj;
class OpTrace;
class OpDet;
class OpInv;
class OpMax;

class OpAdd;
class OpSub;
class OpMul;
class OpDiv;

class OpRevMul;

}; // namespace etLAP

#endif // _ETLAP_TAGS_H_
