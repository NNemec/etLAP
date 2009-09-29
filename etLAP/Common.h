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

#ifndef _ETLAP_COMMON_H_
#define _ETLAP_COMMON_H_

#ifndef ASSIGN_POLICY
#define ASSIGN_POLICY AP_manual
#endif

#define AP_manual 1
#define AP_automatic 2
#define AP_conservative 3
#define AP_norefs 4

#ifndef etLAP_DEFAULT_ENGINE
#define etLAP_DEFAULT_ENGINE Smart
#endif

#include "Tags.h"

namespace etLAP {

template <class _SAME>
class Common;

template <class _SAME,typename T>
class Common_Smart;

template <class _SAME,typename T>
class Common_Packed;

template <typename T=double,int N=0,class EngineTag=etLAP_DEFAULT_ENGINE>
class Vector;

template <typename T=double,int R=0,int C=0,class EngineTag=etLAP_DEFAULT_ENGINE>
class Matrix;

} // namespace etLAP

#include "Util.h"
#include "Types.h"
#include "Traits.h"
#include "Storage.h"

#ifdef DEBUG
#include <typeinfo>
#endif

#include <cassert>


/*****************************************************************************
 *  Common classes
 */


namespace etLAP {

template <class SAME>
class Common {
  protected:
  public:
};

template <typename T> inline bool _is_locked(const T &) { return false; }
template <typename T,int N,class E> inline bool _is_locked(const Vector<T,N,E> &m) { return m.is_locked(); }
template <typename T,int R,int C,class E> inline bool _is_locked(const Matrix<T,R,C,E> &m) { return m.is_locked(); }

template <typename T> inline T buf(const T &t) { return t; }

template <class SAME,typename T>
class Common_Smart
: public Common<SAME> {
  protected:
#if ASSIGN_POLICY==AP_automatic || (ASSIGN_POLICY==AP_manual && !defined(NDEBUG))
    bool l;
    void lock() { l=true; }
    void unlock() { l=false; }
#else
    bool x;  // somehow this does make the program faster ? ? ? ! ! ! ! !
    void lock() { }
    void unlock() { }
#endif
  public:
    Common_Smart<SAME,T>() { unlock(); }

    template<typename X>
    void assign_from(const X &src) {
#if (ASSIGN_POLICY == AP_manual)
        lock();
        assert(!_is_locked(src));
        unlock();
        assign(*(SAME *)this,src,typename TypeCast<T,X>::Cast_t());
#elif (ASSIGN_POLICY == AP_automatic)
        lock();
        if(_is_locked(src)) {
            unlock();
            assign(*(SAME *)this,buf(src),typename TypeCast<T,X>::Cast_t());
        } else {
            unlock();
            assign(*(SAME *)this,src,typename TypeCast<T,X>::Cast_t());
        }
#elif ASSIGN_POLICY == AP_conservative
        assign(*(SAME *)this,buf(src),typename TypeCast<T,X>::Cast_t());
#elif ASSIGN_POLICY == AP_norefs
        assign(*(SAME *)this,src,typename TypeCast<T,X>::Cast_t());
#else
#error ASSIGN_POLICY not defined!!
#endif
    }

#if ASSIGN_POLICY==AP_automatic || (ASSIGN_POLICY==AP_manual && !defined(NDEBUG))
    bool is_locked() const { return l; }
#else
    bool is_locked() const { return false; }
#endif
};

template <class SAME,typename T>
class Common_Packed
: public Common<SAME> {
  public:
    template<typename X>
    void assign_from(const X &src) {
        assign(*(SAME *)this,src,typename TypeCast<T,X>::Cast_t());
    }

    bool is_locked() const { return false; }
};

} // namespace etLAP

#endif // _ETLAP_COMMON_H_
