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

#include "Util.h"

namespace etLAP {

template <typename T,int N>
class Tuple {
    T t[N];

  public:
    Tuple<T,N>() {}
    explicit Tuple<T,N>(const T &t_) { for(int i=N;i-->0;) t[i]=t_; }

    Tuple<T,N>(const Tuple<T,N-1> &tp_,const T &t_) { for(int i=N-1;i-->0;) t[i]=tp_[i]; t[N-1] = t_; }
    Tuple<T,N>(const T &t_,const Tuple<T,N-1> &tp_) { t[0] = t_; for(int i=N-1;i-->0;) t[i+1]=tp_[i]; }
    template <int N1>
    Tuple<T,N>(const Tuple<T,N1> &tp1_,const Tuple<T,N-N1> &tp2_) {
        for(int i=N-N1;i-->0;) t[i+N1]=tp2_[i];
        for(int i=N1;i-->0;)   t[i]   =tp1_[i];
    }

    Tuple<T,N-1> tail() { Tuple<T,N-1> res; for(int i=N-1;i-->0;) res[i]=t[i+1]; return res; }
    Tuple<T,N-1> head() { Tuple<T,N-1> res; for(int i=N-1;i-->0;) res[i]=t[i]; return res; }

    const T operator[](int i) const { return t[i]; }
    T &operator[](int i) { return t[i]; }

    Tuple<T,N> operator+(const Tuple<T,N> &other) const { Tuple<T,N> res; for(int i=N;i-->0;) res[i]=t[i] + other[i]; return res; }
    Tuple<T,N> operator-(const Tuple<T,N> &other) const { Tuple<T,N> res; for(int i=N;i-->0;) res[i]=t[i] - other[i]; return res; }
};

template <typename T>
class Tuple<T,0> {
  public:
    Tuple<T,0>() {}
    explicit Tuple<T,0>(const T &t_) {}

    Tuple<T,0>(const Tuple<T,0> &tp1_,const Tuple<T,0> &tp2_) {}

    const T operator[](int i) const { throw 0; }
    T &operator[](int i) { throw 0; }
};

template <typename T,int N1,int N2>
inline Tuple<T,N1+N2> operator,(const Tuple<T,N1> &t1,const Tuple<T,N2> &t2) {
    return Tuple<T,N1+N2>(t1,t2);
}

template <typename T,int N>
inline Tuple<T,1+N> operator,(const T &t1,const Tuple<T,N> &t2) {
    return Tuple<T,1+N>(t1,t2);
}

template <typename T,int N>
inline Tuple<T,N+1> operator,(const Tuple<T,N> &t1,const T &t2) {
    return Tuple<T,N+1>(t1,t2);
}

template <typename T>
inline Tuple<T,1> _tuple_(const T &t) {
    return Tuple<T,1>(t);
}

#if __GNUC__ == 2
#define tuple(t,args...) (::etLAP::_tuple_(t) , ##args)
#else
#define tuple(t,...) (::etLAP::_tuple_(t) , ##__VA_ARGS__)
#endif

} // namespace etLAP
