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

#ifndef _ETLAP_TUPEL_H_
#define _ETLAP_TUPEL_H_

#include "Util.h"

namespace etLAP {

template <typename T,int N>
class Tupel {
    T t[N];

  public:
    Tupel<T,N>() {};
    explicit Tupel<T,N>(const T &t_) { for(int i=N;i-->0;) t[i]=t_; }

    Tupel<T,N>(const Tupel<T,N-1> &tp_,const T &t_) { for(int i=N-1;i-->0;) t[i]=tp_[i]; t[N-1] = t_; }
    Tupel<T,N>(const T &t_,const Tupel<T,N-1> &tp_) { t[0] = t_; for(int i=N;i-->1;) t[i]=tp_[i-1]; }
    template <int N1>
    Tupel<T,N>(const Tupel<T,N1> &tp1_,const Tupel<T,N-N1> &tp2_) { int i=N; while(i-- > N1) t[i]=tp2_[i-N1]; while(i-->0) t[i]=tp1_[i]; }

    Tupel<T,N-1> tail() { Tupel<T,N-1> res; for(int i=N-1;i-->0;) res[i]=t[i+1]; return res; }
    Tupel<T,N-1> head() { Tupel<T,N-1> res; for(int i=N-1;i-->0;) res[i]=t[i]; return res; }

    const T operator[](int i) const { return t[i]; }
    T &operator[](int i) { return t[i]; }
};

template <typename T,int N1,int N2>
inline Tupel<T,N1+N2> operator,(const Tupel<T,N1> &t1,const Tupel<T,N2> &t2) {
    return Tupel<T,N1+N2>(t1,t2);
};

template <typename T,int N>
inline Tupel<T,1+N> operator,(const T &t1,const Tupel<T,N> &t2) {
    return Tupel<T,1+N>(t1,t2);
};

template <typename T,int N>
inline Tupel<T,N+1> operator,(const Tupel<T,N> &t1,const T &t2) {
    return Tupel<T,N+1>(t1,t2);
};

template <typename T>
inline Tupel<T,1> _tupel_(const T &t) {
    return Tupel<T,1>(t);
};

#if __GNUC__ == 2
#define tupel(t,args...) (::etLAP::_tupel_(t) , ##args)
#else
#define tupel(t,...) (::etLAP::_tupel_(t) , ##__VA_ARGS__)
#endif

}; // namespace etLAP

#endif // _ETLAP_TUPEL_H_
