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

#ifndef _ETLAP_STORAGE_H_
#define _ETLAP_STORAGE_H_

#include <memory>

namespace etLAP {

template <typename T /*,class Allocator = std::allocator<T>*/ >
class Storage {
    struct Rep {
        unsigned int refcount;
        unsigned int size;

        T* data() { return reinterpret_cast<T *>(this + 1); };
        T& operator[] (unsigned int s) { return data() [s]; };
        Rep* grab() { ++refcount; return this; };
        void release() { if (--refcount == 0) delete this; };

        static void *operator new(size_t s, unsigned int extra) {
//            return Allocator::allocate(s + extra * sizeof (T));
            return ::operator new(s + extra * sizeof (T));
        };

        static void operator delete(void *ptr) {
            ::operator delete(ptr);
//            Allocator::deallocate(ptr, sizeof(Rep) +
//                                  reinterpret_cast<Rep *>(ptr)->size * sizeof (T));
        };

        static Rep *create (unsigned int size_) {
            Rep *p = new (size_) Rep;
            p->refcount = 1;
            p->size = size_;
            return p;
        };

        Rep *clone() {
            Rep *p = create(size);
            p->copy_from(data());
            return p;
        };

        void copy_from(const T *src) {
            for(unsigned n=0;n<size;n++) data()[n] = src[n];
        };

        void clear() {
            for(unsigned int n=0;n<size;n++)
                data()[n] = (T)0;
        };

    private:
        Rep &operator= (const Rep &);
    };

    Rep *rep;
    mutable bool souvereign;

  public:
    Storage<T>() { rep = 0; };
    ~Storage<T>() { if(rep) rep->release(); }

    Storage<T>(const Storage<T> &other) {
        if(other.rep) {
            rep = other.rep->grab();
            other.souvereign = false;
            souvereign = false;
        } else {
            rep = 0;
        };
    };

    Storage<T>(unsigned int size) {
        rep = Rep::create(size);
        souvereign = true;
    };

    Storage<T>(T *ptr,unsigned int size) {
        rep = Rep::create(size);
        rep->copy_from(ptr);
        souvereign = true;
    };

    Storage<T> &operator=(const Storage<T> &other) {
        if(rep)
            rep->release();
        if(other.rep) {
            rep = other.rep->grab();
            other.souvereign = false;
            souvereign = false;
        } else {
            rep = 0;
        };
        return *this;
    };

    T operator[](unsigned int i) const { assert(rep); return (*rep)[i]; }
    T &operator[](unsigned int i) { assert(rep); assert(souvereign); return (*rep)[i]; }

    void resize(unsigned int size) {
        if(rep && size==rep->size) return;
        if(rep) rep->release();
        rep = Rep::create(size);
        souvereign = true;
    };

    void clear() {
        prepare_write_noclone();
        rep->clear();
    };

    void prepare_write_clone() {
        assert(rep);
        if(souvereign)
            return;
        if(rep->refcount > 1) {
            Rep *r = rep->clone();
            rep->release();
            rep = r;
        };
        souvereign = true;
    };

    void prepare_write_noclone() {
        assert(rep);
        if(souvereign)
            return;
        if(rep->refcount > 1) {
            unsigned int size = rep->size;
            rep->release();
            rep = Rep::create(size);
        };
        souvereign = true;
    };
    
    T *rawdata() { return rep->data(); };
    // be careful with this one -- only use it if you know exactly what you are doing!
};

}; // namespace etLAP

#endif // _ETLAP_STORAGE_H_
