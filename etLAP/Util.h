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

#ifndef _ETLAP_UTIL_H_
#define _ETLAP_UTIL_H_

namespace etLAP {

template <bool COND,typename IF_TRUE,typename IF_FALSE>
struct IF;
template <typename IF_TRUE,typename IF_FALSE> struct IF<true,IF_TRUE,IF_FALSE> { typedef IF_TRUE T; };
template <typename IF_TRUE,typename IF_FALSE> struct IF<false,IF_TRUE,IF_FALSE> { typedef IF_FALSE T; };

template<bool B> struct CTAssertClass {};
template<> struct CTAssertClass<true> { static void test() {} };
#define CTAssert(c) etLAP::CTAssertClass<(c)>::test()

}; // namespace etLAP

#endif // _ETLAP_UTIL_H_