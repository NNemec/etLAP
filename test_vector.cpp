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

#include <etLAP/Vector.h>
#include <etLAP/Operator.h>
#include <etLAP/Output.h>
#include <iostream>

using namespace etLAP;

int main() {
    Vector<double,3> a;
    for(int i=0;i<3;i++)
        a(i) = i;
    std::cout << "[0,1,2]: " << a << "\n";
    Vector<double,3> b = 2*a;
    b += 3*a;
    std::cout << "[0,5,10]:" << b << "\n";
    
    Vector<double> c;
    c = a-b;
    std::cout << "[0,5,10]:" << c << "\n";
    
    return 0;
}
