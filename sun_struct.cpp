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

#include <etLAP/SU_N.h>
#include <etLAP/Operator.h>
#include <etLAP/Output.h>

using namespace etLAP;

enum { DIM = 2 };

int main() {
/*
    for(int a=1;a<DIM*DIM-2;a++)
    for(int b=a+1;b<DIM*DIM-1;b++)
    for(int c=b+1;c<DIM*DIM;c++)
*/

    for(int a=1;a<DIM*DIM;a++)
    for(int b=1;b<DIM*DIM;b++)
    for(int c=1;c<DIM*DIM;c++)

        std::cout << "f(" << a << "," << b << "," << c << ") = "
                  << SUNstructure<DIM,double>(a,b,c) << "\n";
};
