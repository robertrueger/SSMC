/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of SSMC.
 *
 * SSMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _UTILS_VEC2_H_INCLUDED
#define _UTILS_VEC2_H_INCLUDED

#include <iostream>
#include <cmath>
using namespace std;


class Vec2
{
    // ----- HELPER CLASS: TWO-DIMENSIONAL VECTOR -----

  public:

    double x, y;

    // constructors
    Vec2();
    Vec2( double const& x_init, double const& y_init );
};

// multiplication with double
Vec2 operator*( double const& a, Vec2 const& myvec );

// multiplication with short
Vec2 operator*( short const& a, Vec2 const& myvec );

// vector addition
Vec2 operator+( Vec2 const& left, Vec2 const& right );
Vec2 operator-( Vec2 const& left, Vec2 const& right );

// dot product
double operator*( Vec2 const& left, Vec2 const& right );

// euclidian norm of a vector
double abs( Vec2 const& v );

// ostream output
ostream& operator<<( ostream& stream, Vec2 const& myvec );

#endif // _UTILS_VEC2_H_INCLUDED
