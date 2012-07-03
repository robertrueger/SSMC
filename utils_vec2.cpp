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

#include "utils_vec2.hpp"


// ----- HELPER CLASS: TWO-DIMENSIONAL VECTOR -----

Vec2::Vec2() : x( 0.0 ), y( 0.0 ) { }


Vec2::Vec2( double const& x_init, double const& y_init ) : x( x_init ), y( y_init ) { }


Vec2 operator*( double const& a, Vec2 const& myvec )
{
    return Vec2( a * myvec.x, a * myvec.y );
}


Vec2 operator*( short const& a, Vec2 const& myvec )
{
    return static_cast<double>( a ) * myvec;
}


Vec2 operator+( Vec2 const& left, Vec2 const& right )
{
    return Vec2( left.x + right.x, left.y + right.y );
}


Vec2 operator-( Vec2 const& left, Vec2 const& right )
{
    return left + ( ( -1.0 ) * right );
}


double operator*( Vec2 const& left, Vec2 const& right )
{
    return left.x * right.x + left.y * right.y;
}


double abs( Vec2 const& v )
{
    return sqrt( v * v );
}


ostream& operator<<( ostream& stream, Vec2 const& myvec )
{
    stream << myvec.x << ' ' << myvec.y << endl;
    return stream;
}
