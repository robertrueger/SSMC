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

#include "isingspin.hpp"


// ----- ISING SPIN CLASS -----

IsingSpin::IsingSpin() : S( 1 ) { }


int IsingSpin::get() const
{
  return S;
}


void IsingSpin::set( const int& newS )
{
  if ( ( newS == +1 ) || ( newS == -1 ) ) {
    S = newS;
  }
}


void IsingSpin::flip()
{
  S = -1 * S;
}


int IsingSpin::operator*( const IsingSpin& S ) const
{
  return get() * S.get();
}
