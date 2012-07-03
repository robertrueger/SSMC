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

#include "systemmodel.hpp"


// ----- SYSTEM MODEL TEMPLATE (DEFINES INTERFACE) -----

SystemModel::SystemModel( const bool& periodic, const double& B, const double& T, const string& cwd )
    : periodic( periodic ), B( B ), T( T ), time( 0 ), cwd( cwd )
{
    // initialize the model's random number generator (Mersenne Twister)
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
    gsl_rng_set( rng, rand() );
}


SystemModel::~SystemModel()
{
    // free rng's memory
    gsl_rng_free( rng );
}


void SystemModel::set_T( const double& newT )
{
    if ( newT >= 0 ) {
        T = newT;
    }
}


void SystemModel::set_B( const double& newB ) {
    B = newB;
}
