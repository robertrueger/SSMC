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

#include "sa_coolingschedules.hpp"


// ----- SIMULATED ANNEALING COOLING SCHEDULES -----

double linear_cooling( const double& T_start, const double& T_end,
                       const unsigned long int& t_end, 
                       const unsigned long int& t )
{
  return T_start - ( T_start - T_end ) * ( double( t ) / double( t_end ) );
}


double parabolic_cooling( const double& T_start, const double& T_end,
                          const unsigned long int& t_end,
                          const unsigned long int& t )
{
  return ( T_start - T_end ) * pow( ( 1 - double( t ) / double( t_end ) ), 2 )
                                                                        + T_end;
}
