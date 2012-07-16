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

/* ------------- SPIN SYSTEM MONTE CARLO SIMULATED ANNEALING ----------------
 *
 * arguments: [uint system_type]
 *            [uint N] [bool periodic (0,1)]
 *            [char init (r,u,d,...)]
 *            [double T_start] [double T_end]
 *            [ulint t_end] [char cooling_schedule]
 *            [uint t_boost]
 *            [bool run_plots (0,1)] [uint take_images (0,1)]
 *            [double J] [double g]
 *            [double B]
 *
 * system_type:
 * -- 1: one-dimensional Ising-Model with single-flip Metropolis-Algorithm
 * -- 2: two-dimensional Ising-Model with single-flip Metropolis-Algorithm
 * -- 3: two-dimensional field-free Ising-Model with Wolff-Algorithm
 * -- 4: 2d Ising-Model with long range dipole interaction (Metropolis-Alg.)
 * -- 5: like 4 but improved
 * -- 6: 2d Ising-Model with dip-dip-int on a honeycomb lattice
 *
 * init:
 * -- r: initialize spins randomly
 * -- u: all spins up at t=0
 * -- d: all spins down at t=0
 * -- e: 50:50 up/down (for B=0 and very low temperatures)
 * -- c: checkerboard (2d)
 * -- s: enery minimizing stripes
 * -- [uint]: stripes of manually defined width
 *
 * cooling_schedule:
 * -- l: linear
 * -- p: parabolic
 *
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
using namespace std;

#include "sarun.hpp"


// ----- MAIN: SIMULATED ANNEALING MANAGER -----

int main( int argc, char* argv[] )
{
  cout << "SPIN SYSTEM MONTE CARLO SIMULATED ANNEALING\n";
  cout << "===========================================\n\n";

  // ----- READ ARGUMENTS FROM COMMAND LINE -----
  cout << "reading from command line and preparing the simulated annealing ...\n";

  simann sa;

  sa.par.system_type = atoi( argv[1] );

  sa.par.N = atoi( argv[2] );
  sa.par.periodic = atoi( argv[3] );

  sa.par.init = *argv[4];

  sa.par.T_start = atof( argv[5] );
  sa.par.T_end = atof( argv[6] );

  sa.par.t_end = atoi( argv[7] );
  sa.par.cooling_schedule = *argv[8];

  sa.par.t_boost = atoi( argv[9] );

  sa.par.run_plot = atoi( argv[10] );

  sa.par.take_images = atoi( argv[11] );

  sa.par.J = atof( argv[12] );
  sa.par.g = atof( argv[13] );

  sa.par.B = atof( argv[14] );

  // generate the name of the output directory ...
  stringstream tmp;
  tmp << "sa" << sa.par.system_type << "_"
      << "N" << sa.par.N << ( sa.par.periodic ? "p" : "" ) << "_"
      << sa.par.init << "_"
      << setfill( '0' )
      << "T-" << setw( 5 ) << int( sa.par.T_start * 1000 ) << "-"
              << setw( 5 ) << int( sa.par.T_end * 1000 ) << "_"
      << sa.par.t_end << sa.par.cooling_schedule << "_"
      << "J" << setw( 5 ) << int( sa.par.J * 1000 ) << "_"
      << "g" << setw( 5 ) << int( sa.par.g * 1000 ) << "_"
      << "B" << setw( 5 ) << int( sa.par.B * 1000 );
  const string wdir = tmp.str();
  tmp.str( "" );

  // make an output directory
  if ( system( ( "test -e " + wdir + " && rm -r " + wdir
                + " ; mkdir " + wdir ).c_str() ) != 0 ) {
    cout << "ERROR while making the output directory " << wdir;
    return 1;
  }

  // ----- RUN THE SIMULATED ANNEALING -----
  cout << "running the simulated annealing ( " << wdir << " ) ...\n";

  srand( time( NULL ) );

  sa.res = sarun( sa.par, wdir );
  if ( !sa.res.success ) {
    cout << "WARNING: simulated annealing terminated unsuccessfully!\n";
  }

  cout << "... done!";
  return 0;
}
