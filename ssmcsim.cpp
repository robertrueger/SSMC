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

/* ---------------- SPIN SYSTEM MONTE CARLO SIMULATION ----------------------
 *
 * arguments:	[uint system_type]
 * 				[uint N] [bool periodic (0,1)]
 *				[char init (r,u,d,...)] [uint drysweeps]
 * 				[uint bins] [uint binwidth] [uint intersweeps]
 *				[bool run_plots (0,1)]
 * 				[uint take_images]
 *				[bool calc_autocorr (0,1)] [bool calc_sscorr (0,1)]
 * 				[uint smode_perbin] [uint smode_permcs]
 * 				[uint finite_size_correction (0,1,2)]
 *				[double J] [double g]
 *				[double B]
 *				[double T]
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
 * finite_size_correction:
 * -- 0: disabled
 * -- 1: use absolute magnetization
 * -- 2: determine automatically
 *
*/

#include <cstdlib>
#include <iostream>
#include <sstream>
using namespace std;

#include "simrun.hpp"


// ----- MAIN: SIMULATION MANAGER -----

int main( int argc, char* argv[] )
{
  cout << "SPIN SYSTEM MONTE CARLO SIMULATION\n";
  cout << "==================================\n\n";

  // ----- READ ARGUMENTS FROM COMMAND LINE -----
  cout << "reading from command line and preparing the simulations ...\n";

  simulation sim;

  sim.par.system_type = atoi( argv[1] );

  sim.par.N = atoi( argv[2] );
  sim.par.periodic = atoi( argv[3] );

  sim.par.init = *argv[4];
  sim.par.drysweeps = atoi( argv[5] );

  sim.par.bins = atoi( argv[6] );
  sim.par.binwidth = atoi( argv[7] );
  sim.par.intersweeps = atoi( argv[8] );

  sim.par.run_plot = atoi( argv[9] );

  sim.par.take_images = atoi( argv[10] );

  sim.par.calc_autocorr = atoi( argv[11] );
  sim.par.calc_sscorr = atoi( argv[12] );

  sim.par.smode_perbin = atoi( argv[13] );
  sim.par.smode_permcs = atoi( argv[14] );

  sim.par.use_fsize_correction = atoi( argv[15] );

  sim.par.J = atof( argv[16] );
  sim.par.g = atof( argv[17] );

  sim.par.B = atof( argv[18] );
  sim.par.T = atof( argv[19] );

  // generate the name of the output directory ...
  stringstream tmp;
  tmp << "sim" << sim.par.system_type << "_"
      << "N" << sim.par.N << ( sim.par.periodic ? "p" : "" ) << "_"
      << sim.par.init << sim.par.drysweeps << "_"
      << sim.par.bins << "-" << sim.par.binwidth
      << "-" << sim.par.intersweeps << "_"
      << "fsc" << sim.par.use_fsize_correction << "_"
      << "smode-" << sim.par.smode_perbin << "-" << sim.par.smode_permcs << "_"
      << setfill( '0' )
      << "J" << setw( 5 ) << int( sim.par.J * 1000 ) << "_"
      << "g" << setw( 5 ) << int( sim.par.g * 1000 ) << "_"
      << "B" << setw( 5 ) << int( sim.par.B * 1000 ) << "_"
      << "T" << setw( 5 ) << int( sim.par.T * 1000 );
  const string wdir = tmp.str();
  tmp.str( "" );

  // make an output directory
  if ( system( ( "test -e " + wdir + " && rm -r " + wdir
                 + " ; mkdir " + wdir ).c_str() ) != 0 ) {
    cout << "ERROR while making the output directory " << wdir;
    return 1;
  }

  // ----- RUN THE SIMULATIONS -----
  cout << "running the simulation ( " << wdir << " ) ...\n";

  srand( time( NULL ) );

  sim.res = simrun( sim.par, wdir );
  if ( !sim.res.success ) {
    cout << "WARNING: simulation terminated unsuccessfully!\n";
  }

  cout << "... done!";
  return 0;
}
