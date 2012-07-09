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

#include "model_ising2dsqrdipolemet.hpp"


// ------- 2D ISING MODEL WITH LONG RANGE DIPOLE INTERACTION ---------

IsingModel2dDipole::IsingModel2dDipole( const unsigned int& N,
                                        const bool& periodic, const double& J,
                                        const double& g, const double& B,
                                        const double& T, string cwd )
  : IsingModel2d( N, false, J, B, T, 0, cwd ), g( g ), H_memory( H() )
{
  if ( periodic ) {
    cout << "WARNING: Periodic boundary conditions are not supported!" << endl;
  }
}


bool IsingModel2dDipole::prepare( const char& mode )
{
  bool success = IsingModel2d::prepare( mode );
  H_memory = H();
  return success;
}


void IsingModel2dDipole::metropolis_singleflip()
{
  // (reimplementation because we can't calculate deltaH)

  // find a random spin to flip
  unsigned int flip_line = gsl_rng_uniform_int( rng, size );
  unsigned int flip_col  = gsl_rng_uniform_int( rng, size );

  // flip it and calculate energy difference!
  spin[flip_line][flip_col].flip();
  double deltaH = H() - H_memory;

  if ( deltaH > 0 ) {
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp( - deltaH / T ) ) {
      // new state rejected ... reverting!
      spin[flip_line][flip_col].flip();
      return;
    }
  }

  H_memory += deltaH;
}


void IsingModel2dDipole::metropolis_mirror()
{
  // vertical or horizontal mirroring?
  bool mode = gsl_rng_uniform_int( rng, 2 );
  unsigned int fmline = gsl_rng_uniform_int( rng, size - 1 ) + 1;
  unsigned int fmcol = gsl_rng_uniform_int( rng, size - 1 ) + 1;

  if ( mode ) {
    // mirroring: horizontal
    for ( unsigned int l = fmline; l < size; l++ ) {
      for ( unsigned int c = 0; c < size; c++ ) {
        spin[l][c].flip();
      }
    }
  } else {
    // mirroring: vertical
    for ( unsigned int c = fmcol; c < size; c++ ) {
      for ( unsigned int l = 0; l < size; l++ ) {
        spin[l][c].flip();
      }
    }
  }

  // energy difference
  double deltaH = H() - H_memory;

  if ( deltaH > 0 ) {
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp( - deltaH / T ) ) {
      // new state rejected ... reverting!
      if ( mode ) {
        // mirroring: horizontal
        for ( unsigned int l = fmline; l < size; l++ ) {
          for ( unsigned int c = 0; c < size; c++ ) {
            spin[l][c].flip();
          }
        }
      } else {
        // mirroring: vertical
        for ( unsigned int c = fmcol; c < size; c++ ) {
          for ( unsigned int l = 0; l < size; l++ ) {
            spin[l][c].flip();
          }
        }
      }
      return;
    }
  }

  H_memory += deltaH;
}


void IsingModel2dDipole::metropolis_blockflip()
{
  // generate coordinates of the block
  unsigned int temp;
  unsigned int l1 = gsl_rng_uniform_int( rng, size );
  unsigned int l2 = gsl_rng_uniform_int( rng, size );
  if ( l1 > l2 ) {
    temp = l1;
    l1 = l2;
    l2 = temp;
  }
  unsigned int c1 = gsl_rng_uniform_int( rng, size );
  unsigned int c2 = gsl_rng_uniform_int( rng, size );
  if ( c1 > c2 ) {
    temp = c1;
    c1 = c2;
    c2 = temp;
  }

  // flip the block
  for ( unsigned int l = l1; l <= l2; l++ ) {
    for ( unsigned int c = c1; c <= c2; c++ ) {
      spin[l][c].flip();
    }
  }

  // energy difference
  double deltaH = H() - H_memory;

  if ( deltaH > 0 ) {
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp( - deltaH / T ) ) {
      // new state rejected ... reverting!
      for ( unsigned int l = l1; l <= l2; l++ ) {
        for ( unsigned int c = c1; c <= c2; c++ ) {
          spin[l][c].flip();
        }
      }
      return;
    }
  }

  H_memory += deltaH;
}


void IsingModel2dDipole::mcstep()
{
  for ( unsigned int n = 1; n <= N; n++ ) {
    metropolis_singleflip();
    metropolis_blockflip();
    metropolis_singleflip();
    metropolis_mirror();
  }
  time++;
}


double IsingModel2dDipole::H() const
{
  double H = IsingModel2d::H();	// energy without dipole interaction

  for ( unsigned int i = 0; i < N - 1; i++ ) {
    for ( unsigned int j = i + 1; j < N; j++ ) {
      unsigned int iline = i / size, icol = i % size;
      unsigned int jline = j / size, jcol = j % size;
      double r = sqrt( ( jline - iline ) * ( jline - iline )
                       + ( jcol - icol ) * ( jcol - icol ) );
      H += - g / ( r * r * r ) * ( spin[iline][icol] * spin[jline][jcol] );
    }
  }
  return H;
}
