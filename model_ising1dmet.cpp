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

#include "model_ising1dmet.hpp"


// ----- 1D ISING MODEL -----

IsingModel1d::IsingModel1d( const unsigned int& N, const bool& periodic,
                            const double& J, const double& B, const double& T,
                            const unsigned int& fsize_correction_mode,
                            const string& cwd )
  : SystemModel( periodic, B, T, cwd ), J( J ),
    fsize_correction_mode( fsize_correction_mode ),
    fsize_ordered_phase( false )
{
  spin.resize( N );

  pal.resize( 2 );
  pal[0] = png::color( 0, 0, 0 );
  pal[1] = png::color( 255, 255, 255 );
}


unsigned int IsingModel1d::spin_count() const
{
  // returns the total number of spins in the system
  return spin.size();
}


png::image< png::index_pixel > IsingModel1d::get_image() const
{
  png::image< png::index_pixel > image( spin.size(), 8 );
  image.set_palette( pal );
  for ( size_t n = 0; n < image.get_width(); ++n ) {
    for ( unsigned int l = 0; l < 8; ++l ) {
      if ( spin[n].get() == 1 ) {
        image[l][n] = png::index_pixel( 1 );
      } else if ( spin[n].get() == -1 ) {
        image[l][n] = png::index_pixel( 0 );
      }
    }
  }
  return image;
}


bool IsingModel1d::prepare( const char& mode )
{
  switch ( mode ) {

  case 'r':	// completely random state ...
    for ( unsigned int n = 0; n < spin.size(); n++ ) {
      if ( gsl_rng_uniform_int( rng, 2 ) == 0 ) {
        spin[n].flip();
      }
    }
    break;

  case 'u': // sets all spins up
    for ( unsigned int n = 0; n < spin.size(); n++ ) {
      spin[n].set( +1 );
    }
    break;

  case 'd': // sets all spins down
    for ( unsigned int n = 0; n < spin.size(); n++ ) {
      spin[n].set( -1 );
    }
    break;

  case 'e': // sets first 50% up, other 50% down
    for ( unsigned int n = 0; n < spin.size() / 2; n++ ) {
      spin[n].set( +1 );
    }
    for ( unsigned int n = spin.size() / 2; n < spin.size(); n++ ) {
      spin[n].set( -1 );
    }
    break;

  default: // unknown mode?
    return false;
  }
  return true;
}


void IsingModel1d::metropolis_singleflip()
{

  // find a random spin to flip
  unsigned long int k = gsl_rng_uniform_int( rng, spin.size() );
  spin[k].flip();

  // calculate energy difference
  double deltaH = - 2 * B * spin[k].get();
  if ( k == 0 ) {
    if ( periodic ) {
      deltaH += -2 * J * ( spin[spin.size() - 1] * spin[0] );
    }
    deltaH += -2 * J * ( spin[0] * spin[1] );
  } else if ( k == spin.size() - 1 ) {
    deltaH += -2 * J * ( spin[k - 1] * spin[k] );
    if ( periodic ) {
      deltaH += -2 * J * ( spin[k] * spin[0] );
    }
  } else {
    deltaH += -2 * J * ( spin[k - 1] * spin[k] );
    deltaH += -2 * J * ( spin[k] * spin[k + 1] );
  }

  if ( deltaH > 0 ) {
    // read or calculate exp(- deltaH / T)
    double exp_deltaHoverT;
    if ( exp_precalc.count( deltaH ) == 1 ) {
      exp_deltaHoverT = exp_precalc[deltaH];
    } else {
      exp_deltaHoverT = exp( - deltaH / T );
      exp_precalc[deltaH] = exp_deltaHoverT;
    }
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp_deltaHoverT ) {
      // new state rejected ... reverting!
      spin[k].flip();
      return;
    }
  }
}


void IsingModel1d::mcstep()
{
  // empty precalculated exponentials (needed for SA)
  exp_precalc.clear();

  for ( unsigned int n = 1; n <= spin.size(); n++ ) {
    metropolis_singleflip();
  }
  time++;
}


void IsingModel1d::mcstep_dry( const unsigned int& k_max )
{
  for ( unsigned int k = 0; k < k_max; k++ ) {
    mcstep();
  }
  time -= k_max;

  if ( fsize_correction_mode == 2 ) {
    // try do determine if the system is in the ordered phase
    unsigned int msmall_count = 0, mlarge_count = 0;
    for ( unsigned int k = 0; k < k_max; k++ ) {
      mcstep();
      if ( abs( M() ) < spin.size() / 2 ) {
        msmall_count++;
      } else {
        mlarge_count++;
      }
    }
    if ( mlarge_count > msmall_count ) {
      fsize_ordered_phase = true;
      cout << "assuming ordered phase @ T = " << T << endl;
    } else {
      cout << "assuming disordered phase @ T = " << T << endl;
    }
    time -= k_max;
  }
}


double IsingModel1d::H() const
{
  // measures the system's energy
  double H = 0;

  for ( unsigned int n = 0; n < spin.size() - 1; n++ ) {
    H += - J * ( spin[n] * spin[n + 1] ) - B * spin[n].get();
  }
  H += - B * spin[spin.size() - 1].get();

  if ( periodic ) {
    H += - J * ( spin[spin.size() - 1] * spin[0] );
  }

  return H;
}


double IsingModel1d::h() const
{
  // measures the system's energy per spin
  return H() / spin.size();
}


unsigned long int IsingModel1d::t() const
{
  // measures the system's time in lattice sweeps (MC time units)
  return time;
}


int IsingModel1d::M() const
{
  // measures the system's magnetization
  int M = 0;
  for ( unsigned int n = 0; n < spin.size(); n++ ) {
    M += spin[n].get();
  }

  // finite size corrections to the magnetization
  if ( ( fsize_correction_mode == 1 ) ||
       ( ( fsize_correction_mode == 2 ) && fsize_ordered_phase ) ) {
    M = abs( M );
  }
  return M;
}


double IsingModel1d::m() const
{
  // measures the system's magnetization per spin
  return double( M() ) / spin.size();
}


vector<double> IsingModel1d::ss_corr() const
{
  // measure spin-spin correlations
  vector<double> result;
  vector<unsigned int> samples;
  result.resize( spin.size(), 0 );
  samples.resize( spin.size(), 0 );
  for ( unsigned int i = 0; i < spin.size(); i++ ) {
    for ( unsigned int j = 0; j < spin.size(); j++ ) {
      result[abs( int( i - j ) )] += spin[i] * spin[j];
      samples[abs( int( i - j ) )]++;
    }
  }
  for ( unsigned int d = 0; d < spin.size(); d++ ) {
    result[d] /= samples[d];
  }
  return result;
}
