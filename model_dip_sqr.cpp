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

#include "model_dip_sqr.hpp"


// ------ 2D ISING MODEL ON A SQUARE LATTICE WITH DIP-DIP-INT -----

Ising2dDipSqr::Ising2dDipSqr( const unsigned int& size, const bool& periodic,
                              const double& J, const double& g, const double& B,
                              const double& T,
                              const unsigned int& fsize_correction_mode,
                              const string& cwd )
  : SystemModel( true, B, T, cwd ), J( J ), g( g ), N( size* size ),
    size( size ), fsize_correction_mode( fsize_correction_mode ),
    fsize_ordered_phase( false ), k_range( 2 * M_PI ), k_points( 201 ),
    sf_measurements( 0 )
{
  // wrong use of the model?
  if ( !periodic ) {
    cout <<	"Model: WARNING - Non-periodic boundary conditions"
         << " are not supported!" << endl;
  }
  if ( size % 2 != 0 ) {
    cout <<	"Model: WARNING - Size must be a multiple of 2!" << endl;
  }

  // set up two dimensional NxN array: spin[line][column]
  spin.resize( size );
  for ( unsigned int line = 0; line < size; line++ ) {
    spin[line].resize( size );
  }

  // set up an array for the precalculated interaction coeffs: effint[dl][dc]
  effint.resize( ( size / 2 ) + 1 );
  //cout << "gate1" << endl;
  for ( unsigned int dl = 0; dl < ( size / 2 ) + 1; ++dl ) {
    effint[dl].resize( ( size / 2 ) + 1 );
    for ( unsigned int dc = 0; dc < ( size / 2 ) + 1; ++dc ) {
      effint[dl][dc] = 0;
    }
  }

  // calculate the effective interaction strengths
  ofstream effint_log;
  effint_log.open( ( cwd + "effint.log" ).c_str() );
  if ( !effint_log.is_open() ) {
    cout << "ERROR while opening effective interaction log file in "
         << cwd << endl;
    exit( 1 );
  }
  ostream_setup( effint_log );

  effint[1][0] = J;// nearest neighbour interaction
  for ( unsigned int dl = 0; dl < ( size / 2 ) + 1; ++dl ) {
    for ( unsigned int dc = 0; dc < dl + 1; ++dc ) {
      effint[dl][dc] += calc_effint_dipdip( dl, dc, 1000 );
      effint[dc][dl] = effint[dl][dc];
      effint_log << "effint[" << dl << "][" << dc << "] = "
                 << effint[dl][dc] << endl;
    }
  }

  effint_log.close();

  // set up an array for the structure factor
  sf_sum.resize( k_points );
  for ( unsigned int line = 0; line < k_points; line++ ) {
    sf_sum[line].resize( k_points );
  }

  // set up a black and white palette
  pal.push_back( png::color( 0, 0, 0 ) );
  pal.push_back( png::color( 255, 255, 255 ) );
}


Ising2dDipSqr::~Ising2dDipSqr()
{
  // custom destructor, because we need to write
  // the structure factor to the logfile first

  // write structure factor to file
  ofstream sf_log;
  sf_log.open( ( cwd + "sf.log" ).c_str() );
  if ( !sf_log.is_open() ) {
    cout << "ERROR while opening structure factor log file in " << cwd << endl;
    exit( 1 );
  }
  ostream_setup( sf_log );

  for ( unsigned int line = 0; line < k_points; line++ ) {
    for ( unsigned col = 0; col < k_points; col++ ) {
      double kx = line * k_range / ( k_points - 1 ) - k_range / 2.0;
      double ky = col * k_range / ( k_points - 1 ) - k_range / 2.0;
      sf_log << kx << ' ' << ky << ' '
             << sf_sum[line][col] / sf_measurements / N << endl;
    }
  }

  sf_log.close();

  // write plotting script
  ofstream sf_plot;
  sf_plot.open( ( cwd + "sf_plot.pyx" ).c_str() );
  if ( !sf_plot.is_open() ) {
    cout << "ERROR while creating structure factor pyxplot file in "
         << cwd << endl;
    exit( 1 );
  }
  sf_plot << "\
    set terminal pdf \n\
    set output 'structure_factor.pdf' \n\
    set title 'Structure factor $S(\\vec k)$' \n\
    set size 15 square \n\
    set samples grid 201x201 \n\
    set colourmap rgb(1-c1):(1-c1):(1-c1) \n\
    set tics out \n\
    set grid x y \n\
    set xlabel \"$k_x$\" \n\
    set xtics (\"$\\pi$\" pi, \"$\\frac{3\\pi}{4}$\" 3*pi/4, \
               \"$\\frac{\\pi}{2}$\" pi/2, \"$\\frac{\\pi}{4}$\" pi/4, \
               \"0\" 0, \\  \n\
               \"$-\\pi$\" -pi, \"$-\\frac{3\\pi}{4}$\" -3*pi/4, \
               \"$-\\frac{\\pi}{2}$\" -pi/2, \"$-\\frac{\\pi}{4}$\" -pi/4)  \n\
    set mxtics pi/8 \n\
    set ylabel \"$k_y$\"  \n\
    set ytics (\"$\\pi$\" pi, \"$\\frac{3\\pi}{4}$\" 3*pi/4, \
               \"$\\frac{\\pi}{2}$\" pi/2, \"$\\frac{\\pi}{4}$\" pi/4, \
               \"0\" 0, \\  \n\
               \"$-\\pi$\" -pi, \"$-\\frac{3\\pi}{4}$\" -3*pi/4, \
               \"$-\\frac{\\pi}{2}$\" -pi/2, \"$-\\frac{\\pi}{4}$\" -pi/4) \n\
    set mytics pi/8 \n\
    plot [-pi:pi][-pi:pi] 'sf.log' with colourmap notitle";
  sf_plot.close();
}


double Ising2dDipSqr::calc_effint_dipdip( const int& dl_seed, const int& dc_seed,
    const int& system_clones )
{
  // calculate dipole-dipole-interaction with a seed spin that is dl lines and dc columns
  // far away + all of its copies on a system_clones large system with periodic boundaries
  // system_clones = 1 --> 3x3; 2 --> 5x5; 3 --> 7x7; ...

  double result = 0;
  //uint count = 0;
  for ( long int dl = dl_seed - system_clones * size;
        dl < ( system_clones + 1 )*size; dl += size ) {
    for ( long int dc = dc_seed - system_clones * size;
          dc < ( system_clones + 1 )*size; dc += size ) {
      double r = sqrt( dl * dl + dc * dc );
      if ( r == 0 ) {
        continue;    // no interaction with itself ...
      }
      result += g / ( r * r * r );
      //cout << dl << ' ' << dc << ' ' << ' ' << dl*dl+dc*dc << ' '
      //     << r << ' ' << result << endl;
      //count++;
      //if (count == 100) exit(0);
      /*if ((dl_seed == 1) && (dc_seed == 0)) {
          cout << g / (r*r*r) << endl;
      }*/
    }
  }
  return result;
}


png::image< png::index_pixel > Ising2dDipSqr::get_image() const
{
  png::image< png::index_pixel > image( size, size );
  image.set_palette( pal );
  for ( size_t line = 0; line < image.get_height(); ++line ) {
    for ( size_t col = 0; col < image.get_width(); ++col ) {
      if ( spin[line][col].get() == 1 ) {
        image[line][col] = png::index_pixel( 1 );
      } else if ( spin[line][col].get() == -1 ) {
        image[line][col] = png::index_pixel( 0 );
      }
    }
  }
  return image;
}


bool Ising2dDipSqr::prepare_striped( const unsigned int& stripe_width )
{
  if ( ( stripe_width == 0 ) || ( stripe_width > spin.size() ) ) {
    return false;
  } else {
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        if ( ( col / stripe_width ) % 2 == 0 ) {
          spin[line][col].set( -1 );
        } else {
          spin[line][col].set( 1 );
        }
      }
    }
    return true;
  }
}


bool Ising2dDipSqr::prepare( const char& mode )
{
  switch ( mode ) {

  case 'r':	// completely random state ...
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        if ( gsl_rng_uniform_int( rng, 2 ) == 0 ) {
          spin[line][col].flip();
        }
      }
    }
    break;

  case 'u': // sets all spins up
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        spin[line][col].set( +1 );
      }
    }
    break;

  case 'd': // sets all spins down
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        spin[line][col].set( -1 );
      }
    }
    break;

  case 'c': // checkerboard (afm ground state)
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        if ( ( line + col ) % 2 == 0 ) {
          spin[line][col].set( -1 );
        } else {
          spin[line][col].set( 1 );
        }
      }
    }
    break;

  case 's': { // striped (automatically try to find ground state)
    unsigned int optimal_width = 1;
    prepare_striped( 1 );
    double optimal_energy = h();
    for ( unsigned int width = 2; width < spin.size(); width++ ) {
      prepare_striped( width );
      if ( optimal_energy > h() ) {
        optimal_width = width;
        optimal_energy = h();
      }
    }
    prepare_striped( optimal_width );
  }
  break;

  default:
    if ( isdigit( mode ) ) {
      // striped with user defined width
      return prepare_striped( atoi( &mode ) );
    } else {
      // unknown mode?
      return false;
    }
  }
  return true;
}


void Ising2dDipSqr::metropolis_singleflip()
{
  // find a random spin to flip
  int flip_line = gsl_rng_uniform_int( rng, size );
  int flip_col  = gsl_rng_uniform_int( rng, size );

  // flip it!
  spin[flip_line][flip_col].flip();

  // calculate energy difference
  double deltaH = - 2.0 * B * spin[flip_line][flip_col].get();
  //cout << "  " << deltaH << endl;
  for ( int line = 0; line < size; line++ ) {
    for ( int col = 0; col < size; col++ ) {
      int dl = min( abs( flip_line - line ), size - abs( flip_line - line ) );
      int dc = min( abs( flip_col - col ), size - abs( flip_col - col ) );
      deltaH += - 2.0 * effint[dl][dc]
                * ( spin[flip_line][flip_col] * spin[line][col] );
    }
  }

  if ( deltaH > 0.0 ) {
    // accept the new state?
    if ( gsl_rng_uniform( rng ) > exp( - deltaH / T ) ) {
      // new state rejected ... reverting!
      spin[flip_line][flip_col].flip();
    }
  }
}


void Ising2dDipSqr::mcstep()
{
  for ( unsigned long int n = 1; n <= N; n++ ) {
    metropolis_singleflip();
  }
  time++;
}


void Ising2dDipSqr::mcstep_dry( const unsigned int& k_max )
{
  for ( unsigned int k = 0; k < k_max; k++ ) {
    mcstep();
  }
  time = 0;

  if ( fsize_correction_mode == 2 ) {
    // try do determine if the system is in the ordered phase
    unsigned int msmall_count = 0, mlarge_count = 0;
    for ( unsigned int k = 0; k < k_max; k++ ) {
      mcstep();
      if ( abs( M() ) < N / 2 ) {
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
    time = 0;
  }
}


double Ising2dDipSqr::H() const
{
  // measures the system's energy

  double H = 0;

  // spin-spin-interactions between S_i and S_j
  for ( int iline = 0; iline < size; iline++ ) {
    for ( int icol = 0; icol < size; icol++ ) {

      for ( int jline = 0; jline < size; jline++ ) {
        for ( int jcol = 0; jcol < size; jcol++ ) {

          int dl = min( abs( jline - iline ), size - abs( jline - iline ) );
          int dc = min( abs( jcol - icol ), size - abs( jcol - icol ) );

          H += - effint[dl][dc] * ( spin[iline][icol] * spin[jline][jcol] );

        }
      }

    }
  }
  // TODO: double counting correction?
  H /= 2;

  // energy in external magnetic field
  if ( B != 0 ) {
    for ( int line = 0; line < size; line++ ) {
      for ( int col = 0; col < size; col++ ) {
        H += - B * spin[line][col].get();
      }
    }
  }

  /*double oldH = H_old();
  if (abs(1 - (H / oldH)) > 0.01) {
      cout << "H = " << H << " <--/--> oldH = " << oldH << endl;
      exit(1);
  }*/

  return H;
}


double Ising2dDipSqr::h() const
{
  // measures the system's energy per spin
  return H() / N;
}


unsigned long int Ising2dDipSqr::t() const
{
  // measures the system's time in lattice sweeps (MC time units)
  return time;
}


int Ising2dDipSqr::M() const
{
  // measures the system's magnetization
  int M = 0;
  for ( int line = 0; line < size; line++ ) {
    for ( int col = 0; col < size; col++ ) {
      M += spin[line][col].get();
    }
  }

  // finite size corrections to the magnetization
  if ( ( fsize_correction_mode == 1 ) ||
       ( ( fsize_correction_mode == 2 ) && fsize_ordered_phase ) ) {
    M = abs( M );
  }
  return M;
}


double Ising2dDipSqr::m() const
{
  // measures the system's magnetization per spin
  return double( M() ) / N;
}


vector<double> Ising2dDipSqr::ss_corr() const
{
  // TODO: think about it ...

  // measure spin-spin correlations
  vector<double> result;
  vector<unsigned int> samples;
  result.resize( spin.size(), 0 );
  samples.resize( spin.size(), 0 );
  for ( int i = 0; i < size; i++ ) {
    for ( int j = 0; j < size; j++ ) {
      result[abs( int( i - j ) )] += spin[i][i] * spin[i][j];
      result[abs( int( i - j ) )] += spin[i][i] * spin[j][i];
      samples[abs( int( i - j ) )] += 2;
    }
  }
  for ( int d = 0; d < size; d++ ) {
    result[d] /= samples[d];
  }
  return result;
}


unsigned int Ising2dDipSqr::spin_count() const
{
  // returns the total number of spins in the system
  return N;
}


void Ising2dDipSqr::special_perbin( unsigned int& mode )
{
  // the only special thing we want to do is to calculate the structure factor
  if ( mode != 0 ) {
    sf_calc();
  }
}


void Ising2dDipSqr::sf_calc()
{
  // calculate the structure factor function
  for ( unsigned int line = 0; line < k_points; line++ ) {
    for ( unsigned int col = 0; col < k_points; col++ ) {
      double kx = line * k_range / ( k_points - 1 ) - k_range / 2.0;
      double ky = col * k_range / ( k_points - 1 ) - k_range / 2.0;
      sf_sum[line][col] += sf_k_calc( kx, ky );
    }
  }
  sf_measurements++;
}


double Ising2dDipSqr::sf_k_calc( double& kx, double& ky )
{
  // calculate the structure factor sf(k) for a specific k

  complex<double> z( 0.0, 0.0 );
  complex<double> i( 0.0, 1.0 );
  complex<double> kr( 0.0, 0.0 );

  for ( int y = 0; y < size; y++ ) {
    for ( int x = 0; x < size; x++ ) {
      z += double( spin[y][x].get() ) * exp( i * ( kx * x + ky * y ) );
    }
  }

  return abs( z ) * abs( z );
}
