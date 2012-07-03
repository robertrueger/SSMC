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

#include "model_ising2dsqrmet.hpp"


// ----- 2D ISING-MODEL ON A SQUARE LATTICE (METROPOLIS-ALGORITHM)

IsingModel2d::IsingModel2d( const unsigned int& N, const bool& periodic,
                            const double& J, const double& B, const double& T,
                            const unsigned int& fsize_correction_mode, const string& cwd )
    : SystemModel( periodic, B, T, cwd ), J( J ),
      N( N* N ), size( N ),
      fsize_correction_mode( fsize_correction_mode ),
      fsize_ordered_phase( false )
{
    // set up two dimensional NxN array: spin[line][column]
    spin.resize( size );
    for ( unsigned int line = 0; line < N; line++ ) {
        spin[line].resize( size );
    }

    pal.push_back( png::color( 0, 0, 0 ) );
    pal.push_back( png::color( 255, 255, 255 ) );
}


unsigned int IsingModel2d::spin_count() const
{
    // returns the total number of spins in the system
    return N;
}


png::image< png::index_pixel > IsingModel2d::get_image() const
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


bool IsingModel2d::prepare_striped( const unsigned int& stripe_width )
{
    if ( ( stripe_width == 0 ) || ( stripe_width > spin.size() ) ) {
        return false;
    } else {
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
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


bool IsingModel2d::prepare( const char& mode )
{
    switch ( mode ) {

    case 'r':	// completely random state ...
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
                if ( gsl_rng_uniform_int( rng, 2 ) == 0 ) {
                    spin[line][col].flip();
                }
            }
        }
        break;

    case 'u': // sets all spins up
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
                spin[line][col].set( +1 );
            }
        }
        break;

    case 'd': // sets all spins down
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
                spin[line][col].set( -1 );
            }
        }
        break;

    case 'c': // checkerboard (afm ground state)
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
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


void IsingModel2d::metropolis_singleflip()
{
    // find a random spin to flip
    unsigned int flip_line = gsl_rng_uniform_int( rng, size );
    unsigned int flip_col  = gsl_rng_uniform_int( rng, size );

    // flip it!
    spin[flip_line][flip_col].flip();

    // calculate energy difference
    double deltaH = - 2 * B * spin[flip_line][flip_col].get();
    if ( periodic ) {
        // top neighbour
        deltaH += -2 * J * ( spin[flip_line][flip_col]
                             * spin[( flip_line + size - 1 ) % size][flip_col] );
        // bottom neighbour
        deltaH += -2 * J * ( spin[flip_line][flip_col]
                             * spin[( flip_line + 1 ) % size][flip_col] );
        // left neighbour
        deltaH += -2 * J * ( spin[flip_line][flip_col]
                             * spin[flip_line][( flip_col + size - 1 ) % size] );
        // right neighbour
        deltaH += -2 * J * ( spin[flip_line][flip_col]
                             * spin[flip_line][( flip_col + 1 ) % size] );
    } else {
        if ( flip_line != 0 ) {
            // top neighbour
            deltaH += -2 * J * ( spin[flip_line][flip_col]
                                 * spin[flip_line - 1][flip_col] );
        }
        if ( flip_line != size - 1 ) {
            // bottom neighbour
            deltaH += -2 * J * ( spin[flip_line][flip_col]
                                 * spin[flip_line + 1][flip_col] );
        }
        if ( flip_col != 0 ) {
            // left neighbour
            deltaH += -2 * J * ( spin[flip_line][flip_col]
                                 * spin[flip_line][flip_col - 1] );
        }
        if ( flip_col != size - 1 ) {
            // right neighbour
            deltaH += -2 * J * ( spin[flip_line][flip_col]
                                 * spin[flip_line][flip_col + 1] );
        }
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
            spin[flip_line][flip_col].flip();
            return;
        }
    }
}


void IsingModel2d::mcstep()
{
    // empty precalculated exponentials (needed for SA)
    exp_precalc.clear();

    for ( unsigned long int n = 1; n <= N; n++ ) {
        metropolis_singleflip();
    }
    time++;
}


void IsingModel2d::mcstep_dry( const unsigned int& k_max )
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
        time -= k_max;
    }
}


double IsingModel2d::H() const
{
    // measures the system's energy

    double H = 0;

    // energy in external magnetic field
    if ( B != 0 ) {
        for ( unsigned int line = 0; line < size; line++ ) {
            for ( unsigned int col = 0; col < size; col++ ) {
                H += - B * spin[line][col].get();
            }
        }
    }

    // energy due to interaction within the lattice
    for ( unsigned int line = 0; line < size - 1; line++ ) {
        for ( unsigned int col = 0; col < size - 1; col++ ) {
            H += - J * ( spin[line][col] * spin[line + 1][col] ); // below
            H += - J * ( spin[line][col] * spin[line][col + 1] ); // right
        }
    }
    for ( unsigned int col = 0; col < size - 1; col++ ) {
        // horizontal neighbours in the last line
        H += - J * ( spin[size - 1][col] * spin[size - 1][col + 1] );
    }
    for ( unsigned int line = 0; line < size - 1; line++ ) {
        // vertical neighbours in the last column
        H += - J * ( spin[line][size - 1] * spin[line + 1][size - 1] );
    }

    if ( periodic ) {
        // interaction over vertical and horizontal borders
        for ( unsigned int col = 0; col < size; col++ ) {
            H += - J * ( spin[0][col] * spin[size - 1][col] );
        }
        for ( unsigned int line = 0; line < size; line++ ) {
            H += - J * ( spin[line][0] * spin[line][size - 1] );
        }
    }

    return H;
}


double IsingModel2d::h() const
{
    // measures the system's energy per spin
    return H() / N;
}


unsigned long int IsingModel2d::t() const
{
    // measures the system's time in lattice sweeps (MC time units)
    return time;
}


int IsingModel2d::M() const
{
    // measures the system's magnetization
    int M = 0;
    for ( unsigned int line = 0; line < size; line++ ) {
        for ( unsigned int col = 0; col < size; col++ ) {
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


double IsingModel2d::m() const
{
    // measures the system's magnetization per spin
    return double( M() ) / N;
}


vector<double> IsingModel2d::ss_corr() const
{
    // measure spin-spin correlations
    vector<double> result;
    vector<unsigned int> samples;
    result.resize( spin.size(), 0 );
    samples.resize( spin.size(), 0 );
    for ( unsigned int i = 0; i < size; i++ ) {
        for ( unsigned int j = 0; j < size; j++ ) {
            result[abs( int( i - j ) )] += spin[i][i] * spin[i][j];
            result[abs( int( i - j ) )] += spin[i][i] * spin[j][i];
            samples[abs( int( i - j ) )] += 2;
        }
    }
    for ( unsigned int d = 0; d < size; d++ ) {
        result[d] /= samples[d];
    }
    return result;
}
