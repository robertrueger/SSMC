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

#include "model_ising2dsqrffwolff.hpp"


// ----- HELPER STRUCT: POSITION ON A 2D SQUARE LATICE -----

lsite_2dsquare::lsite_2dsquare() { }


lsite_2dsquare::lsite_2dsquare( const unsigned int& line, const unsigned int& col ) : line( line ), col( col ) { }




// ----- 2D FIELD-FREE ISING MODEL (WOLFF-ALGORITHM) -----

IsingModel2dWolff::IsingModel2dWolff( const unsigned int& N, const bool& periodic,
                                      const double& J, const double& T,
                                      const unsigned int& fsize_correction, const string& cwd )
    : IsingModel2d( N, periodic, J, 0, T, fsize_correction, cwd ),
      add_prob( 1 - exp( -2 * J / T ) )
{
    mask.resize( size );
    for ( unsigned int line = 0; line < N; line++ ) {
        mask[line].resize( size );
    }

    cluster_size.resize( 100 );
    for ( unsigned int k = 0; k < 100; k++ ) {
        cluster_size[k] = N * N / 2;
    }

    pal.push_back( png::color( 255, 0, 0 ) );
    pal.push_back( png::color( 0, 255, 0 ) );
}


png::image< png::index_pixel > IsingModel2dWolff::get_image() const
{
    png::image< png::index_pixel > image( size, size );
    image.set_palette( pal );
    for ( size_t line = 0; line < image.get_height(); ++line ) {
        for ( size_t col = 0; col < image.get_width(); ++col ) {
            if ( spin[line][col].get() == 1 ) {
                if ( mask[line][col] ) {
                    image[line][col] = png::index_pixel( 3 );
                } else {
                    image[line][col] = png::index_pixel( 1 );
                }
            } else if ( spin[line][col].get() == -1 ) {
                if ( mask[line][col] ) {
                    image[line][col] = png::index_pixel( 2 );
                } else {
                    image[line][col] = png::index_pixel( 0 );
                }
            }
        }
    }
    return image;
}


void IsingModel2dWolff::wolff_clusterflip()
{
    // flips a cluster according to the Wolff-Algorithm

    // reset the mask in which we build the cluster
    for ( unsigned int k = 0; k < mask_items.size(); k++ ) {
        mask[mask_items[k].line][mask_items[k].col] = false;
    }
    mask_items.clear();

    // find a seed spin and add it to the cluster
    lsite_2dsquare seed;
    seed.line = gsl_rng_uniform_int( rng, size );
    seed.col  = gsl_rng_uniform_int( rng, size );
    mask[seed.line][seed.col] = true;
    mask_items.push_back( seed );

    // add the seed's neighbours to the candidate list if they are pointing
    // in the same direction ...
    if ( spin[( seed.line + size - 1 ) % size][seed.col].get()
            == spin[seed.line][seed.col].get() ) {
        // top neighbour
        mask_candidates.push_back(
            lsite_2dsquare( ( seed.line + size - 1 ) % size, seed.col )
        );
    }
    if ( spin[( seed.line + 1 ) % size][seed.col].get()
            == spin[seed.line][seed.col].get() ) {
        // bottom neighbour
        mask_candidates.push_back(
            lsite_2dsquare( ( seed.line + 1 ) % size, seed.col )
        );
    }
    if ( spin[seed.line][( seed.col + size - 1 ) % size].get()
            == spin[seed.line][seed.col].get() ) {
        // left neighbour
        mask_candidates.push_back(
            lsite_2dsquare( seed.line, ( seed.col + size - 1 ) % size )
        );
    }
    if ( spin[seed.line][( seed.col + 1 ) % size].get()
            == spin[seed.line][seed.col].get() ) {
        // right neighbour
        mask_candidates.push_back(
            lsite_2dsquare( seed.line, ( seed.col + 1 ) % size )
        );
    }

    // build the cluster ...
    while ( !mask_candidates.empty() ) {

        // read a candidate from the list
        lsite_2dsquare cand = mask_candidates.back();
        mask_candidates.pop_back();

        // candidate has already been added to the cluster?
        if ( mask[cand.line][cand.col] ) {
            continue;
        }

        // add the candidate?
        if ( gsl_rng_uniform( rng ) < add_prob ) {
            mask[cand.line][cand.col] = true;
            mask_items.push_back( lsite_2dsquare( cand.line, cand.col ) );

            // add its neighbours to the candidate list
            if ( ( spin[( cand.line + size - 1 ) % size][cand.col].get()
                    == spin[cand.line][cand.col].get() )
                    && !mask[( cand.line + size - 1 ) % size][cand.col] ) {
                // top neighbour
                mask_candidates.push_back(
                    lsite_2dsquare( ( cand.line + size - 1 ) % size, cand.col )
                );
            }
            if ( ( spin[( cand.line + 1 ) % size][cand.col].get()
                    == spin[cand.line][cand.col].get() )
                    && !mask[( cand.line + 1 ) % size][cand.col] ) {
                // bottom neighbour
                mask_candidates.push_back(
                    lsite_2dsquare( ( cand.line + 1 ) % size, cand.col )
                );
            }
            if ( ( spin[cand.line][( cand.col + size - 1 ) % size].get()
                    == spin[cand.line][cand.col].get() )
                    && !mask[cand.line][( cand.col + size - 1 ) % size] ) {
                // left neighbour
                mask_candidates.push_back(
                    lsite_2dsquare( cand.line, ( cand.col + size - 1 ) % size )
                );
            }
            if ( ( spin[cand.line][( cand.col + 1 ) % size].get()
                    == spin[cand.line][cand.col].get() )
                    && !mask[cand.line][( cand.col + 1 ) % size] ) {
                // right neighbour
                mask_candidates.push_back(
                    lsite_2dsquare( cand.line, ( cand.col + 1 ) % size )
                );
            }
        }
    }

    // flip all the spins in the cluster
    for ( unsigned int k = 0; k < mask_items.size(); k++ ) {
        spin[mask_items[k].line][mask_items[k].col].flip();
    }

    // update cluster size
    cluster_size[gsl_rng_uniform_int( rng, 100 )] = mask_items.size();
}


void IsingModel2dWolff::mcstep()
{
    // recalculate adding probability (needed for SA)
    add_prob = 1 - exp( -2 * J / T );

    // VIDEO MODE
    //wolff_clusterflip();
    //time++;
    //return;

    // calculate the mean cluster size of the last mcsteps
    double mean_cluster_size = 0;
    for ( unsigned int k = 0; k < 100; k++ ) {
        mean_cluster_size += cluster_size[k];
    }
    mean_cluster_size = mean_cluster_size / 100;
    // add a litte bit of additional randomness to the cluster size
    mean_cluster_size *= ( gsl_rng_uniform( rng ) + 0.5 );

    // flip as many clusters as needed so that N spins have been flipped
    for ( unsigned int k = 0; k < ( uint ) ceil( N / mean_cluster_size ); k++ ) {
        wolff_clusterflip();
    }
    time++;
}
