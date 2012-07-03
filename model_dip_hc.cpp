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

#include "model_dip_hc.hpp"


// ----- HELPER CLASS: EFFECTIVE INTERACTION COEFFICIENTS -----

template <typename T>
Ising2dDipHC_effint<T>::Ising2dDipHC_effint( const unsigned int& N_init ) : N( N_init )
{
    elem.resize( N * N, 0.0 );
}

template <typename T>
T& Ising2dDipHC_effint<T>::operator()( const unsigned int& s1, const unsigned int& s2 )
{
    return elem[N * s1 + s2];
}

template <typename T>
T Ising2dDipHC_effint<T>::operator()( const unsigned int& s1, const unsigned int& s2 ) const
{
    return elem[N * s1 + s2];
}




// ----- HELPER STRUCT: SET OF BASIS VECTOR COEFFICIENTS -----

Ising2dDipHC_spinpos::Ising2dDipHC_spinpos( const short& nb1_init, const short& nb2_init, const short& nc_init )
    : nb1( nb1_init ), nb2( nb2_init ), nc( nc_init ) { }




// ------ 2D ISING MODEL ON A HONEYCOMB LATTICE WITH DIP-DIP-INT -----

Ising2dDipHC::Ising2dDipHC( const unsigned int& size_init, const double& J, const double& g, const double& B, const double& T,
                            const unsigned int& fcm_init, const string& cwd )
    : SystemModel( true, B, T, cwd ), J( J ), g( g ), size( size_init ), N( 0 ),
      fsize_correction_mode( fcm_init ), fsize_ordered_phase( false ),
      effint( Ising2dDipHC_effint<float>( 0 ) ), sf_measurements( 0 )
{
    // calculate all 3 kinds of lattice vectors
    // unit cell
    b1_vec.x = sqrt( 3.0 ) * cos( 30.0 * M_PI / 180.0 );
    b1_vec.y = sqrt( 3.0 ) * sin( 30.0 * M_PI / 180.0 );
    b2_vec.x = sqrt( 3.0 ) * cos( -30.0 * M_PI / 180.0 );
    b2_vec.y = sqrt( 3.0 ) * sin( -30.0 * M_PI / 180.0 );
    // basis
    c_vec.x = 1.0;
    c_vec.y = 0.0;
    // system basis (for periodic boundary conditions)
    a1_vec = short( 2 * size - 1 ) * b1_vec - short( size - 1 ) * b2_vec;
    a2_vec = short( size - 1 )   * b1_vec + short( size )    * b2_vec;

    // calculate the spin positions
    for ( short nb1 = -( size - 1 ); nb1 <= size - 1; ++nb1 ) {
        if ( nb1 < 0 ) {
            for ( short nb2 = -( size - 1 ) - nb1; nb2 <= size - 1; ++nb2 ) {
                for ( short nc = -1; nc <= 1; nc += 2 ) {
                    spin_pos.push_back( Ising2dDipHC_spinpos( nb1, nb2, nc ) );
                    ++N;
                }
            }
        } else {
            for ( short nb2 = -( size - 1 ); nb2 <= ( size - 1 ) - nb1; ++nb2 ) {
                for ( short nc = -1; nc <= 1; nc += 2 ) {
                    spin_pos.push_back( Ising2dDipHC_spinpos( nb1, nb2, nc ) );
                    ++N;
                }
            }
        }
    }
    // ... and create new spins for all of these positions
    spin.resize( N );

    // set up an object that holds the effective interations strengths:
    effint = Ising2dDipHC_effint<float>( N );

    // calculate the effective interaction strengths!

    if ( J != 0.0 ) {
        // how many nearest neighbours have we found for each spin?
        vector<float> spin_neighbours;
        spin_neighbours.resize( N, 0 );

        // step 1: nearest neighbour interaction within the system
        for ( unsigned int i = 0; i < N; ++i ) {
            Vec2 si_pos = spin_pos[i].nb1 * b1_vec + spin_pos[i].nb2 * b2_vec + spin_pos[i].nc * c_vec;
            for ( unsigned int j = 0; j < N; ++j ) {
                Vec2 sj_pos = spin_pos[j].nb1 * b1_vec + spin_pos[j].nb2 * b2_vec + spin_pos[j].nc * c_vec;

                // nearest neighbours?
                double distance = abs( sj_pos - si_pos );
                if ( ( distance < 1.5 ) && ( distance > 0.5 ) ) {
                    spin_neighbours[i] += 1;
                    effint( i, j ) = J;
                }
            }
        }

        // step 2: nearest neighbour interaction across the system's boundaries
        short nsystems_na1[6] = {1, 1, 0, 0, -1, -1};
        short nsystems_na2[6] = {0, -1, 1, -1, 1, 0};
        for ( unsigned int i = 0; i < N; ++i ) {
            // iterate over all spins

            // calculate position of S_i
            Vec2 si_pos = spin_pos[i].nb1 * b1_vec + spin_pos[i].nb2 * b2_vec + spin_pos[i].nc * c_vec;

            for ( short nsystem = 0; nsystem < 6; nsystem++ ) {
                // iterate over all 6 neighbouring systems

                // calculate offset into the neighbouring system
                Vec2 nsystem_offset = nsystems_na1[nsystem] * a1_vec + nsystems_na2[nsystem] * a2_vec;

                for ( unsigned int j = 0; j < N; ++j ) {
                    // iterate over all spins in the neighbouring system

                    // calculate position of S_j in the central system
                    Vec2 sj_pos = spin_pos[j].nb1 * b1_vec + spin_pos[j].nb2 * b2_vec + spin_pos[j].nc * c_vec;

                    // add the offset into the neighbouring system
                    sj_pos = sj_pos + nsystem_offset;

                    // S_i and S_j in the neighbouring system are nearest neighbours?
                    double distance = abs( sj_pos - si_pos );
                    if ( ( distance < 1.5 ) && ( distance > 0.5 ) ) {
                        spin_neighbours[i] += 1;
                        effint( i, j ) = J;
                    }
                }
            }
        }

        // check if everyone has 3 nearest neigbours
        for ( unsigned int i = 0; i < N; ++i ) {
            if ( spin_neighbours[i] != 3 ) {
                cout << "WARNING: spin " << i << " has " << spin_neighbours[i] << " nearest neighbours?" << endl;
            }
        }
    }

    if ( g != 0.0 ) {
        // step 3: dip-dip interaction for spins that are not too far away
        // double startTime = current_time();
        for ( unsigned int i = 0; i < N; ++i ) {
            // iterate over S_i

            // calculate position of S_i
            Vec2 si_pos = spin_pos[i].nb1 * b1_vec + spin_pos[i].nb2 * b2_vec + spin_pos[i].nc * c_vec;

            for ( unsigned int j = i; j < N; ++j ) {
                // iterate over S_j

                // calculate position of S_j in the central system
                Vec2 sj_pos = spin_pos[j].nb1 * b1_vec + spin_pos[j].nb2 * b2_vec + spin_pos[j].nc * c_vec;

                // initialize dipole-dipole interaction to 0 (we will add to it later)
                double dipint = 0;

                for ( unsigned int sys = 0; sys < N; ++sys ) {
                    // iterate over system copies
                    // (we also use the spin_pos coefficient as the coefficients to basis a)

                    // calculate offset into the copied system
                    Vec2 sys_offset = spin_pos[sys].nb1 * a1_vec + spin_pos[sys].nb2 * a2_vec;

                    // distance of S_i and S_j (in the system's copy)
                    double distance = abs( ( sj_pos + sys_offset ) - si_pos );

                    // make sure S_i and S_j are not the same spin (can only happen of sys_offset = vec 0)
                    if ( !( distance < 0.1 ) ) {
                        // ... different spins -> there is a dip-dip interaction
                        dipint += g / ( distance * distance * distance );
                    }
                }

                // iteration over all system copies is finished -> add accumulated effint
                effint( i, j ) += dipint;
                effint( j, i ) += dipint;

                // cout << effint(i,j) << endl;
            }
            //cout << "effint calc finished: " << i+1 << "/" << N << endl;
        }
        // cout << "ddCalcTime: " << current_time() - startTime << "ms" << endl;

        // step 4: add a mean field approximation for the dip-dip-int with spins that are very far away!
        double mf_adder = 0.0;
        const short mf_size = 10000;
        for ( short na1 = -( mf_size - 1 ); na1 <= mf_size - 1; ++na1 ) {
            if ( na1 < 0 ) {
                for ( short na2 = -( mf_size - 1 ) - na1; na2 <= mf_size - 1; ++na2 ) {
                    if ( !( na1 >= -( size - 1 ) && na1 <= size - 1 &&
                            na2 >= -( size - 1 ) - na1 && na2 <= size - 1 ) ) {
                        // we are outside of the central region

                        // cout << na1 << ' ' << na2 << endl;
                        Vec2 position = na1 * a1_vec + na2 * a2_vec;
                        double distance = abs( position );
                        mf_adder += g / ( distance * distance * distance );
                    }
                }
            } else {
                for ( short na2 = -( mf_size - 1 ); na2 <= ( mf_size - 1 ) - na1; ++na2 ) {
                    if ( !( na1 >= -( size - 1 ) && na1 <= size - 1 &&
                            na2 >= -( size - 1 ) && na2 <= ( size - 1 ) - na1 ) ) {
                        // we are outside of the central region

                        // cout << na1 << ' ' << na2 << endl;
                        Vec2 position = na1 * a1_vec + na2 * a2_vec;
                        double distance = abs( position );
                        mf_adder += g / ( distance * distance * distance );
                    }
                }
            }
        }
        cout << "approx mf_effint due to long range dip-dip-int: " << mf_adder << endl;
        // add it to all effective interaction coefficients
        for ( unsigned int i = 0; i < N; ++i ) {
            for ( unsigned int j = 0; j < N; ++j ) {
                effint( i, j ) += mf_adder;
            }
        }
    }

    // check if any effective interaction coefficients are asymmetric
    for ( unsigned int i = 0; i < N; ++i ) {
        for ( unsigned int j = 0; j < N; ++j ) {
            if ( abs( effint( i, j ) - effint( j, i ) ) > 0.001 ) {
                cout << "WARNING: effint(" << i << "," << j << ") != effint(" << j << "," << i << ")"
                     << "... diff = " << abs( effint( i, j ) - effint( j, i ) ) << endl;
            }
        }
    }

    // set up an array for the structure factor
    sf_sum.resize( 100 );
    for ( unsigned int nr = 0; nr < 100; ++nr ) {
        sf_sum[nr].resize( 360 );
    }

    // set up a red, green, black palette
    pal.push_back( png::color( 255, 0, 0 ) );
    pal.push_back( png::color( 0, 255, 0 ) );
    pal.push_back( png::color( 0, 0, 0 ) );
}


Ising2dDipHC::~Ising2dDipHC()
{
    // write structure factor to file
    ofstream sf_log;
    sf_log.open( ( cwd + "sf.log" ).c_str() );
    if ( !sf_log.is_open() ) {
        cout << "ERROR while opening structure factor log file in " << cwd << endl;
        exit( 1 );
    }
    ostream_setup( sf_log );

    for ( unsigned short nr = 0; nr < 100; ++nr ) {
        for ( unsigned short nphi = 0; nphi < 360; ++nphi ) {
            double r = ( nr + 1 ) * M_PI / 100.0;
            double phi = nphi * M_PI / 180.0;
            Vec2 k_vec( r * cos( phi ), r * sin( phi ) );
            sf_log << k_vec.x << ' ' << k_vec.y << ' ' << sf_sum[nr][nphi] / sf_measurements / N << endl;
        }
    }

    sf_log.close();

    // write structure factor plotting script
    ofstream sf_plot;
    sf_plot.open( ( cwd + "sf_plot.pyx" ).c_str() );
    if ( !sf_plot.is_open() ) {
        cout << "ERROR while creating structure factor pyxplot file in " << cwd << endl;
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
    set xtics (\"$\\pi$\" pi, \"$\\frac{3\\pi}{4}$\" 3*pi/4, \"$\\frac{\\pi}{2}$\" pi/2, \"$\\frac{\\pi}{4}$\" pi/4, \"0\" 0, \\  \n\
               \"$-\\pi$\" -pi, \"$-\\frac{3\\pi}{4}$\" -3*pi/4, \"$-\\frac{\\pi}{2}$\" -pi/2, \"$-\\frac{\\pi}{4}$\" -pi/4)  \n\
    set mxtics pi/8 \n\
    set ylabel \"$k_y$\"  \n\
    set ytics (\"$\\pi$\" pi, \"$\\frac{3\\pi}{4}$\" 3*pi/4, \"$\\frac{\\pi}{2}$\" pi/2, \"$\\frac{\\pi}{4}$\" pi/4, \"0\" 0, \\  \n\
               \"$-\\pi$\" -pi, \"$-\\frac{3\\pi}{4}$\" -3*pi/4, \"$-\\frac{\\pi}{2}$\" -pi/2, \"$-\\frac{\\pi}{4}$\" -pi/4) \n\
    set mytics pi/8 \n\
    plot [-pi:pi][-pi:pi] 'sf.log' with colourmap notitle";
    sf_plot.close();

    // write the last microstate to a file
    ofstream last_mstate_file;
    last_mstate_file.open( ( cwd + "last_mstate.log" ).c_str() );
    if ( !last_mstate_file.is_open() ) {
        cout << "ERROR while opening output file for the last microstate in " << cwd << endl;
        exit( 1 );
    }
    for ( unsigned int i = 0; i < N; ++i ) {
        last_mstate_file << spin[i].get() << endl;
    }
    last_mstate_file.close();
}


png::image< png::index_pixel > Ising2dDipHC::get_image() const
{
    png::image< png::index_pixel > image( 6 * size + 8, 8 * size );
    image.set_palette( pal );

    // paint the entire picture black
    for ( size_t line = 0; line < image.get_height(); ++line ) {
        for ( size_t col = 0; col < image.get_width(); ++col ) {
            image[line][col] = png::index_pixel( 2 );
        }
    }

    // basis vectors for drawing
    Vec2 b1_drawvec( 3.0, 2.0 );
    Vec2 b2_drawvec( 3.0, -2.0 );
    Vec2  c_drawvec( 2.0, 0.0 );

    // find center of the picture
    Vec2 center( ( image.get_width() - 1 ) / 2.0, ( image.get_height() - 1 ) / 2.0 );

    for ( unsigned int i = 0; i < N; ++i ) {
        Vec2 imgpos = center + spin_pos[i].nb1 * b1_drawvec
                      + spin_pos[i].nb2 * b2_drawvec + spin_pos[i].nc * c_drawvec;
        double offset_x[4] = {0.5, 0.5, -0.5, -0.5};
        double offset_y[4] = {0.5, -0.5, -0.5, 0.5};
        for ( short opix = 0; opix < 4; ++opix ) {
            image[static_cast<unsigned int>( round( imgpos.y + offset_y[opix] ) )]
            [static_cast<unsigned int>( round( imgpos.x + offset_x[opix] ) )]
                = ( spin[i].get() == 1 ) ? png::index_pixel( 0 ) : png::index_pixel( 1 );
        }
    }

    return image;
}


bool Ising2dDipHC::prepare_striped( const int& stripe_width )
{
    if ( ( stripe_width == 0 ) || ( stripe_width > size ) ) {
        return false;
    } else {
        for ( unsigned int i = 0; i < N; ++i ) {
            short new_s;
            if ( ( ( spin_pos[i].nb1 + size - 1 ) / stripe_width ) % 2 == 0 ) {
                new_s = +1;
            } else {
                new_s = -1;
            }
            if ( ( ( spin_pos[i].nb1 + size - 1 ) / stripe_width ) % 2 != ( ( spin_pos[i].nb1 + size ) / stripe_width ) % 2
                    && spin_pos[i].nc == 1 ) {
                new_s *= -1;
            }
            spin[i].set( new_s );
        }
        return true;
    }
}


bool Ising2dDipHC::prepare( const char& mode )
{
    switch ( mode ) {

    case 'r': // completely random state ...
        for ( unsigned int i = 0; i < N; ++i ) {
            if ( gsl_rng_uniform_int( rng, 2 ) == 0 ) {
                spin[i].flip();
            }
        }
        break;

    case 'u': // sets all spins up
        for ( unsigned int i = 0; i < N; ++i ) {
            spin[i].set( +1 );
        }
        break;

    case 'd': // sets all spins down
        for ( unsigned int i = 0; i < N; ++i ) {
            spin[i].set( -1 );
        }
        break;

    case 's': { // striped (automatically try to find ground state)
        unsigned int optimal_width = 1;
        prepare_striped( 1 );
        double optimal_energy = h();
        for ( int width = 2; width < size; width++ ) {
            prepare_striped( width );
            cout << width << ' ' << h() << endl;
            if ( optimal_energy > h() ) {
                optimal_width = width;
                optimal_energy = h();
            }
        }
        prepare_striped( optimal_width );
    }
    break;

    case 'f': { // read microstate from last_mstate.log
        ifstream spin_init_file;
        spin_init_file.open( "spin_init.log" );
        if ( !spin_init_file.is_open() ) {
            cout << "ERROR while opening spin init file in " << cwd << endl;
            exit( 1 );
        }
        short s = 0;
        unsigned int i = 0;
        while ( ( spin_init_file >> s ) && ( i < N ) ) {
            spin[i].set( s );
            ++i;
        }
        spin_init_file.close();
    }
    break;

    default: {
        if ( isdigit( mode ) ) {
            // striped with user defined width
            return prepare_striped( atoi( &mode ) );
        } else {
            // unknown mode?
            return false;
        }
    }

    }
    return true;
}


void Ising2dDipHC::metropolis_singleflip()
{
    // find a random spin to flip
    unsigned int flipper = gsl_rng_uniform_int( rng, N );

    // flip it!
    spin[flipper].flip();

    // calculate energy difference
    double deltaH = - 2.0 * B * spin[flipper].get();
    //cout << "  " << deltaH << endl;
    for ( unsigned int i = 0; i < N; ++i ) {
        deltaH += - 2.0 * effint( flipper, i ) * ( spin[flipper] * spin[i] );
    }

    if ( deltaH > 0 ) {
        // accept the new state?
        if ( gsl_rng_uniform( rng ) > exp( - deltaH / T ) ) {
            // new state rejected ... reverting!
            spin[flipper].flip();
        }
    }
}


void Ising2dDipHC::mcstep()
{
    for ( unsigned long int n = 1; n <= N; n++ ) {
        metropolis_singleflip();
    }
    time++;
}


void Ising2dDipHC::mcstep_dry( const unsigned int& k_max )
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


double Ising2dDipHC::H() const
{
    // measures the system's energy

    double H = 0;

    // spin-spin-interactions between S_i and S_j
    for ( unsigned int i = 0; i < N; ++i ) {
        for ( unsigned int j = i; j < N; ++j ) {
            H += - effint( i, j ) * ( spin[i] * spin[j] );
        }
    }

    // energy in external magnetic field
    if ( B != 0 ) {
        for ( unsigned int i = 0; i < N; ++i ) {
            H += - B * spin[i].get();
        }
    }

    return H;
}


double Ising2dDipHC::h() const
{
    // measures the system's energy per spin
    return H() / N;
}


unsigned long int Ising2dDipHC::t() const
{
    // measures the system's time in lattice sweeps (MC time units)
    return time;
}


int Ising2dDipHC::M() const
{
    // measures the system's magnetization
    int M = 0;
    for ( unsigned int i = 0; i < N; ++i ) {
        M += spin[i].get();
    }

    // finite size corrections to the magnetization
    if ( ( fsize_correction_mode == 1 ) ||	( ( fsize_correction_mode == 2 ) && fsize_ordered_phase ) ) {
        M = abs( M );
    }
    return M;
}


double Ising2dDipHC::m() const
{
    // measures the system's magnetization per spin
    return double( M() ) / N;
}


vector<double> Ising2dDipHC::ss_corr() const
{
    // TODO: implement!

    vector<double> result;
    result.resize( size, 0.0 );
    return result;
}


unsigned int Ising2dDipHC::spin_count() const
{
    // returns the total number of spins in the system
    return N;
}


void Ising2dDipHC::special_perbin( const unsigned int& mode )
{
    // the only special thing we want to do is to calculate the structure factor
    if ( mode != 0 ) {
        sf_calc();
    }
}


void Ising2dDipHC::sf_calc()
{
    // calculate the structure factor function
    for ( unsigned short r = 0; r < 100; ++r ) {
        for ( unsigned short phi = 0; phi < 360; ++phi ) {
            sf_sum[r][phi] += sf_k_calc( r, phi );
        }
    }
    ++sf_measurements;
}


double Ising2dDipHC::sf_k_calc( const unsigned int& nr, const unsigned int& nphi )
{
    // calculate the structure factor sf(k) for a specific k

    double r = ( nr + 1 ) * M_PI / 100.0;
    double phi = nphi * M_PI / 180.0;
    Vec2 k_vec( r * cos( phi ), r * sin( phi ) );

    complex<double> z( 0.0, 0.0 );
    for ( unsigned int i = 0; i < N; ++i ) {
        Vec2 r_vec = spin_pos[i].nb1 * b1_vec + spin_pos[i].nb2 * b2_vec + spin_pos[i].nc * c_vec;
        z += double( spin[i].get() ) * exp( complex<double>( 0.0, 1.0 ) * ( k_vec * r_vec ) );
    }

    return abs( z ) * abs( z );
}
