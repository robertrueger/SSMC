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

#ifndef _MODEL_DIP_HC_H_INCLUDED
#define _MODEL_DIP_HC_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "isingspin.hpp"
#include "utils.hpp"
#include "utils_vec2.hpp"


template <typename T>
class Ising2dDipHC_effint
{
    // ----- HELPER CLASS: EFFECTIVE INTERACTION COEFFICIENTS -----

    protected:

        unsigned int N;
        vector<T> elem;

    public:

        Ising2dDipHC_effint( const unsigned int& N_init );

        T& operator()( const unsigned int& s1, const unsigned int& s2 );
        T operator()( const unsigned int& s1, const unsigned int& s2 ) const;
};


class  Ising2dDipHC_spinpos
{
    // ----- HELPER STRUCT: SET OF BASIS VECTOR COEFFICIENTS -----

    public:

        short nb1, nb2, nc;

        Ising2dDipHC_spinpos( const short& nb1_init, const short& nb2_init, const short& nc_init );
};


class Ising2dDipHC: public SystemModel
{
    // ----- 2D ISING MODEL ON A HONEYCOMB LATTICE WITH DIP-DIP-INT -----

    protected:

        double J, g;
        short size;
        unsigned int N;

        // finite size correction mode
        unsigned int fsize_correction_mode;
        bool fsize_ordered_phase;

        // basis vectors
        Vec2 a1_vec, a2_vec;
        Vec2 b1_vec, b2_vec;
        Vec2 c_vec;

        // spins
        vector<IsingSpin> spin;
        // ... and their positions
        vector<Ising2dDipHC_spinpos> spin_pos;

        // precalculated interaction coefficients
        Ising2dDipHC_effint<float> effint;

        // structure factor calculation parameters
        vector< vector<double> > sf_sum;
        unsigned int sf_measurements;

    public:

        Ising2dDipHC( const unsigned int& size_init, const double& J, const double& g, const double& B, const double& T,
                      const unsigned int& fcm_init, const string& cwd );

        ~Ising2dDipHC();

        png::image< png::index_pixel > get_image() const;

        bool prepare_striped( const int& stripe_width );
        bool prepare( const char& mode );


        // ----- METROPOLIS ALGORITHM WITH SINGLE SPIN FLIP -----

        void metropolis_singleflip();

        void mcstep();
        void mcstep_dry( const unsigned int& k_max );


        // ----- MEASUREMENT FUNCTIONS -----

        double H() const;
        double h() const;

        unsigned long int t() const;

        int M() const;
        double m() const;

        vector<double> ss_corr() const;

        unsigned int spin_count() const;

        void special_perbin( const unsigned int& mode );

        void sf_calc();
        double sf_k_calc( const unsigned int& nr, const unsigned int& nphi );
};

#endif // _MODEL_DIP_HC_H_INCLUDED
