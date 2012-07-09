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

#ifndef _MODEL_DIP_SQR_H_INCLUDED
#define _MODEL_DIP_SQR_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "isingspin.hpp"
#include "utils.hpp"


class Ising2dDipSqr: public SystemModel
{
    // ------ 2D ISING MODEL ON A SQUARE LATTICE WITH DIP-DIP-INT -----

  protected:

    double J, g;
    unsigned long int N;
    int size;
    vector< vector<IsingSpin> > spin;

    // finite size correction mode
    unsigned int fsize_correction_mode;
    bool fsize_ordered_phase;

    // precalculated interaction coefficients
    vector< vector<double> > effint;

    // structure factor calculation parameters
    double k_range;
    unsigned int k_points;
    vector< vector<double> > sf_sum;
    unsigned int sf_measurements;

  public:

    Ising2dDipSqr( const unsigned int& size, const bool& periodic,
                   const double& J, const double& g,
                   const double& B, const double& T,
                   const unsigned int& fsize_correction_mode,
                   const string& cwd );
    ~Ising2dDipSqr();

    double calc_effint_dipdip( const int& dl_seed, const int& dc_seed,
                               const int& system_clones );

    png::image< png::index_pixel > get_image() const;

    bool prepare_striped( const unsigned int& stripe_width );
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

    void special_perbin( unsigned int& mode );

    void sf_calc();
    double sf_k_calc( double& kx, double& ky );
};

#endif // _MODEL_DIP_SQR_H_INCLUDED
