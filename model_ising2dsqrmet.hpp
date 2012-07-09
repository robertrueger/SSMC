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

#ifndef _MODEL_ISING2DSQRMET_H_INCLUDED
#define _MODEL_ISING2DSQRMET_H_INCLUDED

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "isingspin.hpp"


class IsingModel2d: public SystemModel
{
    // ----- 2D ISING-MODEL ON A SQUARE LATTICE (METROPOLIS-ALGORITHM)

  protected:

    double J;
    unsigned long int N;
    unsigned int size;
    vector< vector<IsingSpin> > spin;

    // finite size correction mode
    unsigned int fsize_correction_mode;
    bool fsize_ordered_phase;

    // precalculated exponentials
    map <double, double> exp_precalc;

  public:

    IsingModel2d( const unsigned int& N, const bool& periodic,
                  const double& J, const double& B, const double& T,
                  const unsigned int& fsize_correction_mode, const string& cwd );

    unsigned int spin_count() const;

    virtual png::image< png::index_pixel > get_image() const;

    virtual bool prepare_striped( const unsigned int& stripe_width );

    virtual bool prepare( const char& mode );


    // ----- METROPOLIS ALGORITHM WITH SINGLE SPIN FLIP -----

    virtual void metropolis_singleflip();

    virtual void mcstep();

    void mcstep_dry( const unsigned int& k_max );


    // ----- MEASUREMENT FUNCTIONS -----

    virtual double H() const;
    double h() const;

    unsigned long int t() const;

    int M() const;
    double m() const;

    vector<double> ss_corr() const;
};

#endif // _MODEL_ISING2DSQRMET_H_INCLUDED
