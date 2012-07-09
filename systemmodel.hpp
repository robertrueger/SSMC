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

#ifndef _SYSTEMMODEL_H_INCLUDED
#define _SYSTEMMODEL_H_INCLUDED

#include <string>
#include <vector>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>


class SystemModel
{
    // ----- SYSTEM MODEL TEMPLATE (DEFINES INTERFACE) -----

  protected:

    // the model's properties and environment
    bool periodic;
    double B, T;

    // the model's time in Monte Carlo steps
    unsigned long int time;

    // random number generator
    gsl_rng* rng;

    // color palette for the png output
    png::palette pal;

    // working directory
    string cwd;

  public:

    SystemModel( const bool& periodic, const double& B,
                 const double& T, const string& cwd );
    virtual ~SystemModel();

    // graphical output
    virtual png::image< png::index_pixel > get_image() const = 0;

    // spin preparation
    virtual bool prepare( const char& mode ) = 0;

    // Monte Carlo step
    virtual void mcstep() = 0;
    virtual void mcstep_dry( const unsigned int& k_max ) = 0;

    // external parameter changing
    void set_T( const double& newT );
    void set_B( const double& newB );

    // observables measurements
    virtual unsigned int spin_count() const = 0;
    virtual unsigned long int t() const = 0;
    virtual double h() const = 0;
    virtual double m() const = 0;
    virtual vector<double> ss_corr() const = 0;

    // model specific function invocations
    virtual void special_permcs( const unsigned int& mode ) { };
    virtual void special_perbin( const unsigned int& mode ) { };
};

#endif // _SYSTEMMODEL_H_INCLUDED
