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

#ifndef _MODEL_ISING2DSQRDIPOLEMET_H_INCLUDED
#define _MODEL_ISING2DSQRDIPOLEMET_H_INCLUDED

#include <iostream>
#include <cmath>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "isingspin.hpp"
#include "model_ising2dsqrmet.hpp"


class IsingModel2dDipole: public IsingModel2d
{
    // ------- 2D ISING MODEL WITH LONG RANGE DIPOLE INTERACTION ---------

  protected:

    // strength of dipole-dipole interaction
    double g;
    double H_memory;

  public:

    IsingModel2dDipole( const unsigned int& N, const bool& periodic,
                        const double& J, const double& g,
                        const double& B, const double& T, string cwd );

    bool prepare( const char& mode );

    void metropolis_singleflip();
    void metropolis_mirror();
    void metropolis_blockflip();

    void mcstep();

    double H() const;
};

#endif // _MODEL_ISING2DSQRDIPOLEMET_H_INCLUDED
