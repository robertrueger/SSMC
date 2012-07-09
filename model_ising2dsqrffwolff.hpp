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

#ifndef _MODEL_ISING2DSQRFFWOLFF_H_INCLUDED
#define _MODEL_ISING2DSQRFFWOLFF_H_INCLUDED

#include <cmath>
#include <vector>
using namespace std;

#include <gsl/gsl_rng.h>
#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "isingspin.hpp"
#include "model_ising2dsqrmet.hpp"


struct lsite_2dsquare {
  // ----- HELPER STRUCT: POSITION ON A 2D SQUARE LATICE -----
  unsigned short int line, col;
  lsite_2dsquare();
  lsite_2dsquare( const unsigned int& line, const unsigned int& col );
};


class IsingModel2dWolff: public IsingModel2d
{
    // ----- 2D FIELD-FREE ISING MODEL (WOLFF-ALGORITHM) -----

  protected:

    vector< vector<bool> > mask;
    vector< lsite_2dsquare > mask_items;
    vector< lsite_2dsquare > mask_candidates;

    double add_prob;

    vector<unsigned long int> cluster_size;

  public:

    IsingModel2dWolff( const unsigned int& N, const bool& periodic,
                       const double& J, const double& T,
                       const unsigned int& fsize_correction, const string& cwd );

    png::image< png::index_pixel > get_image() const;

    void wolff_clusterflip();

    void mcstep();
};

#endif // _MODEL_ISING2DSQRFFWOLFF_H_INCLUDED
