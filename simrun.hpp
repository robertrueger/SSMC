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

#ifndef _SIMRUN_H_INCLUDED
#define _SIMRUN_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

#include <png++/png.hpp>

#include "systemmodel.hpp"
#include "model_ising1dmet.hpp"
#include "model_ising2dsqrmet.hpp"
#include "model_ising2dsqrffwolff.hpp"
#include "model_ising2dsqrdipolemet.hpp"
#include "model_dip_sqr.hpp"
#include "model_dip_hc.hpp"
#include "sim_datastruct.hpp"
#include "utils.hpp"


// ----- SIMRUN: PERFORMS A SINGLE SIMULATION -----

sim_results simrun( const sim_parameters par, const string dir_init );

#endif // _SIMRUN_H_INCLUDED
