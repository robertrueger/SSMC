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


// ----- SIMULATED ANNEALING DATA STRUCTURES -----

struct sa_parameters {
  // ----- INPUT STRUCT FOR SARUN() -----

  unsigned int system_type;

  unsigned int N;
  bool periodic;

  char init;

  double T_start, T_end;

  unsigned long int t_end;
  char cooling_schedule;

  unsigned int t_boost;

  bool run_plot;
  unsigned int take_images;

  double J;
  double g;
  double B;
};


struct sa_results {
  // ----- OUTPUT STRUCT FOR SARUN() -----

  bool success;
};


struct simann {
  // ----- COMBINED SIMULATED ANNEALING PARAMS AND RESULTS  -----

  sa_parameters par;
  sa_results res;
};
