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


// ----- SIMULATION DATA STRUCTURES -----

struct sim_parameters {
  // ----- INPUT STRUCT FOR SIMRUN() -----

  unsigned int system_type;

  unsigned int N;
  bool periodic;

  char init;
  unsigned int drysweeps;

  unsigned int bins;
  unsigned int binwidth;
  unsigned int intersweeps;

  bool run_plot;
  unsigned int take_images;

  bool calc_autocorr;
  bool calc_sscorr;

  unsigned int smode_perbin;
  unsigned int smode_permcs;

  unsigned int use_fsize_correction;

  double J;
  double g;
  double B;
  double T;
};


struct sim_results {
  // ----- OUTPUT STRUCT FOR SIMRUN() -----

  bool success;

  double h;
  double sigma3_h;

  double c;
  double sigma3_c;

  double m;
  double sigma3_m;

  double x;
  double sigma3_x;

  double tau;
};


struct simulation {
  // ----- COMBINED SIMULATION PARAMS AND RESULTS -----

  sim_parameters par;
  sim_results res;
};


struct bin_results {
  // ----- RESULTS OF A SINGLE BIN -----

  double h, h2;
  double m, m2;
};
