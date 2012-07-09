#! /bin/bash

# Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
#
# This file is part of SSMC.
#
# SSMC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SSMC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SSMC.  If not, see <http://www.gnu.org/licenses/>.



# ========================= CLI MODE =========================
#         (read parameters from the command line)

if [ $# -eq 14 ]; then
	echo "ssmcsa called in cli mode ..."
	time ./ssmcsa $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14}
	exit
fi



# ======================= SCRIPT MODE ========================
# 				(read parameters from this script)
echo "ssmcssa called in script mode ..."

# ----------------- SIMULATION PARAMETERS --------------------

# type of system to be simulated?
system_type=2
# -- 1: one-dimensional Ising-Model with single-flip Metropolis-Algorithm
# -- 2: two-dimensional Ising-Model with single-flip Metropolis-Algorithm
# -- 3: two-dimensional field-free Ising-Model with Wolff-Algorithm
# -- 4: 2d Ising-Model with long range dipole interaction (Metropolis-Alg.)
# -- 5: like 4 but improved
# -- 6: 2d Ising-Model with dip-dip-int on a honeycomb lattice

# size of the system?
N=50

# periodic boundary conditions?
periodic=1

# initial conditions?
init=r
# -- r: initialize spins randomly
# -- u: all spins up at t=0
# -- d: all spins down at t=0
# -- e: 50:50 up/down (for B=0 and very low temperatures)
# -- c: checkerboard (2d)
# -- s: enery minimizing stripes
# -- [uint]: stripes of manually defined width

# temperature rage to be sweeped
T_start=10
T_end=0

# cooling time and schedule
t_end=100000
cooling_schedule=l

# do measurements every ... steps?
t_boost=10

# strength of nearest-neighbour-coupling?
J=1

# strength of dipole-dipole-coupling (only relevant for system_type 4 and 5)
g=0

# external magnetic field
B=0



# ------------------ DATA ANALYSIS OPTIONS --------------------

# do you want SSMCSA to call the plotting scripts automatically? 
# (requires gnuplot and pyxplot)
run_plots=0

# do you want to take images of the system?
images=10
# -- 0: dont take images
# -- n: take an image every n mcsteps



# ----------------------- EXECUTION --------------------------

time ./ssmcsa	$system_type $N $periodic \
              $init \
              $T_start $T_end \
              $t_end $cooling_schedule \
              $t_boost \
              $run_plots $images \
              $J $g $B
