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


set -e

# tidy up ...
rm -rf $(ls . | grep -v ssmcsim_example)

# create links ...
ln -s ../ssmcsim	ssmcsim 	|| true
ln -s ../ssmcsim.sh	ssmcsim.sh 	|| true

# set fixed parameters ...
system_type=3
N=50
periodic=1
init=u
drysweeps=1000
bins=100
binwidth=100
intersweeps=1
J=1
g=0
B=0
#T=variable
run_plots=0
images=1
calc_autocorr=1
calc_sscorr=1
finite_size_correction=2
smode_perbin=0
smode_permcs=0

# run the simulations
for T in 0 0.5 1 1.2 1.4 1.6 1.7 1.8 1.9 \
    		 2.0 2.1 2.15 2.175 2.2 2.225 2.25 2.275 \
         2.3 2.325 2.35 2.375 2.4 2.425 2.45 2.5 \
		     2.6 2.7 2.8 2.9 3.0 3.2 3.4 3.6 4 4.5 5 
do
  ./ssmcsim.sh $system_type $N $periodic $init $drysweeps $bins $binwidth \
               $intersweeps $run_plots $images $calc_autocorr $calc_sscorr \
               $smode_perbin $smode_permcs $finite_size_correction $J $g $B $T
  echo
done


# collect results and make run_plots
rm results.dat results.raw || true
for dir in `ls -l | grep ^d | awk '{print $9}'`; do

  cat $dir/results.dat >> ./results.raw
  cd $dir

  [ -e ac_plot.png ]  	|| (gnuplot ac_plot.gnu  	|| true)
  [ -e run_plot.png ] 	|| (gnuplot run_plot.gnu 	|| true)
  [ -e sscorr_plot.png ] 	|| (gnuplot sscorr_plot.gnu || true)
  [ -e mhisto_plot.pdf ] 	|| (pyxplot histo_plot.pyx 	|| true)

  if [ -d images ];
  then
    cd images
    mencoder mf://*.png -mf fps=30:type=png -ovc copy -oac copy -o ../video.avi
    cd ..
    tar -caf images.tar.xz images
    rm -r images
  fi

  cd ..
done

# remove infs and nans from results.dat
sed -e 's/inf/0.0000000e+00/g' results.raw \
| sed -e 's/nan/0.0000000e+00/g' > results.dat

# plot results.dat with pyxplot
echo -e "[latex]\nPreamble = \usepackage[ngerman]{babel}" >> .pyxplotrc
pyxplot ssmcsim_example_calcsf.pyx
pyxplot ssmcsim_example_plot.pyx || true
