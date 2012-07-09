#SSMC

SSMC (Spin System Monte Carlo) is a free Monte Carlo simulation code for
classical spin systems like the Ising model. Here are some of its features:

0. __Models:__ nearest-neighbor Ising model in one and two dimensions, dipolar
 Ising model on two-dimensional square and honeycomb lattices, all models with
user defined sizes, support for open and periodic boundary conditions

1. __Algorithms:__ Metropolis, Wolff (only for 2d nearest-neighbor Ising model
without external magnetic field), simulated annealing with various cooling
schedules

2. __Statistics:__ binning, bootstrapping (for the calculation of susceptibility
 errorbars), calculation of the magnetization's autocorrelation

3. __Observables:__ internal energy, magnetization, specific heat, magnetic
susceptibility, spin-spin-correlation (nearest-neighbor models only), structure
factor (dipolar models only), entropy and free energy

4. __Visualization:__ energy and magnetization as a function of MC time, images
of the system's spin configurations (with colored Wolff-clusters), energy and
magnetization histograms

5. __Expandability:__ generic spin system interface (write a new model code
which implements an interface to SSMC; only trivial changes to the main code
are needed)

I originally started to write the SSMC code as a part of my bachlor thesis. It
was never intended to be used by anybody but myself and while there are some
source code comments, there is no formal documentation. I nevertheless wanted to
release it as free software in the hope that a least parts of it will be useful.
I hope this readme gives you a head start in using, understanding and modifying
SSMC, even though the documentation is lacking.

Feel free to contact me with any questions (or bugs!) that come up!


##0. Building

__Note:__ This section assumes that you have the GNU software development tools
like GCC, make and ld already installed on your machine.

Before building SSMC make sure that you have installed its dependencies.
SSMC distinguishes required and optional dependencies: While you don't need the
optional dependencies to build SSMC, you will not be able to automatically
visualize your results without them.

__required:__

* [GNU scientific library](http://www.gnu.org/software/gsl/)
* [png++](http://savannah.nongnu.org/projects/pngpp/)

__optional:__

* [gnuplot](http://www.gnuplot.info)
* [pyxplot](http://www.pyxplot.org.uk)
* [mencoder](http://www.mplayerhq.hu)

Once you have the dependencies installed, building SSMC is really simple and
straightforward:

    git clone git://github.com/robertrueger/SSMC.git
    cd SSMC
    make

You should now have two executables: ssmcsim and ssmcsa


##1. Quick start guide

ssmcsim performs a simulation at constant temperature, while ssmcsa allows you
to perform simulated annealings. Both are command line programs, that read all
simulation parameters directly as arguments from the command line when they are
called. As there are a lot of available options, the program call is quite
complicated and long. See below for a full list of all command line arguments.

There are two scripts, that make using SSMC easier: ssmcsim.sh and ssmcsa.sh.
Instead of writing all arguments to ssmcsim/ssmcsa on the command line, you can
just modify the respective *.sh script and run it. Note that you can use the
*.sh scripts in the same way as the real executables and call them with
arguments on the command line. I therefore suggest, that you _never_ call the
executables directly, but always go through the respective *.sh script.

Go ahead and run ssmcsim.sh now! If you have not made any changes to it, it
should run a simulation of a 100x100 Ising model close to the critical
temperature using the Wolff algorithm. This should not take longer than a
minute. All results are written to an output folder whose name is determined by
the simulations parameters. Your output folder should look somewhat like this:

(output: http://itp.uni-frankfurt.de/~rueger/SSMC/output_example0/)

    sim3_N100p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T02269
    ├── ac.log
    ├── ac_plot.gnu
    ├── ac_plot.png
    ├── hhisto_plot.pdf
    ├── histo_plot.pyx
    ├── h.log
    ├── images
    │   ├── 000000000.png
    │   ├── 000000010.png
    │   ├── 000000020.png
    │   ├── 000000030.png
    │   ├── ...
    │   └── 000009990.png
    ├── mhisto_plot.pdf
    ├── m.log
    ├── results.dat
    ├── run.log
    ├── run_plot.gnu
    ├── run_plot.png
    ├── sscorr.log
    ├── sscorr_plot.gnu
    └── sscorr_plot.png

More complex use of SSMC is demonstrated in the example folder. Have a look at
the ssmcsim_example_runme.sh file and run it! It will perform many simulations
of the 2d Ising model at different temperatures using the Wolff algorithm and
combines the results to calculate entropy and free energy. Output on
stdout/stderr is rather verbose, but within half an hour (at most!) you should
have 36 simulation output folders and most notably a plot of the Ising model's
famous thermodynamics!

(output: http://itp.uni-frankfurt.de/~rueger/SSMC/output_example1/)

    example
    ├── results2.dat
    ├── results.dat
    ├── results.raw
    ├── sim3_N50p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T00000
    │   ├── ac.log
    │   ├── ac_plot.gnu
    │   ├── ac_plot.png
    │   ├── hhisto_plot.pdf
    │   ├── histo_plot.pyx
    │   ├── h.log
    │   ├── images.tar.xz
    │   ├── mhisto_plot.pdf
    │   ├── m.log
    │   ├── results.dat
    │   ├── run.log
    │   ├── run_plot.gnu
    │   ├── run_plot.png
    │   ├── sscorr.log
    │   ├── sscorr_plot.gnu
    │   ├── sscorr_plot.png
    │   └── video.avi
    ├── sim3_N50p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T00500
    │   ├── ac.log
    │   ├── ...
    │   └── video.avi
    ├── sim3_N50p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T01000
    ├── sim3_N50p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T01200
    ├── ...
    ├── sim3_N50p_u1000_100-100-1_fsc2_smode-0-0_J01000_g00000_B00000_T05000
    └── thermodynamics.pdf


##2. Command line arguments

... coming soon ...


##3. Expanding SSMC

... coming soon ...


##License

Copyright (c) 2012, Robert Rüger (rueger@itp.uni-frankfurt.de)

SSMC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SSMC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SSMC.  If not, see <http://www.gnu.org/licenses/>.
