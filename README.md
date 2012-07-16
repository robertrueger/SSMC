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
   susceptibility, spin-spin-correlation (nearest-neighbor models only),
   structure factor (dipolar models only), entropy and free energy

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

Both SSMC executables (ssmcsim and ssmcsa) read all parameters directly from the
command line when you execute them. SSMC is _ridiculously_ stupid when reading
arguments from the command line: It expects exactly the right number of
arguments in exactly the right order and will crash or behave undefined
otherwise! It is therefore absolutely crucial that your call to ssmcsim or
ssmcsa is 100% correct! If you just want to run a single simulation it will be
easier if you modify the respective shell scripts (ssmcsim.sh and ssmcsa.sh) to
suit your needs.

###ssmcsim arguments

    ssmcsim.sh [uint system_type] \
               [uint N] [bool periodic] \
               [char init] [uint drysweeps] \
               [uint bins] [uint binwidth] [uint intersweeps] \
               [bool run_plots] \
               [uint take_images] \
               [bool calc_autocorr] [bool calc_sscorr] \
               [uint smode_perbin] [uint smode_permcs] \
               [uint finite_size_correction] \
               [double J] [double g] \
               [double B] \
               [double T]

__uint systemtype__
Sets the type of the system (Hamiltonian + lattice) to be simulated.
== 1: one-dimensional Ising-Model with single-flip Metropolis-Algorithm
== 2: two-dimensional Ising-Model with single-flip Metropolis-Algorithm
== 3: two-dimensional field-free Ising-Model with Wolff-Algorithm
== 4: DEPRECATED, 2d Ising-Model with long range dipole interaction (Met.alg.)
== 5: proper implementation of 4
== 6: 2d Ising-Model with long range dipole interaction on a honeycomb lattice

__uint N__
Sets the size of the system. The way N is interpreted depends on which system
you choose to simulate: For one-dimensional models (systemtype 1) it is just the
total number of spins in the chain. For two-dimensional square models (2-5) it
is the length of a side, so that the system has N^2 spins in total. For the
honeycomb lattice model (6) it  is the size of the original system as defined in 
[arXiv:1207.1864 [cond-mat.stat-mech]](http://arxiv.org/abs/1207.1864).

__bool periodic__
Sets the boundary conditions.
== 0: open boundary conditions
== 1: periodic boundary conditions
Note: Not all systems support both open and periodic boundary conditions.
Open boundary conditions are only supported for the nearest-neighbor models
only, as it doesn't make much sense to simulate systems with long-range
interactions using open boundaries. The deprecated systemtype 4 does not support
periodic boundary conditions (and this is one of the reasons as to why it is
considered deprecated!).

__char init__
Sets the spin configuration used to initialize the equilibration period.
Options available for all models:
== r: random spin configuration
== u: all spins up
== d: all spins down
Model specific options:
== c: checkerboard configuration (all but systemtype 6)
== f: read initial state from file spin_init.log (only systemtype 6)
== [uint w]: stripes of width w (systemtype 2-6)
== s: stripe width is selected automatically to minimize energy (systemtype 2-6)
== e: first half up, second half down (only systemtype 1)

__uint drysweeps__
Sets the number of Monte Carlo steps used to equilibrate the system.

__uint bins__
Sets the number of bins used for measuring observables.

__uint binwidth__
Sets the number of observable measurements per bin.

__uint intersweeps__
Sets the number of Monte Carlo steps inbetween two observable measurements.
This is usually set to 1 since SSMC uses binning to deal with Markov chain
correlation, so there is no problem if two adjacent measurements are correlated.
A possible reason to increase this value is an _extremely_ long autocorrelation
time: If your AC time is very large you probably don't want to be measuring
after _every_ MCS, since the measurements will be virtually identical and you
waste time measuring very little new information. I would generally recommend
adjusting intersweeps so that you get a binwidth < 1000.

__bool run_plots__
Tells SSMC to automatically call (GNU/PyX)Plot to visualize its results.
== 0: don't plot (you can still execute the scripts later manually)
== 1: execute the plotting scripts after the simulation is complete

__uint take_images__
Tells SSMC to save png images of the spin configurations
== 0: don't take any pictures
== [uint n]: take a picture every n-th time the observables are measured

__bool calc_autocorr__
Tells SSMC to measure the magnetization's autocorrelation.
== 0: do nothing
== 1: measure it and calculate the autocorrelation time

__bool calc_sscorr__
Tells SSMC to measure spin-spin correlation functions.
== 0: do nothing
== 1: measure spin-spin correlation
Note: not implemented for systemtype 6!

__uint smode_perbin__
Sets if any 'special' model specific things should be done once every bin.
The idea behind this option and smode_permcs (see below) is to incorporate model
specific things (like measurements) into SMMC's general simulation procedure:
The structure factor is for example only interesting for some models (those with
dipolar interaction). It would therefore not be a good idea to have structure
factor measurements in the general system model defined in systemmodel.hpp.
Instead, there is only a general function special_perbin() that
is called once in every bin with smode_perbin as the input. For all models but
those with dipolar interaction special_perbin() does absolutely nothing, so it
doesn't matter what you set as smode_perbin! The new dipolar models (5,6) will
measure the structure factor if you set smode_perbin to anything but 0.

__uint smode_permcs__
Same as smode_perbin, but at the level of individual Monte Carlo steps.
Note: Not used at all in _any_ of the models, that come with SSMC!

__uint finite_size_correction__
Affects how magnetization measurements are interpreted.
== 0: use the magnetization as M = sum_i S_i
== 1: use the absolute value as the magnetization (M = abs(sum_i S_i))
== 2: determine automatically which interpretation to use (histogram method)
This setting is supposed to solve the following problem: Imagine simulating the
well known 2d Ising model with a ferromagnetic nearest-neighbor interaction at
temperatures slightly below the phase transition using the Metropolis algorithm.
Lets assume that just after equilibration most of the spins are 'up', so the
magnetization defined as M = sum_i S_i is positive. We know that below the 
critical temperature there should be a non-zero magnetization, so this is fine
for the moment. What can happen during the simulation is the following: With a
non-zero probability the Metropolis algorithm will reverse the the predominant
spin direction to 'down'. This will of course change the sign of the
magnetization. (This will only happen for finite systems, hence the name of the
option ...) If we average the systems magnetization over the course of the
entire simulation we get something close to zero because we have equal
contributions for both signs of the magnetization. This is not surprising, the
Hamiltonian is invariant under reversal of all spins anyway! What went wrong in
our simulation? Our problem is, that the magnetization M = sum_i S_i is _not_
the right observable to look at if one wants to see the Ising model's phase
transition. Analytically one should calculate the _spontaneous_ magnetization
defined as M_s = lim_{B->o^+} lim_{n->infty} M, but how could we do that in a
Monte Carlo simulation? One solution might be to always use the absolute value
of the magnetization. But there is a problem associated with that: For
temperatures above the critical temperature the average will be >0 because we
are only averaging over positive contributions. And >0 is wrong ...
To my knowledge the best way out of this mess is the following: Run the
simulation for some time after the equilibration and record a magnetization
histogram. If there are more states with a magnetization per spin between -1 to
-0.5 and 0.5 to 1 than there are between -0.5 and 0.5, use the absolute value of
the magnetization. On other words the system has to favor states with a high
absolute magnetization over those with low magnetization. The argument here is
that if states with M approx 0 are essentially forbidden, then the real system
(that has a time dependence governed by Hamilton's equations of motion) would
not reverse its magnetization: It's path through configuration space would no
longer be ergodic! I'm not sure if there are other/better methods to measure
spontaneous magnetizations in a Monte Carlo simulation. SSMC's method works
really well for the Ising model (especially with Wolff-algorithm!), but if
anyone has any better ideas, I would be happy to hear them!

__double J__
Sets the nearest-neighbor coupling constant J.

__double g__
Sets the strength of the dipolar intaraction.
Note: Is ignored for all models without dipolar interaction.

__double B__
Sets the external magnetic field B.

__double T__
Sets the temperature T.

###ssmcsa arguments

    ssmcsa.sh  [uint system_type] \
               [uint N] [bool periodic] \
               [char init] \
               [double T_start] [double T_end] \
               [ulint t_end] [char cooling_schedule] \
               [uint t_boost] \
               [bool run_plots] [uint take_images] \
               [double J] [double g] \
               [double B]

__double T_start__
Sets the temperature to start the simulated annealing at.

__double T_end__
Sets the temperature to anneal the system to.

__ulint t_end__
The length of the annealing in Monte Carlo steps.

__char cooling_schedule__
Determines how the temperature is decreased.
== l: linear decrease
== p: parabolic decrease (faster then linear in the beginning, slower at end)

__uint t_boost__
Sets the number of Monte Carlo steps performed before the temperature is
recalculated.


##3. Expanding SSMC

SSMC was designed with some expandability in mind in the sense that the it can
easily be adapted to other classical spin system models (geometries,
dimensionality and even non-Ising spins!). This is possible because all the data
analysis and decorrelation procedure is entirely independent of the actual model
that is simulated. In SSMC this is reflected in the fact that there is a mostly
pure virtual class SystemModel defined in systemmodel.hpp that defines the
interface the main SSMC code uses to communicate with the model systems. It is
not entirely pure but implements a few things that you definitely need for every
model, for example it takes are of initializing a random number generator for
you.

If you want to implement your own model for use with SSMC you have to do two
things:

1. Write a model class that inherits from SystemModel and implements all of its
   pure functions.

2. Add the call to the model's constructor to simrun.cpp and sarun.cpp. (It's
   rather self explanatory, just have a look at how the models that come with
   SSMC do it!)

        // ___ ADD CUSTOM SYSTEM MODELS HERE ___
        } else if ( par.system_type == MY_MODEL_NUMBER ) {
          model = new MyNewModel( MY_MODEL_PARAMETER_LIST );

That's it! You should be able to use SSMC to simulate your own classical spin
system models now!


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
