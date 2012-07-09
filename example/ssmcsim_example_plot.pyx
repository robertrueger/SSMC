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


set terminal pdf

# ----- PAGE LAYOUT -----
scale=2
plotwidth = scale*8
plotheight = plotwidth/2
hspace = scale*2.5
vspace = scale*0.45
set width plotwidth
set size ratio plotheight/plotwidth


# ----- MACROS & STUFF -----
style_points = "with points pointsize 0.666"
style_errorbars = "with errorbars color black"


# ----- PLOT -----
set xrange [0:5]
#set xtics (0, 1, 2, "$T_c$" 2.269, 3, 4, 5)
#set mxtics 0.1
set grid x y

set multiplot
set output 'thermodynamics.pdf'
set nodisplay

set ylabel 'internal energy~$u(T)$'
set origin 0*(plotwidth + hspace), 5*(plotheight + vspace)
plot './results.dat' using 5:6 @style_points notitle, \
	 './results.dat' using 5:6:7 @style_errorbars notitle

set ylabel 'specific heat capacity~$c(T)$'
set origin 0*(plotwidth + hspace), 4*(plotheight + vspace)
plot './results.dat' using 5:8 @style_points notitle, \
	 './results.dat' using 5:8:9 @style_errorbars notitle

set ylabel 'entropy~$s(T)$'
set origin 0*(plotwidth + hspace), 3*(plotheight + vspace)
plot './results2.dat' using 1:2 with lines notitle


set ylabel 'free energy~$f(T)$'
set origin 0*(plotwidth + hspace), 2*(plotheight + vspace)
plot './results2.dat' using 1:3 with lines notitle

set ylabel 'magnetization~$m(T)$'
set origin 0*(plotwidth + hspace), 1*(plotheight + vspace)
plot './results.dat' using 5:10 @style_points notitle, \
	 './results.dat' using 5:10:11 @style_errorbars notitle

set xlabel 'temperature~$T$'
set ylabel 'susceptibility~$\chi(T)$'
set origin 0*(plotwidth + hspace), 0*(plotheight + vspace)
plot './results.dat' using 5:12 @style_points notitle, \
	 './results.dat' using 5:12:13 @style_errorbars notitle


set display
refresh
clear
