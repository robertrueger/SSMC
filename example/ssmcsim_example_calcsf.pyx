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


interpolate linear u() 'results.dat' using 5:6
interpolate akima  c() 'results.dat' using 5:8
s(T) = int_dTprime('c(Tprime)/Tprime', 0, T)
set samples 501
set xrange [0:5]
set output 'results2.dat'
tabulate s(x):u(x)-x*s(x) using 1:2:3
