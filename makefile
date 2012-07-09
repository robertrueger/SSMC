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


CXX				= g++
CXXFLAGS	= -Wall -march=native -O3 -flto -fuse-linker-plugin
DEFINES		=
LDFLAGS		= -lgsl -lgslcblas `libpng-config --ldflags`

EXECTUABLES	= ssmcsim ssmcsa

OBJECTS_SIM	= ssmcsim.o simrun.o
OBJECTS_SA	= ssmcsa.o sarun.o sa_coolingschedules.o
OBJECTS_ALL	= systemmodel.o isingspin.o utils.o utils_vec2.o \
 							model_dip_sqr.o model_dip_hc.o model_ising1dmet.o \
							model_ising2dsqrmet.o model_ising2dsqrffwolff.o \
							model_ising2dsqrdipolemet.o

HEADERS_MODELS = model_dip_sqr.hpp model_dip_hc.hpp model_ising1dmet.hpp \
								 model_ising2dsqrmet.hpp model_ising2dsqrffwolff.hpp \
								 model_ising2dsqrdipolemet.hpp


ssmc all : ssmcsim ssmcsa

ssmcsim : $(OBJECTS_SIM) $(OBJECTS_ALL)
	$(CXX) $(CXXFLAGS) $(OBJECTS_SIM) $(OBJECTS_ALL) $(LDFLAGS) -o ssmcsim

ssmcsa : $(OBJECTS_SA) $(OBJECTS_ALL)
	$(CXX) $(CXXFLAGS) $(OBJECTS_SA) $(OBJECTS_ALL) $(LDFLAGS) -o ssmcsa

clean:
	rm -f $(EXECTUABLES) $(OBJECTS_SIM) $(OBJECTS_SA) $(OBJECTS_ALL)


# ----- OBJECT FILES -----

ssmcsim.o : ssmcsim.cpp simrun.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c ssmcsim.cpp -o ssmcsim.o

ssmcsa.o : ssmcsa.cpp sarun.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c ssmcsa.cpp -o ssmcsa.o

simrun.o : simrun.hpp simrun.cpp systemmodel.hpp utils.hpp sim_datastruct.hpp \
					 $(HEADERS_MODELS)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c simrun.cpp -o simrun.o

sarun.o : sarun.hpp sarun.cpp systemmodel.hpp utils.hpp sa_datastruct.hpp \
					sa_coolingschedules.hpp $(HEADERS_MODELS)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c sarun.cpp -o sarun.o

systemmodel.o : systemmodel.hpp systemmodel.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c systemmodel.cpp -o systemmodel.o

isingspin.o : isingspin.hpp isingspin.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c isingspin.cpp -o isingspin.o

utils.o : utils.hpp utils.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c utils.cpp -o utils.o

utils_vec2.o : utils_vec2.hpp utils_vec2.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c utils_vec2.cpp -o utils_vec2.o

sa_coolingschedules.o : sa_coolingschedules.hpp sa_coolingschedules.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c sa_coolingschedules.cpp \
																-o sa_coolingschedules.o


# ----- SPIN SYSTEM MODEL OBJECT FILES -----

model_ising1dmet.o : model_ising1dmet.hpp model_ising1dmet.cpp systemmodel.hpp \
										 isingspin.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_ising1dmet.cpp -o model_ising1dmet.o

model_ising2dsqrmet.o : model_ising2dsqrmet.hpp model_ising2dsqrmet.cpp \
												systemmodel.hpp isingspin.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_ising2dsqrmet.cpp \
														    -o model_ising2dsqrmet.o

model_ising2dsqrffwolff.o : model_ising2dsqrffwolff.hpp \
														model_ising2dsqrffwolff.cpp model_ising2dsqrmet.hpp \
														systemmodel.hpp isingspin.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_ising2dsqrffwolff.cpp \
																-o model_ising2dsqrffwolff.o

model_ising2dsqrdipolemet.o : model_ising2dsqrdipolemet.hpp \
															model_ising2dsqrdipolemet.cpp \
															model_ising2dsqrmet.hpp systemmodel.hpp \
															isingspin.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_ising2dsqrdipolemet.cpp \
															  -o model_ising2dsqrdipolemet.o

model_dip_sqr.o : model_dip_sqr.hpp model_dip_sqr.cpp systemmodel.hpp \
									isingspin.hpp utils.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_dip_sqr.cpp -o model_dip_sqr.o

model_dip_hc.o : model_dip_hc.hpp model_dip_hc.cpp systemmodel.hpp \
								 isingspin.hpp utils.hpp utils_vec2.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c model_dip_hc.cpp -o model_dip_hc.o
