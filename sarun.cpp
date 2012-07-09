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


#include "sarun.hpp"


// ----- SARUN: PERFORMS A SINGLE SIMULATED ANNEALING -----

sa_results sarun( const sa_parameters par, const string dir_init )
{
  // ----- PREPARE SIMULATED ANNEALING  -----

  // assume something went wrong until we are sure it didn't
  sa_results res;
  res.success = false;

  // make a folder to work in
  stringstream tmp;
  tmp << setfill( '0' );
  tmp << "./" << dir_init << '/';
  const string dir = tmp.str();
  tmp.str( "" );
  // folder for images of the model
  if ( par.take_images != 0 ) {
    tmp << dir << "images/" ;
    const string image_dir = tmp.str();
    tmp.str( "" );
    if ( system( ( "mkdir " + image_dir ).c_str() ) != 0 ) {
      cout << "ERROR while making the image directory " << dir << endl;
      return res;
    }
  }

  // start a logfile
  ofstream run_log( ( dir + "run.log" ).c_str() );
  if ( !run_log.is_open() ) {
    cout << "ERROR while opening the run log file in " << dir << endl;
    return res;
  }
  run_log.precision( numeric_limits<double>::digits10 + 1 );

  // write the simulations parameters to the logfile
  run_log << "Simulated annealing running in " << dir
          << endl << endl
          << "--- PARAMETERS ---" << endl
          << "system          = " << par.system_type << endl
          << "N               = " << par.N << endl
          << "periodic        = " << par.periodic << endl
          << "init            = " << par.init << endl
          << "T               = " << par.T_start << "->" << par.T_end << endl
          << "t_end           = " << par.t_end << endl
          << "coolingschedule = " << par.cooling_schedule << endl
          << "J               = " << par.J << endl
          << "g               = " << par.g << endl
          << "B               = " << par.B << endl << endl;
  run_log.flush();


  // ----- RUN SIMULATED ANNEALING -----

  time_t rawtime;
  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> creating the system\n\n";
  run_log.flush();

  // create a new model
  SystemModel* model;
  if ( par.system_type == 1 ) {
    model = new IsingModel1d( par.N, par.periodic, par.J, par.B,
                              par.T_start, 0, dir );
  } else if ( par.system_type == 2 ) {
    model = new IsingModel2d( par.N, par.periodic, par.J, par.B,
                              par.T_start, 0, dir );
  } else if ( par.system_type == 3 ) {
    model = new IsingModel2dWolff( par.N, par.periodic, par.J,
                                   par.T_start, 0, dir );
  } else if ( par.system_type == 4 ) {
    model = new IsingModel2dDipole( par.N, par.periodic, par.J, par.g, par.B,
                                    par.T_start, dir );
  } else if ( par.system_type == 5 ) {
    model = new Ising2dDipSqr( par.N, par.periodic, par.J, par.g, par.B,
                               par.T_start, 0, dir );
  } else if ( par.system_type == 6 ) {
    model = new Ising2dDipHC( par.N, par.J, par.g, par.B,
                              par.T_start, 0, dir );
  } else {
    cout << "ERROR creating the model system in " << dir << endl;
    return res;
  }
  if ( model->prepare( par.init ) == false ) {
    cout << "ERROR preparing the models spins in " << dir << endl;
    delete model;
    return res;
  }

  double ( *Toft )( double const&, double const&, unsigned long int const&,
                    unsigned long int const& );
  // select a cooling schedule
  if ( par.cooling_schedule == 'l' ) {
    Toft = &linear_cooling;
  } else if ( par.cooling_schedule == 'p' ) {
    Toft = &parabolic_cooling;
  } else {
    cout << "ERROR selecting the cooling schedule in " << dir << endl;
    delete model;
    return res;
  }

  // open measurement logfiles
  ofstream h_log( ( dir + "h.log" ).c_str() );
  ofstream m_log( ( dir + "m.log" ).c_str() );
  if ( !( h_log.is_open() && m_log.is_open() ) ) {
    cout << "ERROR while opening measurement log files in " << dir << endl;
    delete model;
    return res;
  }

  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> simulation started\n\n";
  run_log.flush();

  // sample loop
  while ( model->t() <= par.t_end ) {
    double T_now = Toft( par.T_start, par.T_end, par.t_end, model->t() );
    model->set_T( T_now );
    // write this sample's properties to the logfile
    h_log << model->t() << ' ' << T_now << ' ' << model->h() << endl;
    m_log << model->t() << ' ' << T_now << ' ' << model->m() << endl;
    // make an image of the system
    if ( ( par.take_images != 0 ) && ( model->t() % par.take_images == 0 ) ) {
      tmp << dir << "images/"	<< setw( 9 ) << model->t() << ".png";
      const string image_file = tmp.str();
      tmp.str( "" );
      model->get_image().write( image_file );
    }
    // do t_boost monte carlo steps
    for ( unsigned int i = 0; i < par.t_boost; ++i ) {
      model->mcstep();
    }
  }

  // all measurements done ... let's tidy things up
  delete model;
  h_log.close();
  m_log.close();

  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> simulated annealing finished";
  run_log.flush();


  // ----- PLOTTING -----

  if ( par.run_plot ) {
    time( &rawtime );
    run_log << ctime( &rawtime ) << "-> starting to plot the results\n\n";
    run_log.flush();
  }

  // run_plot
  ofstream run_plot( ( dir + "run_plot.gnu" ).c_str() );
  if ( !run_plot.is_open() ) {
    cout << "ERROR while opening run_plot.gnu in " << dir << endl;
    return res;
  }
  tmp << "system_type = " << par.system_type
      << ", N = " << par.N << ", periodic = " << par.periodic
      << ", init = " << par.init << ", T = " << par.T_start << "->" << par.T_end
      << ", t_end	= " << par.t_end << "cooling_schedule = "
                                                         << par.cooling_schedule
      << ", J = " << par.J << ", g = " << par.g << ", B = " << par.B;
  int plot_size = ( par.t_end < 1600 ) ? 
                  1600 : min( ( unsigned long int )10000, par.t_end );
  run_plot <<
  "  set terminal pngcairo size " << plot_size << ",600 \n\
	set output 'run_plot.png' \n\n\
	set multiplot layout 2,1 title '" << tmp.str() << "'\n\
	set grid x y \n\
	set mxtics 10 \n\
	set key top left \n\
	plot 'm.log' using 1:3 with lines ls 2 title 'magnetization per spin' \n\
	plot 'h.log' using 1:3 with lines ls 1 title 'energy per spin' \n";
  tmp.str( "" );
  run_plot.close();

  if ( par.run_plot && par.t_end <= 1e5 ) {
    if ( system( ( "cd " + dir + " ; gnuplot run_plot.gnu" ).c_str() ) != 0 ) {
      cout << "ERROR while running gnuplot run_plot.gnu in " << dir << endl;
      return res;
    }
  }

  // everything is fine ... return the results!
  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> everything finished!\n\n";
  run_log.close();
  res.success = true;
  return res;
}
