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

#include "simrun.hpp"


// ----- SIMRUN: PERFORMS A SINGLE SIMULATION -----

sim_results simrun( const sim_parameters par, const string dir_init )
{
  // ----- PREPARE SIMULATION -----

  // assume something went wrong until we are sure it didn't
  sim_results res;
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
  run_log << "Simulation running in " << dir
          << endl << endl
          << "--- PARAMETERS ---" << endl
          << "system          = " << par.system_type << endl
          << "N               = " << par.N << endl
          << "periodic        = " << par.periodic << endl
          << "init            = " << par.init << endl
          << "drysweeps       = " << par.drysweeps << endl
          << "bins            = " << par.bins << endl
          << "binwidth        = " << par.binwidth << endl
          << "intersweeps     = " << par.intersweeps << endl
          << "smode_perbin    = " << par.smode_perbin << endl
          << "smode_permcs    = " << par.smode_permcs << endl
          << "J               = " << par.J << endl
          << "g               = " << par.g << endl
          << "B               = " << par.B << endl
          << "T               = " << par.T << endl << endl;
  run_log.flush();


  // ----- RUN SIMULATION -----

  time_t rawtime;
  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> creating the system\n\n";
  run_log.flush();

  // create a new model
  SystemModel* model;
  if ( par.system_type == 1 ) {
    model = new IsingModel1d( par.N, par.periodic, par.J, par.B,
                              par.T, par.use_fsize_correction, dir );
  } else if ( par.system_type == 2 ) {
    model = new IsingModel2d( par.N, par.periodic, par.J, par.B, par.T,
                              par.use_fsize_correction, dir );
  } else if ( par.system_type == 3 ) {
    model = new IsingModel2dWolff( par.N, par.periodic, par.J, par.T,
                                   par.use_fsize_correction, dir );
  } else if ( par.system_type == 4 ) {
    model = new IsingModel2dDipole( par.N, par.periodic, par.J, par.g, par.B,
                                    par.T, dir );
  } else if ( par.system_type == 5 ) {
    model = new Ising2dDipSqr( par.N, par.periodic, par.J, par.g, par.B,
                               par.T, par.use_fsize_correction, dir );
  } else if ( par.system_type == 6 ) {
    model = new Ising2dDipHC( par.N, par.J, par.g, par.B,
                              par.T, par.use_fsize_correction, dir );
  } else {
    cout << "ERROR creating the model system in " << dir << endl;
    return res;
  }
  if ( model->prepare( par.init ) == false ) {
    cout << "ERROR preparing the models spins in " << dir << endl;
    delete model;
    return res;
  }
  unsigned int spin_count = model->spin_count();

  // open measurement logfiles
  ofstream h_log( ( dir + "h.log" ).c_str() );
  ofstream m_log( ( dir + "m.log" ).c_str() );
  if ( !( h_log.is_open() && m_log.is_open() ) ) {
    cout << "ERROR while opening measurement log files in " << dir << endl;
    delete model;
    return res;
  }

  // initialize binning array
  vector <bin_results> binres;
  // initialize the magnetization data log (needed for autocorr calc only)
  vector <float> m_memlog;
  // initialize the spin-spin correlation
  vector <double> ss_corr;
  try {
    // try to allocate enough memory ...
    binres.reserve( par.bins );
    unsigned int ss_corr_size = model->ss_corr().size();
    ss_corr.reserve( ss_corr_size );
    for ( unsigned int i = 0; i < ss_corr_size; i++ ) {
      ss_corr.push_back( 0 );
    }
    if ( par.calc_autocorr ) {
      m_memlog.reserve( par.bins * par.binwidth );
    }
  } catch ( bad_alloc ) {
    cout << "ERROR while allocating memory in " << dir << endl;
    delete model;
    return res;
  }
  run_log
      << 1 + ( ( sizeof( m_memlog ) + m_memlog.capacity() * sizeof( float )
         + sizeof( binres ) + binres.capacity() * sizeof( bin_results ) ) / 1024 )
      << "KiB of memory reserved\n\n";


  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> relaxing the system\n\n";
  run_log.flush();

  // perform dry runs to reach thermal equilibrium
  model->mcstep_dry( par.drysweeps );

  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> simulation started\n\n";
  run_log.flush();


  // binning loop
  for ( unsigned int bin = 0; bin < par.bins; bin++ ) {
    //double startTime = current_time();

    // initialize variables to measure the systems properties
    double h = 0, h2 = 0;
    double m = 0, m2 = 0;
    // sample loop
    for ( unsigned int sample = 0; sample < par.binwidth; sample++ ) {
      double thissample_h = model->h();
      double thissample_m = model->m();
      // write this sample's properties to the logfile
      h_log << model->t() << ' ' << thissample_h << endl;
      m_log << model->t() << ' ' << thissample_m << endl;
      if ( par.calc_autocorr ) {
        m_memlog.push_back( float( thissample_m ) );
      }
      // remember the sample's properties to calculate their mean value
      h 	+= thissample_h;
      h2 	+= thissample_h * thissample_h;
      m 	+= thissample_m;
      m2 	+= thissample_m * thissample_m;
      // make an image of the system
      if ( ( par.take_images != 0 ) &&
           ( ( bin * par.binwidth + sample ) % par.take_images == 0 ) ) {
        tmp << dir << "images/"	<< setw( 9 ) << model->t() << ".png";
        const string image_file = tmp.str();
        tmp.str( "" );
        model->get_image().write( image_file );
      }
      // flip the spins!
      for ( unsigned int step = 0; step < par.intersweeps; step++ ) {
        model->special_permcs( par.smode_permcs );
        model->mcstep();
      }
    }

    if ( par.calc_sscorr ) {
      // spin-spin correlation calculation
      vector <double> ss_corr_thisbin = model->ss_corr();
      for ( unsigned int i = 0; i < par.N; i++ ) {
        ss_corr[i] += ss_corr_thisbin[i];
      }
    }

    // invoke the systems special function
    model->special_perbin( par.smode_perbin );

    // calculate mean
    h = h / par.binwidth;
    h2 = h2 / par.binwidth;
    m = m / par.binwidth;
    m2 = m2 / par.binwidth;

    // write the bin's results to binres
    bin_results this_binres;
    this_binres.h			= h;
    this_binres.h2		= h2;
    this_binres.m 		= m;
    this_binres.m2		= m2;
    binres.push_back( this_binres );

    //cout << "Bin: " << current_time() - startTime <<  "ms" << endl;
  }

  // all measurements done ... let's tidy things up
  delete model;
  h_log.close();
  m_log.close();

  // calculate simulation results from the individual bins
  double h  = 0, 	sigma3_h  = 0;
  double h2 = 0, 	sigma3_h2 = 0;
  double m  = 0, 	sigma3_m  = 0;
  double m2 = 0, 	sigma3_m2 = 0;
  // average values
  for ( unsigned int bin = 0; bin < par.bins; bin++ ) {
    h  += binres[bin].h  / par.bins;
    h2 += binres[bin].h2 / par.bins;
    m  += binres[bin].m  / par.bins;
    m2 += binres[bin].m2 / par.bins;
  }
  if ( par.calc_sscorr ) {
    for ( unsigned int i = 0; i < par.N; i++ ) {
      ss_corr[i] /= par.bins;
    }
  }
  // calculate susceptibilities
  double c = spin_count / par.T / par.T * ( h2 - h * h );
  double x = spin_count / par.T * ( m2 - m * m );

  // calculate variance of the results from the bins first
  for ( unsigned int bin = 0; bin < par.bins; bin++ ) {
    sigma3_h 	+= pow( ( binres[bin].h	 - h ),  2 ) / par.bins;
    sigma3_h2	+= pow( ( binres[bin].h2 - h2 ), 2 ) / par.bins;
    sigma3_m 	+= pow( ( binres[bin].m  - m ),  2 ) / par.bins;
    sigma3_m2	+= pow( ( binres[bin].m2 - m2 ), 2 ) / par.bins;
  }
  // use variances to calculate the error of the average
  sigma3_h 	= 3 * sqrt( sigma3_h  / par.bins );
  sigma3_h2	= 3 * sqrt( sigma3_h2 / par.bins );
  sigma3_m 	= 3 * sqrt( sigma3_m  / par.bins );
  sigma3_m2	= 3 * sqrt( sigma3_m2 / par.bins );

  // calculate errors of the susceptibilities (bootstrapping)
  double sigma3_c = 0, sigma3_x = 0;
  gsl_rng* rng;  // make a new random number generator ...
  rng = gsl_rng_alloc( gsl_rng_mt19937 );
  gsl_rng_set( rng, rand() );
  for ( unsigned int bsset = 0; bsset < par.bins; bsset++ ) {
    double bsset_h = 0, bsset_h2 = 0, bsset_m = 0, bsset_m2 = 0;
    for ( unsigned int bssample = 0; bssample < par.bins; bssample++ ) {
      unsigned int bssample_this = gsl_rng_uniform_int( rng, par.bins );
      bsset_h  += binres[bssample_this].h;
      bsset_h2 += binres[bssample_this].h2;
      bsset_m  += binres[bssample_this].m;
      bsset_m2 += binres[bssample_this].m2;
    }
    bsset_h /= par.bins;
    bsset_h2 /= par.bins;
    bsset_m /= par.bins;
    bsset_m2 /= par.bins;

    // calculate the c and x for the selected set ...
    double bsset_c = spin_count / par.T / par.T
                     * ( bsset_h2 - bsset_h * bsset_h );
    double bsset_x = spin_count / par.T * ( bsset_m2 - bsset_m * bsset_m );

    sigma3_c += pow( ( bsset_c - c ), 2 ) / par.bins;
    sigma3_x += pow( ( bsset_x - x ), 2 ) / par.bins;
  }
  sigma3_c = 3 * sqrt( sigma3_c );
  sigma3_x = 3 * sqrt( sigma3_x );
  gsl_rng_free( rng );

  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> simulation finished";
  if ( par.calc_autocorr ) {
    run_log << " ... starting autocorrelation calculation";
  }
  run_log << "\n\n";
  run_log.flush();


  // ----- AUTOCORRELATION CALCULATION -----

  double tau = 0;
  if ( par.calc_autocorr ) {
    // open output file
    ofstream acout_log( ( dir + "ac.log" ).c_str() );
    if ( !acout_log.is_open() ) {
      cout << "ERROR while opening the acout_log file in " << dir << endl;
      return res;
    }
    // loop over different delta_t
    double norm = 1;
    for ( unsigned int d = 0; d < par.binwidth; d++ ) {
      double product = 0;
      unsigned int acsamples = m_memlog.size() - d;
      for ( unsigned int i = 0; i < acsamples; i++ ) {
        unsigned int sample = i;
        product += m_memlog[sample] * m_memlog[sample + d];
      }
      product = product / acsamples;
      if ( d == 0 ) {
        norm = ( product - m * m );
      }
      tau += ( product - m * m ) / norm * par.intersweeps;
      acout_log << d* par.intersweeps << ' '
                << ( product - m * m ) / norm << endl;
    }
    acout_log.close();

    // estimate if bins are correlated
    if ( tau > ( par.binwidth * par.intersweeps ) / 4 ) {
      cout << "WARNING: bins correlated in " << dir << endl;
      run_log << "!!!!! WARNING !!!!!\n"
              << "Bins correlated: Errors and autocorr time are underestimated!"
              << endl << endl;
      run_log.flush();
    }
  }


  // ----- WRITE SPIN-SPIN CORRELATIONS TO DISK -----

  if ( par.calc_sscorr ) {
    // correction for wrap-around errors with periodic boundaries
    if ( par.periodic ) {
      for ( unsigned int i = 1; i < ss_corr.size(); i++ ) {
        ss_corr[i] = ( ss_corr[i] + ss_corr[ss_corr.size() - 1] ) / 2;
        ss_corr.pop_back();
      }
    }

    // open output file
    ofstream sscorr_log( ( dir + "sscorr.log" ).c_str() );
    if ( !sscorr_log.is_open() ) {
      cout << "ERROR while opening the sscorr_log file in " << dir << endl;
      return res;
    }
    double norm = ss_corr[0] - m * m;
    for ( unsigned int i = 0; i < ss_corr.size(); i++ ) {
      sscorr_log << i << ' ' << ( ss_corr[i] - m * m ) / norm
                 << ' ' << ss_corr[i] / ss_corr[0] << endl;
    }
    sscorr_log.close();
  }


  // ---- RESULT OUTPUT -----

  // write simulation results to the output struct
  res.h = h;
  res.sigma3_h = sigma3_h;
  res.c = c;
  res.sigma3_c = sigma3_c;
  res.m = m;
  res.sigma3_m = sigma3_m;
  res.x = x;
  res.sigma3_x = sigma3_x;
  res.tau = tau;

  // write simulation results to the logfile
  run_log.precision( numeric_limits<float>::digits10 + 1 ); 
  run_log.setf( ios::scientific );
  run_log.setf( ios::showpos );
  run_log << "--- RESULTS ---\n"
          << "h	= " << h << "  +-  " << sigma3_h << endl
          << "c	= " << c << "  +-  " << sigma3_c << endl
          << "m	= " << m << "  +-  " << sigma3_m << endl
          << "chi	= " << x << "  +-  " << sigma3_x << endl
          << "tau = " << tau 		<< "\n\n";
  run_log.flush();

  // write results to a pyxplot readable output file
  ofstream results_file( ( dir + "results.dat" ).c_str() );
  if ( !results_file.is_open() ) {
    cout << "FATAL ERROR: unable to create results.dat in " << dir << endl;
    return res;
  }
  results_file << setiosflags( ios::scientific );
  results_file.setf( ios::showpos );
  results_file.precision( numeric_limits<float>::digits10 + 1 );
  results_file << par.N << ' '
               << par.J << ' '
               << par.g << ' '
               << par.B << ' '
               << par.T << ' '
               << res.h << ' ' << res.sigma3_h << ' '
               << res.c << ' ' << res.sigma3_c << ' '
               << res.m << ' ' << res.sigma3_m << ' '
               << res.x << ' ' << res.sigma3_x << ' '
               << res.tau		<< endl;
  results_file.close();


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
      << ", init = " << par.init << ", drysweeps = " << par.drysweeps
      << ", bins = " << par.bins << ", binwidth = " << par.binwidth
      << ", intersweeps = " << par.intersweeps
      << ", J = " << par.J << ", g = " << par.g
      << ", B = " << par.B << ", T = " << par.T;
  int plot_size = ( par.bins * par.binwidth < 1600 ) ?
                  1600 : min( ( uint )10000, par.bins * par.binwidth );
  run_plot <<
           "	set terminal pngcairo size " << plot_size << ",600 \n\
	set output 'run_plot.png' \n\n\
	set multiplot layout 2,1 title '" << tmp.str() << "'\n\
	set grid x y \n\
	set mxtics 10 \n\
	set key top left \n\
	set arrow from graph 0,first " << m << " to graph 1,first " << m
           << " ls -1 lw 2 nohead" << endl;
  for ( unsigned int bin = 0; bin < par.bins; bin++ ) {
    double sep = double( bin ) / double( par.bins );
    run_plot <<"\
    set arrow from graph " << sep << ",0 to graph " << sep << ",1 nohead\n\
		set arrow from graph " << sep << ", first " << binres[bin].m << "\
		to graph " << sep + 1 / double( par.bins ) << ", first " << binres[bin].m
             << " ls -1 lw 1 nohead\n";
  }
  run_plot << "\
  plot 'm.log' with lines ls 2 title 'magnetization ("
  << m << " +- "	<< sigma3_m << ")' \n\
	unset arrow\n\
	set arrow from graph 0,first " << h << " to graph 1,first "
  << h << " ls -1 lw 2 nohead" << endl;
  for ( unsigned int bin = 0; bin < par.bins; bin++ ) {
    double sep = double( bin ) / double( par.bins );
    run_plot << "\
    set arrow from graph " << sep << ",0 to graph " << sep << ",1 nohead\n\
		set arrow from graph " << sep << ", first " << binres[bin].h << "\
		to graph " << sep + 1 / double( par.bins ) << ", first " << binres[bin].h
    << " ls -1 lw 1 nohead\n";
  }
  run_plot << "\
  plot 'h.log' with lines ls 1 title 'energy per spin ("
  << h << " +- " << sigma3_h << ")' \n";
  tmp.str( "" );
  run_plot.close();

  if ( par.run_plot && par.bins * par.binwidth <= 1e5 ) {
    if ( system( ( "cd " + dir + " ; gnuplot run_plot.gnu" ).c_str() ) != 0 ) {
      cout << "ERROR while running gnuplot run_plot.gnu in " << dir << endl;
      return res;
    }
  }

  // magnetization and energy histogram plots
  ofstream histo_plot( ( dir + "histo_plot.pyx" ).c_str() );
  if ( !histo_plot.is_open() ) {
    cout << "ERROR while opening histo_plot.gnu in " << dir << endl;
    return res;
  }
  histo_plot <<
             "	set terminal pdf\n\
	set output 'mhisto_plot.pdf'\n\
	set title 'magnetization histogram'\n\
	histogram m() 'm.log' using $2 binorigin 0.005 binwidth 0.01\n\
	plot m(x) with boxes notitle \n\
	set output 'hhisto_plot.pdf'\n\
	set title 'energy histogram'\n\
	histogram h() 'h.log' using $2 binorigin 0.005 binwidth 0.01\n\
	plot h(x) with boxes notitle";
  histo_plot.close();
  if ( par.run_plot ) {
    if ( system( ( "cd " + dir + " ; pyxplot histo_plot.pyx" ).c_str() ) != 0 ) {
      cout << "ERROR running pyxplot histo_plot.pyx in " << dir << endl;
      return res;
    }
  }

  // ac_plot
  if ( par.calc_autocorr ) {
    ofstream ac_plot( ( dir + "ac_plot.gnu" ).c_str() );
    if ( !ac_plot.is_open() ) {
      cout << "ERROR while opening ac_plot.gnu in " << dir << endl;
      return res;
    }
    ac_plot <<
            "		set terminal pngcairo size 1000,600\n\
		set output 'ac_plot.png'\n\
		set grid x y\n\
		set arrow from graph 0, first 0 to graph 1, first 0 nohead lw 2\n\
		set title 'magnetization autocorrelation'\n\
		set arrow from first " << tau << ", graph 0 to first "
    << tau << ", graph 1 nohead\n\
		plot 'ac.log' with linespoints notitle";
    ac_plot.close();
    if ( par.run_plot ) {
      if ( system( ( "cd " + dir + " ; gnuplot ac_plot.gnu" ).c_str() ) != 0 ) {
        cout << "ERROR while running gnuplot ac_plot.gnu in " << dir << endl;
        return res;
      }
    }
  }

  // sscorr_plot
  if ( par.calc_sscorr ) {
    ofstream sscorr_plot( ( dir + "sscorr_plot.gnu" ).c_str() );
    if ( !sscorr_plot.is_open() ) {
      cout << "ERROR while opening sscorr_plot.gnu in " << dir << endl;
      return res;
    }
    sscorr_plot <<
                " 		set terminal pngcairo size 1000,600\n\
		set output 'sscorr_plot.png'\n\
		set grid x y\n\
		set arrow from graph 0, first 0 to graph 1, first 0 nohead lw 2\n\
		set title 'spin-spin correlation'\n\
		plot 'sscorr.log' u 1:2 with linespoints title '<s1s2> - <m>^2', \\\n\
			 'sscorr.log' u 1:3 with linespoints title '<s1s2> '";
    sscorr_plot.close();
    if ( par.run_plot ) {
      if ( system( ( "cd " + dir + " ; gnuplot sscorr_plot.gnu" ).c_str() )
                                                                        != 0 ) {
        cout << "ERROR running gnuplot sscorr_plot.gnu in " << dir << endl;
        return res;
      }
    }
  }

  // everything is fine ... return the results!
  time( &rawtime );
  run_log << ctime( &rawtime ) << "-> everything finished!\n\n";
  run_log.close();
  res.success = true;
  return res;
}
