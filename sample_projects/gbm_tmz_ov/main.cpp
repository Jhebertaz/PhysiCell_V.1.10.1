/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>



#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"

// put custom code modules here!
#include "./custom_modules/custom.h"


using namespace BioFVM;
using namespace PhysiCell;


double one_column_bin_reader(int idx)
{
	// open file
	std::ifstream file("./config/CSF_TMZ.bin", std::ios::in|std::ios::binary|std::ios::ate);

	if (file.is_open())
	{

		// stream position
		std::streampos size = file.tellg();
		char * memblock = new char [size];

		file.seekg (0, std::ios::beg);
		file.read (memblock, size);

		//reinterpret as doubles
		double* double_values = (double*)memblock;

		std::cout<<"conc returning: "<<	double_values[idx]<<std::endl;

		delete[] memblock;
		file.close();
		return double_values[idx];
	}
	else
	{
		std::cout<<"failed"<<std::endl;
	}
	return 0.0;
}

double CSF_conc_to_density(double time)
{
	int indexing_CSF_Vec = time-1;
	double CSF_value = one_column_bin_reader(time);
	double density_in_voxels = CSF_value*1e3/1e6; // converting miligram to microgram and liter to microliter

	return density_in_voxels;
}


int main( int argc, char* argv[] )
{
	// load and parse settings file(s)

	bool XML_status = false;
	char copy_command [1024];
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] );
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() );
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() );
	}
	if( !XML_status )
	{ exit(-1); }

	// copy config file to output directry
	system( copy_command );

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// time setup
	std::string time_units = "min";

	/* Microenvironment setup */

	setup_microenvironment(); // modify this in the custom code


	/* PhysiCell setup */

	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30;
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	/* Users typically start modifying here. START USERMODS */

	create_cell_types();
	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */

	// set MultiCellDS save options

	set_save_biofvm_mesh_as_matlab( true );
	set_save_biofvm_data_as_matlab( true );
	set_save_biofvm_cell_data( true );
	set_save_biofvm_cell_data_as_custom_matlab( true );

	// save a simulation snapshot

	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

	// save a quick SVG cross section through z = 0, after setting its
	// length bar to 200 microns

	PhysiCell_SVG_options.length_bar = 200;

	// for simplicity, set a pathology coloring function
	std::vector<std::string> (*cell_coloring_function)(Cell*) = colouring_by_intracellular_virus_amount;

	// std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;

	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() );
	create_plot_legend( filename , cell_coloring_function );

	display_citations();

	// set the performance timers

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() );

		report_file.open(filename); 	// create the data log file
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}

	// main loop
	int number_of_TMZ_updates = 1;
	double last_time_updated = 0; // keep tabs of the last time the radius was increased
	double number_of_times_updated = 1; // the number of times the radius has been increased

	std::vector<double> time_to_radius(62);
	time_to_radius = time_to_radius_increase();

	try
	{
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time.
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout );
				if( PhysiCell_settings.enable_legacy_saves == true )
				{
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}

				if( PhysiCell_settings.enable_full_saves == true )
				{
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );
				}

				PhysiCell_globals.full_output_index++;
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}

			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

					PhysiCell_globals.SVG_output_index++;
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			// static int virus_index = microenvironment.find_density_index( "virus");
			static int wall_index = microenvironment.find_density_index( "wall");
		  static double tumour_radius_initial = parameters.doubles("R");

			// update migratory domain radius
			if( PhysiCell_globals.current_time > last_time_updated + time_to_radius[number_of_times_updated] )
			{
				for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
				{
					std::vector<double> Walldense = microenvironment.mesh.voxels[n].center;
					if( Walldense[0]*Walldense[0]+Walldense[1]*Walldense[1]<(tumour_radius_initial+10*number_of_times_updated)*(tumour_radius_initial+10*number_of_times_updated)
						&& Walldense[0]*Walldense[0]+Walldense[1]*Walldense[1]>(tumour_radius_initial+10*(number_of_times_updated-1))*(tumour_radius_initial+10*(number_of_times_updated-1)) )	//1270*1270)
					{
						microenvironment(n)[wall_index] = 3.5;
					}
				}
				last_time_updated = PhysiCell_globals.current_time;
				number_of_times_updated += 1;
			}

			// CODE TO INCREASE TMZ in the brain at periphery
			// Adding PK simulation data to microenvironment
			static int tmz_index = microenvironment.find_density_index( "tmz");
			int PKPD_time_grid = 10;

			if( PhysiCell_globals.current_time > PKPD_time_grid*number_of_TMZ_updates && PhysiCell_globals.current_time<=7200)
			{
				double density_in_voxels = CSF_conc_to_density( number_of_TMZ_updates );
				for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )//update voxels
				{
					std::vector<double> loc_vector = microenvironment.mesh.voxels[n].center;
					if( loc_vector[0]*loc_vector[0]+loc_vector[1]*loc_vector[1]<(tumour_radius_initial+10)*(tumour_radius_initial+10)
						&& loc_vector[0]*loc_vector[0]+loc_vector[1]*loc_vector[1]>(tumour_radius_initial-10)*(tumour_radius_initial-10))	//1270*1270)
					{
						microenvironment(n)[tmz_index] = microenvironment(n)[tmz_index]+density_in_voxels;
					}
				}
				number_of_TMZ_updates +=1;
			}
			 // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			// run PhysiCell
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );

			/*
			  Custom add-ons could potentially go here.
			*/

			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}

	// save a final simulation snapshot

	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() );
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() );
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

	// timer

	std::cout << std::endl << "Total simulation runtime: " << std::endl;
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );

	return 0;
}
