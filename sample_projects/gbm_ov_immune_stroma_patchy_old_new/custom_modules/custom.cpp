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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>

void create_cell_types( void ) 																									// done
{
	// set the random seed
	SeedRandom( parameters.ints("random_seed") );

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );

	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This intializes cell signal and response dictionaries
	*/

	setup_signal_behavior_dictionaries();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	static int wall_index = microenvironment.find_density_index("wall");
	static int virus_index = microenvironment.find_density_index("virus");
	Cell_Definition* pCD;

	// th cell creation
	pCD = find_cell_definition("th");

	// Make sure we're ready for 2D
	// set functions
	pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = th_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_CD4");

	pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.fraction_released_at_death[wall_index] = 1;

	// motility
	pCD->phenotype.motility.is_motile = true;
	pCD->phenotype.motility.restrict_to_2D = true;
	pCD->phenotype.motility.migration_speed = parameters.doubles("nu");

	pCD->functions.cycle_model = Ki67_basic;
	pCD->phenotype.cycle.sync_to_cycle_model( Ki67_basic);
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );
	pCD->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4");
	pCD->phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("r_10");

	// cell actions
	pCD->phenotype.secretion.uptake_rates[virus_index] = 0.0;
	pCD->phenotype.secretion.secretion_rates[virus_index] = 0.0;
	pCD->phenotype.motility.migration_speed = parameters.doubles("nu");

	// cell morphology
	// pCD->phenotype.geometry.radius = parameters.doubles("R_CD4");
	// pCD->phenotype.volume.total = 4/3*3.141*parameters.doubles("R_CD4")*parameters.doubles("R_CD4")*parameters.doubles("R_CD4");;
	// pCD->phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	// pCD->phenotype.volume.fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.total;
	// pCD->phenotype.volume.solid = pCD->phenotype.volume.total-pCD->phenotype.volume.fluid;
	// pCD->phenotype.volume.nuclear = parameters.doubles("V_N_CD4");
	// pCD->phenotype.volume.nuclear_solid = 23.8;
	// pCD->phenotype.volume.nuclear_fluid = pCD->phenotype.volume.nuclear - pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.cytoplasmic = pCD->phenotype.volume.total - pCD->phenotype.volume.nuclear;
	// pCD->phenotype.volume.cytoplasmic_fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.cytoplasmic;
	// pCD->phenotype.volume.cytoplasmic_solid = pCD->phenotype.volume.cytoplasmic-pCD->phenotype.volume.cytoplasmic_fluid;
	// pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.05;
	// pCD->phenotype.volume.target_solid_cytoplasmic = pCD->phenotype.volume.cytoplasmic_solid;
	// pCD->phenotype.volume.target_solid_nuclear = pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.target_fluid_fraction = pCD->phenotype.volume.fluid_fraction;
	// pCD->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio;

	// update function
	// pCD->functions.update_phenotype = th_phenotype;



	// cancer cell creation
	pCD = find_cell_definition("cancer");

	// Make sure we're ready for 2D
	// set functions
	pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = cancer_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_GBM");

	pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.fraction_released_at_death[wall_index] = 1;

	// motility
	pCD->phenotype.motility.is_motile = true;
	pCD->phenotype.motility.restrict_to_2D = true;
	pCD->phenotype.motility.migration_speed = 0;

	// mechanics
	// pCD->phenotype.mechanics.cell_cell_repulsion_strength = 0.35*pCD->phenotype.mechanics.cell_cell_repulsion_strength;
	//setting cycle model to live
	pCD->phenotype.cycle.sync_to_cycle_model( live );
	cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	std::cout<<"cycle_start_index"<<cycle_start_index<<std::endl;
	std::cout<<"cycle_end_index"<<cycle_end_index<<std::endl;

	pCD->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0; //parameters.doubles("beta");

	// cell actions
	pCD->phenotype.secretion.uptake_rates[virus_index] = 0.0;
	pCD->phenotype.secretion.secretion_rates[virus_index] = 0.0;
	pCD->phenotype.motility.migration_speed = 0.0;

	// cell morphology
	// pCD->phenotype.geometry.radius = parameters.doubles("R_GBM");
	// pCD->phenotype.volume.total = 4/3*3.1416*parameters.doubles("R_GBM")*parameters.doubles("R_GBM")*parameters.doubles("R_GBM");
	// pCD->phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	// pCD->phenotype.volume.fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.total;
	// pCD->phenotype.volume.solid = pCD->phenotype.volume.total-pCD->phenotype.volume.fluid;
	// pCD->phenotype.volume.nuclear = parameters.doubles("V_N_GBM");
	// pCD->phenotype.volume.nuclear_solid = 185;
	// pCD->phenotype.volume.nuclear_fluid = pCD->phenotype.volume.nuclear - pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.cytoplasmic = pCD->phenotype.volume.total - pCD->phenotype.volume.nuclear;
	// pCD->phenotype.volume.cytoplasmic_fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.cytoplasmic;
	// pCD->phenotype.volume.cytoplasmic_solid = pCD->phenotype.volume.cytoplasmic-pCD->phenotype.volume.cytoplasmic_fluid;
	// pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio = pCD->phenotype.volume.cytoplasmic/pCD->phenotype.volume.nuclear;//6.0321;
	// pCD->phenotype.volume.target_solid_cytoplasmic = pCD->phenotype.volume.cytoplasmic_solid;
	// pCD->phenotype.volume.target_solid_nuclear = pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.target_fluid_fraction = pCD->phenotype.volume.fluid_fraction;
	// pCD->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio;
	// pCD->phenotype.volume.calcified_fraction = 0;
	// pCD->phenotype.volume.calcification_rate = 0;
	// pCD->phenotype.mechanics.cell_cell_repulsion_strength = 0.35*pCD->phenotype.mechanics.cell_cell_repulsion_strength;
	//
	// //update phenoypt
	// pCD->functions.update_phenotype = cancer_phenotype;



	// ctl cell creation
	pCD = find_cell_definition("ctl");

	// Make sure we're ready for 2D
	// set functions
	pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = NULL;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_CD8");

	pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.fraction_released_at_death[wall_index] = 1;

	// motility
	pCD->phenotype.motility.is_motile = true;
	pCD->phenotype.motility.restrict_to_2D = true;
	pCD->phenotype.motility.migration_speed = parameters.doubles("nu");


	// proliferation
	pCD->functions.cycle_model = Ki67_basic;
	pCD->phenotype.cycle.sync_to_cycle_model( Ki67_basic);
	cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive);
	pCD->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");//7.2206*1e-5;
	pCD->phenotype.cycle.data.transition_rate(cycle_end_index,cycle_start_index) = parameters.doubles("r_10");//0.00143;

	// cell actions
	pCD->phenotype.secretion.uptake_rates[virus_index] = 0.0;
	pCD->phenotype.secretion.secretion_rates[virus_index] = 0.0;
	pCD->phenotype.motility.migration_speed = parameters.doubles("nu");

	//cell morphology
	// pCD->phenotype.geometry.radius = parameters.doubles("R_CD8");
	// pCD->phenotype.volume.total = 4/3*3.141*parameters.doubles("R_CD8")*parameters.doubles("R_CD8")*parameters.doubles("R_CD8");
	// pCD->phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	// pCD->phenotype.volume.fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.total;
	// pCD->phenotype.volume.solid = pCD->phenotype.volume.total-pCD->phenotype.volume.fluid;
	// pCD->phenotype.volume.nuclear = parameters.doubles("V_N_CD8");
	// pCD->phenotype.volume.nuclear_solid = 24.1;
	// pCD->phenotype.volume.nuclear_fluid = pCD->phenotype.volume.nuclear - pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.cytoplasmic = pCD->phenotype.volume.total - pCD->phenotype.volume.nuclear;
	// pCD->phenotype.volume.cytoplasmic_fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.cytoplasmic;
	// pCD->phenotype.volume.cytoplasmic_solid = pCD->phenotype.volume.cytoplasmic-pCD->phenotype.volume.cytoplasmic_fluid;
	// pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio = 1.03;
	// pCD->phenotype.volume.target_solid_cytoplasmic = pCD->phenotype.volume.cytoplasmic_solid;
	// pCD->phenotype.volume.target_solid_nuclear = pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.target_fluid_fraction = pCD->phenotype.volume.fluid_fraction;
	// pCD->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio;

	// cell update phenotype
	// pCD->functions.update_phenotype = ctl_phenotype;


	// stroma cell creation
	pCD = find_cell_definition("stromal");


	// Make sure we're ready for 2D
	// set functions
	pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = stromal_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_stroma");

	pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	pCD->phenotype.molecular.fraction_released_at_death[wall_index] = 1;

	// motility
	pCD->phenotype.motility.is_motile = false;
	pCD->phenotype.motility.restrict_to_2D = true;
	pCD->phenotype.motility.migration_speed = 0;


	// turn off proliferation
	cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	cycle_end_index = live.find_phase_index( PhysiCell_constants::live);
	std::cout << "stroma_cell cycle_start_index" << cycle_start_index << std::endl;
	pCD->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0;

	// cell actions
	pCD->phenotype.secretion.uptake_rates[virus_index] = parameters.doubles("u_s");
	pCD->phenotype.secretion.secretion_rates[virus_index] = 0.0;
	pCD->phenotype.motility.migration_speed = 0;
	pCD->phenotype.motility.is_motile = false;
	pCD->phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;

	// cell morphology
	// pCD->phenotype.geometry.radius = parameters.doubles("R_stromal");//7.5;
	// pCD->phenotype.volume.total = 4/3*3.141*parameters.doubles("R_stromal")*parameters.doubles("R_stromal")*parameters.doubles("R_stromal");//1767;
	// pCD->phenotype.volume.fluid_fraction = parameters.doubles("f_F");
	// pCD->phenotype.volume.fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.total;
	// pCD->phenotype.volume.solid = pCD->phenotype.volume.total-pCD->phenotype.volume.fluid;
	// pCD->phenotype.volume.nuclear = parameters.doubles("V_N_stroma");
	// pCD->phenotype.volume.nuclear_solid = 125;
	// pCD->phenotype.volume.nuclear_fluid = pCD->phenotype.volume.nuclear - pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.cytoplasmic = pCD->phenotype.volume.total - pCD->phenotype.volume.nuclear;
	// pCD->phenotype.volume.cytoplasmic_fluid = pCD->phenotype.volume.fluid_fraction*pCD->phenotype.volume.cytoplasmic;
	// pCD->phenotype.volume.cytoplasmic_solid = pCD->phenotype.volume.cytoplasmic-pCD->phenotype.volume.cytoplasmic_fluid;
	// pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio = 2.53;
	// pCD->phenotype.volume.target_solid_cytoplasmic = pCD->phenotype.volume.cytoplasmic_solid;
	// pCD->phenotype.volume.target_solid_nuclear = pCD->phenotype.volume.nuclear_solid;
	// pCD->phenotype.volume.target_fluid_fraction = pCD->phenotype.volume.fluid_fraction;
	// pCD->phenotype.volume.target_cytoplasmic_to_nuclear_ratio = pCD->phenotype.volume.cytoplasmic_to_nuclear_ratio;

	// update phenotype
	// pCD->functions.update_phenotype = stromal_phenotype;


	/*
		 This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}
void setup_microenvironment(void) 																							// done
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.


	// initialize BioFVM

	initialize_microenvironment();
	testing();

	return;
}
void setup_tissue( void ) 																											// done
{
	// double Xmin = microenvironment.mesh.bounding_box[0];
	// double Ymin = microenvironment.mesh.bounding_box[1];
	// double Zmin = microenvironment.mesh.bounding_box[2];
	//
	// double Xmax = microenvironment.mesh.bounding_box[3];
	// double Ymax = microenvironment.mesh.bounding_box[4];
	// double Zmax = microenvironment.mesh.bounding_box[5];
	//
	// if( default_microenvironment_options.simulate_2D == true )
	// {
	// 	Zmin = 0.0;
	// 	Zmax = 0.0;
	// }
	//
	// double Xrange = Xmax - Xmin;
	// double Yrange = Ymax - Ymin;
	// double Zrange = Zmax - Zmin;
	//
	// // create some of each type of cell
	//
	// Cell* pC;
	//
	// for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	// {
	// 	Cell_Definition* pCD = cell_definitions_by_index[k];
	// 	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	// 	for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
	// 	{
	// 		std::vector<double> position = {0,0,0};
	// 		position[0] = Xmin + UniformRandom()*Xrange;
	// 		position[1] = Ymin + UniformRandom()*Yrange;
	// 		position[2] = Zmin + UniformRandom()*Zrange;
	//
	// 		pC = create_cell( *pCD );
	// 		pC->assign_position( position );
	// 	}
	// }
	// std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();

	return;
}

// color function

// default
std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }
// custom
std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" );

	static int v_index = microenvironment.find_density_index( "virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );

	double p_min = 1;
	double p_max = parameters.doubles("alpha");

	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	if(n_I>10000 && pCell->type ==2)
	{
		//std::cout<<" High viral infection m_i = "<<n_I<<std::endl;
	}

	if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 );
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u

		if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{
					output[0] = "orange";
					output[1] = "orange";
					output[2] = "coral";
					output[3] = "coral";
					return output;
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{
					output[0] = "darkred";
					output[1] = "darkred";
					output[2] = "firebrick";
					output[3] = "firebrick";
			}
		}
		if(pCell->type == 3)
		{
			if(pCell-> phenotype.cycle.data.current_phase_index==0)
			{
			  output[0] = "aquamarine";
				output[1] = "lightsteelblue";
				output[2] = "lightskyblue";
				output[3] = "aquamarine";
			}
			else if(pCell-> phenotype.cycle.data.current_phase_index==1)
			{

				output[0] = "darkslateblue";
				output[1] = "darkblue";
				output[2] = "darkblue";
				output[3] = "aquamarine";
			}
		}
		if(pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if(n_I>1)
			{
				double p_min = 1;
				int oncoprotein1 = (int) round((210-139)*(1.0/(p_max-p_min)) * (n_I-1));
				int oncoprotein2 = (int) round((180-69)*(1.0/(p_max-p_min)) * (n_I-1));
				int oncoprotein3 = (int) round((140-19)*(1.0/(p_max-p_min)) * (n_I-1));

				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", 210-oncoprotein1, 180-oncoprotein2, 140-oncoprotein3); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u

				output[0].assign( szTempString );
				output[1]="brown";
				output[2].assign( szTempString );
				output[3]="brown";

				return output;
			}
			else
			{
					output[0] = "rgb(104, 55, 99)";//"orchid";//"rgb(255,230,230)";
					output[1] = "rgb(104, 55, 99)";
					output[2] = "rgb(85, 50, 70)";//"plum";//"rgb(255,230,230)";
					output[3] = "rgb(85, 50, 70)";

					return output;
			}
		}
	if( pCell->phenotype.death.dead == true )
	{
			output[0] = "rgb(255, 255, 224)";
			output[1] = "rgb(255, 255, 224)";
			output[2] = "rgb(255, 228, 181)";
			output[3] = "rgb(255, 228, 181)";

		return output;
	}
	if( pCell->type == 4)
	{
		output[0] = "rgb(234, 172, 199)";
		output[1] = "rgb(234, 172, 199)";
		output[2] = "rgb(243, 186, 211)";
		output[3] = "rgb(243, 186, 211)";
	}
	return output;
}


// Custom distribution
std::vector<double> time_to_radius_increase()
{
	std::vector<double> time_to_radius(32);
	time_to_radius[0] = 460;
	time_to_radius[1] = 457;
	time_to_radius[2] = 455;
	time_to_radius[3] = 452;
	time_to_radius[4] = 450;
	time_to_radius[5] = 448;
	time_to_radius[6] = 445;
	time_to_radius[7] = 443;
	time_to_radius[8] = 441;
	time_to_radius[9] = 439;
	time_to_radius[10] = 436;
	time_to_radius[11] = 434;
	time_to_radius[12] = 432;
	time_to_radius[13] = 430;
	time_to_radius[14] = 428;
	time_to_radius[15] = 426;
	time_to_radius[16] = 424;
	time_to_radius[17] = 422;
	time_to_radius[18] = 420;
	time_to_radius[19] = 418;
	time_to_radius[20] = 416;
	time_to_radius[21] = 414;
	time_to_radius[22] = 412;
	time_to_radius[23] = 410;
	time_to_radius[24] = 409;
	time_to_radius[25] = 407;
	time_to_radius[26] = 405;
	time_to_radius[27] = 403;
	time_to_radius[28] = 401;
	time_to_radius[29] = 400;
	time_to_radius[30] = 398;
	time_to_radius[31] = 10000000;
	return time_to_radius;
}
std::vector<double> go_times_cumulative()
{
	// setting up distributions for movement and persistance of cells
	std::vector<double> go_times_cumul(8);
  go_times_cumul[0] = 0.01;
	go_times_cumul[1] = 0.962;
	go_times_cumul[2] = 0.9735;
	go_times_cumul[3] = 0.9835;
	go_times_cumul[4] = 0.9935;
	go_times_cumul[5] = 0.9955;
	go_times_cumul[6] = 0.9975;
	go_times_cumul[7] = 1;

	return go_times_cumul;
}
std::vector<double> persistance_distribution()
{
	std::vector<double> persistance_times_vec(8);
  persistance_times_vec[0] = 0;
	persistance_times_vec[1] = 30;
	persistance_times_vec[2] = 60;
	persistance_times_vec[3] = 90;
	persistance_times_vec[4] = 120;
	persistance_times_vec[5] = 150;
	persistance_times_vec[6] = 180;
	persistance_times_vec[7] = 240;

	return  persistance_times_vec;
}
std::vector<double> speed_cumulative()
{
	std::vector<double> speed_cumul(12);
  speed_cumul[0] = 0.0014;
	speed_cumul[1] = 0.0317;
	speed_cumul[2] = 0.2441;
	speed_cumul[3] = 0.5137;
	speed_cumul[4] = 0.7598;
	speed_cumul[5] = 0.8822;
	speed_cumul[6] = 0.9453;
	speed_cumul[7] = 0.9787;
	speed_cumul[8] = 0.9882;
	speed_cumul[9] = 0.9937;
	speed_cumul[10] = 0.9963;
	speed_cumul[11] = 1;

	return speed_cumul;
}
std::vector<double> speed_distribution()
{
	std::vector<double> speed_vec(12);
  speed_vec[0] = 0.0833;
	speed_vec[1] = 0.1667;
	speed_vec[2] = 0.25;
	speed_vec[3] = 0.333;
	speed_vec[4] = 0.4167;
	speed_vec[5] = 0.5;
	speed_vec[6] = 0.5833;
	speed_vec[7] = 0.667;
	speed_vec[8] = 0.75;
	speed_vec[9] = 0.833;
	speed_vec[10] = 0.9167;
	speed_vec[11] = 1;

	return speed_vec;
}

// usefull functions
bool am_i_dead(Cell* pCell)
{
	bool dead = pCell->phenotype.death.dead;
	return dead;
}
bool am_i_infected(Cell* pCell)
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	int cell_intern_virions = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	bool infected = false;

	if (cell_intern_virions>=1)
	{
		infected = true;
	}
	return infected;
}
bool am_i_an_detectable_infected_cell(Cell* pCell)
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	double infection_threshold_half_effect = parameters.doubles("m_half");
	int cell_intern_virions = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	bool infected = false;

	if (cell_intern_virions>infection_threshold_half_effect)
	{
		infected = true;
	}
	return infected;
}
bool is_there_an_infected_cell_around(Cell* pCell)
{
	// neighborhood
	std::vector<Cell*> nearby = pCell->cells_in_my_container();
	// infection
	double m_half = parameters.doubles("m_half");

	// loop variables
	Cell* pC = NULL;
	bool stop = false;
	int i=0;

	while (!stop && i<nearby.size())
	{
		pC = nearby[i]; // ith cell object

		bool dead = pC->phenotype.death.dead;
		bool itself = (pC == pCell);
		bool infected_cell = am_i_an_detectable_infected_cell(pC);

		// infected cell around
		if (infected_cell && !dead && !itself)
		{
			stop = true;
		}
		i++;
	}
	return stop;
}

// th related functions
void th_movement(Cell* pCell, Phenotype& phenotype)
{
	static int wall_index = microenvironment.find_density_index("wall");
	double wall_amount = pCell->nearest_density_vector()[wall_index];

	// cell movement
	double nu = parameters.doubles("nu");

	// migration bias direction
	std::vector<double> dbias(3);

	if (wall_amount<2)
	{
		// turn around
		dbias = -1*pCell->position;
		normalize( &( dbias ) );

		phenotype.motility.migration_bias_direction = dbias;
		phenotype.motility.migration_bias = 1;
		phenotype.motility.migration_speed = nu;
	}
	else
	{
		phenotype.motility.migration_bias = 0;
		phenotype.motility.migration_speed = nu;
	}
	return;
}
void th_infection_dynamics(Cell* pCell, Phenotype& phenotype)
{
	static int chemokine_index = microenvironment.find_density_index("chemokine");
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );

	// Recurrent variables
	double S_chemokine_CD4 = parameters.doubles("S_chemokine_CD4");
	double r_01_CD4 = parameters.doubles("r_01_CD4");
	double prolif_increase = parameters.doubles("TH_prolif_increase_due_to_stimulus");

	// search for infected cell around
	bool infected_cell_around = is_there_an_infected_cell_around(pCell);

	if (infected_cell_around)
	{
		phenotype.secretion.secretion_rates[chemokine_index] = S_chemokine_CD4;
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD4*prolif_increase;
		phenotype.motility.migration_speed = 0.1;
	}
	else
	{
		phenotype.secretion.secretion_rates[chemokine_index] = 0;
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD4;
	}
	return;
}


// cancer related functions
void resample_movement_type(Cell* pCell, Phenotype& phenotype)
{
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
	double new_type_rand = UniformRandom();
	double speed;

	std::vector<double> speed_cumul(12);
	std::vector<double> speed_vec(12);

	speed_cumul = speed_cumulative();
	speed_vec = speed_distribution();

	if(new_type_rand<=0.5)
	{
		// assign GO type
		// pCell->custom_data["cell_motility_type"] = 1;
		double speed_var = UniformRandom();

		for (int k=0; k<12;)
		{
			if (speed_var> speed_cumul[k])
			{
				k++;
			}
			else
			{
				// assign migration speed
				// phenotype.motility.migration_speed = speed_vec[k];
				pCell->phenotype.motility.migration_speed = speed_vec[k];
				speed = speed_vec[k];
				k = 12;
			}
		}
	}
	else
	{
		// assign STOP type
		// pCell->custom_data["cell_motility_type"] = 2;
		// phenotype.motility.migration_speed = 0;
		pCell->custom_data.variables[cell_motility_type_index].value = 2;
		pCell->phenotype.motility.migration_speed = 0;;
	}

	// phenotype.motility.migration_speed = speed;
	return;
}
void resample_persistence_time(Cell* pCell, Phenotype& phenotype)
{
	double go_stop_var = UniformRandom();

	std::vector<double> go_times_cumul(8);
	std::vector<double> persistence_times_vec(12);

	go_times_cumul = go_times_cumulative();
	persistence_times_vec = persistance_distribution();

	// assign persistence time - needs to be a real time!
	for (int j=0; j<8;)
	{
		if (go_stop_var> go_times_cumul[j])
		{
			j++;
		}
		else
		{
			// assign persist time
			// pCell->custom_data["persistence_time"] = persistence_times_vec[j]+PhysiCell_globals.current_time;
			phenotype.motility.persistence_time = persistence_times_vec[j]+PhysiCell_globals.current_time;
			j = 8;
		}
	}
	return;
}
void update_virus_uptake_rate(Cell* pCell,  Phenotype& phenotype)
{
	static int virus_signal_index = microenvironment.find_density_index("virus");

	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];

	double p = pCell->nearest_density_vector()[virus_signal_index];
	double pmax = parameters.doubles("rho_max");

	double u = parameters.doubles("u_g");
	double Vvoxel = microenvironment.mesh.voxels[0].volume;
	double m_half = parameters.doubles("m_half");


	pCell->custom_data["intracellular_virus_amount"] = cell_intern_virions;
	if (p < pmax)
	{
		phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(cell_intern_virions/Vvoxel + m_half/Vvoxel);
	}
	else
	{
		phenotype.secretion.uptake_rates[virus_signal_index] = (u*pmax*pmax/p)*(1/(cell_intern_virions/Vvoxel + m_half/Vvoxel));
	}
	return;
}
void virus_replication(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	static int apoptosis_model_index = phenotype.death.find_death_model_index("apoptosis");

	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];

	double alpha = parameters.doubles("alpha");
	double gamma = parameters.doubles("gamma");
	double delta_V = parameters.doubles("delta_V");

	if (cell_intern_virions > 1 && cell_intern_virions <= alpha)
	{
		pCell->custom_data["intracellular_virus_amount"] = cell_intern_virions+dt*(gamma*cell_intern_virions);

		phenotype.molecular.internalized_total_substrates[virus_signal_index] = cell_intern_virions+dt*(gamma*cell_intern_virions);
		phenotype.secretion.secretion_rates[virus_signal_index] = 0;
	}
	else if (cell_intern_virions > alpha-1)
	{
		pCell->custom_data["intracellular_virus_amount"] = cell_intern_virions;

		phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
		phenotype.secretion.uptake_rates[virus_signal_index] = 0;
		// phenotype.secretion.secretion_rates[virus_signal_index] = delta_V;
		pCell->start_death(apoptosis_model_index);
	}
	else if (cell_intern_virions<0)
	{
		std::cout<<"Negative intracellular virus!! m_i = "<<cell_intern_virions<<std::endl;
	}
	return;
}
void virus_induced_lysis(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int virus_signal_index = microenvironment.find_density_index( "virus");

	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];
	// Virus saturation concentration
	double pstar = phenotype.secretion.saturation_densities[virus_signal_index];
  //
	double p = pCell->nearest_density_vector()[virus_signal_index];
	// volume of voxel
	double Vvoxel = microenvironment.mesh.voxels[0].volume;
	// release rate of virus
	double delta_V = parameters.doubles("delta_V");


	if (cell_intern_virions > 1)
	{
			double amount_to_add = (cell_intern_virions-cell_intern_virions*exp(-delta_V*dt))/Vvoxel;
			if (amount_to_add > pstar-p)
			{
				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
				phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
				pCell->custom_data["intracellular_virus_amount"] -= (pstar-p)*Vvoxel*dt;
			}
			else
			{
				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
				phenotype.molecular.internalized_total_substrates[virus_signal_index] = cell_intern_virions*exp(-delta_V*dt);
				pCell->custom_data["intracellular_virus_amount"] = cell_intern_virions*exp(-delta_V*dt);
			}
	}

	return;
}

// void cancer_movement(Cell* pCell, Phenotype& phenotype)
// {
// 	//cell movement
// 	static int wall_index = microenvironment.find_density_index( "wall" );
// 	double wall_amount = pCell->nearest_density_vector()[wall_index];
//
// 	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
// 	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );
//
// 	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
// 	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; // 1 = go, 2 = stop
// 	// double persistence_time = phenotype.motility.persistence_time;
// 	std::vector<double> go_times_cumul(8);
// 	std::vector<double> persistence_times_vec(12);
// 	std::vector<double> speed_cumul(12);
// 	std::vector<double> speed_vec(12);
// 	go_times_cumul = go_times_cumulative();
// 	persistence_times_vec = persistance_distribution();
// 	speed_cumul = speed_cumulative();
// 	speed_vec = speed_distribution();
//
// 	if( wall_amount<2 & pCell->type == 2 )
// 	{
// 		pCell->phenotype.motility.migration_speed = 0;
// 	}
// 	else
// 	{
// 		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
// 		{
// 			// assign new type (stop = 2, or go = 1)
// 			double new_type_rand = UniformRandom();
// 			if(new_type_rand<=0.5)// GO
// 			{
// 				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type
//
// 				double speed_var = UniformRandom();
//
// 				for( int k=0; k<12; )
// 				{
// 					if( speed_var> speed_cumul[k] )
// 					{k++;}
// 						else
// 					{
// 						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
// 						k = 12;
// 					}
// 				}
// 			}
// 			else
// 			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
// 			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type
//
// 			// assign persistence time - needs to be a real time!
// 			double go_stop_var = UniformRandom();
// 			for( int j=0; j<8; )
// 			{
// 				if( go_stop_var> go_times_cumul[j] )
// 				{j++;}
// 					else
// 				{
// 					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
// 				  phenotype.motility.persistence_time = persistence_times_vec[j]+PhysiCell_globals.current_time;
// 					j = 8;
// 				}
// 			}
// 		}
//
// 	}
// 	return;
// }
// 	//cell movement
// 	static int wall_index = microenvironment.find_density_index("wall");
// 	double wall_amount = pCell->nearest_density_vector()[wall_index];
// 	if (wall_amount<2)
// 	{
// 		phenotype.motility.migration_speed = 0;
// 	}
// 	else
// 	{
// 		// Has my motility persistence time expired
// 		// double persistence_time = pCell->custom_data["persistence_time"];
// 		double persistence_time = phenotype.motility.persistence_time;
// 		double wall_amount = pCell->nearest_density_vector()[wall_index];
//
// 		if (persistence_time <= PhysiCell_globals.current_time)
// 		{
// 			resample_movement_type(pCell, phenotype);
// 			resample_persistence_time(pCell, phenotype);
// 		}
//
// 	}
// 	return;
// }
// void cancer_infection_dynamics(Cell* pCell, Phenotype& phenotype, double dt)
// {
// 	if (!am_i_dead(pCell))// cell not dead
// 	{
// 		update_virus_uptake_rate(pCell, phenotype);
// 		if (am_i_infected(pCell))
// 		{
// 			virus_replication(pCell, phenotype, dt);
// 		}
// 	}
// 	else
// 	{
// 		virus_induced_lysis(pCell, phenotype, dt);
// 	}
// 	return;
// }






// ctl related functions
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype)
{
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["c_s"] );
	}

	return;
}
void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position;
	axpy( &(pActingOn->velocity) , elastic_constant , displacement );

	return;
}
Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker)
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while( i < nearby.size() )
	{
		// don't try to kill yourself
		if( nearby[i] != pAttacker )
		{
			if( immune_cell_attempt_attachment( pAttacker, nearby[i]) )
			{ return nearby[i]; }
		}
		i++;
	}

	return NULL;
}
bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget)
{
	static double max_attachment_distance = parameters.doubles("d_attach"); // 18.0;

	static int virus_signal_index = microenvironment.find_density_index( "virus");
	double internal_virus = pTarget->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;

	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment_lifetime" );

	double kill_time = parameters.doubles("tau"); // how long the cell needs to attach for  the infected cell to be killed
	double nstar = parameters.doubles("m_half"); // how long the cell needs to attach for  the infected cell to be killed

	if( internal_virus > nstar && pTarget->phenotype.death.dead == false && pTarget->type!=4)
	{
		std::vector<double> displacement = pTarget->position - pAttacker->position;
		double distance_scale = norm( displacement );
		if( distance_scale > max_attachment_distance )
		{ return false; }

		attach_cells( pAttacker, pTarget );
		pAttacker->custom_data[attach_lifetime_i] = PhysiCell_globals.current_time + kill_time;

		return true;
	}

	return false;
}
bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget)
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );
	static int attach_lifetime_i = pAttacker->custom_data.find_variable_index( "attachment_lifetime" );

	// CTL kills cell if it has been attached for enough time
	if( pAttacker->custom_data[attach_lifetime_i] < PhysiCell_globals.current_time )
	{
		return true;
	}
	return false;
}
bool immune_cell_trigger_apoptosis( Cell* pAttacker, Cell* pTarget )
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index( "apoptosis" );

	static int virus_index = microenvironment.find_density_index( "virus" );

	// if the Target cell is already dead, don't bother!
	if( pTarget->phenotype.death.dead == true )
	{ return false; }

	pTarget->start_death( apoptosis_model_index );
	pTarget->phenotype.molecular.fraction_released_at_death[virus_index] = 0;//1;
	return true;
}
void ctl_movement(Cell* pCell, Phenotype& phenotype)
{
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );

	static int chemokine_index = microenvironment.find_density_index("chemokine");
	static int wall_index = microenvironment.find_density_index("wall");

	// Recurrent values
	double b_CD8 = parameters.doubles("b_CD8");
	double nu = parameters.doubles("nu");
	double nu_max = parameters.doubles("nu_max");
	double nu_star = parameters.doubles("nu_star");
	double r_01_CD8 = parameters.doubles("r_01_CD8");
	double prolif_increase = parameters.doubles("CTL_prolif_increase_due_to_stimulus");

	// substrate amount
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	// gradient vectors
	std::vector<double> chemokine_gradient = pCell->nearest_gradient(chemokine_index);
	std::vector<double> wall_gradient = pCell->nearest_gradient(wall_index);


	std::vector<double> ae_ini(3);

	if (wall_amount<2)
	{
		ae_ini = -1*pCell->position;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;
		pCell->phenotype.motility.migration_bias = 1;
	}
	else if (chemokine_amount>1e-8)
	{
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD8;
		phenotype.motility.migration_bias_direction = chemokine_gradient;
		phenotype.motility.migration_bias = b_CD8;
		phenotype.motility.migration_speed = nu+(nu_max-nu)*(chemokine_amount/(nu_star+chemokine_amount));
	}
	else if (chemokine_amount>1e-3)
	{
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD8*prolif_increase;
		phenotype.motility.migration_bias_direction = chemokine_gradient;
		phenotype.motility.migration_bias = b_CD8;
		phenotype.motility.migration_speed = nu+(nu_max-nu)*(chemokine_amount/(nu_star+chemokine_amount));
	}
	else
	{
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD8;
		pCell->phenotype.motility.migration_bias = 0;
		pCell->phenotype.motility.migration_speed = nu;
	}
	return;
}
void ctl_infection_dynamics(Cell* pCell, Phenotype& phenotype)
{
	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );

	double r_01_CD8 = parameters.doubles("r_01_CD8");
	double prolif_increase = parameters.doubles("CTL_prolif_increase_due_to_stimulus");

	// infection dynamics
	if (pCell->state.neighbors.size() > 0)
	{
		extra_elastic_attachment_mechanics( pCell, phenotype);
		bool detach_me = false;

		if (immune_cell_attempt_apoptosis(pCell, pCell->state.neighbors[0]))
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] );
			detach_me = true;
		}

		if( detach_me )
		{
			detach_cells(pCell, pCell->state.neighbors[0]);
			phenotype.motility.is_motile = true;

			ctl_movement(pCell, phenotype);
			phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = r_01_CD8*prolif_increase;
		}
		return;
	}

	if( immune_cell_check_neighbors_for_attachment(pCell))
	{
		phenotype.motility.is_motile = false;
		return;
	}
	phenotype.motility.is_motile = true;
	return;
}


// custom phenotype functions
void th_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	//cell movement
	th_movement(pCell, phenotype);
	// infection_dynamics
	th_infection_dynamics(pCell, phenotype);
	return;
}
void cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	// static int wall_index = microenvironment.find_density_index("wall");
	//
	// int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	// int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
	//
	// // pressure around
	// double radius = parameters.doubles("R_GBM");
	// double s = parameters.doubles("xhi");
	// double pressure = 6*(1-1/(2*radius)*s)*(1-1/(2*radius)*s);
	// double pressure_scale = 0.027288820670331;
	// double maximal_pressure = pressure/pressure_scale;
	//
	// double beta = parameters.doubles("beta");
	// double cell_pressure = pCell->state.simple_pressure;
	//
	// // cell attemps proliferation
	// if (cell_pressure*cell_pressure > maximal_pressure)
	// {
	// 	phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0;
	// }
	// else
	// {
	// 	phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = beta;
	// }
	//
	// // cell movement
	//
	// cancer_movement(pCell, phenotype);
	// // infection dynamics
	// cancer_infection_dynamics(pCell, phenotype, dt);
	// return;


	double R = parameters.doubles("R_cell_GBM");
	double SA = 4*3.1416*R*R;
	double s = parameters.doubles("xhi");
	double pressure = 6*(1-1/(2*R)*s)*(1-1/(2*R)*s);
	double pressure_scale = 0.027288820670331;
	double max_pressure = pressure/pressure_scale;

	//tumour cell proliferation
	if(pCell->type ==2)
	{
		int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
		int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = parameters.doubles("beta");

		if( pCell->state.simple_pressure*pCell->state.simple_pressure>max_pressure) // if cell under too much pressure -> no proliferation
		{pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0;}

	}

	cell_movement( pCell, phenotype, dt);
	infection_dynamics( pCell, phenotype, dt );
	return;

}
void ctl_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}
	// infection dynamics
	// ctl_infection_dynamics(pCell, phenotype);
	// ctl_movement(pCell, phenotype);
	// return;

	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attachment_lifetime" );

	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics( pCell, phenotype);
		bool dettach_me = false;

		if( immune_cell_attempt_apoptosis( pCell, pCell->state.neighbors[0]) )
		{
			immune_cell_trigger_apoptosis( pCell, pCell->state.neighbors[0] );
			dettach_me = true;
		}

		if( dettach_me )
		{

			detach_cells( pCell, pCell->state.neighbors[0] );
			phenotype.motility.is_motile = true;
			CTL_movement( pCell, phenotype, dt);
			pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8")*parameters.doubles("CTL_prolif_increase_due_to_stimulus");
		}
		return;
	}

	if( immune_cell_check_neighbors_for_attachment( pCell) )
	{
		phenotype.motility.is_motile = false;
		return;
	}

	phenotype.motility.is_motile = true;
	CTL_movement( pCell, phenotype, dt);

	return;
}
void stromal_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	if(pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]>1e5)
	{
		pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
		pCell->custom_data.variables[intracellular_virus_index].value = 0;
	}
	else
	{
		pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
	}

	return;
}

// default
void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }
void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }
void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }




// // custom version
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}

	pCell->state.simple_pressure = 0.0;
	pCell->state.neighbors.clear(); // new 1.8.0

	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->custom_add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index =
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end;
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->custom_add_potentials(*neighbor);
		}
	}

	pCell->update_motility_vector(dt);
	pCell->velocity += phenotype.motility.motility_vector;

	return;
}







void cell_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	//cell movement
	static int wall_index = microenvironment.find_density_index( "wall" );
	double wall_amount = pCell->nearest_density_vector()[wall_index];

	static int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
	static int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );

	double persistence_time = pCell->custom_data.variables[persistence_time_index].value;
	double cell_motility_type = pCell->custom_data.variables[cell_motility_type_index].value; // 1 = go, 2 = stop

	std::vector<double> go_times_cumul(8);
	std::vector<double> persistence_times_vec(12);
	std::vector<double> speed_cumul(12);
	std::vector<double> speed_vec(12);
	go_times_cumul = go_times_cumulative();
	persistence_times_vec = persistance_distribution();
	speed_cumul = speed_cumulative();
	speed_vec = speed_distribution();

	if( wall_amount<2 & pCell->type == 2 )
	{
		pCell->phenotype.motility.migration_speed = 0;
	}
	else if(pCell->type != 2)
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
	}
	else
	{
		if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
		{
			// assign new type (stop = 2, or go = 1)
			double new_type_rand = UniformRandom();
			if(new_type_rand<=0.5)// GO
			{
				pCell->custom_data.variables[cell_motility_type_index].value = 1; // assign go type

				double speed_var = UniformRandom();

				for( int k=0; k<12; )
				{
					if( speed_var> speed_cumul[k] )
					{k++;}
					else
					{
						pCell->phenotype.motility.migration_speed = speed_vec[k]; // assign migration speed
						k = 12;
					}
				}
			}
			else
			{pCell->custom_data.variables[cell_motility_type_index].value = 2;
			pCell->phenotype.motility.migration_speed = 0;} // assign STOP type

			// assign persistence time - needs to be a real time!
			double go_stop_var = UniformRandom();
			for( int j=0; j<8; )
			{
				if( go_stop_var> go_times_cumul[j] )
				{j++;}
					else
				{
					pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j]+PhysiCell_globals.current_time; // assign persist time
					j = 8;
				}
			}
		}

	}
	return;
}

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{

	//cell infection
	static int virus_signal_index = microenvironment.find_density_index( "virus");
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );

	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];//custom_data.variables[intracellular_virus_index].value;
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double u = parameters.doubles("u_g");
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double nstar = parameters.doubles("m_half");//10;//infection threshol
	double alp = parameters.doubles("alpha");//1000;//virus burst number
	double nu = parameters.doubles("gamma");
	double pmax = parameters.doubles("rho_max");//0.0125;

	if( n>1)
	{
	}
	pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];

	if( pCell->phenotype.death.dead == false )// cell not dead
	{
		if(p<pmax)
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(n/Vvoxel+nstar/Vvoxel);}
		else
		{pCell->phenotype.secretion.uptake_rates[virus_signal_index] = u*pmax*pmax/(n/Vvoxel+nstar/Vvoxel)/p;}

		if( n > 1 && n <= alp) // update amount inside due to replication
		{
			pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n+dt*(nu*n);
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];

		}
		else if( n > alp-1)
		{
			pCell->custom_data.variables[intracellular_virus_index].value = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
			pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->start_death( apoptosis_model_index );
		}
		else if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<0 )
		{std::cout<<"Negative intracellular virus!! m_i = "<<pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]<<std::endl;}
	}
	if( pCell->phenotype.death.dead == true && pCell->phenotype.molecular.fraction_released_at_death[virus_signal_index]>0)
	{
		virus_induced_lysis(pCell, phenotype, dt );
	}

	return;
}

// void virus_induced_lysis( Cell* pCell, Phenotype& phenotype, double dt )
// {
// 	static int virus_signal_index = microenvironment.find_density_index( "virus");
// 	double pstar = pCell->phenotype.secretion.saturation_densities[virus_signal_index];
// 	double delta_V = parameters.doubles("delta_V");//0.1466;
// 	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
// 	double p = pCell->nearest_density_vector()[virus_signal_index];
// 	double n = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];
//
// 	if( pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index]> 1 )
// 	{
// 			double amount_to_add = (n-n*exp(-delta_V*dt))/Vvoxel;
// 			if( amount_to_add > pstar-p )
// 			{
// 				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
// 				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
// 			}
// 			else
// 			{
// 				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
// 				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = n*exp(-delta_V*dt);
// 			}
// 	}
//
// 	return;
// }








void CTL_movement( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int wall_index = microenvironment.find_density_index( "wall" );
	static int chemokine_index = microenvironment.find_density_index( "chemokine" );

	double wall_amount = pCell->nearest_density_vector()[wall_index];
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];



	int cycle_start_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int cycle_end_index = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );


	std::vector<double> ae_ini(3);

	// TH movement
	if (wall_amount<2)
	{
		std::cout<<"outside domain "<<pCell->phenotype.motility.migration_speed<<std::endl;
		pCell->phenotype.motility.migration_bias = 1;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(wall_index);
		ae_ini = -1*pCell->position;

		pCell->phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		pCell->phenotype.motility.migration_bias_direction = ae_ini;

		return;
	}
	else if(chemokine_amount>1e-8)// sample chemotaxis gradient and random walk in that direction
	{
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu")+(parameters.doubles("nu_max")-parameters.doubles("nu"))*(chemokine_amount/(parameters.doubles("nu_star")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("b_CD8");//0.85;
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");

		return;
	}
	else if(chemokine_amount>1e-3)
	{
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8")*parameters.doubles("CTL_prolif_increase_due_to_stimulus");//7.9026*1e-5*1e1;
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu")+(parameters.doubles("nu_max")-parameters.doubles("nu"))*(chemokine_amount/(parameters.doubles("nu_star")+chemokine_amount));
		pCell->phenotype.motility.migration_bias = parameters.doubles("b_CD8");
		pCell->phenotype.motility.migration_bias_direction = pCell->nearest_gradient(chemokine_index);

		return;
	}
	else
	{
		pCell->phenotype.motility.migration_bias = 0;
		pCell->phenotype.motility.migration_speed = parameters.doubles("nu");
		pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD8");

		return;
	}
}



// simulation related functions
void testing()
{
	static int virus_index = microenvironment.find_density_index( "virus");
	static int wall_index = microenvironment.find_density_index( "wall");
	static int chemokine_index = microenvironment.find_density_index( "chemokine");

	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")+10)*(parameters.doubles("R")+10))
		{
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 1;
			microenvironment(n)[chemokine_index] = 0;
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")-10)*(parameters.doubles("R")-10))
		{
			microenvironment(n)[virus_index] = 0;//parameters.doubles("initial_virus_density");
			microenvironment(n)[wall_index] = 3.5;
			microenvironment(n)[chemokine_index] = 0;
		}
		else if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>250*250 )
		{
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 5;
			microenvironment(n)[chemokine_index] = 0;

		}
		else
		{
			microenvironment(n)[virus_index] = 0;
			microenvironment(n)[wall_index] = 10;
			microenvironment(n)[chemokine_index] = 0;
		}
	}

	// locations of dense patches
	double x_offset_patch_1 = 173.38;
	double y_offset_patch_1 = -270.5;
	double x_offset_patch_2 = 41.32;
	double y_offset_patch_2 = 1015.1;
	double x_offset_patch_3 = -757.84;
	double y_offset_patch_3 = 128.4;
	double x_offset_patch_4 = 899.79;
	double y_offset_patch_4 = -391.6;
	double x_offset_patch_5 = -9.99;
	double y_offset_patch_5 = -859.9;

	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;
		if( (x_offset_patch_1-ECMdense[0])*(x_offset_patch_1-ECMdense[0])+(y_offset_patch_1-ECMdense[1])*(y_offset_patch_1-ECMdense[1])<20*20)// Centre of patch 1
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_2-ECMdense[0])*(x_offset_patch_2-ECMdense[0])+(y_offset_patch_2-ECMdense[1])*(y_offset_patch_2-ECMdense[1])<20*20)// Centre of patch 2
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_3-ECMdense[0])*(x_offset_patch_3-ECMdense[0])+(y_offset_patch_3-ECMdense[1])*(y_offset_patch_3-ECMdense[1])<20*20) //Centre of patch 3
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_4-ECMdense[0])*(x_offset_patch_4-ECMdense[0])+(y_offset_patch_4-ECMdense[1])*(y_offset_patch_4-ECMdense[1])<20*20) // Centre of patch 4
		{
			microenvironment(n)[virus_index] = 7.12;
		}
		else if((x_offset_patch_5-ECMdense[0])*(x_offset_patch_5-ECMdense[0])+(y_offset_patch_5-ECMdense[1])*(y_offset_patch_5-ECMdense[1])<20*20) // Centre of patch 5
		{
			microenvironment(n)[virus_index] = 7.12;
		}

	}
	return;
}
void immune_cell_placement(void)
{
	double x = 0.0;
	double y = 0.0;

	Cell* pCell = NULL;
	Cell_Definition* pCD;

	pCD = find_cell_definition("th");

	//UNPROLIFERATIVE TH CELLS
	x = -18.0725;
	y = 922.2948;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 954.1439;
	y = 56.7601;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	// /
  x = -115.9499;
	y = -349.9367;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
  x = 773.5733;
	y = -613.7738;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

  x = 753.4041;
	y = -219.3846;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///
	x = -292.5135;
	y = -91.2428;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
	x = 449.1283;
	y = 407.6992;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 452.4130;
	y = 616.9196;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  60.4108;
	y =  -1021.8;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 534.5535;
	y =  -882.0731;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  15.5367;
	y = 400.0095;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 412.5616;
	y = -921.4520;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  64.9108;
	y = 1021.7;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 669.1877;
	y = -318.7248;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -483.3989;
	y = 759.6808;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 134.3949;
	y = 347.2437;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -47.5182;
	y = 611.5124;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -738.2632;
	y = -659.5292;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -937.1317;
	y = 156.6883;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -614.4497;
	y = 789.3966;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  437.2258;
	y = -793.7856;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -170.1507;
	y = -132.5528;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -885.0803;
	y = -277.3879;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  838.1983;
	y = -501.9633;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -260.7010;
	y =  604.0757;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  103.6960;
	y = -899.8739;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  -23.0385;
	y = -922.1671;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -458.5337;
	y = 474.7504;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -765.6844;
	y = -326.6257;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  360.4622;
	y = 228.2000;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  771.4487;
	y = 319.7186;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -195.5453;
	y = -5.9414;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  114.0970;
	y = -561.2747;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///
	x =  203.4536;
	y = -89.5478;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  235.8782;
	y = 265.3679;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -868.0050;
	y = -405.4009;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -811.6261;
	y = 187.1180;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	//PROLIFERATIVE TH CELLS
	x = 454.1934;
	y = -46.7240;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	x = -427.8415;
	y = 616.7429;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;


	pCD = find_cell_definition("ctl");
	//UNPROLIFERATIVE CTL CELLS
	x = 27.9074;
	y = 119.1812;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 225.7304;
	y = -728.5426;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = -258.8621;
	y =  265.2323;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -809.2029;
	y =  -60.6306;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 473.7135;
	y =  887.8817;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = -366.2289;
	y = -332.2614;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -197.1932;
	y =  446.5596;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -393.1067;
	y =  -597.4458;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  -264.8455;
	y =  -840.0492;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x =  -62.2700;
	y = -690.4383;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 968.1414;
	y =  -363.0428;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 592.2715;
	y =   125.8000;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 231.3138;
	y =  809.8071;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 547.2984;
	y =  -545.0700;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 215.9077;
	y =  557.6151;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	//*
	x = 141.7866;
	y = -321.7799;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;
	///*
	x = -557.2198;
	y = -155.5629;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 809.8675;
	y =  128.2654;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 682.2515;
	y =  547.2600;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = -717.2219;
	y =  432.8312;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	x = 588.0764;
	y =  345.4055;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 0;

	// PROLIFERATIVE CTL cell
	x = 736.2712;
	y = -90.0571;
	pCell = create_cell( *pCD );
	pCell->assign_position( x , y , 0.0 );
	pCell-> phenotype.cycle.data.current_phase_index = 1;

	return;
}
void setup_tissue_circle_immune( void )
{
	double Radius = parameters.doubles("R");

	std::vector<double> go_times_cumul(8);
	std::vector<double> persistence_times_vec(12);
	std::vector<double> speed_cumul(12);
	std::vector<double> speed_vec(12);
	go_times_cumul = go_times_cumulative();
	persistence_times_vec = persistance_distribution();
	speed_cumul = speed_cumulative();
	speed_vec = speed_distribution();

	immune_cell_placement();

	Cell* pCell = NULL;

	double x = 0.0;
	double y = 0.0;

	double GBM_NO_dense = parameters.ints("N_GBM_dense");
	double stroma_NO_dense = parameters.ints("N_stroma_dense");
	double GBM_NO_sparse = parameters.ints("N_GBM_sparse");
	double stroma_NO_sparse = parameters.ints("N_stroma_sparse");

	double fragment_radius = 1267.7;
	double area_of_fragment = 3.141*fragment_radius*fragment_radius;

	double patch_radius_1 = 317.93*parameters.doubles("kappa");
	double x_offset_patch_1 = 173.38;
	double y_offset_patch_1 = -270.5;
	double area_of_patch_1 = 3.141*patch_radius_1*patch_radius_1;
	double no_GBM_cells_patch_1 = GBM_NO_dense/area_of_fragment*area_of_patch_1;

	double patch_radius_2 = 638.2*parameters.doubles("kappa");
	double x_offset_patch_2 = 41.32;
	double y_offset_patch_2 = 1015.1;
	double area_of_patch_2 = 3.141*patch_radius_2*patch_radius_2;
	double no_GBM_cells_patch_2 = GBM_NO_dense/area_of_fragment*area_of_patch_2;

	double patch_radius_3 = (407.26+42.5)*parameters.doubles("kappa");
	double x_offset_patch_3 = -757.84;
	double y_offset_patch_3 = 128.4;
	double area_of_patch_3 = 3.141*patch_radius_3*patch_radius_3;
	double no_GBM_cells_patch_3 = GBM_NO_dense/area_of_fragment*area_of_patch_3;

	double patch_radius_4 = 307*parameters.doubles("kappa");
	double x_offset_patch_4 = 899.79;
	double y_offset_patch_4 = -391.6;
	double area_of_patch_4 = 3.141*patch_radius_4*patch_radius_4;
	double no_GBM_cells_patch_4 = GBM_NO_dense/area_of_fragment*area_of_patch_4;

	double patch_radius_5 = 407.67*parameters.doubles("kappa");
	double x_offset_patch_5 = -9.99;
	double y_offset_patch_5 = -859.9;
	double area_of_patch_5 = 3.141*patch_radius_5*patch_radius_5;
	double no_GBM_cells_patch_5 = GBM_NO_dense/area_of_fragment*area_of_patch_5;

	double no_GBM_cells_remaining = GBM_NO_sparse*(1-parameters.doubles("proportion"));
	double no_stroma_cells_remaining = stroma_NO_sparse*(1-parameters.doubles("proportion"));
	double no_GBM_cells_inpatch = GBM_NO_dense*(parameters.doubles("proportion"));
	double no_stroma_cells_inpatch = stroma_NO_dense*(parameters.doubles("proportion"));


	Cell_Definition* pCD;

	//GBM cells
	pCD = find_cell_definition("cancer");

	for( int i=0; i<no_GBM_cells_remaining; i++ )
	{
		double r_cell = sqrt(UniformRandom())*patch_radius_1;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp)+x_offset_patch_1;
		y = r_cell*sin(alp)+y_offset_patch_1;

		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))<patch_radius_1 ||
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))<patch_radius_2 ||
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))<patch_radius_3 ||
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))<patch_radius_4 ||
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))<patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}

		pCell = create_cell( *pCD );
		pCell->assign_position( x , y , 0.0 );

		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_g");//distribution2(generator);//distribution2(generator);

		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );

		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();

			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pCell->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		}
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP

		double go_stop_var = UniformRandom();

		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}

		}
	}



	// Stromal cells
	pCD = find_cell_definition("stromal");

	for( int l=0; l<no_stroma_cells_remaining; l++ )
	{
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);

		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))<patch_radius_1 ||
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))<patch_radius_2 ||
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))<patch_radius_3 ||
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))<patch_radius_4 ||
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))<patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		pCell = create_cell( *pCD );
		pCell->assign_position( x , y , 0.0 );

		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_s");//distribution2(generator);//distribution2(generator);

	}


	//GBM cells
	pCD = find_cell_definition("cancer");

	for( int i=0; i<no_GBM_cells_inpatch; i++ )
	{
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);

		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))>patch_radius_1 &&
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))>patch_radius_2 &&
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))>patch_radius_3 &&
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))>patch_radius_4 &&
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))>patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}

		pCell = create_cell( *pCD );
		pCell->assign_position( x , y , 0.0 );

		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_g");//distribution2(generator);//distribution2(generator);

		int persistence_time_index = pCell->custom_data.find_variable_index( "persistence_time" );
		int cell_motility_type_index = pCell->custom_data.find_variable_index( "cell_motility_type" );

		double p = UniformRandom();
		if(p<=0.5)// GO
		{
			pCell->custom_data.variables[cell_motility_type_index].value = 1;
			double speed_var = UniformRandom();

			for( int k=0; k<12; )
			{
				if( speed_var> speed_cumul[k] )
				{k++;}
				else
				{
					pCell->phenotype.motility.migration_speed = speed_vec[k];
					k = 12;
				}
			}
		}
		else
		{pCell->custom_data.variables[cell_motility_type_index].value = 2;} // STOP

		double go_stop_var = UniformRandom();

		for( int j=0; j<8; )
		{
			if( go_stop_var> go_times_cumul[j] )
			{j++;}
			else
			{
				pCell->custom_data.variables[persistence_time_index].value = persistence_times_vec[j];
				j = 8;
			}

		}
	}

	// Stromal cells
	pCD = find_cell_definition("stromal");

	for( int l=0; l<no_stroma_cells_inpatch; l++ )
	{
		double r_cell = sqrt(UniformRandom())*fragment_radius;
		double alp = UniformRandom()*2*3.141;
		x = r_cell*cos(alp);
		y = r_cell*sin(alp);

		while( sqrt((x-x_offset_patch_1)*(x-x_offset_patch_1)+(y-y_offset_patch_1)*(y-y_offset_patch_1))>patch_radius_1 &&
		sqrt((x-x_offset_patch_2)*(x-x_offset_patch_2)+(y-y_offset_patch_2)*(y-y_offset_patch_2))>patch_radius_2 &&
		sqrt((x-x_offset_patch_3)*(x-x_offset_patch_3)+(y-y_offset_patch_3)*(y-y_offset_patch_3))>patch_radius_3 &&
		sqrt((x-x_offset_patch_4)*(x-x_offset_patch_4)+(y-y_offset_patch_4)*(y-y_offset_patch_4))>patch_radius_4 &&
		sqrt((x-x_offset_patch_5)*(x-x_offset_patch_5)+(y-y_offset_patch_5)*(y-y_offset_patch_5))>patch_radius_5 ) //cell is located within the patch keep trying to get new position
		{
			r_cell = sqrt(UniformRandom())*fragment_radius;
			alp = UniformRandom()*2*3.141;
			x = r_cell*cos(alp);
			y = r_cell*sin(alp);
		}
		pCell = create_cell( *pCD );
		pCell->assign_position( x , y , 0.0 );

		static int virus_signal_index = microenvironment.find_density_index("virus");
		pCell->phenotype.secretion.uptake_rates[virus_signal_index] = parameters.doubles("u_s");//distribution2(generator);//distribution2(generator);

	}

	return;
}