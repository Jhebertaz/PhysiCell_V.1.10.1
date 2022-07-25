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
	// cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );

	// Make sure we're ready for 2D
	// cell_defaults.functions.set_orientation = up_orientation;
	// cell_defaults.phenotype.geometry.polarity = 1.0;
	// cell_defaults.phenotype.motility.restrict_to_2D = true;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;//standard_update_cell_velocity;

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

	Cell_Definition* pCD; // to use less memory

	// th cell creation
	pCD = find_cell_definition("th");

	// pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	// Make sure we're ready for 2D
	// pCD->functions.set_orientation = up_orientation;
	// pCD->phenotype.geometry.polarity = 1.0;
	// pCD->phenotype.motility.restrict_to_2D = true;
	// set functions
	pCD->functions.update_phenotype = th_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// cancer cell creation
	pCD = find_cell_definition("cancer");

	// pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	// Make sure we're ready for 2D
	// pCD->functions.set_orientation = up_orientation;
	// pCD->phenotype.geometry.polarity = 1.0;
	// pCD->phenotype.motility.restrict_to_2D = true;
	// set functions
	pCD->functions.update_phenotype = cancer_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// ctl cell creation
	pCD = find_cell_definition("ctl");

	// pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	// Make sure we're ready for 2D
	// pCD->functions.set_orientation = up_orientation;
	// pCD->phenotype.geometry.polarity = 1.0;
	// pCD->phenotype.motility.restrict_to_2D = true;
	// set functions
	pCD->functions.update_phenotype = ctl_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

	// stroma cell creation
	pCD = find_cell_definition("stromal");

	// pCD->phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// pCD->phenotype.molecular.sync_to_microenvironment( &microenvironment );
	// Make sure we're ready for 2D
	// pCD->functions.set_orientation = up_orientation;
	// pCD->phenotype.geometry.polarity = 1.0;
	// pCD->phenotype.motility.restrict_to_2D = true;
	// set functions
	pCD->functions.update_phenotype = stromal_phenotype;
	pCD->functions.update_velocity = custom_update_cell_velocity;

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
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	// create some of each type of cell

	Cell* pC;

	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k];
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0};
			position[0] = Xmin + UniformRandom()*Xrange;
			position[1] = Ymin + UniformRandom()*Yrange;
			position[2] = Zmin + UniformRandom()*Zrange;

			pC = create_cell( *pCD );
			pC->assign_position( position );
		}
	}
	std::cout << std::endl;

	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();

	return;
}
std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }


void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }



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
		if( pCell->type == 3)
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
		if( pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if( n_I>1)
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


// custom functions
bool am_i_dead(Cell* pCell) 																										// done
{
	bool dead = pCell->phenotype.death.dead;
	if (dead)
	{
		pCell->functions.update_phenotype = NULL;
	}
	return dead;
}
bool am_i_infected(Cell* pCell) 																								// done
{
	static int virus_signal_index = microenvironment.find_density_index("virus");

	double infection_threshold_half_effect = parameters.doubles("m_half");
	int cell_intern_virions = pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index];

	bool infected = false;

	if (cell_intern_virions > infection_threshold_half_effect)
	{
		infected = true;
	}
	return infected;
}
bool is_there_an_infected_cell_around(Cell* pCell) 															// done
{
	// neighborhood
	std::vector<Cell*> nearby = pCell->cells_in_my_container();

	// loop variables
	Cell* pC = NULL;
	bool stop = false;
	int i=0;

	while (!stop && i<nearby.size())
	{
		pC = nearby[i]; // ith cell object

		// cancer cell
		if (pC->ID==2)
		{
			bool dead = pC->phenotype.death.dead;
			bool itself = (pC == pCell);
			bool infected_cell = am_i_infected(pC);

			// infected cell around
			if (infected_cell && !dead && !itself)
			{
				stop = true;
			}
		}
		i++;
	}
	return stop;
}
bool am_i_attached(Cell* pCell)
{
	return pCell->state.attached_cells.size()>0;
}

// for cancer cells
void resample_movement_type(Cell* pCell) 																				// done
{
	double new_type_rand = UniformRandom();

	std::vector<double> speed_cumul(12);
	speed_cumul = speed_cumulative();

	std::vector<double> speed_vec(12);
	speed_vec = speed_distribution();

	if(new_type_rand<=0.5)
	{
		// assign GO type
		pCell->custom_data["cell_motility_type"] = 1;
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
				pCell->phenotype.motility.migration_speed = speed_vec[k];
				k = 12;
			}
		}
	}
	else
	{
		// assign STOP type
		pCell->custom_data["cell_motility_type"] = 2;
		pCell->phenotype.motility.migration_speed = 0;
	}
	return;
}
void resample_persistence_time(Cell* pCell) 																		// done
{
	std::vector<double> go_times_cumul(8);
	go_times_cumul = go_times_cumulative();

	std::vector<double> persistence_times_vec(12);
	persistence_times_vec = persistance_distribution();

	// assign persistence time - needs to be a real time!
	double go_stop_var = UniformRandom();
	for (int j=0; j<8;)
	{
		if (go_stop_var> go_times_cumul[j])
		{
			j++;
		}
		else
		{
			// assign persist time
			pCell->phenotype.motility.persistence_time = persistence_times_vec[j]+PhysiCell_globals.current_time;
			j = 8;
		}
	}
	return;
}

// for ctl cells
bool immune_cell_trigger_apoptosis(Cell* pAttacker, Cell* pTarget)
{
	static int apoptosis_model_index = pTarget->phenotype.death.find_death_model_index("apoptosis");
	static int virus_index = microenvironment.find_density_index("virus");

	// if the Target cell is already dead, don't bother!
	if (am_i_dead(pTarget))
	{
		return false;
	}
	else
	{
		std::cout<< "CTL successfully trigger apoptosis"<<std::endl;
		pTarget->start_death(apoptosis_model_index);
		pTarget->phenotype.molecular.fraction_released_at_death[virus_index] = 0;
		return true;
	}
}
Cell* immune_cell_check_neighbors_for_attachment(Cell* pAttacker)
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while (i < nearby.size())
	{
		// don't try to kill yourself
		if (nearby[i] != pAttacker)
		{
			// conditions to be a candidate
			// be infected
			bool is_candidate_infected = am_i_infected(nearby[i]);
			// be alive
			bool is_candidate_dead = am_i_dead(nearby[i]);
			// be a cancer cell
			bool is_candidate_gbm = (nearby[i]->ID == 2);
			// be close
			std::vector<double> displacement = (nearby[i]->position) - (pAttacker->position);
			bool is_candidate_close = norm(displacement) > parameters.doubles("d_attach");

			if (is_candidate_infected && is_candidate_dead && is_candidate_gbm &&is_candidate_close)
			{
				return nearby[i];
			}
		}
		i++;
	}
	return NULL;
}
void CTL_movement(Cell* pCell, Phenotype& phenotype) 														// done
{
	static int chemokine_index = microenvironment.find_density_index("chemokine");
	static int wall_index = microenvironment.find_density_index("wall");

	double wall_amount = pCell->nearest_density_vector()[wall_index];
	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];

	std::vector<double> chemokine_gradient = pCell->nearest_gradient(chemokine_index);
	std::vector<double> wall_gradient = pCell->nearest_gradient(wall_index);


	int cycle_start_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_negative);
	int cycle_end_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_positive);

	double nu = parameters.doubles("nu");
	double nu_max = parameters.doubles("nu_max");
	double nu_star = parameters.doubles("nu_star");
	double b_CD8 = parameters.doubles("b_CD8");
	double r_01_CD8 = parameters.doubles("r_01_CD8");
	double prolif_increase = parameters.doubles("CTL_prolif_increase_due_to_stimulus");

	std::vector<double> ae_ini(3);// ??
	// TH movement
	if (wall_amount < 2)
	{
		//std::cout<<"outside domain "<<phenotype.motility.migration_speed<<std::endl;

		// Outside

		ae_ini = -1*pCell->position;
		normalize( &( ae_ini ) );

		phenotype.motility.migration_bias = 1;
		phenotype.motility.migration_bias_direction = wall_gradient;
		phenotype.motility.migration_bias_direction = ae_ini;

		return;
	}
	else if (chemokine_amount > 1e-8) // sample chemotaxis gradient and random walk in that direction
	{
		// Presence of chemokine (low)
		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = r_01_CD8;
		phenotype.motility.migration_speed = nu+(nu_max-nu)*(chemokine_amount/(nu_star+chemokine_amount));
		phenotype.motility.migration_bias = b_CD8;
		phenotype.motility.migration_bias_direction = chemokine_gradient;

		return;
	}
	else if (chemokine_amount > 1e-3)
	{
		// Presence of chemokine (high)
		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = r_01_CD8*prolif_increase;
		phenotype.motility.migration_speed = nu+(nu_max-nu)*(chemokine_amount/(nu_star+chemokine_amount));
		phenotype.motility.migration_bias = b_CD8;
		phenotype.motility.migration_bias_direction = chemokine_gradient;

		return;
	}
	else
	{
		// normal
		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = r_01_CD8;
		phenotype.motility.migration_speed = nu;
		phenotype.motility.migration_bias = 0;

		return;
	}
}
void add_elastic_velocity(Cell* pActingOn, Cell* pAttachedTo , double elastic_constant)
{
	std::vector<double> displacement = (pAttachedTo->position) - (pActingOn->position);
	axpy( &(pActingOn->velocity), elastic_constant , displacement );

	return;
}
void extra_elastic_attachment_mechanics(Cell* pCell)
{
	for (int i=0; i < pCell->state.neighbors.size() ; i++)
	{
		add_elastic_velocity(pCell, pCell->state.neighbors[i], parameters.doubles("c_s"));
	}
	return;
}


// infection dynamics
void update_virus_uptake_rate(Cell* pCell,  Phenotype& phenotype) 							// done
{
	static int virus_signal_index = microenvironment.find_density_index("virus");

	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];

	double p = pCell->nearest_density_vector()[virus_signal_index];
	double pmax = parameters.doubles("rho_max");

	double u = parameters.doubles("u_g");
	double Vvoxel = microenvironment.mesh.voxels[0].volume;
	double m_half = parameters.doubles("m_half");

	if (p < pmax)
	{
		phenotype.secretion.uptake_rates[virus_signal_index] = u*p/(cell_intern_virions/Vvoxel+m_half/Vvoxel);
	}
	else
	{
		phenotype.secretion.uptake_rates[virus_signal_index] = u*pmax*pmax/(cell_intern_virions/Vvoxel+m_half/Vvoxel)/p;
	}
	return;
}
void virus_replication(Cell* pCell, Phenotype& phenotype, double dt) 						// done
{
	static int virus_signal_index = microenvironment.find_density_index("virus");
	static int apoptosis_model_index = phenotype.death.find_death_model_index("apoptosis");

	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];
	double alpha = parameters.doubles("alpha");
	double gamma = parameters.doubles("gamma");

	if (cell_intern_virions > 1 && cell_intern_virions <= alpha)
	{
		phenotype.molecular.internalized_total_substrates[virus_signal_index] = cell_intern_virions+dt*(gamma*cell_intern_virions);
	}
	else if (cell_intern_virions > alpha-1)
	{
		phenotype.molecular.fraction_released_at_death[virus_signal_index] = 1;
		phenotype.secretion.uptake_rates[virus_signal_index] = 0;
		pCell->start_death(apoptosis_model_index);
	}
	else if (cell_intern_virions<0)
	{
		std::cout<<"Negative intracellular virus!! m_i = "<<cell_intern_virions<<std::endl;
	}
	return;
}
void virus_induced_lysis(Cell* pCell, Phenotype& phenotype, double dt) 					// done
{
	static int virus_signal_index = microenvironment.find_density_index( "virus");

	double pstar = phenotype.secretion.saturation_densities[virus_signal_index];
	double delta_V = parameters.doubles("delta_V");//0.1466;
	double Vvoxel = microenvironment.mesh.voxels[0].volume;//volume of voxel
	double p = pCell->nearest_density_vector()[virus_signal_index];
	double cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];

	if (cell_intern_virions > 1)
	{
			double amount_to_add = (cell_intern_virions-cell_intern_virions*exp(-delta_V*dt))/Vvoxel;
			if (amount_to_add > pstar-p)
			{
				pCell->nearest_density_vector()[virus_signal_index] += (pstar-p)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] -= (pstar-p)*Vvoxel*dt;
			}
			else
			{
				pCell->nearest_density_vector()[virus_signal_index] += (amount_to_add)*dt;
				pCell->phenotype.molecular.internalized_total_substrates[virus_signal_index] = cell_intern_virions*exp(-delta_V*dt);
			}
	}

	return;
}
void checked_for_lysis(Cell* pCell, Phenotype& phenotype, double dt) 						// done
{
	static int virus_signal_index = microenvironment.find_density_index("virus");

	if (phenotype.molecular.fraction_released_at_death[virus_signal_index]>0)
	{
		virus_induced_lysis(pCell, phenotype, dt);
	}
	return;
}



// custom from Source
// custom version
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

// Custom phenotype update fonctions
void th_phenotype(Cell* pCell, Phenotype& phenotype, double dt) 								// done
{
	if (am_i_dead(pCell))
	{
		return;
	}

	static int chemokine_index = microenvironment.find_density_index("chemokine");

	int cycle_start_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_negative);
	int cycle_end_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_positive);

	// Is there an infected cell around
	bool infected_neighbor = is_there_an_infected_cell_around(pCell);
	bool activated = (phenotype.secretion.secretion_rates[chemokine_index] == parameters.doubles("S_chemokine_CD4"));

	if (infected_neighbor == true || activated)
	{
		// there is an infected cells
		// need to be activated
		std::cout<<"TH secretion increased"<<std::endl;
		phenotype.secretion.secretion_rates[chemokine_index] = parameters.doubles("S_chemokine_CD4");
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4")*parameters.doubles("TH_prolif_increase_due_to_stimulus");
		phenotype.motility.migration_speed = 0.1;
	}
	else
	{
		// there is no infected cells and not already activated
		// keep it desactivated
		phenotype.secretion.secretion_rates[chemokine_index] = 0;
		phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = parameters.doubles("r_01_CD4");
	}
	return;
}

void cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	// Has my motility persistence time expired
	double persistence_time = phenotype.motility.persistence_time;

	if (persistence_time <= PhysiCell_globals.current_time)
	{
		resample_movement_type(pCell);
		resample_persistence_time(pCell);
	}

	if (!am_i_dead(pCell))
	{
		update_virus_uptake_rate(pCell, phenotype);
		if (am_i_infected(pCell))
		{
			virus_replication(pCell, phenotype, dt);
		}
		else
		{
			// cell attemps proliferation
			int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
			int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);

			// pressure around
			double radius = parameters.doubles("R_GBM");
			double s = parameters.doubles("xhi");

			double pressure = 6*(1-1/(2*radius)*s)*(1-1/(2*radius)*s);
			double pressure_scale = 0.027288820670331;
			double maximal_pressure = pressure/pressure_scale;

			if (pCell->state.simple_pressure*pCell->state.simple_pressure >= maximal_pressure)
			{
				phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;
			}
			else
			{
				phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = parameters.doubles("beta");
			}
		}
	}
	else
	{
		checked_for_lysis(pCell, phenotype, dt);
		return;
	}


	//
	// if (am_i_infected(pCell))
	// {
	// 	// Intracellelar virus replicates
	// 	if (!am_i_dead(pCell))
	// 	{
	// 		update_virus_uptake_rate(pCell, phenotype);
	// 		virus_replication(pCell, phenotype, dt);
	// 	}
	// }
	// // cell is checked for lysis
	// if (am_i_dead(pCell))
	// {
	// 	checked_for_lysis(pCell, phenotype, dt);
	// 	return;
	// }
	// else
	// {
	// 	// cell attemps proliferation
	// 	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	// 	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
	//
	// 	// pressure around
	// 	double radius = parameters.doubles("R_GBM");
	// 	double s = parameters.doubles("xhi");
	//
	// 	double pressure = 6*(1-1/(2*radius)*s)*(1-1/(2*radius)*s);
	// 	double pressure_scale = 0.027288820670331;
	// 	double maximal_pressure = pressure/pressure_scale;
	//
	// 	if (pCell->state.simple_pressure*pCell->state.simple_pressure >= maximal_pressure)
	// 	{
	// 		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;
	// 	}
	// 	else
	// 	{
	// 		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = parameters.doubles("beta");
	// 	}
	//
	// 	//checks whether it becomes infected
	//
	// }
	return;
}

void ctl_phenotype(Cell* pCell, Phenotype& phenotype, double dt) 								// done
{
	if (am_i_dead(pCell))
	{
		return;
	}

	int cycle_start_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_negative);
	int cycle_end_index = Ki67_basic.find_phase_index(PhysiCell_constants::Ki67_positive);

	// Verify if it is attached to an infected cell.
	bool is_attached = am_i_attached(pCell);
	double tau = parameters.doubles("tau");

	// ???
	if( pCell->state.neighbors.size() > 0 )
	{
		extra_elastic_attachment_mechanics(pCell);
	}

	if (is_attached)
	{
		// Attached

		//How long it was attached
		double attach_time = pCell->custom_data["attachment_lifetime"];
		Cell* attached_cell = pCell->state.neighbors[0]; // ok why ?

		if (attach_time >= tau)
		{
			//CTL kills attached cell
			bool success = immune_cell_trigger_apoptosis(pCell, attached_cell);

			if (success)
			{
				detach_cells(pCell, attached_cell);

				phenotype.motility.is_motile = true;
				CTL_movement(pCell, phenotype);

				//increases proliferation
				double proliferation_rate = parameters.doubles("r_01_CD8")*parameters.doubles("CTL_prolif_increase_due_to_stimulus");
				phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = proliferation_rate;
				return;
			}
		}
		else
		{
			// There is nothing else to do
			return;
		}
	}
	else
	{
		// Not attached
		// search for other infected cells around
		bool infected_neighbor = is_there_an_infected_cell_around(pCell);

		if (infected_neighbor)
		{
			//Infected cell around
			//Attemp attachment
			Cell* candidate = immune_cell_check_neighbors_for_attachment(pCell);

			if (candidate != NULL)
			{
				std::cout<<"CTL attached"<<std::endl;
				attach_cells(pCell, candidate);
				pCell->custom_data["attachment_lifetime"] = PhysiCell_globals.current_time + tau;
				phenotype.motility.is_motile = false;
				return;
			}
			else
			{
				// There is nothing else to do
				phenotype.motility.is_motile = true;
				CTL_movement(pCell, phenotype);
				return;
			}
		}
		else
		{
			// There is nothing else to do
			phenotype.motility.is_motile = true;
			CTL_movement(pCell, phenotype);
			return;
		}
	}
	return;
}

void stromal_phenotype(Cell* pCell, Phenotype& phenotype, double dt) 						// done
{
	if (am_i_dead(pCell))
	{
		return;
	}

	static int virus_signal_index = microenvironment.find_density_index("virus");

	int cell_intern_virions = phenotype.molecular.internalized_total_substrates[virus_signal_index];
	// double infection_threshold_half_effect = parameters.doubles("m_half");

	// weird but ok
	if (cell_intern_virions > 1e5)
	{
		phenotype.molecular.internalized_total_substrates[virus_signal_index] = 0;
	}

	return;
	// if infected
	// if (cell_intern_virions > infection_threshold_half_effect)
	// {
	// 	// Do nothing
	// 	// Question4: virus_induced_lysis ?
	// 	//
	// }
	// else
	// {
	// 	// Not enought to induced cell lysis
	// 	// Call infection dynamics ?
	// }
}

void testing()
{
	// setup initial concentration of "wall", virus

	// for testing only
	static int virus_index = microenvironment.find_density_index("virus");
	static int wall_index = microenvironment.find_density_index("wall");
	static int chemokine_index = microenvironment.find_density_index("chemokine");

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

	// location of virus
	double x_virus = 58;
	double y_virus = 17;
	double zone_radius = 20;
	// double tumour_radius = parameters.doubles("R");


	for (int n = 0 ; n < microenvironment.mesh.voxels.size(); n++)
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if((x_virus-ECMdense[0])*(x_virus-ECMdense[0])+(y_virus-ECMdense[1])*(y_virus-ECMdense[1])<zone_radius*zone_radius)
		{
			microenvironment(n)[virus_index] = 7.12;
		}
	}


	// location of virus
	x_virus = -58;
	y_virus = -17;
	zone_radius = 20;
	for (int n = 0 ; n < microenvironment.mesh.voxels.size(); n++)
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if((x_virus-ECMdense[0])*(x_virus-ECMdense[0])+(y_virus-ECMdense[1])*(y_virus-ECMdense[1])<zone_radius*zone_radius)
		{
			microenvironment(n)[virus_index] = 7.12;
		}
	}

	// location of virus
	x_virus = -100;
	y_virus = -100;
	zone_radius = 50;
	for (int n = 0 ; n < microenvironment.mesh.voxels.size(); n++)
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if((x_virus-ECMdense[0])*(x_virus-ECMdense[0])+(y_virus-ECMdense[1])*(y_virus-ECMdense[1])<zone_radius*zone_radius)
		{
			microenvironment(n)[virus_index] = 7.12;
		}
	}


	// location of virus
	x_virus = -200;
	y_virus = 100;
	zone_radius = 50;
	for (int n = 0 ; n < microenvironment.mesh.voxels.size(); n++)
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if((x_virus-ECMdense[0])*(x_virus-ECMdense[0])+(y_virus-ECMdense[1])*(y_virus-ECMdense[1])<zone_radius*zone_radius)
		{
			microenvironment(n)[virus_index] = 7.12;
		}
	}

	return;
}
