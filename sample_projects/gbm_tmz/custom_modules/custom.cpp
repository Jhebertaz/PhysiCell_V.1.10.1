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

	static int wall_index = microenvironment.find_density_index("wall");
	static int chemokine_index = microenvironment.find_density_index("chemokine");

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );

	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	// cell_defaults.phenotype.motility.restrict_to_2D = true;

	cell_defaults.phenotype.molecular.fraction_released_at_death[wall_index] = 1;
	cell_defaults.phenotype.molecular.fraction_released_at_death[chemokine_index] = 1;

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

	Cell_Definition* pCD;

	// th cell creation
	pCD = find_cell_definition("th");

	// Make sure we're ready for 2D
	// set functions
	// pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = th_phenotype;
	// pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	// pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_CD4");


	// cancer cell creation
	pCD = find_cell_definition("cancer");

	// Make sure we're ready for 2D
	// set functions
	// pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = cancer_phenotype;
	// pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	// pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_GBM");


	// ctl cell creation
	pCD = find_cell_definition("ctl");

	// Make sure we're ready for 2D
	// set functions
	// pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = ctl_phenotype;
	// pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	// pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_CD8");


	// stroma cell creation
	pCD = find_cell_definition("stromal");


	// Make sure we're ready for 2D
	// set functions
	// pCD->functions.set_orientation = up_orientation;
	pCD->functions.update_phenotype = stromal_phenotype;
	// pCD->functions.update_velocity = custom_update_cell_velocity;

	// custom phenotype
	// pCD->phenotype.geometry.polarity = 1.0;
	pCD->phenotype.geometry.radius = parameters.doubles("R_stromal");

	// motility
	pCD->phenotype.motility.is_motile = false;


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

	// make sure not override and go back to 2D
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl;
		default_microenvironment_options.simulate_2D = true;
	}

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.


	// initialize BioFVM

	initialize_microenvironment();
	set_up_wall();

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
	// init_persistence_time for cancer_cell
	// init_speed for cancel_cell

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


// Custom distribution
double time_to_radius_increase_formula(int i)
{
  // fabulous pi number
  double pi = 3.141592653589793238462;
  // fix beta, not really appropriate with tmz treatment
  // initial tumour radius (not dynamic radius...)
  double i_r = parameters.doubles("R");
  // tumour carrying capacity
  double k_V = parameters.doubles("k_V");
  // GBM cell proliferation rate
  double beta = parameters.doubles("beta");
  double const_nume = log((1/k_V)*(4/3)*pi*i_r*i_r*i_r);
  double const_deno = log((1/k_V)*(4/3)*pi*(i_r+i*10)*(i_r+i*10)*(i_r+i*10));

  return (1/beta)*log(const_nume/const_deno);
}
std::vector<double> time_to_radius_increase() // vector for the times the migratory domain radius can increase
{
	std::vector< double > time_to_radius(62);
  for (int i=0; i<63;i++)
  {
    time_to_radius[i] = time_to_radius_increase_formula(i);
  }
  return time_to_radius;
	// time_to_radius[0] = 460;
	// time_to_radius[1] = 457;
	// time_to_radius[2] = 455;
	// time_to_radius[3] = 452;
	// time_to_radius[4] = 450;
	// time_to_radius[5] = 448;
	// time_to_radius[6] = 445;
	// time_to_radius[7] = 443;
	// time_to_radius[8] = 441;
	// time_to_radius[9] = 439;
	// time_to_radius[10] = 436;
	// time_to_radius[11] = 434;
	// time_to_radius[12] = 432;
	// time_to_radius[13] = 430;
	// time_to_radius[14] = 428;
	// time_to_radius[15] = 426;
	// time_to_radius[16] = 424;
	// time_to_radius[17] = 422;
	// time_to_radius[18] = 420;
	// time_to_radius[19] = 418;
	// time_to_radius[20] = 416;
	// time_to_radius[21] = 414;
	// time_to_radius[22] = 412;
	// time_to_radius[23] = 410;
	// time_to_radius[24] = 409;
	// time_to_radius[25] = 407;
	// time_to_radius[26] = 405;
	// time_to_radius[27] = 403;
	// time_to_radius[28] = 401;
	// time_to_radius[29] = 400;
	// time_to_radius[30] = 398;
	// time_to_radius[31] = 396;
	// time_to_radius[32] = 394;
	// time_to_radius[33] = 393;
	// time_to_radius[34] = 391;
	// time_to_radius[35] = 390;
	// time_to_radius[36] = 388;
	// time_to_radius[37] = 386;
	// time_to_radius[38] = 385;
	// time_to_radius[39] = 383;
	// time_to_radius[40] = 382;
	// time_to_radius[41] = 380;
	// time_to_radius[42] = 379;
	// time_to_radius[43] = 377;
	// time_to_radius[44] = 376;
	// time_to_radius[45] = 374;
	// time_to_radius[46] = 373;
	// time_to_radius[47] = 371;
	// time_to_radius[48] = 370;
	// time_to_radius[49] = 369;
	// time_to_radius[50] = 367;
	// time_to_radius[51] = 366;
	// time_to_radius[52] = 364;
	// time_to_radius[53] = 363;
	// time_to_radius[54] = 362;
	// time_to_radius[55] = 360;
	// time_to_radius[56] = 359;
	// time_to_radius[57] = 358;
	// time_to_radius[58] = 357;
	// time_to_radius[59] = 355;
	// time_to_radius[60] = 354;
	// time_to_radius[61] = 10000000;
	// return time_to_radius;
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



// custom functions
bool am_i_dead(Cell* pCell) 																										// done
{
	bool dead = pCell->phenotype.death.dead;
	return dead;
}
bool is_someone_around(Cell* pCell) 																						// done
{
	// return pCell->nearby_cells().size()>0;
	return pCell->cells_in_my_container().size()>0;
}

// pressure functions
double equilibrium_spacing()
{
	// tumour volume carrying capacity
	double k_C = parameters.doubles("K_C");

	// tumour surface area
	double a_frag = parameters.doubles("A_frag");

	// equilibrium spacing between cells
	double s = sqrt( (2*a_frag) / (sqrt(3)*k_C) );
	return s;
}
double pressure_scale(Cell* pCell)
{
	//cell-cell repulsion force coefficient
	double c_j_ccr = pCell->phenotype.mechanics.cell_cell_repulsion_strength;//10

	// normal pi
	double pi_ = 3.141592653589793238462;

	// cell radius
	double radius =  pCell->phenotype.geometry.radius;

	// cell surface area
	double a_cell = pi_*radius*radius;

	// pressure scale
	double p_s = c_j_ccr/a_cell;
	return p_s;
}
double maximal_pressure(Cell* pCell)
{
	// equilibrium spacing between cells
	static double s = equilibrium_spacing();

	// cell radius
	double radius =  pCell->phenotype.geometry.radius;

	// pressure scale
	double p_s = pressure_scale(pCell);

	// maximal pressure
	double m_p = 6*p_s*(1-(1/(2*radius))*s)*(1-(1/(2*radius))*s);
	return m_p;
}

// tmz effect
void cancer_tmz_effect(Cell* pCell, Phenotype& phenotype)
{
	static int tmz_index = microenvironment.find_density_index("tmz");
	static int apoptosis_model_index = phenotype.death.find_death_model_index("apoptosis");

	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);

	double ec50 = parameters.doubles("IC50");

	// GBM cell default proliferation rate
	double beta = parameters.doubles("beta");
	double tmz_amount = pCell->nearest_density_vector()[tmz_index];
	double tmz_effect = tmz_amount/(tmz_amount+ec50);

  phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = beta*(1-tmz_effect);

	//need to add it when including TMZ
	phenotype.death.rates[apoptosis_model_index] = beta*(tmz_effect);
	return;
}

// for cancer cells
double resample_movement_type() 									                              // done
{
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
				speed = speed_vec[k];
				return speed;
			}
		}
	}
	else
	{
		// assign STOP type
		// pCell->custom_data["cell_motility_type"] = 2;
		// phenotype.motility.migration_speed = 0;
		speed = 0.0;
		return speed;
	}
	return speed;
}
double resample_persistence_time() 																							// done
{
	double go_stop_var = UniformRandom();
	double time;
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
			time = persistence_times_vec[j]+PhysiCell_globals.current_time;
			return time;
			// j = 8;
		}
	}
	return time;
}



bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt )
{
	double attached_time = 	pAttacker->custom_data["attachment_lifetime"];
	double real_time = PhysiCell_globals.current_time;

	// CTL kills cell if it has been attached for enough time
	if( attached_time < real_time && pTarget->type==2)
	{
		return true;
	}
	return false;
}
void add_elastic_velocity(Cell* pActingOn, Cell* pAttachedTo , double elastic_constant) // done
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position;
	axpy( &(pActingOn->velocity), elastic_constant , displacement );

	return;
}
void extra_elastic_attachment_mechanics(Cell* pCell) 														// done
{
	for (int i=0; i < pCell->state.neighbors.size() ; i++)
	{
		add_elastic_velocity(pCell, pCell->state.neighbors[i], parameters.doubles("c_s"));
	}
	return;
}



// Custom cell movement
void th_movement(Cell* pCell, Phenotype& phenotype)
{
	static int wall_index = microenvironment.find_density_index("wall");

	double wall_amount = pCell->nearest_density_vector()[wall_index];
	double nu = parameters.doubles("nu");
	std::vector<double> ae_ini(3);

	if( wall_amount<2 )// Make TH cells that have left the diameter of the tumour turn around
	{
		phenotype.motility.migration_speed = nu;
		ae_ini = -1*pCell->position;
		phenotype.motility.migration_bias = 1;
		normalize( &( ae_ini ) );
		phenotype.motility.migration_bias_direction = ae_ini;
	}
	else
	{
		phenotype.motility.migration_speed = nu;
		phenotype.motility.migration_bias = 0;
	}
	return;
}
void cancer_movement(Cell* pCell, Phenotype& phenotype) 												// done
{
	//cell movement
	static int wall_index = microenvironment.find_density_index("wall");
	double wall_amount = pCell->nearest_density_vector()[wall_index];
	double persistence_time = phenotype.motility.persistence_time;

	if( persistence_time <= PhysiCell_globals.current_time ) // if the cell's persistence time is up
	{
		phenotype.motility.migration_speed = resample_movement_type(); // assign migration speed
		// std::cout<<phenotype.motility.migration_speed<<std::endl;
		if (phenotype.motility.migration_speed > 0.0)
		{
			phenotype.motility.is_motile = true;
		}
		else
		{
			phenotype.motility.is_motile = false;
		}
		phenotype.motility.persistence_time = resample_persistence_time();// assign persistence time - needs to be a real time!
	}
	// to be sure that the cell don't go outside
	if(wall_amount<2)
	{
		phenotype.motility.migration_speed = 0.0;
		phenotype.motility.is_motile = false;
	}
	return;
}
void ctl_movement(Cell* pCell, Phenotype& phenotype)
{
	static int chemokine_index = microenvironment.find_density_index("chemokine");
	static int wall_index = microenvironment.find_density_index("wall");

	double b;
	double psi;

	double nu = parameters.doubles("nu");
	double nu_max = parameters.doubles("nu_max");
	double nu_star = parameters.doubles("nu_star");
	double b_CD8 = parameters.doubles("b_CD8");

	double chemokine_amount = pCell->nearest_density_vector()[chemokine_index];
	std::vector<double> chemokine_gradient = pCell->nearest_gradient(chemokine_index);

	double wall_amount = pCell->nearest_density_vector()[wall_index];
	std::vector<double> dbias(3);


	// verify if cell go outside of the tumour_radius

	// Extra condition (chemokine) for ctl only

	if(chemokine_amount>1e-8)// sample chemotaxis gradient and random walk in that direction
	{
	  dbias = chemokine_gradient;
	  b = b_CD8;
	  psi = nu+(nu_max-nu)*(chemokine_amount/(nu_star/2+chemokine_amount));
		phenotype.motility.migration_bias_direction = dbias;
	}
	else if(chemokine_amount>1e-3)
	{
	  dbias = chemokine_gradient;
	  b = b_CD8;
	  psi = nu+(nu_max-nu)*(chemokine_amount/(nu_star/2+chemokine_amount));
		phenotype.motility.migration_bias_direction = dbias;
	}
	else
	{
	  b = 0.0;
	  psi = nu;
	}

	if (wall_amount<2)
	{
		// Make TH/CTL cells that have left the diameter of the tumour turn around
		dbias = -1*pCell->position;
		normalize( &( dbias ) );

		// dbias = ae_ini;
		b = b_CD8;
		psi = nu;
		phenotype.motility.migration_bias_direction = dbias;
	}
	//
	phenotype.motility.migration_bias = b;
	phenotype.motility.migration_speed = psi;
	return;
}
void stromal_movement(Cell* pCell, Phenotype& phenotype) 												// done
{
	//stromal doesn't move
	return;
}

// Custom phenotype update fonctions
void th_phenotype(Cell* pCell, Phenotype& phenotype, double dt) 								// done
{
	if (am_i_dead(pCell))
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	th_movement(pCell, phenotype);
	return;
}
void cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	// TMZ effect before pressure
	cancer_tmz_effect(pCell, phenotype);

	static int wall_index = microenvironment.find_density_index("wall");

	int cycle_start_index = live.find_phase_index(PhysiCell_constants::live);
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);

	// maximal pressure
	double m_p = maximal_pressure(pCell);
	// gbm proliferation rate
	double beta = parameters.doubles("beta");
	// cell pressure
	double cell_pressure = pCell->state.simple_pressure;

	// cell attemps proliferation
	if (cell_pressure*cell_pressure > m_p)
	{
		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;
	}
	else
	{
		phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = beta;
	}

	cancer_movement(pCell, phenotype);
	return;
}
void ctl_phenotype(Cell* pCell, Phenotype& phenotype, double dt)
{
	if (am_i_dead(pCell))
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	ctl_movement(pCell, phenotype);
	return;
}
void stromal_phenotype(Cell* pCell, Phenotype& phenotype, double dt) 						// done
{
	if (am_i_dead(pCell))
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	stromal_movement(pCell, phenotype);
	return;
}




void set_up_wall()
{
	static int wall_index = microenvironment.find_density_index("wall");

	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center;

		if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")+10)*(parameters.doubles("R")+10))
		{
			microenvironment(n)[wall_index] = 1;
		}
		else if(ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>(parameters.doubles("R")-10)*(parameters.doubles("R")-10))
		{
			microenvironment(n)[wall_index] = 3.5;
		}
		else if( ECMdense[0]*ECMdense[0]+ECMdense[1]*ECMdense[1]>250*250 )
		{
			microenvironment(n)[wall_index] = 5;
		}
		else
		{
			microenvironment(n)[wall_index] = 10;
		}
	}
	return;
}
