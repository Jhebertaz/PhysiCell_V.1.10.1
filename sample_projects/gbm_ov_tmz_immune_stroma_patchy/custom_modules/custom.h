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

#include "../core/PhysiCell.h"
#include "../core/PhysiCell_cell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

// setup functions to help us along

void create_cell_types( void );
void setup_tissue( void );

// set up the BioFVM microenvironment
void setup_microenvironment( void );

// custom pathology coloring function
std::vector<std::string> my_coloring_function(Cell*);
std::vector<std::string> colouring_by_intracellular_virus_amount(Cell* pCell);

// custom from Source
void custom_update_cell_velocity(Cell* pCell, Phenotype& phenotype, double dt);

// custom distribution
std::vector<double> time_to_radius_increase(void);
std::vector<double> speed_distribution(void);
std::vector<double> go_times_cumulative(void);
std::vector<double> persistance_distribution(void);
std::vector<double> speed_cumulative(void);

// for main
double CSF_vals(int indexing_CSF_Vec);
double CSF_conc_to_density(double time);

void testing();

// old
void old_testing();
void old_immune_cell_placement();
void old_setup_tissue_circle_immune();
bool immune_cell_attempt_apoptosis( Cell* pAttacker, Cell* pTarget, double dt );// old one


// custom functions can go here
bool am_i_dead(Cell* pCell);
bool am_i_infected(Cell* pCell);
bool am_i_an_detectable_infected_cell(Cell* pCell);
Cell* is_there_an_infected_cell_around(Cell* pCell);
Cell* am_i_attached_to_an_infected_cell(Cell* pCell);
bool is_someone_around(Cell* pCell);

double resample_movement_type();
double resample_persistence_time();

// ctl extra functions
bool immune_cell_trigger_apoptosis(Cell* pAttacker, Cell* pTarget);
Cell* immune_cell_check_neighbors_for_attachment(Cell* pAttacker);
bool immune_cell_attempt_attachment( Cell* pAttacker, Cell* pTarget);
void ctl_virus(Cell* pCell, Phenotype& phenotype);
void add_elastic_velocity(Cell* pActingOn, Cell* pAttachedTo, double elastic_constant);
void extra_elastic_attachment_mechanics(Cell* pCell);

// pressure functions
double equilibrium_spacing();
double pressure_scale(Cell* pCell);
double maximal_pressure(Cell* pCell);

// tmz effect
void cancer_tmz_effect(Cell* pCell, Phenotype& phenotype);



// virus infection dynamics
double update_virus_uptake_rate(Cell* pCell, Phenotype& phenotype);
void virus_replication(Cell* pCell, Phenotype& phenotype, double dt);
void virus_induced_lysis(Cell* pCell, Phenotype& phenotype, double dt);
void checked_for_lysis(Cell* pCell, Phenotype& phenotype, double dt);


// custom movement functions
void th_movement(Cell* pCell, Phenotype& phenotype);
void cancer_movement(Cell* pCell, Phenotype& phenotype);
void ctl_movement(Cell* pCell, Phenotype& phenotype);
void stromal_movement(Cell* pCell, Phenotype& phenotype);

// custom infection dynamics functions
void th_virus_infection_dynamics(Cell* pCell, Phenotype& phenotype, double dt);
void cancer_virus_infection_dynamics(Cell* pCell, Phenotype& phenotype, double dt);
void ctl_virus_infection_dynamics(Cell* pCell, Phenotype& phenotype);
void stromal_virus_infection_dynamics(Cell* pCell, Phenotype& phenotype, double dt);

// custom phenotype functions
void th_phenotype(Cell* pCell, Phenotype& phenotype, double dt);
void cancer_phenotype(Cell* pCell, Phenotype& phenotype, double dt);
void ctl_phenotype(Cell* pCell, Phenotype& phenotype, double dt);
void stromal_phenotype(Cell* pCell, Phenotype& phenotype, double dt);

// default phenotype functions
void phenotype_function(Cell* pCell, Phenotype& phenotype, double dt);
void custom_function(Cell* pCell, Phenotype& phenotype , double dt);
void contact_function(Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt);
