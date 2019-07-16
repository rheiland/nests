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
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	cell_defaults.functions.update_phenotype = NULL; 
	
	// needed for a 2-D simulation: 
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 
	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	
	int signal_index = microenvironment.find_density_index( "signal" ); 

	int start_index = live.find_phase_index( PhysiCell_constants::live );
	int end_index = live.find_phase_index( PhysiCell_constants::live );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[signal_index] = 0; 
	cell_defaults.phenotype.secretion.secretion_rates[signal_index] = 
		parameters.doubles("signal_secretion"); 
	cell_defaults.phenotype.secretion.saturation_densities[signal_index] = 1; 
	
	// turn off birth and death 
	cell_defaults.phenotype.cycle.data.transition_rate(start_index,end_index) = 0.0; 
	cell_defaults.phenotype.death.rates[ apoptosis_model_index ] = 0.0; 
	cell_defaults.phenotype.death.rates[ necrosis_model_index ] = 0.0; 
	
	// now, set motilty parameters 
	cell_defaults.phenotype.motility.is_motile = true; 
	cell_defaults.phenotype.motility.persistence_time = parameters.doubles("persistence_time"); 
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("migration_speed"); 
	cell_defaults.phenotype.motility.migration_bias_direction ={ 1.0, 0.0, 0.0 };  
	cell_defaults.phenotype.motility.migration_bias = parameters.doubles("migration_bias" ); 
	
	cell_defaults.phenotype.motility.fixed_persistence = parameters.bools("fixed_persistence"); 
	
	// turn off adhesion and repulsoin os tracks are independent 
	cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0; 
	cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0; 

	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// no gradients need for this example 
	default_microenvironment_options.calculate_gradients = false; 
	default_microenvironment_options.use_oxygen_as_first_field = false;
	
	microenvironment.set_density( 0, "signal" , "dimensionless" , 0.0 , 0.0 ); 
	
	// set Dirichlet conditions 
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector = { 0.0 }; // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	int nXnodes = microenvironment.mesh.x_coordinates.size(); 
	int nYnodes = microenvironment.mesh.y_coordinates.size(); 
	
	double y_min = microenvironment.mesh.y_coordinates[0] ; 
	double y_max = microenvironment.mesh.y_coordinates[nYnodes-1] ; 
	
	int number_of_cells = parameters.ints("number_of_cells" ); 
	double spacing = (y_max - y_min) / ( number_of_cells + 1.0 );
		
	double x_position = microenvironment.mesh.x_coordinates[0] + 100.0;  
	double y_position = y_min + spacing;
	
	Cell* pC;
	
	std::cout << "Seeding cell number " ; 
	for( int n = 0 ; n < number_of_cells ; n++ )
	{
		std::cout << n << " "; 
		pC = create_cell(); 
		pC->assign_position( x_position , y_position , 0.0 ); 
		
		y_position += spacing; 
	}
	std::cout << " and done!" << std::endl; 
		
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output( 4 , "red" ); 
	output[0] = "black"; 
	
	return output; 
}
