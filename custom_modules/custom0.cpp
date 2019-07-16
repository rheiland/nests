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
Cell_Definition cell0;   // dividing
Cell_Definition cell1;   // differentiated


void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	
	// Force 2D sim BEFORE using cell_defaults to define our custom cell types
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 

	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );
	

	// Define the 2 cell types:

	// --- dividing cells
	cell0 = cell_defaults;
	cell0.type = 0;
	cell0.phenotype.motility.is_motile = true;
	cell0.phenotype.differentiation.differentiation_possible = true;
	
	std::vector<double> probabilities;
	probabilities.push_back(0.5);
	probabilities.push_back(0.5);
	cell0.phenotype.differentiation.probabilities = probabilities;
	
	std::vector<Differentiation_Outcome> outcomes;
	Differentiation_Outcome* symmetric_0 = new Differentiation_Outcome(&cell0, &cell0);
	Differentiation_Outcome* asymmetric_0 = new Differentiation_Outcome(&cell0, &cell1);
	outcomes.push_back(*symmetric_0);
	outcomes.push_back(*asymmetric_0);
	cell0.phenotype.differentiation.outcomes = outcomes;
	
	// --- differentiated cells
	cell1 = cell_defaults;
	cell1.type = 1;
	cell1.phenotype.motility.is_motile = false;
	cell1.phenotype.differentiation.differentiation_possible = false;
	
	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* 
	// no gradients need for this example 
	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
*/

/* 
sets the custom rule (custom_cell_rule) to NULL and update phenotype function
(update_phenotype) to update_cell_and_death_parameters_O2_based, so the cell changes its cycle entry rate and
necrosis rate according to its local oxygenation conditions. (See Section 17.6.)
*/


	// cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 10; 
	// cell_defaults.phenotype.secretion.uptake_rates[0] = 10; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 

	default_microenvironment_options.calculate_gradients = false;
	default_microenvironment_options.track_internalized_substrates_in_each_agent = true;
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	// default_microenvironment_options.outer_Dirichlet_conditions = true;

	// default_microenvironment_options.use_oxygen_as_first_field = true;
	default_microenvironment_options.use_oxygen_as_first_field = false;
	microenvironment.set_density(0, "signal", "dimensionless" );
	microenvironment.decay_rates[0] = 0; 	 // no decay
	microenvironment.diffusion_coefficients[0] = 0.5; 	

	std::vector<double> bc_vector( 1 , 38  ); 
	
	// default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // physioxic conditions 
	// default_microenvironment_options.Dirichlet_activation_vector[0] = false;  
	
	// initialize BioFVM 
	initialize_microenvironment(); 	

	double factor = 1.0;
	for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		// factor = 1.0 - ((n%20)/20.);
		// bc_vector[0] = factor*38.;
		bc_vector[0] = 0.0;
		microenvironment(n) = bc_vector; 
	}	
	
	return; 
}

void setup_tissue( void )
{
	Cell* pC;
	pC = create_cell( cell0 );
	pC->assign_position( 0.0, 0.0, 0.0 );

	// create some cells near the origin
/* 	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double nest_radius = 60; 
	for(int i = -nest_radius; i < nest_radius; i+=cell_spacing)
	{
		for(int j = -nest_radius; j < nest_radius; j+=15){
			if(pow(i,2) + pow(j,2) <=  pow(nest_radius,2))
			{
				Cell* bC;
				bC = create_cell( stem_cell );
				bC->assign_position( j+25, i, 0.0 );
			}
		}
	}
*/	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
		
	if(pCell->type == 0 )
	{
		 output[0] = "green"; 
		 output[2] = "green"; 
	} else if(pCell->type == 1 ){
		 output[0] = "red"; 
		 output[2] = "red"; 
	} else if(pCell->type == 2 ){
		 output[0] = "blue"; 
		 output[2] = "blue"; 
	}
	
	return output; 
}