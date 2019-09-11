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

#include "./custom_nests.h"

// declare cell definitions here
Cell_Definition dividing_cell;   // dividing
Cell_Definition differentiated_cell;   // differentiated


void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	//SeedRandom( parameters.ints("random_seed") ); // or specify a seed here


	// Force 2D sim BEFORE using cell_defaults to define our custom cell types
	cell_defaults.functions.set_orientation = up_orientation;
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true;

	// housekeeping
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions );

	// no death
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );

	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0;
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	// Turn off Cell Cycling
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	cell_defaults.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;

	// cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	// cell_defaults.parameters.o2_proliferation_saturation = 30.0;  
	// cell_defaults.parameters.o2_reference = 38.0; 

	// int oxygen_ID = microenvironment.find_density_index( "oxygen" ); // 0
	int inhibitor_ID = microenvironment.find_density_index( "inhibitor" ); // 1
	// std::cout << "----- create_cell_types: oxygen_ID, inhibitor_ID =" << oxygen_ID << ", "<<inhibitor_ID <<std::endl;
	std::cout << "----- create_cell_types: inhibitor_ID =" << inhibitor_ID <<std::endl;



	// cell_defaults.phenotype.secretion.secretion_rates[oxygen_ID] = 0;
	// cell_defaults.phenotype.secretion.uptake_rates[oxygen_ID] = 1000;
	// cell_defaults.phenotype.secretion.saturation_densities[oxygen_ID] = 38;

	// cell_defaults.phenotype.secretion.saturation_densities[oxygen_ID] = 38;

// Move these down *after* they are created/copied from cell_defaults
	// dividing_cell.phenotype.secretion.secretion_rates[inhibitor_ID] = 0;
	// differentiated_cell.phenotype.secretion.secretion_rates[inhibitor_ID] = 100;
	// dividing_cell.phenotype.secretion.uptake_rates[inhibitor_ID] = 100;
	// differentiated_cell.phenotype.secretion.uptake_rates[inhibitor_ID] = 0;
	// cell_defaults.phenotype.secretion.saturation_densities[inhibitor_ID] = 100;

	// Define the 2 cell types:

	// --- dividing cells
	dividing_cell = cell_defaults;
	dividing_cell.type = 0;
	dividing_cell.phenotype.motility.is_motile = true;
	dividing_cell.phenotype.differentiation.differentiation_possible = true;

	// dividing_cell.phenotype.secretion.secretion_rates[inhibitor_ID] = parameters.doubles("dividing_cell_inhibitor_secretion");
	// dividing_cell.phenotype.secretion.uptake_rates[inhibitor_ID] = parameters.doubles("dividing_cell_inhibitor_uptake");

	//Furkan
//	int glucose_substrate_index = microenvironment.find_density_index( "inhibitor" );
    // cell_defaults.phenotype.secretion.uptake_rates[inhibitor_ID] =0.0;
    // cell_defaults.phenotype.secretion.secretion_rates[inhibitor_ID] = 0.1;
    // cell_defaults.phenotype.secretion.saturation_densities[inhibitor_ID] = 100.0;
    dividing_cell.phenotype.secretion.uptake_rates[inhibitor_ID] =0.0;
    dividing_cell.phenotype.secretion.secretion_rates[inhibitor_ID] = 0.1;
    dividing_cell.phenotype.secretion.saturation_densities[inhibitor_ID] = 100.0;

	dividing_cell.functions.update_phenotype = dividing_cell_phenotype_rule; 

	// dividing_cell.functions.update_velocity = dividing_cell_velocity_rule; 

	// dividing_cell.functions.custom_cell_rule = cargo_cell_rule; 
	// dividing_cell.functions.update_migration_bias = NULL;	

	// dividing_cell.phenotype.secretion.secretion_rates[oxygen_ID] = 0;
	// dividing_cell.phenotype.secretion.uptake_rates[oxygen_ID] = 10;
	// dividing_cell.phenotype.secretion.saturation_densities[oxygen_ID] = 38;

	// --- differentiated cells
	differentiated_cell = cell_defaults;
	differentiated_cell.type = 1;
	differentiated_cell.phenotype.motility.is_motile = false;
	differentiated_cell.phenotype.differentiation.differentiation_possible = false;
	// differentiated_cell.functions.update_phenotype = nondividing_probability_update;
	differentiated_cell.functions.update_phenotype = nondividing_cell_phenotype_rule; 


	// no birth  (turn off proliferation)
	// int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	// int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	differentiated_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0; 



	// rwh: now set secretion and uptake rates - get them from the config file (e.g., config/PhysiCell_settings.xml)
	//
	// "non-dividing cells create some sort of inhibitor molecule and subsequently secrete this molecule. 
	// Then, the dividing cells would uptake this inhibitor and the probability that these dividing cells differentiate 
	// into non-dividing cells would vary with the concentration/absolute number of “inhibitor” molecules consumed.
	differentiated_cell.phenotype.secretion.secretion_rates[inhibitor_ID] = parameters.doubles("differentiated_cell_inhibitor_secretion");

	differentiated_cell.phenotype.secretion.uptake_rates[inhibitor_ID] = parameters.doubles("differentiated_cell_inhibitor_uptake");
	// cell_defaults.phenotype.secretion.saturation_densities[inhibitor_ID] = 100;
	std::cout << "\n------- sanity check inhibitor_ID rates:" <<std::endl;
	std::cout << "       dividing secretion=" << dividing_cell.phenotype.secretion.secretion_rates[inhibitor_ID] <<std::endl;
	std::cout << "       dividing uptake=" << dividing_cell.phenotype.secretion.uptake_rates[inhibitor_ID] <<std::endl;
	std::cout << "       differentiated secretion=" << differentiated_cell.phenotype.secretion.secretion_rates[inhibitor_ID] <<std::endl;
	std::cout << "       differentiated uptake=" << differentiated_cell.phenotype.secretion.uptake_rates[inhibitor_ID] <<std::endl;
	std::cout << "----------------------------\n" <<std::endl;


	// Do this AFTER defining the cell types??
	std::vector<Differentiation_Outcome> outcomes;
	Differentiation_Outcome* symmetric_0 = new Differentiation_Outcome(&dividing_cell, &dividing_cell);
	Differentiation_Outcome* asymmetric_0 = new Differentiation_Outcome(&dividing_cell, &differentiated_cell);
	Differentiation_Outcome* symmetricD_0 = new Differentiation_Outcome(&differentiated_cell, &differentiated_cell);
	outcomes.push_back(*symmetric_0);
	outcomes.push_back(*asymmetric_0);
	outcomes.push_back(*symmetricD_0);
	dividing_cell.phenotype.differentiation.outcomes = outcomes;


	std::vector<double> probabilities;  // pcts, percentages
	probabilities.push_back(0);
	probabilities.push_back(1);
	probabilities.push_back(0);

	// probabilities.push_back(.33);
	// probabilities.push_back(.34);
	// probabilities.push_back(.33);
	dividing_cell.phenotype.differentiation.probabilities = probabilities;

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

	/*
	cell1.phenotype.secretion.secretion_rates[Dfac_ID] = 0;
 	cell0.phenotype.secretion.secretion_rates[Dfac_ID] = 0;
	cell1.phenotype.secretion.uptake_rates[Dfac_ID] = 0;
	cell0.phenotype.secretion.uptake_rates[Dfac_ID] = 20;
	cell_defaults.phenotype.secretion.saturation_densities[Dfac_ID] = 100;

	default_microenvironment_options.calculate_gradients = true;
	default_microenvironment_options.track_internalized_substrates_in_each_agent = true;
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	*/
	// default_microenvironment_options.outer_Dirichlet_conditions = true;

	// default_microenvironment_options.use_oxygen_as_first_field = true;
	/*
	default_microenvironment_options.use_oxygen_as_first_field = false;
	microenvironment.set_density(0, "signal", "dimensionless" );
	microenvironment.decay_rates[0] = 0; 	 // no decay
	microenvironment.diffusion_coefficients[0] = 0.5;

	default_microenvironment_options.initial_condition_vector = { 100.0, 0.0, 0.0 };
	*/
	//std::vector<double> bc_vector( 1 , 38  );

	// default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // physioxic conditions
	// default_microenvironment_options.Dirichlet_activation_vector[0] = false;


	// parameters.doubles("dividing_cell_inhibitor_secretion");
	// parameters.bools("dividing_cell_inhibitor_secretion");
	// default_microenvironment_options.track_internalized_substrates_in_each_agent = parameters.bools("track_internalized_substrates_in_each_agent");
	// std::cout << "default_microenvironment_options.track_internalized_substrates_in_each_agent " << default_microenvironment_options.track_internalized_substrates_in_each_agent << std::endl;
	default_microenvironment_options.track_internalized_substrates_in_each_agent = true;

	// initialize BioFVM
	initialize_microenvironment();

	/*
	double factor = 1.0;
	for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		// factor = 1.0 - ((n%20)/20.);
		// bc_vector[0] = factor*38.;
		bc_vector[0] = 0.0;
		microenvironment(n) = bc_vector;
	}
	*/
	return;
}

void setup_tissue( void )
{

	Cell* pC;
	// create a single cell, at the origin
	pC = create_cell( dividing_cell );
	pC->assign_position( 0.0, 0.0, 0.0 );

	static int inhibitor_i = microenvironment.find_density_index("inhibitor");
	pC->phenotype.molecular.internalized_total_substrates[ inhibitor_i ] = 100; 

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

	if(pCell->type == 0 ) //dividing_cell
	{
		 output[0] = "green";
		 output[1] = "green";
		 output[2] = "green";
	} else if(pCell->type == 1 ) { //non dividing_cell (differentiated)
		 output[0] = "red";
		 output[1] = "red";
		 output[2] = "red";
	}
	return output;
}


void dividing_probability_update(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int inhibitor_i = microenvironment.find_density_index("inhibitor");
	double inhibitor_Con = phenotype.molecular.internalized_total_substrates[inhibitor_i];
	//pCell -> advance_bundled_phenotype_functions( dt);
	//phenotype.secretion.advance(pCell, phenotype, dt);
	double SymR = (8)/(10 + 10*inhibitor_Con);
	double Asym = 1 - SymR;
	std::cout << "------- dividing_probability_update: SymR, Asym = "<< SymR << ",  " << Asym << std::endl;

	std::vector<double> probabilities{ SymR, 0 , Asym };
	// phenotype.differentiation.probabilities = probabilities;

//	pCell->phenotype.volume *= 1.1; 
	std::cout << "------- : total vol = "<< pCell->get_total_volume() << std::endl;
	if (pCell->get_total_volume() > 10000) {
		// pCell->flag_for_division();
	}
	else {
		// pCell->set_total_volume(1.02 * pCell->get_total_volume());
	}
	return;
}

void nondividing_probability_update(Cell* pCell, Phenotype& phenotype, double dt)
{
	static int inhibitor_i = microenvironment.find_density_index("inhibitor");
	phenotype.molecular.internalized_total_substrates[inhibitor_i] += 1.0;

  	double inhibitor_con2 = phenotype.molecular.internalized_total_substrates[inhibitor_i];
	//phenotype.secretion.advance(pCell, phenotype, dt);
	//pCell -> advance_bundled_phenotype_functions( dt);
	std::cout << "nondividing_probability_update: inhibitor_con2= " << inhibitor_con2 << " ";
  return;
}

void nondividing_cell_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int inhibitor_i = microenvironment.find_density_index( "inhibitor" ); 

	// std::cout << "nondividing_cell_phenotype_rule: cell ID, vol= " << pCell->ID<<", " <<pCell->get_total_volume() << std::endl;
  	// double inhibitor_con2 = phenotype.molecular.internalized_total_substrates[inhibitor_i];
	// std::cout << "					: inhibitor_con2= " << inhibitor_con2 << " ";


	// std::cout << "------- : total vol = "<< pCell->get_total_volume() << std::endl;
	// pCell->set_total_volume(1.05 * pCell->get_total_volume());
	// if (pCell->get_total_volume() > 5000) {
	// 	pCell->flag_for_division();
	// }

	// if inhibitor high, secrete nothing, receptor off 
	// if( pCell->nearest_density_vector()[inhibitor_i] > pCell->custom_data[drop_index] )
	// {
	// 	phenotype.secretion.secretion_rates[inhibitor_i] = 0.0; 
	// 	// pCell->custom_data[receptor_index] = 0.0; 
	// }
	// else
	// {
	// 	phenotype.secretion.secretion_rates[inhibitor_i] = 0.0; 
	// 	// pCell->custom_data[receptor_index] = 0.0; 		
	// }
	
	return; 
}

void dividing_cell_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int inhibitor_i = microenvironment.find_density_index( "inhibitor" ); 

	std::cout << "dividing_cell_phenotype_rule: cell ID, vol= " << pCell->ID<<", " <<pCell->get_total_volume() << std::endl;
	// std::cout << "------- : total vol = "<< pCell->get_total_volume() << std::endl;
	pCell->set_total_volume(1.05 * pCell->get_total_volume());
	if (pCell->get_total_volume() > 5000) {
		std::cout << "         -- flag_for_division(): \n";
		pCell->flag_for_division();
	}

	// if inhibitor high, secrete nothing, receptor off 
	// if( pCell->nearest_density_vector()[inhibitor_i] > pCell->custom_data[drop_index] )
	// {
	// 	phenotype.secretion.secretion_rates[inhibitor_i] = 0.0; 
	// 	// pCell->custom_data[receptor_index] = 0.0; 
	// }
	// else
	// {
	// 	phenotype.secretion.secretion_rates[inhibitor_i] = 0.0; 
	// 	// pCell->custom_data[receptor_index] = 0.0; 		
	// }
	
	return; 
}