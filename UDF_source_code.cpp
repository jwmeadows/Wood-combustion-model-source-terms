#include "udf.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

/* The following code sets the Navier-Stokes source terms reflecting the pyrolysis volatile production, pyrolysis mass loss and char oxidation energy and volatile source terms. The code also defines the transient change in wood density and thermal conductivity. */
/* READ INSTRUCTIONS CAREFULLY BEFORE HOOKING UDFs */
/*IMPORTANT: Check fluid and solid zone, boundary IDs before proceeding. Also check to make sure if the n_cells variable is set to the correct number of cells on the surface of the wood logs */

/* This code is comprised of four different types of macros:
DEFINE_INIT macro sets the initialization conditions. In the current code, this macro sets a flag based on density for the pyrolysis source terms.
DEFINE_PROPERTY macro sets the thermal properties of a material. In the current code, this macro sets the wood density and thermal conductivity. This macro executes at the beginning of each time step.
DEFINE_EXECUTE_AT_END does not alter the N-S equation source terms and executes at the end of each time step. However, this macro is important for setting macro user defined memories (UDMs) which can be later used in the property or source term macros.
DEFINE_SOURCE sets the energy and species source terms in the current code. This macro executes after the DEFINE_PROPERTY macro at the beginning of each time step.*/

static real total_mass_loss; /* define total mass loss variable */
static real total_density_loss; /* define total density loss variable */
static real initial_density = 660; /* initial density of wood */ /* change according to mass */
static real final_density = 0.265 * initial_density; /* final density of wood */
static real ambient_temp = 300; /* ambient temperature */
static real n_cells = 173965; /* Only accounting for pyrolysis on the surface, change this when refining or coarsening the mesh */
static real cutoff_temperature = 873; /* cutoff temperature for high speed flow cases */
static real time_loss;
static real total_loss;
static real e = 2.71828; /* e number */
static real E_char = 124000; /* char oxidation activation energy kJ/mol*/
static real A_char = pow(10, 6.55); /* pre-exponential factor char oxidation */
static real n_char = 0.56; /* order of reaction */
static real oxygen_exp_char = 0.68; /* oxygen exponent */
static real R = 8.314; /* universal gas constant kJ/mol*/
static real E_hemicellulose = 102 * pow(10, 3); /* hemicellulose activation energy */
static real E_cellulose = 193.7 * pow(10, 3); /* cellulose activation energy */
static real E_lignin = 78.23 * pow(10, 3); /* lignin activiation energy */
static real A_hemicellulose = pow(10, 7.29); /* hemicellulose pre-exponential factor */
static real A_cellulose = pow(10, 14.06); /* cellulose pre-exponential factor */
static real A_lignin = pow(10, 3.042); /* lignin pre-exponential factor */
static real f_hemicellulose = 0.25; /* hemicellulose fraction */
static real f_cellulose = 0.33; /* cellulose fraction */
static real f_lignin = 0.155; /* lignin fraction */
static real f_res = 0.265; /* residual char fraction */
static real f_hemicellulose_prev; /* hemicellulose fraction previous*/
static real f_cellulose_prev; /* cellulose fraction previous*/
static real f_lignin_prev; /* lignin fraction previous*/
//static float solid_ID_begin = 4247; /* ID of solid zone */
static int fluid_ID = 17453; /* ID of fluid zone */
static int solid_ID = 17474; /* ID of solid zone */
static int fluid_interface_ID = 17693; /* ID of fluid interface */
static int solid_interface_ID = 17623; /* ID of solid interface */
//static float solid_interface_ID_begin = 2.0; /* begin ID of wood surfaces for multiple solid zones */
//static float solid_interface_ID_end = 81.0; /* end ID of wood surfaces for multiple solid zones */
//static int step_size_solid_bounds = (solid_interface_ID_end - solid_interface_ID_begin) + 1; /* step size of wood boundaries for multiple solid zones */
//static float fluid_interface_ID_begin = 504.0; /* begin ID of wood surface fluid for multiple solid zones */
//static float fluid_interface_ID_end = 583.0; /* end ID of wood surface fluid for multiple solid zones */
//static int step_size_fluid_bounds = (fluid_interface_ID_end - fluid_interface_ID_begin) + 1; /* step size of wood boundaries for multiple solid zones */
//static real total_char_loss_initial = 0.002995231684578447; /* total char loss from previous session */
//static real factor_initial = 0.6425676019291235; /* initial shrinkage factor from previous run */ 

/* populate and search arrays for multiple solid zones */

float* populate_linear_space(int size, float start, float end) {
	float* numbers = (float*)malloc(size * sizeof(float)); // Allocate memory 
	if (numbers == NULL) {
		printf("Memory allocation failed.\n");
		exit(1); // Exit if memory allocation fails
	}

	float step = (end - start) / (size - 1); // Calculate the step size

	for (int i = 0; i < size; i++) {
		numbers[i] = start + i * step; // Populate the array
	}

	return numbers; // Return the pointer to the array
}

/* result of index search in solid zone array for multiple solid zones */
int search_array(float* numbers, int size, float search_value) {
	for (int i = 0; i < size; i++) {
		if (fabs(numbers[i] - search_value) < 1e-6) { // Tolerance for floating-point comparison
			return i; // Return the index if the value is found
		}
	}
	return -1; // Return -1 if the value is not found
}

/* calculate mole fraction of species from mass fractions */

float get_mole_fractions(float mass_fracs[], int index)
{
	float oxygen_molecular_weight = 0.032; /* oxygen molecular weight */
	float water_molecular_weight = 0.018; /* water vapor molecular weight */
	float methane_molecular_weight = 0.016; /* methane molecular weight */
	float co_molecular_weight = 0.028; /* CO molecular weight */
	float co2_molecular_weight = 0.044; /* co2 molecular weight */
	float h2_molecular_weight = 0.002; /* h2 molecular weight */
	float n2_molecular_weight = 0.028; /* n2 molecular weight */
	float Mol_weights[] = { oxygen_molecular_weight, water_molecular_weight, methane_molecular_weight, co_molecular_weight, co2_molecular_weight, h2_molecular_weight, n2_molecular_weight }; /* molecular weights array. The species order in the array is based on the species list in your gas phase combustion mechanism. Look at the species names in Fluent to determine the order. */
	int length_weights = sizeof(Mol_weights) / sizeof(Mol_weights[0]); /* length of the species array */
	int i;
	float part_mole_frac;
	float mole_frac; /* species mole fraction */
	float mixture_mol_weight = 0.0; /* mixture molecular weight */
	for (i = 0; i < length_weights; i++)
	{
		mixture_mol_weight += mass_fracs[i] / Mol_weights[i]; /* sum up to get the mixture molecular weight */
	}
	part_mole_frac = mass_fracs[index] / Mol_weights[index];
	mole_frac = part_mole_frac / mixture_mol_weight;

	return mole_frac;
}


/* Initialize the density value */

DEFINE_INIT(density_initial, domain)
{
	
	cell_t c;
	Thread* t = Lookup_Thread(domain, solid_ID); /*look up the thread pointer to the solid zone */
	real rho_crit = final_density;
	total_mass_loss = 0; /*initialize global variables */
	total_density_loss = 0;
	real temp;
	real xc[ND_ND];
	real Tref = 288.16;
	
	thread_loop_c(t, domain) /*loop over cells in the zone*/
	{

		begin_c_loop(c, t)
		{

			/*Define density at previous time step and current time step*/
			temp = C_T(c, t);
			C_UDMI(c, t, 0) = C_R(c, t); // no apparent density. uncomment this stuff
			C_UDMI(c, t, 1) = C_R(c, t);
	
			if (C_R(c, t) <= rho_crit)
			{
				C_UDMI(c, t, 2) = 1; /*flag to indicate if the density variable is within the critical region*/
			}
			else
			{
				C_UDMI(c, t, 2) = 0;

			}
			C_CENTROID(xc, c, t);
			/*Message0("The cell center is %.5f, %.5f, %.5f\n", xc[0], xc[1], xc[2]);*/
		}
		end_c_loop(c, t)
		
	}

}

/* set wood properties */

DEFINE_PROPERTY(rho_wood, c, t)   /* actual wood density profile computation*/
{
	float rho_w;
	float rho_apparent = 0;
	float rho_w_prev = 0;
	float temp;
	float temp_prev = 0;
	float f_total;
	float time_step;
	real rho_crit;
	real celsius_temp;
	real flow_time;
	temp = C_T(c, t);
	rho_crit = final_density;
	celsius_temp = temp - 273.15;
	int n = 0; /* loop over boundary faces */
	//int zone_ID = 2; /* solid wall ID */
	flow_time = CURRENT_TIME;
	time_step = CURRENT_TIMESTEP;
	real Cp_cellulose;
	real Cp_hemicellulose;
	real Cp_lignin;
	real Cp_char;
	real Cp_total;
	real sensible_total;
	real sensible_cellulose;
	real sensible_hemicellulose;
	real sensible_lignin;
	real sensible_char;
	real massfrac_hemicellulose;
	real massfrac_cellulose;
	real massfrac_char;
	real massfrac_lignin;
	real Tref = 288.16;

	if (C_UDMI(c, t, 2) == 0)   /*UDM indicating if the density falls within the critical density cut-off region*/
	{


		if (flow_time == 0)
		{
			rho_w = initial_density;
			C_UDMI(c, t, 13) = f_hemicellulose;
			C_UDMI(c, t, 14) = f_cellulose;
			C_UDMI(c, t, 15) = f_lignin;
			initial_density = rho_w; /* initial density update after setting it */
		}
		else
		{
			f_total = C_UDMI(c, t, 13) + C_UDMI(c, t, 14) + C_UDMI(c, t, 15) + f_res;
			rho_w = f_total * initial_density;
		}

	}
	else
	{
		rho_w = final_density;
	}

	/* uncomment above loop for constant specific heat */
	C_UDMI(c, t, 40) = rho_w; /* store density to use in specific heat profile */

	
	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
		face_t f = C_FACE(c, t, n);
		//float* solid_interface_IDs = populate_linear_space(step_size_solid_bounds, solid_interface_ID_begin, solid_interface_ID_end); /* array containing solid interface IDs */
		//int index = search_array(solid_interface_IDs, step_size_solid_bounds, THREAD_ID(tf)); /* search whether the curret solid thread belongs to the solid interface */
		//if (index != -1)
		if (THREAD_ID(tf) == solid_interface_ID) /* check if the face belongs to the solid interface */
		{
			F_UDMI(f, tf, 19) = rho_w; /* store density of the boundary cells in F_UDMI. This value is not accessible through report plots */
		}

	}

	//rho_apparent = rho_w * sensible_total; /* apparent density. Uncomment to use */
	/*printf("This worked\n");*/ /*print statement indicating that the solver passed through the property macro*/
	return rho_w;
	//return rho_apparent;
}


/* set thermal conductivity property */

DEFINE_PROPERTY(thermal_conductivity, c, t)
{
	real k_wood;
	real k_char;
	real k_eff;
	float temp;
	real progress_variable;
	temp = C_T(c, t);

	k_wood = 0.1327 + (0.000004813 * temp);
	k_char = 0.05 + (0.0000002396 * temp);

	progress_variable = C_UDMI(c, t, 48); /* store the evaluated progress variable */
	k_eff = (progress_variable * k_char) + ((1 - progress_variable) * k_wood);
	return k_eff;
}

/* functions to execute at the end of the timestep */

/* main execute at end function */
DEFINE_EXECUTE_AT_END(execute_at_end) /*print density difference at the end of each time step */
{

#if !RP_HOST
	int i = 1; /*integer index corresponding to oxygen */
	real mass_loss = 0.;
	real density_loss;
	real vol_tot = 0.;
	real global_volume = 0.;
	real global_mass_loss = 0.;
	real temp;
	real char_mass_loss;
	cell_t c;
	cell_t c0;
	face_t f;
	int n = 0;
	int m = 0;
	Domain* domain;
	real time_step;
	real flow_time;
	real mlr;
	real cellulose_loss_rate = 0.0;
	real hemicellulose_loss_rate = 0.0;
	real lignin_loss_rate = 0.0;
	real hemicellulose_loss = 0.0;
	real cellulose_loss = 0.0;
	real lignin_loss = 0.0;
	real global_hemicellulose_loss = 0.0;
	real global_cellulose_loss = 0.0;
	real global_lignin_loss = 0.0;
	real rho_crit = final_density;
	domain = Get_Domain(1);
	//Thread* t = Lookup_Thread(domain, solid_ID); /* solid domain thread */
	Thread* t; /* initialize the solid domain threads */
	Thread* t0 = Lookup_Thread(domain, fluid_ID); /* fluid domain thread */
	real loss_cells = 0.0;
	real loss_boundary = 0.0;
	real progress_variable; /* define charring progress variable representing char fraction */
	real hemicellulose_time_scale;
	real cellulose_time_scale;
	real lignin_time_scale;
	real pyrolysis_time_scale;
	real thermal_diffusivity;
	real char_thickness;
	real thermal_diffusion_timescale;
	real damkohler_number;
	time_step = CURRENT_TIMESTEP;
	flow_time = CURRENT_TIME;
	int thread_id = 0;

	/* gather only cells in the solid domain and define UDF source terms for them */
	thread_loop_c(t, domain) /* loop through all the threads in the domain */
	{
		thread_id = THREAD_ID(t); /* gather the ID for the thread */

		if (thread_id == solid_ID) /* identify whether the zone is solid or not*/
		{
			begin_c_loop(c, t)
			{
				if (C_UDMI(c, t, 2) == 0)
				{
					C_UDMI(c, t, 1) = C_R(c, t); /*updated constant specific heat density at current time step*/
					C_UDMI(c, t, 0) = C_R_M1(c, t); /* constant specific heat density at previous time step is given by C_R_M1(c,t) */

					if (C_UDMI(c, t, 1) > C_UDMI(c, t, 0))  /*find if density is increasing with each time step */
					{
						C_UDMI(c, t, 1) = C_UDMI(c, t, 0);
						/*entered this loop */
						/*Message0("Entered this loop.");*/
					}
					C_UDMI(c, t, 3) = C_UDMI(c, t, 1) - C_UDMI(c, t, 0); /*store density difference*/
					C_UDMI(c, t, 11) = (C_UDMI(c, t, 3) * C_VOLUME(c, t)) / time_step; /* cell mass loss rate */
					temp = C_T(c, t); /* cell temperature */
					f_hemicellulose_prev = C_UDMI(c, t, 13) * pow(e, -A_hemicellulose * pow(e, -(E_hemicellulose / (R * temp))) * time_step);
					C_UDMI(c, t, 53) = -C_UDMI(c, t, 13) * A_hemicellulose * pow(e, -(E_hemicellulose / (R * temp))); /* derivative of hemicellulose */
					C_UDMI(c, t, 66) = f_hemicellulose_prev - C_UDMI(c, t, 13); /* change in hemicellulose fraction from previous time step. To be used in supplying species source terms. */
					C_UDMI(c, t, 13) = f_hemicellulose_prev;
					f_cellulose_prev = C_UDMI(c, t, 14) * pow(e, -A_cellulose * pow(e, -(E_cellulose / (R * temp))) * time_step);
					C_UDMI(c, t, 54) = -C_UDMI(c, t, 14) * A_cellulose * pow(e, -(E_cellulose / (R * temp))); /* derivative of cellulose */
					C_UDMI(c, t, 67) = f_cellulose_prev - C_UDMI(c, t, 14); /* change in cellulose fraction from previous time step. To be used in supplying species source terms. */
					C_UDMI(c, t, 14) = f_cellulose_prev;
					f_lignin_prev = C_UDMI(c, t, 15) * pow(e, -A_lignin * pow(e, -(E_lignin / (R * temp))) * time_step);
					C_UDMI(c, t, 55) = -C_UDMI(c, t, 15) * A_lignin * pow(e, -(E_lignin / (R * temp))); /* derivative of lignin */
					C_UDMI(c, t, 68) = f_lignin_prev - C_UDMI(c, t, 15); /* change in lignin fraction from previous time step. To be used in supplying species source terms. */
					C_UDMI(c, t, 15) = f_lignin_prev;
				}

				if (C_R(c, t) <= rho_crit)
				{
					C_UDMI(c, t, 2) = 1; /*update cut-off criteria for constant specific heat */

				}
				else
				{
					C_UDMI(c, t, 2) = 0;

				}

			}
			end_c_loop(c, t)

			begin_c_loop_int(c, t) /* loop to perform global summing operations over interior cells of a node */
			{
				if (C_UDMI(c, t, 2) == 0)
				{
					mass_loss += C_UDMI(c, t, 3) * C_VOLUME(c, t);
					vol_tot += C_VOLUME(c, t);
					hemicellulose_loss += (C_UDMI(c, t, 66) * initial_density * C_VOLUME(c, t));
					cellulose_loss += (C_UDMI(c, t, 67) * initial_density * C_VOLUME(c, t));
					lignin_loss += (C_UDMI(c, t, 68) * initial_density * C_VOLUME(c, t));
				}

			}
			end_c_loop_int(c, t) /* loop over interior cells only */

			begin_c_loop(c, t) /* loop to compute progress variable */
			{
				progress_variable = (initial_density - C_UDMI(c, t, 1)) / (initial_density - final_density);
				C_UDMI(c, t, 48) = progress_variable; /* store progress variable in UDMI */

				/* calculate pyrolysis timescales */
				if (f_hemicellulose_prev / f_hemicellulose <= 0.368 && C_UDMI(c, t, 72) != 1)
				{
					C_UDMI(c, t, 72) = 1; /* flag for calculating time constant */
					hemicellulose_time_scale = flow_time;
					C_UDMI(c, t, 75) = hemicellulose_time_scale;
					//C_UDMI(c, t, 72) = 0; /* reset the flag */
				}

				if (f_cellulose_prev / f_cellulose <= 0.368 && C_UDMI(c, t, 73) != 1)
				{
					C_UDMI(c, t, 73) = 1; /* flag for calculating time constant */
					cellulose_time_scale = flow_time;
					C_UDMI(c, t, 76) = cellulose_time_scale;
					//C_UDMI(c, t, 73) = 0; /* reset the flag */
				}

				if (f_lignin_prev / f_lignin <= 0.368 && C_UDMI(c, t, 74) != 1)
				{
					C_UDMI(c, t, 74) = 1; /* flag for calculating time constant */
					lignin_time_scale = flow_time;
					C_UDMI(c, t, 77) = lignin_time_scale;
					//C_UDMI(c, t, 74) = 0; /* reset the flag */
				}

				/* calculate max pyrolysis timescale */
				if (lignin_time_scale >= hemicellulose_time_scale)
				{
					if (lignin_time_scale >= cellulose_time_scale)
						pyrolysis_time_scale = lignin_time_scale;
					else
						pyrolysis_time_scale = cellulose_time_scale;
				}
				else
				{
					if (hemicellulose_time_scale >= cellulose_time_scale)
						pyrolysis_time_scale = hemicellulose_time_scale;
					else
						pyrolysis_time_scale = cellulose_time_scale;
				}

				thermal_diffusivity = C_K_L(c, t) / (C_R(c, t) * C_CP(c, t)); /* thermal diffusivity */
				C_UDMI(c, t, 80) = thermal_diffusivity;
				char_thickness = progress_variable * pow(C_VOLUME(c, t), 0.333);
				thermal_diffusion_timescale = pow(char_thickness, 2) / thermal_diffusivity; /* thermal diffusion timescale */
				C_UDMI(c, t, 78) = thermal_diffusion_timescale;

				if (pyrolysis_time_scale != 0)
					damkohler_number = thermal_diffusion_timescale / pyrolysis_time_scale; /* Damkohler number calculation */
				C_UDMI(c, t, 79) = damkohler_number;
			}
			end_c_loop(c, t)


#if RP_NODE /* perform global summing operation and sum over all parallel processors to calculate mass loss rate */
			global_volume = PRF_GRSUM1(vol_tot);
			global_mass_loss = PRF_GRSUM1(mass_loss);
			global_hemicellulose_loss = PRF_GRSUM1(hemicellulose_loss);
			global_cellulose_loss = PRF_GRSUM1(cellulose_loss);
			global_lignin_loss = PRF_GRSUM1(lignin_loss);
#endif /* RP NODE */

			mlr = global_mass_loss / time_step;
			density_loss = global_mass_loss / global_volume;
			hemicellulose_loss_rate = global_hemicellulose_loss / time_step;
			cellulose_loss_rate = global_cellulose_loss / time_step;
			lignin_loss_rate = global_lignin_loss / time_step;

			begin_c_loop(c, t)
			{
				C_UDMI(c, t, 4) = density_loss;
				C_UDMI(c, t, 5) = mlr;
				C_UDMI(c, t, 43) = global_volume;
			}
			end_c_loop(c, t)

			begin_c_loop(c, t) /* loop over solid domain */
			{
				c_face_loop(c, t, n)
				{
					Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
					face_t f = C_FACE(c, t, n);
					//float* solid_interface_IDs = populate_linear_space(step_size_solid_bounds, solid_interface_ID_begin, solid_interface_ID_end); /* array containing solid interface IDs */
					//int index = search_array(solid_interface_IDs, step_size_solid_bounds, THREAD_ID(tf)); /* search whether the curret solid thread belongs to the solid interface */
					//if (index != -1)
					if (THREAD_ID(tf) == solid_interface_ID)
					{
						F_UDMI(f, tf, 41) = mlr; /* store density of the boundary cells */
						F_UDMI(f, tf, 9) = C_UDMI(c, t, 48); /* store progress variable of boundary cells */

					}

				}
			}
			end_c_loop(c, t)
			
			/* update global total mass and density loss variables only on compute node-0. Summing over other nodes will count extra entries. */

			if (I_AM_NODE_ZERO_P)
			{
				total_mass_loss += global_mass_loss;
				total_density_loss += density_loss;
			}

			begin_c_loop(c, t)
			{
				C_UDMI(c, t, 6) = mass_loss;
				C_UDMI(c, t, 7) = total_mass_loss;
				C_UDMI(c, t, 8) = total_density_loss;
			}
			end_c_loop(c, t)
		}
	}


	begin_c_loop(c0, t0) /* loop over fluid domain */
	{
		c_face_loop(c0, t0, m)
		{
			Thread* tf0 = C_FACE_THREAD(c0, t0, m); /* boundary face thread */
			face_t f0 = C_FACE(c0, t0, m); /* get the boundary face index */

			if PRINCIPAL_FACE_P(f0, tf0)
			{
				//float* fluid_interface_IDs = populate_linear_space(step_size_fluid_bounds, fluid_interface_ID_begin, fluid_interface_ID_end); /* populate fluid interface ID array */
				//int index = search_array(fluid_interface_IDs, step_size_fluid_bounds, THREAD_ID(tf0));

				//if (index != -1) /* check if the cell face thread is equivalent to the boundary thread */
				if (THREAD_ID(tf0) == fluid_interface_ID)
				{

					face_t f_shadow = F_SHADOW(f0, tf0); /* get the shadow face index on the fluid side  */
					Thread* t_shadow = THREAD_SHADOW(tf0); /* get the shadow face thread */
					total_loss += C_UDMI(c0, t0, 81) * time_step; /* get total mass loss over the given node */
					loss_cells += C_UDMI(c0, t0, 81); /* get the total mass loss rate at the current time step over the given node */
				}
			}

		}
	}
	end_c_loop(c0, t0)

#if RP_NODE /* perform global summing operation */

	time_loss = PRF_GRSUM1(total_loss); /* sum over all nodes to get the total mass lost until the current time step */
	loss_boundary = PRF_GRSUM1(loss_cells); /* sum over all nodes to get the total mass loss rate */
#endif /* RP NODE */


	begin_c_loop(c0, t0) /* loop over fluid domain */
	{
		c_face_loop(c0, t0, m)
		{
			face_t f0 = C_FACE(c0, t0, m); /* get the boundary face index */
			Thread* tf0 = C_FACE_THREAD(c0, t0, m); /* boundary face thread */
			//float* fluid_interface_IDs = populate_linear_space(step_size_fluid_bounds, fluid_interface_ID_begin, fluid_interface_ID_end); /* populate fluid interface ID array */
			//int index = search_array(fluid_interface_IDs, step_size_fluid_bounds, THREAD_ID(tf0));

			//if (index != -1) /* check if the cell face thread is equivalent to the boundary thread */
			if (THREAD_ID(tf0) == fluid_interface_ID)
			{
				face_t f_shadow = F_SHADOW(f0, tf0); /* get the shadow face index */
				Thread* t_shadow = THREAD_SHADOW(tf0); /* get the shadow face thread */
				C_UDMI(c0, t0, 83) = time_loss; /* store total mass lost until current time step */
				C_UDMI(c0, t0, 33) = C_UDMI(c0, t0, 29); /* store the char loss distribution in UDM */
				C_UDMI(c0, t0, 42) = F_UDMI(f_shadow, t_shadow, 41); /* store mass loss rate in pyrolysis */
				C_UDMI(c0, t0, 10) = F_UDMI(f_shadow, t_shadow, 9); /* store the char progress variable of boundary cells in the fluid domain */
				C_UDMI(c0, t0, 20) = F_UDMI(f_shadow, t_shadow, 19); /* assign the density value */
				C_UDMI(c0, t0, 63) = hemicellulose_loss_rate / n_cells; /* hemicellulose loss rate */
				C_UDMI(c0, t0, 64) = cellulose_loss_rate / n_cells; /* cellulose loss rate */
				C_UDMI(c0, t0, 65) = lignin_loss_rate / n_cells; /* hemicellulose loss rate */

			}
		}
	}
	end_c_loop(c0, t0)

#endif /* Compute nodes block */

#if !RP_NODE /* print in the host node only */
		Message("Entered the first loop.\n");
#endif /* Host node block */
}

DEFINE_EXECUTE_AT_END(execute_at_end_solid) /* loop over solid to assign char distribution values */
{
#if !RP_HOST
	cell_t c;
	int i = 1; /*integer index corresponding to oxygen */
	Domain* domain;
	domain = Get_Domain(1);
	int n = 0;
	//Thread* t = Lookup_Thread(domain, solid_ID); /* solid domain thread */
	Thread* t;
	int thread_id;

	/* gather only cells in the solid domain and define UDF source terms for them */
	thread_loop_c(t, domain) /* loop through all the threads in the domain */
	{
		thread_id = THREAD_ID(t); /* gather the ID for the thread */

		if (thread_id == solid_ID) /* identify whether the zone is solid or not */
		{
			begin_c_loop(c, t) /* loop over solid domain */
			{
				c_face_loop(c, t, n)
				{
					Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
					face_t f = C_FACE(c, t, n); /* get the boundary face index */
					//float* solid_interface_IDs = populate_linear_space(step_size_solid_bounds, solid_interface_ID_begin, solid_interface_ID_end); /* array containing solid interface IDs */
					//int index = search_array(solid_interface_IDs, step_size_solid_bounds, THREAD_ID(tf)); /* search whether the current solid thread belongs to the solid interface */
					//if (index != -1)
					if (THREAD_ID(tf) == solid_interface_ID)
					{

						face_t f_shadow = F_SHADOW(f, tf); /* get the shadow face index */
						Thread* t_shadow = THREAD_SHADOW(tf); /* get the shadow face thread */
						C_UDMI(c, t, 35) = F_UDMI(f_shadow, t_shadow, 82); /* store the local char loss distribution (i.e., C_UDMI(c, t, 29)) in the solid domain */
						C_UDMI(c, t, 36) = time_loss; /* store the total mass lost until current time step in the solid domain */
					}
				}
			}
			end_c_loop(c, t)
		}

	}
	

#endif /* Compute nodes block */
#if !RP_NODE /* print in the host node only */
		Message("Entered the second loop.\n");
#endif /* Host node block */
}



DEFINE_EXECUTE_ON_LOADING(load, libname)
{
	Message("Loading UDF library %s\n", libname);
	Message("Warning: Hook any UDFs using the 'function hooks' option before initializing the flow !!!\n");

	int udf_mem = 80; /*defined amount of udf memory allocations*/

	if (n_udm < udf_mem)
	{
		Message("Error: This fluent session requires %i user-defined memory allocations but no allocations have been made !!\n", udf_mem);
		return;
	}
}

/* Pyrolysis species source terms */

DEFINE_SOURCE(CO_kinetics, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real volume; /* define volume variable */
	real hemicellulose_fraction_CO = 0.0217; /* fraction of evolved gases */
	real cellulose_fraction_CO = 0.0336;
	real lignin_fraction_CO = 0.0988;
	real hemicellulose_fraction_CO2 = 0.144; /* fraction of evolved gases */
	real cellulose_fraction_CO2 = 0.0238;
	real lignin_fraction_CO2 = 0.12;
	real hemicellulose_fraction_CH4 = 0.0142; /* fraction of evolved gases */
	real cellulose_fraction_CH4 = 0.0107;
	real lignin_fraction_CH4 = 0.072;
	real total_hemicellulose_fraction = 0.1799;
	real total_cellulose_fraction = 0.0681;
	real total_lignin_fraction = 0.2908;
	float CO_mf;
	float CO2_mf;
	float CH4_mf;
	real total_flow_rate;
	int CO_index = 3;
	real mole_frac_co;
	int n = 0;
	volume = C_VOLUME(c, t);

	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
		//float* fluid_interface_IDs = populate_linear_space(step_size_fluid_bounds, fluid_interface_ID_begin, fluid_interface_ID_end); /* populate fluid interface ID array */
		//int index = search_array(fluid_interface_IDs, step_size_fluid_bounds, THREAD_ID(tf));

		//if (index != -1) /* check if the cell face thread is equivalent to the boundary thread */
		if (THREAD_ID(tf) == fluid_interface_ID)
		{
			total_flow_rate = ((total_hemicellulose_fraction * -C_UDMI(c, t, 63)) + (total_cellulose_fraction * -C_UDMI(c, t, 64)) + (total_lignin_fraction * -C_UDMI(c, t, 65)));			

			if (total_flow_rate != 0)
			{
				CO_mf = ((hemicellulose_fraction_CO * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO * -C_UDMI(c, t, 64)) + (lignin_fraction_CO * -C_UDMI(c, t, 65))) / total_flow_rate;
				CO2_mf = ((hemicellulose_fraction_CO2 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO2 * -C_UDMI(c, t, 64)) + (lignin_fraction_CO2 * -C_UDMI(c, t, 65))) / total_flow_rate;
				CH4_mf = ((hemicellulose_fraction_CH4 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CH4 * -C_UDMI(c, t, 64)) + (lignin_fraction_CH4 * -C_UDMI(c, t, 65))) / total_flow_rate;
				float mass_fracs[] = { 0, 0, CH4_mf, CO_mf, CO2_mf, 0, 0 }; /* species mass fraction array */
				mole_frac_co = get_mole_fractions(mass_fracs, CO_index);
				source = (mole_frac_co / volume) * (-C_UDMI(c, t, 42) / n_cells);
			}
			C_UDMI(c, t, 56) = source;
			C_UDMI(c, t, 69) = mole_frac_co;
		}
	}
	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;

}

DEFINE_SOURCE(CO2_kinetics, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real volume; /* define volume variable */
	real hemicellulose_fraction_CO = 0.0217; /* fraction of evolved gases */
	real cellulose_fraction_CO = 0.0336;
	real lignin_fraction_CO = 0.0988;
	real hemicellulose_fraction_CO2 = 0.144; /* fraction of evolved gases */
	real cellulose_fraction_CO2 = 0.0238;
	real lignin_fraction_CO2 = 0.12;
	real hemicellulose_fraction_CH4 = 0.0142; /* fraction of evolved gases */
	real cellulose_fraction_CH4 = 0.0107;
	real lignin_fraction_CH4 = 0.072;
	real total_hemicellulose_fraction = 0.1799;
	real total_cellulose_fraction = 0.0681;
	real total_lignin_fraction = 0.2908;
	float CO_mf;
	float CO2_mf;
	float CH4_mf;
	real total_flow_rate;
	int CO2_index = 4;
	real mole_frac_co2;
	int n = 0;
	volume = C_VOLUME(c, t);

	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
		//float* fluid_interface_IDs = populate_linear_space(step_size_fluid_bounds, fluid_interface_ID_begin, fluid_interface_ID_end); /* populate fluid interface ID array */
		//int index = search_array(fluid_interface_IDs, step_size_fluid_bounds, THREAD_ID(tf));

		//if (index != -1) /* check if the cell face thread is equivalent to the boundary thread */
		if (THREAD_ID(tf) == fluid_interface_ID)
		{
			total_flow_rate = ((total_hemicellulose_fraction * -C_UDMI(c, t, 63)) + (total_cellulose_fraction * -C_UDMI(c, t, 64)) + (total_lignin_fraction * -C_UDMI(c, t, 65)));
			//C_UDMI(c, t, 69) = ((total_hemicellulose_fraction * -C_UDMI(c, t, 63)) + (total_cellulose_fraction * -C_UDMI(c, t, 64)) + (total_lignin_fraction * -C_UDMI(c, t, 65))) / volume;
			//source = ((hemicellulose_fraction * -C_UDMI(c, t, 63)) + (cellulose_fraction * -C_UDMI(c, t, 64)) + (lignin_fraction * -C_UDMI(c, t, 65))) / volume;

			if (total_flow_rate != 0)
			{
				CO_mf = ((hemicellulose_fraction_CO * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO * -C_UDMI(c, t, 64)) + (lignin_fraction_CO * -C_UDMI(c, t, 65))) / total_flow_rate;
				CO2_mf = ((hemicellulose_fraction_CO2 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO2 * -C_UDMI(c, t, 64)) + (lignin_fraction_CO2 * -C_UDMI(c, t, 65))) / total_flow_rate;
				CH4_mf = ((hemicellulose_fraction_CH4 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CH4 * -C_UDMI(c, t, 64)) + (lignin_fraction_CH4 * -C_UDMI(c, t, 65))) / total_flow_rate;
				float mass_fracs[] = { 0, 0, CH4_mf, CO_mf, CO2_mf, 0, 0 }; /* species mass fraction array */
				mole_frac_co2 = get_mole_fractions(mass_fracs, CO2_index);
				//source = (((hemicellulose_fraction * -C_UDMI(c, t, 63)) + (cellulose_fraction * -C_UDMI(c, t, 64)) + (lignin_fraction * -C_UDMI(c, t, 65))) / volume) * (-C_UDMI(c, t, 42) / (n_cells * total_flow_rate));
				source = (mole_frac_co2 / volume) * (-C_UDMI(c, t, 42) / n_cells);
			}
			C_UDMI(c, t, 57) = source;
			C_UDMI(c, t, 70) = mole_frac_co2;
		}
	}

	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;
}

DEFINE_SOURCE(CH4_kinetics, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real volume; /* define volume variable */
	real hemicellulose_fraction_CO = 0.0217; /* fraction of evolved gases */
	real cellulose_fraction_CO = 0.0336;
	real lignin_fraction_CO = 0.0988;
	real hemicellulose_fraction_CO2 = 0.144; /* fraction of evolved gases */
	real cellulose_fraction_CO2 = 0.0238;
	real lignin_fraction_CO2 = 0.12;
	real hemicellulose_fraction_CH4 = 0.0142; /* fraction of evolved gases */
	real cellulose_fraction_CH4 = 0.0107;
	real lignin_fraction_CH4 = 0.072;
	real total_hemicellulose_fraction = 0.1799;
	real total_cellulose_fraction = 0.0681;
	real total_lignin_fraction = 0.2908;
	float CO_mf;
	float CO2_mf;
	float CH4_mf;
	real total_flow_rate;
	int CH4_index = 2;
	real mole_frac_ch4;
	int n = 0;
	volume = C_VOLUME(c, t);

	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */
		//float* fluid_interface_IDs = populate_linear_space(step_size_fluid_bounds, fluid_interface_ID_begin, fluid_interface_ID_end); /* populate fluid interface ID array */
		//int index = search_array(fluid_interface_IDs, step_size_fluid_bounds, THREAD_ID(tf));

		//if (index != -1) /* check if the cell face thread is equivalent to the boundary thread */
		if (THREAD_ID(tf) == fluid_interface_ID)
		{
			total_flow_rate = ((total_hemicellulose_fraction * -C_UDMI(c, t, 63)) + (total_cellulose_fraction * -C_UDMI(c, t, 64)) + (total_lignin_fraction * -C_UDMI(c, t, 65)));
			

			if (total_flow_rate != 0)
			{
				CO_mf = ((hemicellulose_fraction_CO * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO * -C_UDMI(c, t, 64)) + (lignin_fraction_CO * -C_UDMI(c, t, 65))) / total_flow_rate;
				CO2_mf = ((hemicellulose_fraction_CO2 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CO2 * -C_UDMI(c, t, 64)) + (lignin_fraction_CO2 * -C_UDMI(c, t, 65))) / total_flow_rate;
				CH4_mf = ((hemicellulose_fraction_CH4 * -C_UDMI(c, t, 63)) + (cellulose_fraction_CH4 * -C_UDMI(c, t, 64)) + (lignin_fraction_CH4 * -C_UDMI(c, t, 65))) / total_flow_rate;
				float mass_fracs[] = { 0, 0, CH4_mf, CO_mf, CO2_mf, 0, 0 }; /* species mass fraction array */
				mole_frac_ch4 = get_mole_fractions(mass_fracs, CH4_index);
				//source = (((hemicellulose_fraction * -C_UDMI(c, t, 63)) + (cellulose_fraction * -C_UDMI(c, t, 64)) + (lignin_fraction * -C_UDMI(c, t, 65))) / volume) * (-C_UDMI(c, t, 42) / (n_cells * total_flow_rate));
				source = (mole_frac_ch4 / volume) * (-C_UDMI(c, t, 42) / n_cells);
			}
			C_UDMI(c, t, 58) = source;
			C_UDMI(c, t, 71) = mole_frac_ch4;
		}
	}


	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;

}

/* Pyrolysis energy source terms */

DEFINE_SOURCE(hemicellulose_hod, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real hod = -39.53 * pow(10, 3); /* kJ/kg. Exothermic */
	source = hod * C_UDMI(c, t, 53) * initial_density;
	C_UDMI(c, t, 59) = source;
	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;
}

DEFINE_SOURCE(cellulose_hod, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real hod = 146.35 * pow(10, 3); /* kJ/kg. Endothermic */
	source = hod * C_UDMI(c, t, 54) * initial_density;
	C_UDMI(c, t, 60) = source;
	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;
}

DEFINE_SOURCE(lignin_hod, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real hod = -579.87 * pow(10, 3); /* kJ/kg. Endothermic */
	source = hod * C_UDMI(c, t, 55) * initial_density;
	C_UDMI(c, t, 61) = source;
	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;
}

/* Char oxidation species source terms */

/* The char oxidation macros include species source terms, i.e., CO, CO2, and O2, which are modeled using the following single-step global reaction.
	C + (1 - (f/2))O2 -> fCO + (1-f)CO2
 In the default state below, partial oxidation to CO is considered. The CO2 source term is commented out. Practically, smoldering during wood combustion only involves partial oxidation. Thus, complete oxidation to CO2 is not likely. Associated with the CO/CO2 split, the energy source term is assigned. For the present code, the energy source is 9 MJ/kg for partial oxidation to CO. Complete oxidation to CO2 produces ~32 MJ/kg.
 The char mass loss rate is computed using Arrhenius rate parameters provided by Anca-Couce et al. https://doi.org/10.1016/j.combustflame.2011.11.015
*/

/* Oxygen sink */
DEFINE_SOURCE(oxygen_sink, c, t, ds, eqn)
{
	real source = 0.0;
	real volume; /*cell volume */
	int n = 0;
	real r = 1.333;
	real loss_rate;
	real temp_surface;
	real oxygen_massfrac;
	real k;
	real density_surface;
	real shrinkage, factor;
	int oxygen_index = 0; /* oxygen index in the species names */
	int water_index = 1; /* water vapor index in the species names */
	int methane_index = 2; /* methane index in the species names */
	int co_index = 3; /* co index in the species names */
	int co2_index = 4; /* co2 index in the species names */
	int h2_index = 5; /* h2 index in the species names */
	int n2_index = 6; /* n2 index in the species names */
	float oxygen_mf; /* species mass fractions */
	float water_mf;
	float methane_mf;
	float co_mf;
	float co2_mf;
	float h2_mf;
	float n2_mf;
	float mole_frac_o2; /* oxygen mole fraction */
	real conversion; /* conversion from TGA */
	real cumulative_loss; /* cumulative char mass loss */
	real dt = CURRENT_TIMESTEP; /* timestep in seconds */
	real exp_factor_char;
	real char_progress;
	real mass_char;
	real mass_ash;

	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */

		if (THREAD_ID(tf) == fluid_zone_ID) /* check if the cell face thread is equivalent to the boundary thread */
		{
			face_t f = C_FACE(c, t, n); /* get face index */
			temp_surface = F_T(f, tf); /* surface temperature */
			char_progress = C_UDMI(c, t, 10);
			volume = C_VOLUME(c, t);
			mass_char = ((char_progress * final_density)) * volume; /* mass of char */
			mass_ash = 0.02 * mass_char; /* mass of ash */
			density_surface = C_R(c, t); /* adjacent cell density at the surface */
			// Retrieve cumulative loss from last time step
			cumulative_loss = C_UDMI(c, t, 88);
			shrinkage = cumulative_loss / (mass_char - mass_ash); /* shrinkage factor */
			if (shrinkage >= 1.0) shrinkage = 1.0; /* check if shrinkage factor is exceeding 1*/
			if (isnan(shrinkage)) shrinkage = 0.0; /* check if shrinkage factor is NaN */
			factor = 1 - shrinkage;
			oxygen_mf = C_YI(c, t, oxygen_index); /* mass fraction of oxidizer */
			water_mf = C_YI(c, t, water_index); /* mass fraction of water */
			methane_mf = C_YI(c, t, methane_index); /* mass fraction of methane */
			co_mf = C_YI(c, t, co_index); /* mass fraction of co */
			co2_mf = C_YI(c, t, co2_index); /* mass fraction of co2 */
			h2_mf = C_YI(c, t, h2_index); /* mass fraction of h2 */
			n2_mf = C_YI(c, t, n2_index); /* mass fraction of n2 */
			float mass_fracs[] = { oxygen_mf, water_mf, methane_mf, co_mf, co2_mf, h2_mf, n2_mf }; /* species mass fraction array */
			mole_frac_o2 = get_mole_fractions(mass_fracs, oxygen_index); /* normalized mole fraction of o2 */
			conversion = ((mass_char * factor) - mass_ash) / (mass_char - mass_ash);
			if (isnan(conversion)) conversion = 1; /* check if conversion is NaN */
			exp_factor_char = -(E_char / (R * temp_surface));
			loss_rate = A_char * pow(e, exp_factor_char) * pow((mole_frac_o2 / 0.205), oxygen_exp_char) * 0.98 * char_progress * final_density * pow(conversion, n_char) * volume; // use an if flag with C_UDMI(c, t, 10) >= 0.95, i.e., progress variable greater than 95% if the oxidation source terms lead to instabilities.
			//loss_rate = A_char * pow(e, exp_factor_char) * pow((mole_frac_o2 / 0.205), oxygen_exp_char) * 0.98 * char_progress * final_density * pow(conversion, n_char) * volume; // use this version if the char oxidation source terms still lead to instabilities. Conversion timescales for char oxidation are much longer than the pyrolytic degradation timescales.
			C_UDMI(c, t, 81) = loss_rate; /* store the char loss rate */
			F_UDMI(f, tf, 82) = loss_rate; /* store the char loss rate in a face thread for shadow face operations */

			// Update cumulative loss
			cumulative_loss += loss_rate * dt;
			C_UDMI(c, t, 88) = cumulative_loss;
			
			source = -(r * loss_rate) / (volume);
		}
	}

	C_UDMI(c, t, 87) = source; /* store the oxygen sink in UDM */
	ds[eqn] = 0.0;
	return source;
}

/* CO2 source */
DEFINE_SOURCE(char_co2_flux, c, t, ds, eqn)
{
	real source = 0.0;
	real volume; /*cell volume */
	int n = 0;
	real r = 1.333;
	real loss_rate;
	real temp_surface;
	real density_surface;
	real shrinkage, factor;
	int oxygen_index = 0; /* oxygen index in the species names */
	int water_index = 1; /* water vapor index in the species names */
	int methane_index = 2; /* methane index in the species names */
	int co_index = 3; /* co index in the species names */
	int co2_index = 4; /* co2 index in the species names */
	int h2_index = 5; /* h2 index in the species names */
	int n2_index = 6; /* n2 index in the species names */
	float oxygen_mf; /* species mass fractions */
	float water_mf;
	float methane_mf;
	float co_mf;
	float co2_mf;
	float h2_mf;
	float n2_mf;
	float mole_frac_o2; /* oxygen mole fraction */
	real conversion; /* conversion from TGA */
	real cumulative_loss; /* cumulative char mass loss */
	real dt = CURRENT_TIMESTEP; /* timestep in seconds */
	real exp_factor_char;
	real flux_num;
	real char_progress; /* progress variable */
	real mass_char;
	real mass_ash;

	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */

		if (THREAD_ID(tf) == fluid_zone_ID) /* check if the cell face thread is equivalent to the boundary thread */
		{
			face_t f = C_FACE(c, t, n); /* get face index */
			temp_surface = F_T(f, tf); /* surface temperature */
			density_surface = C_R(c, t); /* adjacent cell density at the surface */
			char_progress = C_UDMI(c, t, 10);
			volume = C_VOLUME(c, t);
			mass_char = ((char_progress * final_density)) * volume;
			mass_ash = 0.02 * mass_char;
			// Retrieve cumulative loss from last time step
			cumulative_loss = C_UDMI(c, t, 89);
			shrinkage = cumulative_loss / (mass_char - mass_ash); /* shrinkage factor */
			if (shrinkage >= 1.0) shrinkage = 1.0; /* check if shrinkage factor is exceeding 1*/
			if (isnan(shrinkage)) shrinkage = 0.0; /* check if shrinkage factor is NaN */
			C_UDMI(c, t, 86) = shrinkage; /* store the shrinkage factor */
			factor = 1 - shrinkage;
			oxygen_mf = C_YI(c, t, oxygen_index); /* mass fraction of oxidizer */
			water_mf = C_YI(c, t, water_index); /* mass fraction of water */
			methane_mf = C_YI(c, t, methane_index); /* mass fraction of methane */
			co_mf = C_YI(c, t, co_index); /* mass fraction of co */
			co2_mf = C_YI(c, t, co2_index); /* mass fraction of co2 */
			h2_mf = C_YI(c, t, h2_index); /* mass fraction of h2 */
			n2_mf = C_YI(c, t, n2_index); /* mass fraction of n2 */
			float mass_fracs[] = { oxygen_mf, water_mf, methane_mf, co_mf, co2_mf, h2_mf, n2_mf }; /* species mass fraction array */
			mole_frac_o2 = get_mole_fractions(mass_fracs, oxygen_index); /* mole fraction of o2 */
			conversion = ((mass_char * factor) - mass_ash) / (mass_char - mass_ash);
			if (isnan(conversion)) conversion = 1; /* check if conversion is NaN */
			exp_factor_char = -(E_char / (R * temp_surface));
			loss_rate = A_char * pow(e, exp_factor_char) * pow((mole_frac_o2 / 0.205), oxygen_exp_char) * 0.98 * char_progress * final_density * pow(conversion, n_char) * volume;

			// Update cumulative loss
			cumulative_loss += loss_rate * dt;
			C_UDMI(c, t, 89) = cumulative_loss;
			
			//source = (r * loss_rate) / (volume);
			
		}
	}

	C_UDMI(c, t, 84) = source; /* store the CO2 source in UDM */
	ds[eqn] = 0.0;
	return source;
}

/* CO source */
DEFINE_SOURCE(char_co_flux, c, t, ds, eqn)
{
	real source = 0.0;
	real volume; /*cell volume */
	int n = 0;
	real r = 0.5;
	real loss_rate;
	real temp_surface;
	real density_surface;
	real shrinkage, factor;
	int oxygen_index = 0; /* oxygen index in the species names */
	int water_index = 1; /* water vapor index in the species names */
	int methane_index = 2; /* methane index in the species names */
	int co_index = 3; /* co index in the species names */
	int co2_index = 4; /* co2 index in the species names */
	int h2_index = 5; /* h2 index in the species names */
	int n2_index = 6; /* n2 index in the species names */
	float oxygen_mf; /* species mass fractions */
	float water_mf;
	float methane_mf;
	float co_mf;
	float co2_mf;
	float h2_mf;
	float n2_mf;
	float mole_frac_o2; /* oxygen mole fraction */
	real conversion; /* conversion from TGA */
	real cumulative_loss; /* cumulative char mass loss */
	real dt = CURRENT_TIMESTEP; /* timestep in seconds */
	real exp_factor_char;
	real char_progress;
	real mass_char;
	real mass_ash;
	
	c_face_loop(c, t, n)
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* boundary face thread */

		if (THREAD_ID(tf) == fluid_zone_ID) /* check if the cell face thread is equivalent to the boundary thread */
		{
			face_t f = C_FACE(c, t, n); /* get face index */
			char_progress = C_UDMI(c, t, 10);
			volume = C_VOLUME(c, t);
			mass_char = ((char_progress * final_density)) * volume;
			mass_ash = 0.02 * mass_char;
			temp_surface = F_T(f, tf); /* surface temperature */
			density_surface = C_R(c, t); /* adjacent cell density at the surface */
			// Retrieve cumulative loss from last time step
			cumulative_loss = C_UDMI(c, t, 90);
			shrinkage = cumulative_loss / (mass_char - mass_ash); /* shrinkage factor */
			if (shrinkage >= 1.0) shrinkage = 1.0; /* check if shrinkage factor is exceeding 1*/
			if (isnan(shrinkage)) shrinkage = 0.0; /* check if shrinkage factor is NaN */
			factor = 1 - shrinkage;
			oxygen_mf = C_YI(c, t, oxygen_index); /* mass fraction of oxidizer */
			water_mf = C_YI(c, t, water_index); /* mass fraction of water */
			methane_mf = C_YI(c, t, methane_index); /* mass fraction of methane */
			co_mf = C_YI(c, t, co_index); /* mass fraction of co */
			co2_mf = C_YI(c, t, co2_index); /* mass fraction of co2 */
			h2_mf = C_YI(c, t, h2_index); /* mass fraction of h2 */
			n2_mf = C_YI(c, t, n2_index); /* mass fraction of n2 */
			float mass_fracs[] = { oxygen_mf, water_mf, methane_mf, co_mf, co2_mf, h2_mf, n2_mf }; /* species mass fraction array */
			mole_frac_o2 = get_mole_fractions(mass_fracs, oxygen_index); /* mole fraction of o2 */
			C_UDMI(c, t, 49) = mole_frac_o2;
			conversion = ((mass_char * factor) - mass_ash) / (mass_char - mass_ash);
			if (isnan(conversion)) conversion = 1; /* check if conversion is NaN */
			C_UDMI(c, t, 50) = conversion;
			exp_factor_char = -(E_char / (R * temp_surface));
			C_UDMI(c, t, 51) = exp_factor_char;
			loss_rate = A_char * pow(e, exp_factor_char) * pow((mole_frac_o2 / 0.205), oxygen_exp_char) * 0.98 * char_progress * final_density * pow(conversion, n_char) * volume * C_UDMI(c, t, 10);

			// Update cumulative loss
			cumulative_loss += loss_rate * dt;
			C_UDMI(c, t, 90) = cumulative_loss;
		
			source = (0.5833 * loss_rate) / (volume * r);
		}
	}

	C_UDMI(c, t, 85) = source; /* store the CO source in UDM */
	ds[eqn] = 0.0;
	return source;
}

/* Char oxidation energy source term */

DEFINE_SOURCE(char_energy_source, c, t, ds, eqn)
{
	real source = 0.0; /* initialize source term to zero */
	real volume; /* define volume variable */
	int n = 0; /* initialize flag for cell face loops */
	//int zone_ID = 2; /* surface ID of wood solid interface */
	real shrinkage, factor;
	real char_progress = 0.95;
	real mass_char = (((1 - char_progress) * initial_density) + (char_progress * final_density)) * surf_vol;
	real mass_ash = 0.02 * mass_char;
	c_face_loop(c, t, n) /* begin looping through faces of a cell */
	{
		Thread* tf = C_FACE_THREAD(c, t, n); /* get global face threads */
		face_t f = C_FACE(c, t, n); /* get global face ID for the boundary cell face */

		if (THREAD_ID(tf) == fluid_zone_ID && C_UDMI(c, t, 10) >= 0.95)
		{
			face_t f = C_FACE(c, t, n); /* get global face ID for the boundary cell face */
			volume = C_VOLUME(c, t); /* get volume of the neighbouring cell*/
			source = (9000 * 1000 * C_UDMI(c, t, 81)) / volume; /* C_UDMI(c,t,16) is the char loss distribution. 9 MJ/kg is the char oxidation energy source term for partial char oxidation to CO */
			C_UDMI(c, t, 38) = source;
		}
	}

	ds[eqn] = 0.0; /* explicit evaluation of source term */
	return source;

}
