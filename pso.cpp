#include "pso.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "rand_gen.h"

using namespace std;

CPso::CPso()
{
	n_parameters = 3;
	strcpy(method_name, "PSO");
}

CPso::~CPso()
{
}

void CPso::allocate_memory()
{
	int i;
	int n_nodes = fcm->n_nodes;

	dimension = n_nodes+1; //an additional dimension is used to store lambda
	
	// variables from parent class
	composed_best_solution = new float[dimension*n_nodes];
	solution_lowerbound = new float[dimension];
	solution_upperbound = new float[dimension];
	solution_range = new float[dimension];
	for (i=0; i<dimension; i++){
		solution_lowerbound[i]=-1.0;
		solution_upperbound[i]=1.0;
		solution_range[i] = 2.0;
	}

	//variables for PSO
	new_population = new float*[pop_size];
	population = new float*[pop_size];
	pbest = new float*[pop_size];
	velocity = new float*[pop_size];
	for (i=0; i<pop_size; i++){
		new_population[i] = new float[dimension];
		population[i] = new float[dimension];
		pbest[i] = new float[dimension];
		velocity[i] = new float[dimension];
	}
	fit = new float[pop_size];
	new_fit = new float[pop_size];
	pbest_fit = new float[pop_size];

	vmax = new float[dimension];
}

void CPso::release_memory()
{
	int i;

	delete[] composed_best_solution;
	delete[] solution_lowerbound;
	delete[] solution_upperbound;
	delete[] solution_range;

	for (i=0; i<pop_size; i++){
		delete[] new_population[i];
		delete[] population[i];
		delete[] pbest[i];
		delete[] velocity[i];
	}
	delete[] new_population;
	delete[] population;
	delete[] pbest;
	delete[] velocity;
	delete[] fit;
	delete[] new_fit;
	delete[] pbest_fit;

	delete[] vmax;
	
}

void CPso::write_parameters(FILE *fp)
{
	// fprintf(fp, "pop_size = %d;\n", pop_size);
	fprintf(fp, "omega = %f;\n", omega);
	fprintf(fp, "phi1 = %f;\n", phi1);
	fprintf(fp, "phi2 = %f;\n", phi2);
}

void CPso::set_parameters(float para[])
{
	omega = para[0];
	phi1 = para[1];
	phi2 = para[2];
}

void CPso::init_pso()
{
	int i;
	
	init_random(seed);
	
	// set vmax
	for (i=0; i<dimension; i++){
		vmax[i] = solution_range[i];
	}
}

void CPso::init_population_rand()
{
	int i, j;
	// generate random initial solutions
	for (i=0; i<pop_size; i++){
		for (j=0; j<dimension; j++){
			population[i][j] = rand0to1() * solution_range[j] + solution_lowerbound[j];
		}
	}
}

void CPso::compose_best_solu(int cur_node)
{
	int i;
	int n_nodes = fcm->n_nodes;
	for (i=0; i<n_nodes; i++){
		composed_best_solution[i*n_nodes+cur_node] = pbest[gbest_id][i];
	}
	composed_best_solution[n_nodes*n_nodes+cur_node] = pbest[gbest_id][i];

	composed_bestfit += pbest_fit[gbest_id];
}

void CPso::init_pbest_fit()
{
	int i;
	for (i=0; i<pop_size; i++){
		pbest_fit[i] = FLT_MAX;
	}
	gbest_id = 0;
}

int CPso::get_iter_best_id()
{
	int i, id;
	float iter_best_fit = fit[0];
	id = 0;
	for (i=1; i<pop_size; i++){
		if (fit[i]<iter_best_fit){
			iter_best_fit = fit[i];
			id = i;
		}
	}

	return id;
}

float CPso::evaluate(float *solu, int cur_node)
{
	float temp_fit;
	temp_fit = fcm->fcm_fitness_partial(solu, cur_node);

	int i;
	if (use_sparse){
		for (i=0; i<dimension-1; i++){
			temp_fit += sparse_penalty * fabs(solu[i]/(dimension-1));
		}
	}

	return temp_fit;
}

void CPso::evaluate_init_population(int cur_node)
{
	int i;
	for (i=0; i<pop_size; i++){
		fit[i] = evaluate(population[i], cur_node);
	}
}

void CPso::update_best_solution()
{
	int i, id_best = -1;
	for (i=0; i<pop_size; i++){
		if (fit[i] < pbest_fit[i]){
			pbest_fit[i] = fit[i];
			memcpy(pbest[i], population[i], dimension*sizeof(float));
			if (fit[i] < pbest_fit[gbest_id])
				gbest_id = i;
		}
	}
}

void CPso::init_velocity_rand()
{
	int i, j;
	for (i=0; i<pop_size; i++){
		for (j=0; j<dimension; j++){
			velocity[i][j] = rand0to1() * solution_range[j] + solution_lowerbound[j];
		}
	}
}

void CPso::particle_move(int pid)
{
	int i;
	float temp_v;

	for (i=0; i<dimension; i++){
		temp_v = omega*velocity[pid][i]
				+ phi1*rand0to1()*(pbest[pid][i] - population[pid][i])
				+ phi2*rand0to1()*(pbest[gbest_id][i] - population[pid][i]);
		if (temp_v < vmax[i] && temp_v > -vmax[i])
			velocity[pid][i] = temp_v;
		else {
			if (temp_v >= 0)
				velocity[pid][i] = vmax[i];
			else
				velocity[pid][i] = -vmax[i];
		}
		population[pid][i] = population[pid][i] + velocity[pid][i];
		if (population[pid][i] > solution_upperbound[i])
			population[pid][i] -= solution_range[i];
		else if (population[pid][i] < solution_lowerbound[i])
			population[pid][i] += solution_range[i];
	}
}

void CPso::run()
{
}

void CPso::run_decomposed()
{
	int cur_iter, cur_particle, cur_node;

	init_pso();
	composed_bestfit = 0.0;

	for (cur_node=0; cur_node<fcm->n_nodes; cur_node++){
		init_pbest_fit();
		init_population_rand();
		evaluate_init_population(cur_node);
		update_best_solution();

		init_velocity_rand();

		for (cur_iter=0; cur_iter<max_iter; cur_iter++){
			for (cur_particle=0; cur_particle<pop_size; cur_particle++){
				// modify solutions
				particle_move(cur_particle);
				fit[cur_particle] = evaluate( population[cur_particle], cur_node );
			}
			update_best_solution();

			if (save_intermediate_results==true){
				save_iter_best_solution(cur_iter, cur_node, population[get_iter_best_id()]);
			}

			if (cur_iter % 500 == 0)
				cout<< "cur_node = "<<cur_node<<" cur_iter = "<<cur_iter<<" best fit = "<<pbest_fit[gbest_id] << endl;
		}
		compose_best_solu(cur_node);
	}
}
