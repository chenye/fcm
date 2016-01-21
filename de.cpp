#include "de.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "rand_gen.h"

using namespace std;

CDe::CDe()
{
	n_parameters = 2;
	strcpy(method_name, "DE");
}

CDe::~CDe()
{
}

void CDe::allocate_memory()
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

	//variables for DE
	new_population = new float*[pop_size];
	population = new float*[pop_size];
	for (i=0; i<pop_size; i++){
		new_population[i] = new float[dimension];
		population[i] = new float[dimension];
	}
	fit = new float[pop_size];
	new_fit = new float[pop_size];

	best_solution = new float[dimension];
	trial_vector = new float[dimension];
}

void CDe::release_memory()
{
	int i;

	delete[] composed_best_solution;
	delete[] solution_lowerbound;
	delete[] solution_upperbound;
	delete[] solution_range;

	for (i=0; i<pop_size; i++){
		delete[] new_population[i];
		delete[] population[i];
	}
	delete[] new_population;
	delete[] population;
	delete[] fit;
	delete[] new_fit;

	delete[] best_solution;
	delete[] trial_vector;
}

void CDe::write_parameters(FILE *fp)
{
	// fprintf(fp, "pop_size = %d;\n", pop_size);
	fprintf(fp, "crossover_rate = %f;\n", crossover_rate);
	fprintf(fp, "scale_factor = %f;\n", scale_factor);
}

void CDe::set_parameters(float para[])
{
	crossover_rate = para[0];
	scale_factor = para[1];
}

void CDe::init_de()
{
	init_random(seed);
}

void CDe::init_population_rand()
{
	int i, j;
	// generate random initial solutions
	for (i=0; i<pop_size; i++){
		for (j=0; j<dimension; j++){
			population[i][j] = rand0to1() * solution_range[j] + solution_lowerbound[j];
		}
	}
}

void CDe::mutation(int cur_individual)
{
	int r1, r2, i;
	r1 = (int)(rand0to1() * pop_size) % pop_size;
	r2 = (int)(rand0to1() * pop_size) % pop_size;
	while(r2==r1)
		r2 = (int)(rand0to1() * pop_size) % pop_size;

	for (i=0; i<dimension; i++){
		new_population[cur_individual][i] = best_solution[i] + scale_factor*(population[r1][i] - population[r2][i]);
		if (new_population[cur_individual][i] > solution_upperbound[i])
			new_population[cur_individual][i] -= solution_range[i];
		else if (new_population[cur_individual][i] < solution_lowerbound[i])
			new_population[cur_individual][i] += solution_range[i];
	}
}

void CDe::crossover(int cur_individual)
{
	int i;
	for (i=0; i<dimension; i++){
		if (rand0to1() > crossover_rate)
			new_population[cur_individual][i] = population[cur_individual][i];
	}
}

float CDe::evaluate(float *solu, int cur_node)
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

void CDe::evaluate_init_population(int cur_node)
{
	int i;
	for (i=0; i<pop_size; i++){
		fit[i] = evaluate(population[i], cur_node);
	}
}

void CDe::selection()
{
	int i;
	for (i=0; i<pop_size; i++){
		if (fit[i] > new_fit[i]){
			memcpy(population[i], new_population[i], dimension*sizeof(float));
			fit[i] = new_fit[i];
		}

	}
}

void CDe::update_best_solution()
{
	int i, id_best=-1;
	for (i=0; i<pop_size; i++) {
		if (fit[i] < best_fit) {
			best_fit = fit[i];
			id_best = i;
		}
	}
	if (id_best!=-1){
		memcpy(best_solution, population[id_best], dimension*sizeof(float));
	}
}

void CDe::compose_best_solu(int cur_node)
{
	int i;
	int n_nodes = fcm->n_nodes;
	for (i=0; i<n_nodes; i++){
		composed_best_solution[i*n_nodes+cur_node] = best_solution[i];
	}
	composed_best_solution[n_nodes*n_nodes+cur_node] = best_solution[n_nodes];

	composed_bestfit += best_fit;
}

int CDe::get_iter_best_id()
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

void CDe::run()
{
}

void CDe::run_decomposed()
{
	int cur_iter, cur_individual, cur_node;

	init_de();
	composed_bestfit = 0.0;
	
	for (cur_node=0; cur_node<fcm->n_nodes; cur_node++){
		best_fit = FLT_MAX;
		init_population_rand();
		evaluate_init_population(cur_node);
		update_best_solution();

		for (cur_iter=0; cur_iter<max_iter; cur_iter++){
			for (cur_individual=0; cur_individual<pop_size; cur_individual++){
				mutation(cur_individual);
				crossover(cur_individual);
				new_fit[cur_individual] = evaluate(new_population[cur_individual], cur_node);
			}
			selection();
			update_best_solution();

			if (save_intermediate_results==true){
				save_iter_best_solution(cur_iter, cur_node, population[get_iter_best_id()]);
			}

			if (cur_iter %500 ==0)
				cout<< "cur_node = "<<cur_node<<" cur_iter = "<<cur_iter<<" best fit = "<<best_fit<<endl;
		}
		compose_best_solu(cur_node);
	}
}
