#include "rcga.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <iostream>
#include "rand_gen.h"

using namespace std;

CRcga::CRcga()
{
	n_parameters = 3;
	strcpy(method_name, "RCGA");
}

CRcga::~CRcga()
{
}

void CRcga::allocate_memory()
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

	//variables for RCGA
	new_population = new float*[pop_size];
	population = new float*[pop_size];
	for (i=0; i<pop_size; i++){
		new_population[i] = new float[dimension];
		population[i] = new float[dimension];
	}
	fit = new float[pop_size];
	new_fit = new float[pop_size];
	trans_fit = new float[pop_size];

	best_solution = new float[dimension];
	// trial_vector = new float[dimension];
}

void CRcga::release_memory()
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
	delete[] trans_fit;

	delete[] best_solution;
	// delete[] trial_vector;
}

int CRcga::get_iter_best_id()
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

void CRcga::write_parameters(FILE *fp)
{
	// fprintf(fp, "pop_size = %d;\n", pop_size);
	fprintf(fp, "crossover_rate = %f;\n", crossover_rate);
	fprintf(fp, "mutation_rate = %f;\n", mutation_rate);
	fprintf(fp, "mutation_b = %f;\n", mutation_b);
	// fprintf(fp, "scale_factor = %f;\n", scale_factor);
}

void CRcga::set_parameters(float para[])
{
	crossover_rate = para[0];
	mutation_rate = para[1];
	mutation_b = para[2];
}

void CRcga::init_rcga()
{
	init_random(seed);
}

void CRcga::init_population_rand()
{
	int i, j;
	// generate random initial solutions
	for (i=0; i<pop_size; i++){
		for (j=0; j<dimension; j++){
			population[i][j] = rand0to1() * solution_range[j] + solution_lowerbound[j];
		}
	}
}

/*
void CRcga::crossover(int cur_individual)
{
	int i;
	for (i=0; i<dimension; i++){
		if (rand0to1() > crossover_rate)
			new_population[cur_individual][i] = population[cur_individual][i];
	}
}
*/

float CRcga::evaluate(float *solu, int cur_node)
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

void CRcga::evaluate_init_population(int cur_node)
{
	int i;
	for (i=0; i<pop_size; i++){
		fit[i] = evaluate(population[i], cur_node);
	}
}

/*
void CRcga::prepare_selection()
{
	int i;
	sum_fit = 0.0;
	for (i=0; i<pop_size; i++){
		trans_fit[i] = 1.0/(fit_function_a*fit[i]+1);
		sum_fit = sum_fit + trans_fit[i];
	}
}
*/
void CRcga::prepare_selection() // 2013-09-26
{
	int i;
	float avg_fit;
	
	sum_fit = 0.0;
	for (i=0; i<pop_size; i++){
		sum_fit += fit[i];
	}
	
	avg_fit = sum_fit / pop_size;
	if (avg_fit == 0)
		avg_fit = 0.000001;
		
	sum_fit = 0.0;
	for (i=0; i<pop_size; i++) {
		trans_fit[i] = avg_fit/(fit[i]+0.000001); // 2013-10-01
		sum_fit += trans_fit[i];
	}
}

int CRcga::select()
{
	float rand_sum;
	float partial_sum_fit;
	int cur_ind = 0;

	partial_sum_fit = trans_fit[0];
	rand_sum = rand0to1() * sum_fit;
	while(rand_sum < partial_sum_fit && cur_ind < pop_size-1){   // don't need to try the last individual
		partial_sum_fit = partial_sum_fit + trans_fit[cur_ind];
		cur_ind ++;
	}

	return cur_ind;
}

int CRcga::select2()
{
	int p1, p2;
	p1 = int(rand0to1()*pop_size) % pop_size;
	while ( (p2 = int(rand0to1()*pop_size) % pop_size) == p1 )
		;

	if (fit[p1] < fit[p2])
		return p1;
	else
		return p2;
	
}

int CRcga::select3()
{
	int p1, p2, p3;
	p1 = int(rand0to1()*pop_size) % pop_size;
	while((p2 = int(rand0to1()*pop_size) % pop_size) == p1)
		;
	while(true){
		p3 = int(rand0to1()*pop_size) % pop_size;
		if (p3!=p1 && p3!=p2)
			break;
	}

	if (fit[p1] > fit[p2])
		p1=p2;
	if (fit[p1] > fit[p3])
		p1=p3;

	return p1;
}

void CRcga::crossover(int p1, int p2, int cur_individual)
{
	// single point crossover
	int i;
	int loc = rand0to1()*(dimension-1)+1;
	int c1, c2;
	c1 = cur_individual;
	c2 = cur_individual+1;
	memcpy(new_population[c1], population[p1], loc*sizeof(float));
	memcpy(&new_population[c1][loc], &population[p2][loc], (dimension-loc)*sizeof(float));
	memcpy(new_population[c2], population[p2], loc*sizeof(float));
	memcpy(&new_population[c2][loc], &population[p1][loc], (dimension-loc)*sizeof(float));
}

void CRcga::copy_individual(int p1, int cur_individual)
{
	memcpy(new_population[cur_individual], population[p1], dimension*sizeof(float));
}

void CRcga::mutation(int cur_indvidual, int cur_iter)
{
	int i;
	float delta;
	for (i=0; i<dimension; i++){
		if (rand0to1() < mutation_rate) { //2013-09-26
			delta = 1 - pow(rand0to1(), pow( 1-(float)cur_iter/max_iter, mutation_b));
			if (rand0to1()<0.5){
				new_population[cur_indvidual][i] += (solution_upperbound[i] - new_population[cur_indvidual][i]) * delta;
			}
			else {
				new_population[cur_indvidual][i] -= (new_population[cur_indvidual][i] - solution_lowerbound[i]) * delta;
			}
		}
	}
}

void CRcga::update_best_solution()
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

void CRcga::compose_best_solu(int cur_node)
{
	int i;
	int n_nodes = fcm->n_nodes;
	for (i=0; i<n_nodes; i++){
		composed_best_solution[i*n_nodes+cur_node] = best_solution[i];
	}
	composed_best_solution[n_nodes*n_nodes+cur_node] = best_solution[n_nodes];

	composed_bestfit += best_fit;
}

void CRcga::update_population()
{
	int i;
	for (i=0; i<pop_size; i++){
		memcpy(population[i], new_population[i], dimension*sizeof(float));
		fit[i] = new_fit[i];
	}
}

void CRcga::run()
{
}

void CRcga::run_decomposed()
{
	int cur_iter, cur_individual, cur_node;
	int p1, p2;

	init_rcga();
	composed_bestfit = 0.0;
	
	for (cur_node=0; cur_node<fcm->n_nodes; cur_node++){
		best_fit = FLT_MAX;
		init_population_rand();
		evaluate_init_population(cur_node);
		update_best_solution();

		for (cur_iter=0; cur_iter<max_iter; cur_iter++){
#if	SELECTION_METHOD == 1
			prepare_selection();
#endif
			for (cur_individual=0; cur_individual<pop_size; cur_individual+=2){
#if SELECTION_METHOD == 1
				p1 = select();
				p2 = select();
#else
#if SELECTION_METHOD == 2
				p1 = select2();
				p2 = select2();
#else
				p1 = select3();
				p2 = select3();
#endif
#endif
				if (rand0to1() < crossover_rate)  // bug fixed 2013-09-29. changed from ">" to "<"
					crossover(p1, p2, cur_individual);
				else {
					// bug-fixed 2013-09-29
					copy_individual(p1, cur_individual);
					copy_individual(p2, cur_individual+1);
				}
				// if (rand0to1() > mutation_rate) //2013-09-26
					mutation(cur_individual, cur_iter);
				// if (rand0to1() > mutation_rate) //2013-09-26
					mutation(cur_individual+1, cur_iter);

				new_fit[cur_individual] = evaluate(new_population[cur_individual], cur_node);
				new_fit[cur_individual+1] = evaluate(new_population[cur_individual+1], cur_node);
			}
			// selection();
			update_population();
			update_best_solution();
			
			if (save_intermediate_results==true){
				save_iter_best_solution(cur_iter, cur_node, population[get_iter_best_id()]);
			}

			if (cur_iter %500 ==0){
				cout<< "cur_node = "<<cur_node<<" cur_iter = "<<cur_iter<<" best fit = "<<best_fit<<endl;
				// TEMP: show average fit
				//float avg_fit = 0.0;
				// float min_fit = fit[0];
				// for (int i=0; i<pop_size; i++){
					//avg_fit += fit[i];
				//	if (min_fit > fit[i])
				//		min_fit = fit[i];
				// }
				//avg_fit /= pop_size;
				// cout<< "cur_node = "<<cur_node<<" cur_iter = "<<cur_iter<<" best fit = "<<best_fit<< " min = "<<min_fit<<endl;
			}
		}
		compose_best_solu(cur_node);
	}
}
