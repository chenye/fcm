#ifndef _BASELINE_H
#define _BASELINE_H

#include "opt.h"

class CBaseline : public COpt
{
public:
	// int pop_size;
	float crossover_rate;
	float scale_factor;

	float **population;
	float *fit;
	float **new_population;
	float *new_fit;
	float *trial_vector;
	
	float *best_solution;
	float best_fit;

public:
	CBaseline();
	~CBaseline();

	void allocate_memory();
	void release_memory();
	void run();
	void run_decomposed();
	void write_parameters(FILE *fp);
	void set_parameters(float para[]);

	void init_baseline();
	void init_population_rand();
	float evaluate(float *solu, int cur_node);
	void evaluate_init_population(int cur_node);
	void update_best_solution();
	void compose_best_solu(int cur_node);

	void mutation(int cur_individual);
	void crossover(int cur_individual);
	void selection();

	int get_iter_best_id();
};

#endif
