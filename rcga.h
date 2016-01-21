#ifndef _RCGA_H
#define _RCGA_H

// selection method == 2----> using select2()
#define	SELECTION_METHOD	2

#include "opt.h"

class CRcga : public COpt
{
public:
	// int pop_size;
	float crossover_rate;
	float mutation_rate;
	// float fit_function_a;
	// float scale_factor;
	float mutation_b;

	float **population;
	float *fit;
	float *trans_fit;  //transformed fit
	float **new_population;
	float *new_fit;
	// float *trial_vector;
	
	float *best_solution;
	float best_fit;

public:
	CRcga();
	~CRcga();

	void allocate_memory();
	void release_memory();
	void run();
	void run_decomposed();
	void write_parameters(FILE *fp);
	void set_parameters(float para[]);

	void init_rcga();
	void init_population_rand();
	float evaluate(float *solu, int cur_node);
	void evaluate_init_population(int cur_node);
	void update_best_solution();
	void compose_best_solu(int cur_node);

	void mutation(int cur_individual, int cur_iter);
	// void crossover(int cur_individual);

	int select();
	int select2(); //2013-10-02
	int select3(); //2013-10-05
	void crossover(int p1, int p2, int cur_individual);
	void copy_individual(int p1, int cur_individual);
	void update_population();

	int get_iter_best_id();

private:
	float sum_fit;
	void prepare_selection();
};

#endif
