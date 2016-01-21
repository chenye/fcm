
#ifndef _PSO_H
#define _PSO_H

#include "opt.h"
#include <stdlib.h>

class CPso : public COpt
{
public:
	// int pop_size;
	float omega;
	float phi1;
	float phi2;

	float **population;
	float *fit;
	float **new_population;
	float *new_fit;

	float *vmax;
	float **velocity;
	float **pbest;
	float *pbest_fit;
	int gbest_id;

public:
	CPso();
	~CPso();

	void allocate_memory();
	void release_memory();
	void run();
	void run_decomposed();
	void write_parameters(FILE *fp);
	void set_parameters(float para[]);

	void init_pso();
	void init_population_rand();
	float evaluate(float *solu, int cur_node);
	void evaluate_init_population(int cur_node);
	void compose_best_solu(int cur_node);

	void init_pbest_fit();
	void update_best_solution();
	void init_velocity_rand();
	void particle_move(int pid);
	
	int get_iter_best_id();
};

#endif
