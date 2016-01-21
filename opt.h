#ifndef _OPT_H
#define _OPT_H

#include "fcm.h"
#include <stdio.h>

class COpt
{
public:
	char method_name[100];
	CFcm *fcm;
	long seed;
	int dimension;
	int pop_size;
	int max_iter;
	float *solution_upperbound;
	float *solution_lowerbound;
	float *solution_range;

	float *composed_best_solution;
	float composed_bestfit;

	bool use_sparse;
	float sparse_penalty;

	int n_parameters; // Number of algorithm-specific parameters

	int iter_best_solution_id;
	float *intermediate_fit;
	float *intermediate_model_error;
	bool save_intermediate_results;
public:
	
	virtual void allocate_memory() = 0;
	virtual void release_memory() = 0;

	virtual void run() = 0;
	virtual void run_decomposed() = 0;	//the decompsed approach

	virtual void write_parameters(FILE *fp) = 0;
	virtual void set_parameters(float para[]) = 0;

	void write_result(char filename[]);
	void write_simple_result(char filename[], bool write_new_line_only);

	void init_intermediate_result();
	void save_iter_best_solution(int cur_iter, int cur_node, float *solu);
	void write_intermediate_result_to_file(char filename[]);
	void free_intermediate_result();
};

#endif
