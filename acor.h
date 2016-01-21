#ifndef _ACOR_H
#define _ACOR_H

#include "fcm.h"
#include "opt.h"
#include <stdio.h>

#ifndef PI
	#define PI 3.1415926536
#endif

class CAcoR : public COpt
{
public:
	// CFcm *fcm;

	// long seed;
	// int num_ants;    ////////////// change to pop_size- --
	// int pop_size;
	int k_best;
	// int dimension;
	// int max_iter;

	float p0;
	float q;
	float xi;

	// float *solution_upperbound;
	// float *solution_lowerbound;

	float **bestkant;	//best k ants
	float *bestkantfit;
	int *bestkantrank;
	float **ant;
	float *antfit;
	int *antrank;
	int *temprank;
	float *best_ant_weight;
	float sum_best_ant_weight;

	int get_iter_best_id();

	// float *composed_bestant;   //   change to composed_best_solution
	// float *composed_best_solution;
	// float composed_bestfit;

	// bool use_sparse;
	// float sparse_penalty;
public:
	CAcoR(long seed, bool use_sparse=false, float sparse_penalty=0.0);//2012-05-25 default value changed to 0.0
	CAcoR();
	~CAcoR();

	void allocate_memory();
	void release_memory();
	void write_parameters(FILE *fp); /////////////////////////////////////////////////////////////////////
	void set_parameters(float para[]); /////////////////////////////////////////////////////////////////////
	
	void init_acor();
	void init_ant_rand(int n);
	void run();
	void run_decomposed();	//the decompsed approach
	float evaluate(float *solu, int cur_node);
	void sort_solutions(float *fit, int *rank, int n);
	void calculate_best_ant_weight();
	void build_tour(int cur_ant);
	void update_archive(/*float *bestkantfit, float *antfit, int *bestkantrank, int *antrank, int n1, int n2*/);

	void compose_bestant(int n);	//put the n-th optimization result in the composed_bestant matrix
};

#endif
