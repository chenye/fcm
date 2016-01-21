#ifndef _FCM_H
#define _FCM_H

#include <stdio.h>

#ifndef PI
	#define PI 3.1415926536
#endif

class CFcm
{
public:
	int n_nodes;
	float map_density;
	float *lambda;	//lambda can be different for different node
	float lambda_in_normalization;

	float **weight;
	float **target_weight;

	//data
	int len_data;
	int num_seq;
private:
	int data_len_per_seq;
public:
	float **data_seq;

	float noise_level;

	//for output
	float spec;
	float sens;
	float ss_mean;
	//structure comparison results
	int true_pos;
	int true_neg;
	int false_pos;
	int false_neg;
	// true positives defined by Stach
	int true_pos_s;
	int true_neg_s;
	int false_pos_s;
	int false_neg_s;

private:
	float *cur_node_value;
	float *new_node_value;

public:

	void read_problem_from_file_old(char *filename, bool normalize=false);	//for DREAM project
	void set_problem(int n_nodes, float map_density, float lambda, int num_seq, int len_data, float noise_level);	//generate target_weight and data_seq

	// data sequence preprocessing
	void normalize_data_seq(int method);

	// Fitness calculation
	void calculate_next_matrix(float *output, float *input, float **weight);
	void calculate_next(float *output, float *input);
	float fcm_fitness(float *x);		//calculate the fitness of an FCM
	float fcm_fitness_partial(float *x, int n);	//calculate only part of the FCM

	// Performance measure calculation
	float get_model_error(float x[]);
	void get_undirected_error(float x[], float threshold);//set true_pos, true_neg, false_pos, false_neg
	float get_auroc(float x[]);
	void get_auc(float x[], float *p_auroc, float *p_aupr);

	// Problem generation and reading
	void generate_random_problem(const char* output_filename, int n_nodes, float map_density, float lambda, int num_seq, int len_data, float noise_level);
	void generate_random_problem(const char* input_filename, const char* output_filename, int num_seq=0, int data_len_per_seq=0);
	void read_problem_from_file(const char* filename);
private:
	void write_problem_to_file(FILE *fp);

	// functions related to reading input file
	void read_n_nodes(char *linebuf, const int len, FILE *fp);
	void read_seq_size(char *linebuf, const int len, FILE *fp);
	void read_lambda(char *linebuf, const int len, FILE *fp);
	void read_target_weight(char *linebuf, const int len, FILE *fp); // number of edges will also be read in this function
	void read_seq_data(char *linebuf, const int len, FILE *fp);
public:

	CFcm();
	CFcm(int n_nodes, float map_density, float lambda, int num_seq, int len_data, float noise_level=0.0);
	CFcm(char filename[], float lambda, bool normalize=false);
	~CFcm();
	void release_fcm_memory();
	void allocate_fcm_memory(); // set len_data and n_nodes before invoking this function

private:
	int get_valid_line(char *line, int len, FILE *fp);
};

inline float calculate_next_and_error(int n_nodes, float **w, float *lambda, float *output, float *input, float *desired_output);
inline float calculate_next_and_error_partial(int n_nodes, int cur_node, float *w, float *input, float desired_output);//lambda is at the end of w
#endif

