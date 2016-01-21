#include "opt.h"
#include "fcm.h"
#include "acor.h"
#include "de.h"
#include "pso.h"
#include "rcga.h"
#include "baseline.h"

COpt* create_new_copt(int method)
{
	COpt *opt;
	switch(method){
	case 1:
		opt = new CDe();
		break;
	case 2:
		opt = new CAcoR();
		break;
	case 3:
		opt = new CPso();
		break;
	case 4:
		opt = new CRcga();
		break;
	case 5:
		opt = new CBaseline();
	}
	return opt;
}

void set_default_parameters(int method, float para[])
{
	float tmp_omega = 0.7298;
	float tmp_phi1 = 1.49618;
	float tmp_phi2 = 1.49618;


	float tmp_crossover_rate = 0.8;
	float tmp_scale_factor = 0.5;

	float tmp_q=0.6;
	float tmp_xi=0.6;
	int tmp_k_best = 50;
	float tmp_p0 = 0.0;

	// parameters for RCGA
	float tmp_rcga_crossover_rate = 0.6; //old: 0.8
	float tmp_rcga_mutation_rate = 0.3; // old: 0.2
	// float tmp_rcga_fit_function_a = 0.1;
	float tmp_rcga_mutation_b = 4;

	switch(method){
	case 1:
		para[0] = tmp_crossover_rate;
		para[1] = tmp_scale_factor;
		para[2] = 0.0;
		para[3] = 0.0;
		break;
	case 2:
		para[0] = tmp_k_best;
		para[1] = tmp_q;
		para[2] = tmp_xi;
		para[3] = tmp_p0;
		break;
	case 3:
		para[0] = tmp_omega;
		para[1] = tmp_phi1;
		para[2] = tmp_phi2;
		break;
	case 4:
		para[0] = tmp_rcga_crossover_rate;
		para[1] = tmp_rcga_mutation_rate;
		para[2] = tmp_rcga_mutation_b;
		break;
	case 5:
		// no parameter
		break;
	}
}

void try_one_parameter(int method, char seq_name[], float sparse_penalty, float lambda_min, float lambda_max, int normalization_method, float lambda_in_normalization, bool save_intermediate_results)
{
	int i, j, k;
	COpt *opt;
	char simple_result_filename[200];
	char result_filename[200];
	char intermediate_result_filename[200];
	char input_data_filename[200];
	
	int tmp_pop_size = 100;

	float para[4];

	for (k=1; k<=20; k++){
		// STEP 1: create new COpt
		opt = create_new_copt(method);
		
		// STEP 2: set up FCM
		opt->fcm = new CFcm();
		sprintf(input_data_filename, "%s.txt", seq_name);
		opt->fcm->lambda_in_normalization = lambda_in_normalization;
		opt->fcm->read_problem_from_file(input_data_filename);
		opt->fcm->normalize_data_seq(normalization_method);

		// STEP 3: set parameters
		set_default_parameters(method, para);
		opt->set_parameters(para);

		opt->pop_size = tmp_pop_size;
		opt->max_iter = 0;//15000*100/(opt->pop_size);
		opt->seed = k*100;
		if (sparse_penalty != 0)
			opt->use_sparse = true;
		else
			opt->use_sparse = false;
		opt->sparse_penalty = sparse_penalty;
		sprintf(simple_result_filename, "result-s-%s-para1-np%d-sp%0.2f-%s-norm%d-lam%1.1f-%1.1f-%1.2f.txt", opt->method_name, tmp_pop_size, sparse_penalty, seq_name, normalization_method, lambda_min, lambda_max, lambda_in_normalization);
		sprintf(result_filename, "result-%s-para1-np%d-sp%0.2f-%s-norm%d-lam%1.1f-%1.1f-%1.2f.txt", opt->method_name, tmp_pop_size, sparse_penalty, seq_name, normalization_method, lambda_min, lambda_max, lambda_in_normalization);

		// STEP 4: allocate memory
		opt->allocate_memory();
		opt->save_intermediate_results=save_intermediate_results;
		if (save_intermediate_results==true)
			opt->init_intermediate_result();

		// STEP 5: modify settings
		opt->solution_lowerbound[ opt->fcm->n_nodes ] = lambda_min;
		opt->solution_upperbound[ opt->fcm->n_nodes ] = lambda_max;
		opt->solution_range[ opt->fcm->n_nodes ] = lambda_max - lambda_min;

		// STEP 6: run the algorithm
		opt->run_decomposed();

		// STEP 7: write the results and release the memory
		opt->write_result(result_filename);
		opt->write_simple_result(simple_result_filename, false);

		if (save_intermediate_results==true){
			sprintf(intermediate_result_filename, "result-intermediate-%s-para1-np%d-sp%0.2f-%s-norm%d-lam%1.1f-%1.1f-%1.2f.txt", opt->method_name, tmp_pop_size, sparse_penalty, seq_name, normalization_method, lambda_min, lambda_max, lambda_in_normalization);
			opt->write_intermediate_result_to_file(intermediate_result_filename);
			opt->free_intermediate_result();
		}

		delete opt->fcm;
		opt->release_memory();
		delete opt;
	}
	opt->write_simple_result(simple_result_filename, true);
}




void try_diff_parameters(int method, char seq_name[], float sparse_penalty, float lambda_min, float lambda_max, int normalization_method)
{
	int i, j, k;
	COpt *opt;
	char simple_result_filename[200];
	char result_filename[200];
	char input_data_filename[200];
	
	int tmp_pop_size = 100;

	float para[4];

	set_default_parameters(method, para);
	for (i=1; i<=9; i++){
		for (j=1; j<=9; j++){
				para[0] = (float)i/10.0;  //////////////////////////// Set parameters
				para[1] = (float)j/10.0;

				for (k=1; k<=20; k++){
				// STEP 1: create new COpt
				opt = create_new_copt(method);
		
				// STEP 2: set up FCM
				opt->fcm = new CFcm();
				sprintf(input_data_filename, "%s.txt", seq_name);
				opt->fcm->read_problem_from_file(input_data_filename);
				opt->fcm->normalize_data_seq(normalization_method);

				// STEP 3: set parameters
				opt->set_parameters(para);

				opt->pop_size = tmp_pop_size;
				opt->max_iter = 15000*100/(opt->pop_size);
				opt->seed = k*100;
				if (sparse_penalty != 0)
					opt->use_sparse = true;
				else
					opt->use_sparse = false;
				opt->sparse_penalty = sparse_penalty;
				sprintf(simple_result_filename, "result-s-%s-diffpara-np%d-sp%0.2f-%s-norm%d-lam%1.1f-%1.1f.txt", opt->method_name, tmp_pop_size, sparse_penalty, seq_name, normalization_method, lambda_min, lambda_max);
				sprintf(result_filename, "result-%s-diffpara-np%d-sp%0.2f-%s-norm%d-lam%1.1f-%1.1f.txt", opt->method_name, tmp_pop_size, sparse_penalty, seq_name, normalization_method, lambda_min, lambda_max);

				// STEP 4: allocate memory
				opt->allocate_memory();

				// STEP 5: modify settings
				opt->solution_lowerbound[ opt->fcm->n_nodes ] = lambda_min;
				opt->solution_upperbound[ opt->fcm->n_nodes ] = lambda_max;
				opt->solution_range[ opt->fcm->n_nodes ] = lambda_max - lambda_min;

				// STEP 6: run the algorithm
				opt->run_decomposed();

				// STEP 7: write the results and release the memory
				opt->write_result(result_filename);
				opt->write_simple_result(simple_result_filename, false);

				delete opt->fcm;
				opt->release_memory();
				delete opt;
			}
			opt->write_simple_result(simple_result_filename, true);

		}
	}
}



void try_diff_parameters(int method, char seq_name[])
{
	int i, j, k;
	COpt *opt;
	char simple_result_filename[200];
	char result_filename[200];
	char input_data_filename[200];

	float tmp_sparse_penalty = 0.0;
	int tmp_pop_size = 100;

	float para[4];
	set_default_parameters(method, para);

	for (i=1; i<=9; i++){
		for (j=1; j<=9; j++){
			para[0] = (float)i/10.0;
			para[1] = (float)j/10.0;

			for (k=1; k<=10; k++){
				// STEP 1: create new COpt
				opt = create_new_copt(method);
		
				// STEP 2: set up FCM
				// srand(k);
				// opt->fcm = new CFcm(n_nodes, density, lambda, n_seq, n_seq*n_len_per_seq, tmp_noise_level);
				opt->fcm = new CFcm();
				sprintf(input_data_filename, "%s.txt", seq_name);
				opt->fcm->read_problem_from_file(input_data_filename);

				// STEP 3: set parameters
				// set_default_parameters(method, para);
				opt->set_parameters(para);

				opt->pop_size = tmp_pop_size;
				opt->max_iter = 15000*100/(opt->pop_size);
				opt->seed = k*100;
				opt->use_sparse = false;
				opt->sparse_penalty = tmp_sparse_penalty;
				sprintf(simple_result_filename, "result-s-%s-diffpara-np%d-sp%0.2f-%s.txt", opt->method_name, tmp_pop_size, tmp_sparse_penalty, seq_name);
				sprintf(result_filename, "result-%s-diffpara-np%d-sp%0.2f-%s.txt", opt->method_name, tmp_pop_size, tmp_sparse_penalty, seq_name);

				// STEP 4: allocate memory
				opt->allocate_memory();

				// STEP 5: modify settings
				opt->solution_lowerbound[ opt->fcm->n_nodes ] = 5.0;
				opt->solution_upperbound[ opt->fcm->n_nodes ] = 5.0;
				opt->solution_range[ opt->fcm->n_nodes ] = 0.0;

				// STEP 6: run the algorithm
				opt->run_decomposed();

				// STEP 7: write the results and release the memory
				opt->write_result(result_filename);
				opt->write_simple_result(simple_result_filename, false);

				delete opt->fcm;
				opt->release_memory();
				delete opt;
			}
			opt->write_simple_result(simple_result_filename, true);
		}
	}
}
