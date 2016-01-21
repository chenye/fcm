#include "main_opt.h"
#include "fcm.h"
#include "rand_gen.h"
#include <stdlib.h>
#include <string.h>

// --- for random number generator ---
// #define	DSFMT_MEXP	216091
#include "twister/dSFMT.h"
// --- for random number generator ---

const int n_net = 15;
const char net_name[n_net][100]={"test-n10-1", "test-n10-2", "test-n10-3", "test-n10-4", "test-n10-5",
									"test-n50-1", "test-n50-2", "test-n50-3", "test-n50-4", "test-n50-5",
									"test-n100-1", "test-n100-2", "test-n100-3", "test-n100-4", "test-n100-5"};

void random_number_generator_test()
{
	int i;

	dsfmt_gv_init_gen_rand(1);
	for (i=0; i<100; i++){
		printf("%f\n", dsfmt_gv_genrand_close_open());
	}
}

void generate_random_problem()
{
	// NOTE: set n_seq and n_data_per_seq in this function affects only the output filename
	// NOTE: change the numbers in the file indicated by "net_filename"
	// NOTE: in order to change the actual number of sequences and number of data per seq
	CFcm *fcm;
	fcm = new CFcm();
	int i, j;
	char seq_filename[100];
	char net_filename[100];
	int n_seq = 5;
	int n_data_per_seq = 10;
	long cur_seed = 1;

	for (i=0; i<n_net; i++){
		cur_seed = (i+1)*100;
		init_random(cur_seed);

		sprintf(seq_filename, "%s-seq-%d-%d.txt", net_name[i], n_seq, n_data_per_seq);
		sprintf(net_filename, "%s.txt", net_name[i]);
		fcm->generate_random_problem(net_filename, seq_filename, n_seq, n_data_per_seq);
	}
}

void generate_random_problem(char net_name[], int n_seq, int n_data_per_seq)
{
	CFcm *fcm;
	fcm = new CFcm();
	char seq_filename[1000];
	char net_filename[1000];
	long cur_seed = 100;

	init_random(cur_seed);
	sprintf(seq_filename, "%s-seq-%d-%d.txt", net_name, n_seq, n_data_per_seq);
	sprintf(net_filename, "%s.txt", net_name);
	fcm->generate_random_problem(net_filename, seq_filename, n_seq, n_data_per_seq);
}

int main(int argc, char* argv[])
{
	// generate_random_problem("test-n8-ecoli-sos", 5,10);
	// generate_random_problem("test-n8-ecoli-sos", 10,10);
	// generate_random_problem("test-n8-ecoli-sos", 20,10);
	// return 0;
	/*
	generate_random_problem("test-n50-1", 100,10);
	generate_random_problem("test-n50-2", 100,10);
	generate_random_problem("test-n50-3", 100,10);
	generate_random_problem("test-n50-4", 100,10);
	generate_random_problem("test-n50-5", 100,10);

	return 0;
	generate_random_problem("test-n200-1", 50,10);
	generate_random_problem("test-n300-1", 50,10);
	return 0;
	*/

	// random_number_generator_test();
	// return 0;
	// generate_random_problem();
	// return 0;
	int alg_id;
	float sparse_penalty;
	int normalization_method;
	float lambda_max, lambda_min, lambda_norm;

	switch (argc) {
	case 3:
		alg_id = atoi(argv[1]);
		if (alg_id==0)
			return -1;
		printf("alg_id = %d, data file: %s\n", alg_id, argv[2]);
		try_one_parameter(alg_id, argv[2]);
		break;

	case 4:
		if (strcmp(argv[1], "dp")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			printf("alg_id = %d, data file: %s\n", alg_id, argv[3]);
			try_diff_parameters(alg_id, argv[3]);
		}
		break;
	case 5:
		if (strcmp(argv[1], "sparse")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			sparse_penalty = atof(argv[4]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f\n", alg_id, argv[3], sparse_penalty);
			try_one_parameter(alg_id, argv[3], sparse_penalty);
		}
		break;
	case 6:
		if (strcmp(argv[1], "sparse")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d\n", alg_id, argv[3], sparse_penalty, normalization_method);
			try_one_parameter(alg_id, argv[3], sparse_penalty, 5.0, 5.0, normalization_method);
		}
		break;
	case 8:
		if (strcmp(argv[1], "sparse")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			lambda_min = atof(argv[6]);
			lambda_max = atof(argv[7]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d, lambda_min=%f, lambda_max=%f\n", alg_id, argv[3], sparse_penalty, normalization_method, lambda_min, lambda_max);
			try_one_parameter(alg_id, argv[3], sparse_penalty, lambda_min, lambda_max, normalization_method);
		}
		else if (strcmp(argv[1], "sparsedp")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			lambda_min = atof(argv[6]);
			lambda_max = atof(argv[7]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d, lambda_min=%f, lambda_max=%f\n", alg_id, argv[3], sparse_penalty, normalization_method, lambda_min, lambda_max);
			try_diff_parameters(alg_id, argv[3], sparse_penalty, lambda_min, lambda_max, normalization_method);

		}
		else if (strcmp(argv[1], "sparseint")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>4)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			lambda_min = atof(argv[6]);
			lambda_max = atof(argv[7]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d, lambda_min=%f, lambda_max=%f\n", alg_id, argv[3], sparse_penalty, normalization_method, lambda_min, lambda_max);
			try_one_parameter(alg_id, argv[3], sparse_penalty, lambda_min, lambda_max, normalization_method, 2.0, true);
		}
		break;
	case 9:
		if (strcmp(argv[1], "sparse")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>5)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			lambda_min = atof(argv[6]);
			lambda_max = atof(argv[7]);
			lambda_norm = atof(argv[8]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d, lambda_min=%f, lambda_max=%f, lambda_norm=%f\n", alg_id, argv[3], sparse_penalty, normalization_method, lambda_min, lambda_max, lambda_norm);
			try_one_parameter(alg_id, argv[3], sparse_penalty, lambda_min, lambda_max, normalization_method, lambda_norm);
		}
		else if (strcmp(argv[1], "sparseint")==0){
			alg_id = atoi(argv[2]);
			if (alg_id<1 || alg_id>5)
				return -1;
			sparse_penalty = atof(argv[4]);
			normalization_method = atoi(argv[5]);
			lambda_min = atof(argv[6]);
			lambda_max = atof(argv[7]);
			lambda_norm = atof(argv[8]);
			// lambda max is set to 5.0, lambda min is set to 0.0 by default for the dream dataset
			printf("alg_id = %d, data file: %s, sparse_penalty=%f, normalization_method=%d, lambda_min=%f, lambda_max=%f, lambda_norm=%f\n", alg_id, argv[3], sparse_penalty, normalization_method, lambda_min, lambda_max, lambda_norm);
			try_one_parameter(alg_id, argv[3], sparse_penalty, lambda_min, lambda_max, normalization_method, lambda_norm, true);
		}
		break;
	default:
		printf("Wrong parameters.\n");
		for (int i=0; i<argc; i++)
			printf("%s\n", argv[i]);
		printf("Examples of correct command lines:\n");
		printf("opt4fcm7 alg_id data_file_name_without_file_name_extension\n");
		printf("opt4fcm7 dp alg_id data_file_name:  experiments with different parameters\n");
		printf("opt4fcm7 sparse alg_id data_file_name sparse_penalty: set sparse_penalty value\n");
		printf("opt4fcm7 sparse alg_id data_file_name sparse_pently norm_method lam_min lam_max: choose normalization method and set lambda min and max\n");
	}

	return 0;
}
