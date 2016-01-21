#include "fcm.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "rand_gen.h"

#define	FCM_ALWAYS_HAS_SELF_ACTIVATION	1	//2013-08-17
const bool renormalize_after_fuzzification=false;
using namespace std;

CFcm::CFcm()
{
	lambda_in_normalization=2.0;
}

CFcm::CFcm(int n_nodes, float map_density, float lambda, int num_seq, int len_data, float noise_level)
{
	set_problem(n_nodes, map_density, lambda, num_seq, len_data, noise_level);
}

CFcm::CFcm(char *filename, float lambda, bool normalize)
{
	int i;

	read_problem_from_file_old(filename, normalize);
	//this->lambda=lambda;
	for (i=0; i<n_nodes; i++){
		this->lambda[i]=lambda;
	}
	map_density=-1.0;
}

CFcm::~CFcm()
{
	release_fcm_memory();
}

int CFcm::get_valid_line(char *line, int len, FILE *fp)
{
	int r;
	while(fgets(line, len, fp)!=NULL) {
		if (line[0] != '%')
			break;
	}
	return feof(fp);
}

void CFcm::read_n_nodes(char *linebuf, const int len, FILE *fp)
{
	if (get_valid_line(linebuf, len, fp))
		exit(-1);
	n_nodes = atoi(linebuf);
}

void CFcm::read_seq_size(char *linebuf, const int len, FILE *fp)
{
	if (get_valid_line(linebuf, len, fp))
		exit (-1);
	sscanf(linebuf, "%d%d", &num_seq, &data_len_per_seq);
	// num_seq = 30; // overide the num_seq in the text input file
	len_data = num_seq * data_len_per_seq;
}

void CFcm::read_lambda(char *linebuf, const int len, FILE *fp)
{
	int i;
	if (get_valid_line(linebuf, len, fp))
		exit(-1);
	lambda[0] = atof(strtok(linebuf, "\t "));
	for (i=1; i<n_nodes; i++){
		lambda[i] = atof(strtok(NULL, "\t "));
	}
}

void CFcm::read_target_weight(char *linebuf, const int len, FILE *fp)
{
	int i, n_edges, n1, n2;
	float w;

	// read number of edges
	if (get_valid_line(linebuf, len, fp))
		exit(-1);
	n_edges = atoi(linebuf);

	// read weights
	for (i=0; i<n_edges; i++){
		if (get_valid_line(linebuf, len, fp))
			exit(-1);
		sscanf(linebuf, "%d%d%f", &n1, &n2, &w);
		target_weight[n1-1][n2-1] = w;
	}
}

void CFcm::read_seq_data(char *linebuf, const int len, FILE *fp)
{
	int i, j;
	float temp_time_stamp;

	for (i=0; i<len_data; i++){
		if (get_valid_line(linebuf, len, fp))
			exit(-1);
		temp_time_stamp = atof(strtok(linebuf, "\t ")); // ignore time stamp in the current version of this code
		for (j=0; j<n_nodes; j++){
			data_seq[i][j] = atof(strtok(NULL, "\t "));
		}
	}
}

void CFcm::generate_random_problem(const char* input_filename, const char* output_filename, int num_seq, int data_len_per_seq)
{
	const int len = 30000;
	char linebuf[len];
	int r, i, j;
	int file_type = 0;
	FILE *fpout, *fpin;

	fpin = fopen(input_filename, "r");
	fpout = fopen(output_filename, "w");

	r = get_valid_line(linebuf, len, fpin);
	file_type = atoi(linebuf);
	switch (file_type)
	{
	case 1:
		read_n_nodes(linebuf, len, fpin);
		read_seq_size(linebuf, len, fpin);
		
		// override the numbers read from the file, if num_seq and data_len_per_seq == 0;
		if (num_seq!=0 && data_len_per_seq!=0){
			this->num_seq = num_seq;
			this->data_len_per_seq = data_len_per_seq;
			len_data = num_seq * data_len_per_seq;
		}

		// allocate memory
		allocate_fcm_memory();

		read_lambda(linebuf, len, fpin);
		read_target_weight(linebuf, len, fpin);

		/////////////////////////////////////////////////////////////////////////////////////////////
		// assign random weights 
		for (i=0; i<n_nodes; i++){
			for (j=0; j<n_nodes; j++){
				if (target_weight[i][j]!=0) {
					// generate values in interval [0.05, 1] U [-1, -0.05]
					// target_weight[i][j] = ((float)rand()/RAND_MAX)*1.9 - 0.95;
					target_weight[i][j] = rand0to1()*1.9 - 0.95;
					if (target_weight[i][j]>=0)
						target_weight[i][j] += 0.05;
					else
						target_weight[i][j] -= 0.05;
				}
			}
		}

		// generate random sequence data
		for (i=0; i<len_data; i++)
		{
			if (i%(len_data/num_seq)==0){
				for (j=0; j<n_nodes; j++) {
					// data_seq[i][j]=(float)rand()/RAND_MAX;
					data_seq[i][j]=rand0to1();
				}
			}
			else {
				calculate_next_matrix(data_seq[i], data_seq[i-1], target_weight);
			}
		}

		// write generated data to file
		write_problem_to_file(fpout);
		break;
	default:
		exit(-1);
	}

	fclose(fpin);
	fclose(fpout);
}

void CFcm::write_problem_to_file(FILE *fp)
{
	int i, j, k;
	char *str1 = "%%  --- Type of this file\n"
				"%% --- 0 - data file: network structure, weights, data sequences are all provided\n"
				"%% --- 1 - network structure is provided; weights are not provided; initial states are not provided\n";
	char *str2 = "%% --- Number of nodes\n";
	char *str3 = "%% --- Size of sequence data: number of sequences, number of data points per sequence\n";
	char *str4 = "%% --- Lambda\n";
	char *str5 = "%% --- Network Structure\n"
				"%% ----- Number of edges in the network\n";
	char *str6 = "%% ----- Edge list\n";
	char *str7 = "%% --- Data Sequence: time stamp, values for every node\n";

	fprintf(fp, str1);
	fprintf(fp, "%d\n", 0);
	fprintf(fp, str2);
	fprintf(fp, "%d\n", n_nodes);
	fprintf(fp, str3);
	fprintf(fp, "%d\t%d\n", num_seq, data_len_per_seq);

	// lambda
	fprintf(fp, str4);
	for (i=0; i<n_nodes; i++)
		fprintf(fp, "%f ", lambda[i]);
	fprintf(fp, "\n");

	// number of edges
	fprintf(fp, str5);
	int n_edges = 0;
	for (i=0; i<n_nodes; i++){
		for (j=0; j<n_nodes; j++){
			if (target_weight[i][j]!=0)
				n_edges++;
		}
	}
	fprintf(fp, "%d\n", n_edges);

	// edge list
	fprintf(fp, str6);
	for (i=0; i<n_nodes; i++){
		for (j=0; j<n_nodes; j++){
			if (target_weight[i][j]!=0)
				fprintf(fp, "%d\t%d\t%f\n", i+1, j+1, target_weight[i][j]);
		}
	}

	// data sequence
	int n=0;
	fprintf(fp, str7);
	for (i=0; i<num_seq; i++){
		for (j=0; j<data_len_per_seq; j++){
			fprintf(fp, "%d\t", j);
			for (k=0; k<n_nodes; k++)
				fprintf(fp, "%f\t", data_seq[n][k]);
			n++;
			fprintf(fp, "\n");
		}
		fprintf(fp, "%% ---\n");
	}
}

void CFcm::read_problem_from_file(const char* filename)
{
	const int len = 30000;
	char linebuf[len];
	int i, j, r;
	int file_type = 0;

	FILE *fp = fopen(filename, "r");
	r = get_valid_line(linebuf, len, fp);
	file_type = atoi(linebuf);
	switch(file_type)
	{
	case 0:
		// data file
		read_n_nodes(linebuf, len, fp);
		read_seq_size(linebuf, len, fp);

		// allocate memory
		allocate_fcm_memory();

		read_lambda(linebuf, len, fp);
		read_target_weight(linebuf, len, fp);
		read_seq_data(linebuf, len, fp);
		break;
	default:
		if (fp!=NULL)
			fclose(fp);
		exit(-1);
	}

	fclose(fp);
}

void CFcm::normalize_data_seq(int method)
{
	int i, j ,k;
	int start_pos;
	float maxval, minval, meanval;

	switch(method) {
	case 0:
		// no normalization
		break;
	case 1:
		for (i=0; i<n_nodes; i++){
			maxval=data_seq[0][i];
			minval=data_seq[0][i];
			for (j=0; j<len_data; j++){
				if (data_seq[j][i] > maxval)
					maxval=data_seq[j][i];
				else if (data_seq[j][i]<minval)
					minval=data_seq[j][i];
			}
			for (j=0; j<len_data; j++){
				data_seq[j][i]=(data_seq[j][i]-minval) / (maxval-minval);
			}
		}
		break;
	case 2:
		float maxval2, minval2;
		minval2 = 1/(1+exp(-lambda_in_normalization*0.5));// possible bug 2013-11-26
		maxval2 = 1/(1+exp( lambda_in_normalization*0.5));
		for (i=0; i<n_nodes; i++){
			maxval=data_seq[0][i];
			minval=data_seq[0][i];
			for (j=0; j<len_data; j++){
				if (data_seq[j][i] > maxval)
					maxval=data_seq[j][i];
				else if (data_seq[j][i]<minval)
					minval=data_seq[j][i];
			}
			for (j=0; j<len_data; j++){
				data_seq[j][i]=(data_seq[j][i]-minval) / (maxval-minval);
				data_seq[j][i]=1/(1+exp(-lambda_in_normalization*(data_seq[j][i]-0.5)));
				if (renormalize_after_fuzzification==true){
					data_seq[j][i] = (data_seq[j][i]-minval2)/(maxval2-minval2);
				}
			}
		}
		break;
	case 3:
		// normalize individual sequences separately
		for (i=0; i<n_nodes; i++){
			for (j=0; j<num_seq; j++){
				start_pos = j*data_len_per_seq;
				maxval=data_seq[start_pos][i];
				minval=data_seq[start_pos][i];
				for (k=0; k<data_len_per_seq; k++){
					if (data_seq[start_pos+k][i] > maxval)
						maxval=data_seq[start_pos+k][i];
					else if (data_seq[start_pos+k][i]<minval)
						minval=data_seq[start_pos+k][i];
				}
				for (k=0; k<data_len_per_seq; k++){
					data_seq[start_pos+k][i]=(data_seq[start_pos+k][i]-minval) / (maxval-minval);
				}
			}
		}
		break;
	case 4:
		for (i=0; i<n_nodes; i++){
			meanval = 0.0;
			maxval=data_seq[0][i];
			minval=data_seq[0][i];
			for (j=0; j<len_data; j++){
				meanval = meanval + data_seq[j][i];
				if (data_seq[j][i] > maxval)
					maxval=data_seq[j][i];
				else if (data_seq[j][i]<minval)
					minval=data_seq[j][i];
			}
			meanval = meanval / len_data;
			for (j=0; j<len_data; j++){
				data_seq[j][i]=(data_seq[j][i]-minval) / (maxval-minval);
				data_seq[j][i]=1/(1+exp(-lambda_in_normalization*(data_seq[j][i]-meanval)));
			}
		}
		break;
	case 5:
		for (i=0; i<n_nodes; i++){
			maxval=data_seq[0][i];
			minval=data_seq[0][i];
			for (j=0; j<len_data; j++){
				if (data_seq[j][i] > maxval)
					maxval=data_seq[j][i];
				else if (data_seq[j][i]<minval)
					minval=data_seq[j][i];
			}
			for (j=0; j<len_data; j++){
				data_seq[j][i]=(data_seq[j][i]-minval) / (maxval-minval) * lambda_in_normalization + (1-lambda_in_normalization)/2.0;
			}
		}
		break;
	default:
		break;
	}
}

void CFcm::read_problem_from_file_old(char *filename, bool normalize)
{
	int i, j;
	float maxval, minval;

	printf("%s\n", filename);
	FILE *fp_input=fopen(filename, "r");
	fscanf(fp_input, "%d", &n_nodes);
	//dimension = n_nodes * n_nodes;
	//num_layers = dimension * layers_per_dimension;

	lambda=new float[n_nodes];
	weight=new float*[n_nodes];
	for (i=0; i<n_nodes; i++) {
		weight[i]=new float[n_nodes];
	}

	fscanf(fp_input, "%d", &len_data);// this is not the final value of len_data
	fscanf(fp_input, "%d", &num_seq);
	len_data=len_data*num_seq;
	data_seq=new float*[len_data];
	for (i=0; i<len_data; i++){
		data_seq[i]=new float[n_nodes];
	}

	cur_node_value=new float[n_nodes];
	new_node_value=new float[n_nodes];

	//read gold standard
	int num_non_zero_weight;
	int node_start, node_end, node_weight;

	target_weight=new float*[n_nodes];
	for (i=0; i<n_nodes; i++){
		target_weight[i]=new float[n_nodes];
		for (j=0; j<n_nodes; j++){
			target_weight[i][j]=0;
		}
	}
	fscanf(fp_input, "%d", &num_non_zero_weight);
	for (i=0; i<num_non_zero_weight; i++){
		fscanf(fp_input, "%d%d%d", &node_start, &node_end, &node_weight);
		target_weight[node_start-1][node_end-1]=node_weight;
	}

	//read data sequence
	float tmp_read;
	for (i=0; i<len_data; i++)
	{
		//if (i%(len_data/num_seq)==0){
		fscanf(fp_input, "%f", &tmp_read);	//this is the time stamp, ignore
		for (j=0; j<n_nodes; j++) {
			//data_seq[i][j]=(float)rand()/RAND_MAX;
			fscanf(fp_input, "%f", &data_seq[i][j]);
		}
	}

	//normalize
	if (normalize==true){
		for (i=0; i<n_nodes; i++){
			maxval=data_seq[0][i];
			minval=data_seq[0][i];
			for (j=0; j<len_data; j++){
				if (data_seq[j][i] > maxval)
					maxval=data_seq[j][i];
				else if (data_seq[j][i]<minval)
					minval=data_seq[j][i];
			}
			for (j=0; j<len_data; j++){
				data_seq[j][i]=(data_seq[j][i]-minval) / (maxval-minval);
				
				if (data_seq[j][i]<0.1)
					data_seq[j][i]=0;
				else if (data_seq[j][i]>0.9)
					data_seq[j][i]=0.9;
				else
					data_seq[j][i]=(data_seq[j][i]-0.1)/0.8;
				
			}
		}
	}
	fclose(fp_input);
}

void CFcm::set_problem(int n_nodes, float map_density, float lambda, int num_seq, int len_data, float noise_level)
{
	int i, j;

	this->n_nodes=n_nodes;
	this->map_density=map_density;
	//this->lambda=lambda;
	this->lambda=new float[n_nodes];
	for (i=0; i<n_nodes; i++){
		this->lambda[i]=lambda;
	}

	this->num_seq=num_seq;
	this->len_data=len_data;
	this->noise_level=noise_level;

	//dimension = fcm_num_node * fcm_num_node;
	//num_layers= dimension * layers_per_dimension;	//layers_per_dimension


	weight=new float*[n_nodes];
	for (i=0; i<n_nodes; i++){
		weight[i]=new float[n_nodes];
	}


	data_seq=new float*[len_data];
	for (i=0; i<len_data; i++){
		data_seq[i]=new float[n_nodes];
	}

	cur_node_value=new float[n_nodes];
	new_node_value=new float[n_nodes];

	//assign memory for the weights
	target_weight=new float*[n_nodes];
	for (i=0; i<n_nodes; i++){
		target_weight[i]=new float[n_nodes];
		for (j=0; j<n_nodes; j++){
			target_weight[i][j]=0;
		}
	}

	//generate random network
	int cur_num_edges=0;
	int n1, n2;
	while(cur_num_edges<=map_density*n_nodes*n_nodes){
		// n1=int(float(rand())/RAND_MAX*n_nodes);
		// n2=int(float(rand())/RAND_MAX*n_nodes);
		n1=int(rand0to1()*n_nodes);
		n2=int(rand0to1()*n_nodes);
		//if (n1==n2)
		//	continue;
		if (n1>n_nodes-1)
			n1=n_nodes-1;
		if (n2>n_nodes-1)
			n2=n_nodes-1;
		if (target_weight[n1][n2]==0){
			// target_weight[n1][n2]=int(((float)rand()/RAND_MAX -0.5)*2 *1000)/1000.0;
			target_weight[n1][n2]=int((rand0to1() -0.5)*2 *1000)/1000.0;
			if (target_weight[n1][n2]<0.05 && target_weight[n1][n2]>0)
				target_weight[n1][n2]=0.05;
			if (target_weight[n1][n2]>-0.05 && target_weight[n1][n2]<0)
				target_weight[n1][n2]=-0.05;
			cur_num_edges++;
		}
	}

	//show target_fcm_weight
	for (i=0; i<n_nodes; i++) {
		for (j=0; j<n_nodes; j++) {
			cout<<target_weight[i][j]<<endl;
		}
	}

	//generate data from the target network
	for (i=0; i<len_data; i++)
	{
		if (i%(len_data/num_seq)==0){
			for (j=0; j<n_nodes; j++) {
				// data_seq[i][j]=(float)rand()/RAND_MAX;
				data_seq[i][j]=rand0to1();
			}
		}
		else {
			calculate_next_matrix(data_seq[i], data_seq[i-1], target_weight);


		}
	}

	//add noise
	float debug_t1, debug_t2;  //variables for generating noise;
	//add Gaussian noise to the data for every node
	if (noise_level!=0.0){
		for (i=0; i<len_data; i++){
			for (j=0; j<n_nodes; j++){
				do{
					// debug_t1=(double)rand()/RAND_MAX;
					debug_t1=rand0to1();
				}while (debug_t1==0.0);
				// debug_t2=(double)rand()/RAND_MAX;
				debug_t2=rand0to1();
				data_seq[i][j]=noise_level*sqrt(-2.0*log(debug_t1))*cos(2*PI*debug_t2)+data_seq[i][j];
			}
		}
	}
}

#if FCM_ALWAYS_HAS_SELF_ACTIVATION == 1
void CFcm::calculate_next_matrix(float *output, float *input, float **weight)
{
	int i, j;
	float sum;	//sum for current node

	for (i=0; i<n_nodes; i++){
		sum=0.0;
		for (j=0; j<n_nodes; j++){
			sum=sum+weight[j][i]*input[j];
		}
		sum = sum + input[i];  // 2013-08-16
		output[i]=1.0/(1.0+exp(-lambda[i]*sum));
	}
}
#else
void CFcm::calculate_next_matrix(float *output, float *input, float **weight)
{
	int i, j;
	float sum;	//sum for current node

	for (i=0; i<n_nodes; i++){
		sum=0.0;
		for (j=0; j<n_nodes; j++){
			sum=sum+weight[j][i]*input[j];
		}
		output[i]=1.0/(1.0+exp(-lambda[i]*sum));
	}
}
#endif

void CFcm::allocate_fcm_memory()
{
	// This function uses class member variable n_nodes and len_data.
	// Set these two variable before invoke this function
	int i, j;
	
	data_seq = new float*[len_data];
	for (i=0; i<len_data; i++)
		data_seq[i] = new float[n_nodes];
	lambda = new float[n_nodes];
	target_weight = new float*[n_nodes];
	weight = new float*[n_nodes];
	for (i=0; i<n_nodes; i++){
		weight[i] = new float[n_nodes];
		target_weight[i] = new float[n_nodes];
		for (j=0; j<n_nodes; j++){
			weight[i][j] = 0.0;
			target_weight[i][j] = 0.0;
		}
	}
	cur_node_value = new float[n_nodes];
	new_node_value = new float[n_nodes];
}

void CFcm::release_fcm_memory()
{	
	int i;
	
	delete[] lambda;

	for (i=0; i<n_nodes; i++)
		delete[] weight[i];
	delete[] weight;

	// possible bug: has not release memory for target_weight ???

	for (i=0; i<len_data; i++)
		delete[] data_seq[i];
	delete[] data_seq;

	delete[] cur_node_value;
	delete[] new_node_value;
}

float CFcm::fcm_fitness(float *x)
{
	int i, j;
	int index = 0;
	float tmp, error = 0.0, error2=0.0;
	int non_zeros=0;
	float tmp2;

	for (i=0; i<n_nodes; i++){
		for (j=0; j<n_nodes; j++){
			if (x[index]<0.05 && x[index]>-0.05)
				x[index]=0;
			if (x[index]!=0) {
				non_zeros++;
			}
			weight[i][j]=x[index++];
		}
	}

	for (i=0; i<len_data; i++){
		if (i%(len_data/num_seq)!=0) {
			//error+=calculate_next_and_error(n_nodes, weight, lambda, cur_node_value, data_seq[i-1], data_seq[i]);
			error+=calculate_next_and_error(n_nodes, weight, &x[n_nodes*n_nodes], cur_node_value, data_seq[i-1], data_seq[i]);
		}
	}

	return error/(len_data-num_seq)/n_nodes;// + non_zeros*0.05* (float)(num_iters-cur_iter)/num_iters; //error/((len_data-1)*fcm_num_node);
}

float CFcm::fcm_fitness_partial(float *x, int n)
{
	int i;
	float error=0.0;

	for (i=0; i<len_data; i++) {
		if (i%(len_data/num_seq)!=0) {
			error+=calculate_next_and_error_partial(n_nodes, n, x, data_seq[i-1], data_seq[i][n]);
		}
	}

	return error/(len_data-num_seq);
}

float CFcm::get_model_error(float *x)
{
	float model_error=0.0;
	int ind=0;
	for (int i=0; i<n_nodes; i++) {
		for (int j=0; j<n_nodes; j++) {
			model_error=model_error + fabs(target_weight[i][j] - x[ind++]);
		}
	}
	return model_error/n_nodes/n_nodes;
}

void CFcm::get_auc(float x[], float *p_auroc, float *p_aupr)
{
	int *sort_ind;
	int i, j, ind;
	float maxval;
	float *predictx;
	int *targetx;
	int maxind, t;
	int fp, tp, p, n;
	float *rocx, *rocy;
	int rocn;
	float prevx;
	float auroc;

	float *prx, *pry;
	int prn;
	float aupr;

	sort_ind=new int[n_nodes*(n_nodes-1)];
	predictx=new float[n_nodes*(n_nodes-1)];
	targetx=new int[n_nodes*(n_nodes-1)];
	rocx=new float[n_nodes*(n_nodes-1)+1];
	rocy=new float[n_nodes*(n_nodes-1)+1];
	rocn=0;

	prx=new float[n_nodes*(n_nodes-1)+1];
	pry=new float[n_nodes*(n_nodes-1)+1];
	prn=0;

	ind=0;
	for (i=0; i<n_nodes; i++){
		for (j=0; j<n_nodes; j++){
			if (i!=j){
				predictx[ind]=fabs(x[i*n_nodes+j]);
				targetx[ind]=fabs(target_weight[i][j])>0.04?1:0;
				ind++;
			}
		}
	}
	// sort decreasing
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		sort_ind[i]=i;
	}
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		maxval=predictx[sort_ind[i]];
		maxind=i;
		for (j=i+1; j<n_nodes*(n_nodes-1); j++){
			if (maxval<predictx[sort_ind[j]]){
				maxval=predictx[sort_ind[j]];
				maxind=j;
			}
		}
		t=sort_ind[i];
		sort_ind[i]=sort_ind[maxind];
		sort_ind[maxind]=t;
	}

	fp=0;tp=0;n=0;p=0;
	prevx=-1.0;
	aupr=0.0;//for details refer to: Stolovitzky et al. Lessons from the DREAM2 Challenges: A Community Effort to Assess Biological Network Inference.
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		ind=sort_ind[i];
		if (prevx!=predictx[ind]){
			rocx[rocn]=fp;
			rocy[rocn]=tp;
			rocn++;
			prevx=predictx[ind];
		}


		if (targetx[ind]==1){
			tp++;
			//use non-linear interpolation appraoch
			if (i==0)
				aupr=aupr+1;
			else
				aupr=aupr+1-fp*log((float)(i+1)/(float)(i));
		}
		else {
			fp++;
		}
		
		//2012-05-25 moved from before the "if" to after, because aupr is calcualted in a different way than auroc
		//precision recall
		prx[prn]=tp;// recall (it will be divided by P in the end)
		pry[prn]=(float)tp/(i+1); //precision
		prn++;

	}
	rocx[rocn]=fp;
	rocy[rocn]=tp;
	rocn++;

	//calculate auc
	auroc=rocx[0]*rocy[0]/2.0;	//starting from (0,0)
	for (i=1; i<rocn; i++){
		auroc=auroc+ (rocx[i]-rocx[i-1])*(rocy[i]+rocy[i-1])/2.0;
	}

	if (fp==0 || tp==0)
		auroc=-1.0;
	else
		auroc=auroc/fp/tp;
	//

	//for (i=1; i<prn; i++){
	//	aupr=aupr+ (prx[i]-prx[i-1])*(pry[i]+pry[i-1])/2.0;
	//}
	if (tp==0)
		aupr=-1.0;
	else
		aupr=aupr/tp;

	*p_auroc = auroc;
	*p_aupr = aupr;

	delete[] sort_ind;
	delete[] predictx;
	delete[] targetx;
	delete[] rocx;
	delete[] rocy;
	delete[] prx;
	delete[] pry;
}

float CFcm::get_auroc(float x[])
{
	int *sort_ind;
	int i, j, ind;
	float maxval;
	float *predictx;
	int *targetx;
	int maxind, t;
	int fp, tp, p, n;
	float *rocx, *rocy, prevx;
	int rocn;
	float auroc;
	
	sort_ind=new int[n_nodes*(n_nodes-1)];
	predictx=new float[n_nodes*(n_nodes-1)];
	targetx=new int[n_nodes*(n_nodes-1)];
	rocx=new float[n_nodes*(n_nodes-1)+1];
	rocy=new float[n_nodes*(n_nodes-1)+1];
	rocn=0;
	
	ind=0;
	for (i=0; i<n_nodes; i++){
		for (j=0; j<n_nodes; j++){
			if (i!=j){
				predictx[ind]=fabs(x[i*n_nodes+j]);
				targetx[ind]=fabs(target_weight[i][j])>0.04?1:0;
				ind++;
			}
		}
	}
	// sort decreasing
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		sort_ind[i]=i;
	}
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		maxval=predictx[sort_ind[i]];
		maxind=i;
		for (j=i+1; j<n_nodes*(n_nodes-1); j++){
			if (maxval<predictx[sort_ind[j]]){
				maxval=predictx[sort_ind[j]];
				maxind=j;
			}
		}
		t=sort_ind[i];
		sort_ind[i]=sort_ind[maxind];
		sort_ind[maxind]=t;
	}

	fp=0;tp=0;n=0;p=0;
	prevx=-1.0;
	for (i=0; i<n_nodes*(n_nodes-1); i++){
		ind=sort_ind[i];
		if (prevx!=predictx[ind]){
			rocx[rocn]=fp;
			rocy[rocn]=tp;
			rocn++;
			prevx=predictx[ind];
		}
		
		if (targetx[ind]==1){
			tp++;
		}
		else {
			fp++;
		}
	}
	rocx[rocn]=fp;
	rocy[rocn]=tp;
	rocn++;

	//calculate auc
	auroc=rocx[0]*rocy[0]/2.0;	//starting from (0,0)
	for (i=1; i<rocn; i++){
		auroc=auroc+ (rocx[i]-rocx[i-1])*(rocy[i]+rocy[i-1])/2.0;
	}

	if (fp==0 || tp==0)
		auroc=-1.0;
	else
		auroc=auroc/fp/tp;
	//

	delete[] sort_ind;
	delete[] predictx;
	delete[] targetx;
	delete[] rocx;
	delete[] rocy;
	
	return auroc;
}

void CFcm::get_undirected_error(float x[], float threshold)
{
	true_pos=0;
	true_neg=0;
	false_pos=0;
	false_neg=0;

	true_pos_s=0;
	true_neg_s=0;
	false_pos_s=0;
	false_neg_s=0;

	int i, j, ind;
	ind=0;
	for (i=0; i<n_nodes; i++) {
		for (j=0; j<n_nodes; j++) {
			if (fabs(x[ind])>threshold && fabs(target_weight[i][j])>threshold){
				true_neg_s++;
				true_pos++;
			}
			if (fabs(x[ind])>threshold && fabs(target_weight[i][j])<=threshold){
				false_neg_s++;
				false_pos++;
			}
			if (fabs(x[ind])<=threshold && fabs(target_weight[i][j])>threshold){
				false_pos_s++;
				false_neg++;
			}
			if (fabs(x[ind])<=threshold && fabs(target_weight[i][j])<=threshold){
				true_pos_s++;
				true_neg++;
			}
			ind++;
		}
	}

	spec=(float)true_pos_s/(true_pos_s+false_neg_s);
	sens=(float)true_neg_s/(true_neg_s+false_pos_s);
	ss_mean=2.0*spec*sens/(spec+sens);
}


#if FCM_ALWAYS_HAS_SELF_ACTIVATION == 1

inline float calculate_next_and_error(int n_nodes, float **w, float *lambda, float *output, float *input, float *desired_output)
{
	int i, j;
	float sum;	//sum for current node
	float error, t1, t2;

	error=0.0;
	switch(n_nodes){
	case 100:
		for (i=0; i<100; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				w[40][i]*input[40] + w[41][i]*input[41] + w[42][i]*input[42] + w[43][i]*input[43] + w[44][i]*input[44] +
				w[45][i]*input[45] + w[46][i]*input[46] + w[47][i]*input[47] + w[48][i]*input[48] + w[49][i]*input[49] +
				w[50][i]*input[50] + w[51][i]*input[51] + w[52][i]*input[52] + w[53][i]*input[53] + w[54][i]*input[54] +
				w[55][i]*input[55] + w[56][i]*input[56] + w[57][i]*input[57] + w[58][i]*input[58] + w[59][i]*input[59] +
				w[60][i]*input[60] + w[61][i]*input[61] + w[62][i]*input[62] + w[63][i]*input[63] + w[64][i]*input[64] +
				w[65][i]*input[65] + w[66][i]*input[66] + w[67][i]*input[67] + w[68][i]*input[68] + w[69][i]*input[69] +
				w[70][i]*input[70] + w[71][i]*input[71] + w[72][i]*input[72] + w[73][i]*input[73] + w[74][i]*input[74] +
				w[75][i]*input[75] + w[76][i]*input[76] + w[77][i]*input[77] + w[78][i]*input[78] + w[79][i]*input[79] +
				w[80][i]*input[80] + w[81][i]*input[81] + w[82][i]*input[82] + w[83][i]*input[83] + w[84][i]*input[84] +
				w[85][i]*input[85] + w[86][i]*input[86] + w[87][i]*input[87] + w[88][i]*input[88] + w[89][i]*input[89] +
				w[90][i]*input[90] + w[91][i]*input[91] + w[92][i]*input[92] + w[93][i]*input[93] + w[94][i]*input[94] +
				w[95][i]*input[95] + w[96][i]*input[96] + w[97][i]*input[97] + w[98][i]*input[98] + w[99][i]*input[99] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 50:
		for (i=0; i<50; i++){
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				w[40][i]*input[40] + w[41][i]*input[41] + w[42][i]*input[42] + w[43][i]*input[43] + w[44][i]*input[44] +
				w[45][i]*input[45] + w[46][i]*input[46] + w[47][i]*input[47] + w[48][i]*input[48] + w[49][i]*input[49] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 40:
		for (i=0; i<40; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 30:
		for (i=0; i<30; i++){
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
	case 20:
		for (i=0; i<20; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 10:
		for (i=0; i<10; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			if (t2>0.05 || t2<-0.05)
				error+=(t2*t2);
			//else
			//	error+=0.5*(t2*t2);
		}
		break;
	case 5:
		for (i=0; i<5; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] + input[i];  // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	default:
		for (i=0; i<n_nodes; i++){
			sum=0.0;
			for (j=0; j<n_nodes; j++){
				sum=sum+w[j][i]*input[j];
			}
			sum = sum + input[i];   // 2013-08-16
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			//output[i]=t1;
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	}

	return error;
}

inline float calculate_next_and_error_partial(int n_nodes, int cur_node, float *w, float *input, float desired_output)
{
	int i;
	float sum, t1;
	float error;

	switch (n_nodes){
	case 100:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  w[40]*input[40] + w[41]*input[41] + w[42]*input[42] + w[43]*input[43] + w[44]*input[44] +
			  w[45]*input[45] + w[46]*input[46] + w[47]*input[47] + w[48]*input[48] + w[49]*input[49] +
			  w[50]*input[50] + w[51]*input[51] + w[52]*input[52] + w[53]*input[53] + w[54]*input[54] +
			  w[55]*input[55] + w[56]*input[56] + w[57]*input[57] + w[58]*input[58] + w[59]*input[59] +
			  w[60]*input[60] + w[61]*input[61] + w[62]*input[62] + w[63]*input[63] + w[64]*input[64] +
			  w[65]*input[65] + w[66]*input[66] + w[67]*input[67] + w[68]*input[68] + w[69]*input[69] +
			  w[70]*input[70] + w[71]*input[71] + w[72]*input[72] + w[73]*input[73] + w[74]*input[74] +
			  w[75]*input[75] + w[76]*input[76] + w[77]*input[77] + w[78]*input[78] + w[79]*input[79] +
			  w[80]*input[80] + w[81]*input[81] + w[82]*input[82] + w[83]*input[83] + w[84]*input[84] +
			  w[85]*input[85] + w[86]*input[86] + w[87]*input[87] + w[88]*input[88] + w[89]*input[89] +
			  w[90]*input[90] + w[91]*input[91] + w[92]*input[92] + w[93]*input[93] + w[94]*input[94] +
			  w[95]*input[95] + w[96]*input[96] + w[97]*input[97] + w[98]*input[98] + w[99]*input[99] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[100]*sum));//bug 2012-05-24 w[10] --> w[100]
		error=desired_output-t1;
		break;
	case 50:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  w[40]*input[40] + w[41]*input[41] + w[42]*input[42] + w[43]*input[43] + w[44]*input[44] +
			  w[45]*input[45] + w[46]*input[46] + w[47]*input[47] + w[48]*input[48] + w[49]*input[49] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[50]*sum));
		error=desired_output-t1;
		break;
	case 40:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[40]*sum));
		error=desired_output-t1;
		break;
	case 20:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[20]*sum));
		error=desired_output-t1;
		break;
	case 10:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[10]*sum));
		error=desired_output-t1;
		break;
	case 5:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[5]*sum));
		error=desired_output-t1;
		break;
	default:
		sum=0.0;
		for (i=0; i<n_nodes; i++){
			sum+=w[i]*input[i];
		}
		sum = sum + input[cur_node];  // 2013-08-16
		t1=1.0/(1.0+exp(-w[n_nodes]*sum));//bug 2012-05-24 w[100] --> w[n_nodes]
		error=desired_output-t1;
		break;
	}

	return (error*error);
}

#else
inline float calculate_next_and_error(int n_nodes, float **w, float *lambda, float *output, float *input, float *desired_output)
{
	int i, j;
	float sum;	//sum for current node
	float error, t1, t2;

	error=0.0;
	switch(n_nodes){
	case 300:
		for (i=0; i<300; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				w[40][i]*input[40] + w[41][i]*input[41] + w[42][i]*input[42] + w[43][i]*input[43] + w[44][i]*input[44] +
				w[45][i]*input[45] + w[46][i]*input[46] + w[47][i]*input[47] + w[48][i]*input[48] + w[49][i]*input[49] +
				w[50][i]*input[50] + w[51][i]*input[51] + w[52][i]*input[52] + w[53][i]*input[53] + w[54][i]*input[54] +
				w[55][i]*input[55] + w[56][i]*input[56] + w[57][i]*input[57] + w[58][i]*input[58] + w[59][i]*input[59] +
				w[60][i]*input[60] + w[61][i]*input[61] + w[62][i]*input[62] + w[63][i]*input[63] + w[64][i]*input[64] +
				w[65][i]*input[65] + w[66][i]*input[66] + w[67][i]*input[67] + w[68][i]*input[68] + w[69][i]*input[69] +
				w[70][i]*input[70] + w[71][i]*input[71] + w[72][i]*input[72] + w[73][i]*input[73] + w[74][i]*input[74] +
				w[75][i]*input[75] + w[76][i]*input[76] + w[77][i]*input[77] + w[78][i]*input[78] + w[79][i]*input[79] +
				w[80][i]*input[80] + w[81][i]*input[81] + w[82][i]*input[82] + w[83][i]*input[83] + w[84][i]*input[84] +
				w[85][i]*input[85] + w[86][i]*input[86] + w[87][i]*input[87] + w[88][i]*input[88] + w[89][i]*input[89] +
				w[90][i]*input[90] + w[91][i]*input[91] + w[92][i]*input[92] + w[93][i]*input[93] + w[94][i]*input[94] +
				w[95][i]*input[95] + w[96][i]*input[96] + w[97][i]*input[97] + w[98][i]*input[98] + w[99][i]*input[99] +
				w[100][i]*input[100] + w[101][i]*input[101] + w[102][i]*input[102] + w[103][i]*input[103] + w[104][i]*input[104] +
				w[105][i]*input[105] + w[106][i]*input[106] + w[107][i]*input[107] + w[108][i]*input[108] + w[109][i]*input[109] +
				w[110][i]*input[110] + w[111][i]*input[111] + w[112][i]*input[112] + w[113][i]*input[113] + w[114][i]*input[114] +
				w[115][i]*input[115] + w[116][i]*input[116] + w[117][i]*input[117] + w[118][i]*input[118] + w[119][i]*input[119] +
				w[120][i]*input[120] + w[121][i]*input[121] + w[122][i]*input[122] + w[123][i]*input[123] + w[124][i]*input[124] +
				w[125][i]*input[125] + w[126][i]*input[126] + w[127][i]*input[127] + w[128][i]*input[128] + w[129][i]*input[129] +
				w[130][i]*input[130] + w[131][i]*input[131] + w[132][i]*input[132] + w[133][i]*input[133] + w[134][i]*input[134] +
				w[135][i]*input[135] + w[136][i]*input[136] + w[137][i]*input[137] + w[138][i]*input[138] + w[139][i]*input[139] +
				w[140][i]*input[140] + w[141][i]*input[141] + w[142][i]*input[142] + w[143][i]*input[143] + w[144][i]*input[144] +
				w[145][i]*input[145] + w[146][i]*input[146] + w[147][i]*input[147] + w[148][i]*input[148] + w[149][i]*input[149] +
				w[150][i]*input[150] + w[151][i]*input[151] + w[152][i]*input[152] + w[153][i]*input[153] + w[154][i]*input[154] +
				w[155][i]*input[155] + w[156][i]*input[156] + w[157][i]*input[157] + w[158][i]*input[158] + w[159][i]*input[159] +
				w[160][i]*input[160] + w[161][i]*input[161] + w[162][i]*input[162] + w[163][i]*input[163] + w[164][i]*input[164] +
				w[165][i]*input[165] + w[166][i]*input[166] + w[167][i]*input[167] + w[168][i]*input[168] + w[169][i]*input[169] +
				w[170][i]*input[170] + w[171][i]*input[171] + w[172][i]*input[172] + w[173][i]*input[173] + w[174][i]*input[174] +
				w[175][i]*input[175] + w[176][i]*input[176] + w[177][i]*input[177] + w[178][i]*input[178] + w[179][i]*input[179] +
				w[180][i]*input[180] + w[181][i]*input[181] + w[182][i]*input[182] + w[183][i]*input[183] + w[184][i]*input[184] +
				w[185][i]*input[185] + w[186][i]*input[186] + w[187][i]*input[187] + w[188][i]*input[188] + w[189][i]*input[189] +
				w[190][i]*input[190] + w[191][i]*input[191] + w[192][i]*input[192] + w[193][i]*input[193] + w[194][i]*input[194] +
				w[195][i]*input[195] + w[196][i]*input[196] + w[197][i]*input[197] + w[198][i]*input[198] + w[199][i]*input[199] +
				w[200][i]*input[200] + w[201][i]*input[201] + w[202][i]*input[202] + w[203][i]*input[203] + w[204][i]*input[204] +
				w[205][i]*input[205] + w[206][i]*input[206] + w[207][i]*input[207] + w[208][i]*input[208] + w[209][i]*input[209] +
				w[210][i]*input[210] + w[211][i]*input[211] + w[212][i]*input[212] + w[213][i]*input[213] + w[214][i]*input[214] +
				w[215][i]*input[215] + w[216][i]*input[216] + w[217][i]*input[217] + w[218][i]*input[218] + w[219][i]*input[219] +
				w[220][i]*input[220] + w[221][i]*input[221] + w[222][i]*input[222] + w[223][i]*input[223] + w[224][i]*input[224] +
				w[225][i]*input[225] + w[226][i]*input[226] + w[227][i]*input[227] + w[228][i]*input[228] + w[229][i]*input[229] +
				w[230][i]*input[230] + w[231][i]*input[231] + w[232][i]*input[232] + w[233][i]*input[233] + w[234][i]*input[234] +
				w[235][i]*input[235] + w[236][i]*input[236] + w[237][i]*input[237] + w[238][i]*input[238] + w[239][i]*input[239] +
				w[240][i]*input[240] + w[241][i]*input[241] + w[242][i]*input[242] + w[243][i]*input[243] + w[244][i]*input[244] +
				w[245][i]*input[245] + w[246][i]*input[246] + w[247][i]*input[247] + w[248][i]*input[248] + w[249][i]*input[249] +
				w[250][i]*input[250] + w[251][i]*input[251] + w[252][i]*input[252] + w[253][i]*input[253] + w[254][i]*input[254] +
				w[255][i]*input[255] + w[256][i]*input[256] + w[257][i]*input[257] + w[258][i]*input[258] + w[259][i]*input[259] +
				w[260][i]*input[260] + w[261][i]*input[261] + w[262][i]*input[262] + w[263][i]*input[263] + w[264][i]*input[264] +
				w[265][i]*input[265] + w[266][i]*input[266] + w[267][i]*input[267] + w[268][i]*input[268] + w[269][i]*input[269] +
				w[270][i]*input[270] + w[271][i]*input[271] + w[272][i]*input[272] + w[273][i]*input[273] + w[274][i]*input[274] +
				w[275][i]*input[275] + w[276][i]*input[276] + w[277][i]*input[277] + w[278][i]*input[278] + w[279][i]*input[279] +
				w[280][i]*input[280] + w[281][i]*input[281] + w[282][i]*input[282] + w[283][i]*input[283] + w[284][i]*input[284] +
				w[285][i]*input[285] + w[286][i]*input[286] + w[287][i]*input[287] + w[288][i]*input[288] + w[289][i]*input[289] +
				w[290][i]*input[290] + w[291][i]*input[291] + w[292][i]*input[292] + w[293][i]*input[293] + w[294][i]*input[294] +
				w[295][i]*input[295] + w[296][i]*input[296] + w[297][i]*input[297] + w[298][i]*input[298] + w[299][i]*input[299];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 100:
		for (i=0; i<100; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				w[40][i]*input[40] + w[41][i]*input[41] + w[42][i]*input[42] + w[43][i]*input[43] + w[44][i]*input[44] +
				w[45][i]*input[45] + w[46][i]*input[46] + w[47][i]*input[47] + w[48][i]*input[48] + w[49][i]*input[49] +
				w[50][i]*input[50] + w[51][i]*input[51] + w[52][i]*input[52] + w[53][i]*input[53] + w[54][i]*input[54] +
				w[55][i]*input[55] + w[56][i]*input[56] + w[57][i]*input[57] + w[58][i]*input[58] + w[59][i]*input[59] +
				w[60][i]*input[60] + w[61][i]*input[61] + w[62][i]*input[62] + w[63][i]*input[63] + w[64][i]*input[64] +
				w[65][i]*input[65] + w[66][i]*input[66] + w[67][i]*input[67] + w[68][i]*input[68] + w[69][i]*input[69] +
				w[70][i]*input[70] + w[71][i]*input[71] + w[72][i]*input[72] + w[73][i]*input[73] + w[74][i]*input[74] +
				w[75][i]*input[75] + w[76][i]*input[76] + w[77][i]*input[77] + w[78][i]*input[78] + w[79][i]*input[79] +
				w[80][i]*input[80] + w[81][i]*input[81] + w[82][i]*input[82] + w[83][i]*input[83] + w[84][i]*input[84] +
				w[85][i]*input[85] + w[86][i]*input[86] + w[87][i]*input[87] + w[88][i]*input[88] + w[89][i]*input[89] +
				w[90][i]*input[90] + w[91][i]*input[91] + w[92][i]*input[92] + w[93][i]*input[93] + w[94][i]*input[94] +
				w[95][i]*input[95] + w[96][i]*input[96] + w[97][i]*input[97] + w[98][i]*input[98] + w[99][i]*input[99];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 50:
		for (i=0; i<50; i++){
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39] +
				w[40][i]*input[40] + w[41][i]*input[41] + w[42][i]*input[42] + w[43][i]*input[43] + w[44][i]*input[44] +
				w[45][i]*input[45] + w[46][i]*input[46] + w[47][i]*input[47] + w[48][i]*input[48] + w[49][i]*input[49];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 40:
		for (i=0; i<40; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29] +
				w[30][i]*input[30] + w[31][i]*input[31] + w[32][i]*input[32] + w[33][i]*input[33] + w[34][i]*input[34] +
				w[35][i]*input[35] + w[36][i]*input[36] + w[37][i]*input[37] + w[38][i]*input[38] + w[39][i]*input[39];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 30:
		for (i=0; i<30; i++){
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19] +
				w[20][i]*input[20] + w[21][i]*input[21] + w[22][i]*input[22] + w[23][i]*input[23] + w[24][i]*input[24] +
				w[25][i]*input[25] + w[26][i]*input[26] + w[27][i]*input[27] + w[28][i]*input[28] + w[29][i]*input[29];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
	case 20:
		for (i=0; i<20; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9] +
				w[10][i]*input[10] + w[11][i]*input[11] + w[12][i]*input[12] + w[13][i]*input[13] + w[14][i]*input[14] +
				w[15][i]*input[15] + w[16][i]*input[16] + w[17][i]*input[17] + w[18][i]*input[18] + w[19][i]*input[19];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	case 10:
		for (i=0; i<10; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7] + w[8][i]*input[8] + w[9][i]*input[9];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			if (t2>0.05 || t2<-0.05)
				error+=(t2*t2);
			//else
			//	error+=0.5*(t2*t2);
		}
		break;
	case 8:
		for (i=0; i<8; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4] +
				w[5][i]*input[5] + w[6][i]*input[6] + w[7][i]*input[7];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			if (t2>0.05 || t2<-0.05)
				error+=(t2*t2);
			//else
			//	error+=0.5*(t2*t2);
		}
		break;
	case 5:
		for (i=0; i<5; i++){
			//sum=0.0;
			sum=w[0][i]*input[0] + w[1][i]*input[1] + w[2][i]*input[2] + w[3][i]*input[3] + w[4][i]*input[4];
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	default:
		for (i=0; i<n_nodes; i++){
			sum=0.0;
			for (j=0; j<n_nodes; j++){
				sum=sum+w[j][i]*input[j];
			}
			t1=1.0/(1.0+exp(-lambda[i]*sum));
			//output[i]=t1;
			t2=desired_output[i]-t1;
			error+=(t2*t2);
		}
		break;
	}

	return error;
}

inline float calculate_next_and_error_partial(int n_nodes, int cur_node, float *w, float *input, float desired_output)
{
	int i;
	float sum, t1;
	float error;

	switch (n_nodes){
	case 300:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  w[40]*input[40] + w[41]*input[41] + w[42]*input[42] + w[43]*input[43] + w[44]*input[44] +
			  w[45]*input[45] + w[46]*input[46] + w[47]*input[47] + w[48]*input[48] + w[49]*input[49] +
			  w[50]*input[50] + w[51]*input[51] + w[52]*input[52] + w[53]*input[53] + w[54]*input[54] +
			  w[55]*input[55] + w[56]*input[56] + w[57]*input[57] + w[58]*input[58] + w[59]*input[59] +
			  w[60]*input[60] + w[61]*input[61] + w[62]*input[62] + w[63]*input[63] + w[64]*input[64] +
			  w[65]*input[65] + w[66]*input[66] + w[67]*input[67] + w[68]*input[68] + w[69]*input[69] +
			  w[70]*input[70] + w[71]*input[71] + w[72]*input[72] + w[73]*input[73] + w[74]*input[74] +
			  w[75]*input[75] + w[76]*input[76] + w[77]*input[77] + w[78]*input[78] + w[79]*input[79] +
			  w[80]*input[80] + w[81]*input[81] + w[82]*input[82] + w[83]*input[83] + w[84]*input[84] +
			  w[85]*input[85] + w[86]*input[86] + w[87]*input[87] + w[88]*input[88] + w[89]*input[89] +
			  w[90]*input[90] + w[91]*input[91] + w[92]*input[92] + w[93]*input[93] + w[94]*input[94] +
			  w[95]*input[95] + w[96]*input[96] + w[97]*input[97] + w[98]*input[98] + w[99]*input[99] +
			  w[100]*input[100] + w[101]*input[101] + w[102]*input[102] + w[103]*input[103] + w[104]*input[104] +
			  w[105]*input[105] + w[106]*input[106] + w[107]*input[107] + w[108]*input[108] + w[109]*input[109] +
			  w[110]*input[110] + w[111]*input[111] + w[112]*input[112] + w[113]*input[113] + w[114]*input[114] +
			  w[115]*input[115] + w[116]*input[116] + w[117]*input[117] + w[118]*input[118] + w[119]*input[119] +
			  w[120]*input[120] + w[121]*input[121] + w[122]*input[122] + w[123]*input[123] + w[124]*input[124] +
			  w[125]*input[125] + w[126]*input[126] + w[127]*input[127] + w[128]*input[128] + w[129]*input[129] +
			  w[130]*input[130] + w[131]*input[131] + w[132]*input[132] + w[133]*input[133] + w[134]*input[134] +
			  w[135]*input[135] + w[136]*input[136] + w[137]*input[137] + w[138]*input[138] + w[139]*input[139] +
			  w[140]*input[140] + w[141]*input[141] + w[142]*input[142] + w[143]*input[143] + w[144]*input[144] +
			  w[145]*input[145] + w[146]*input[146] + w[147]*input[147] + w[148]*input[148] + w[149]*input[149] +
			  w[150]*input[150] + w[151]*input[151] + w[152]*input[152] + w[153]*input[153] + w[154]*input[154] +
			  w[155]*input[155] + w[156]*input[156] + w[157]*input[157] + w[158]*input[158] + w[159]*input[159] +
			  w[160]*input[160] + w[161]*input[161] + w[162]*input[162] + w[163]*input[163] + w[164]*input[164] +
			  w[165]*input[165] + w[166]*input[166] + w[167]*input[167] + w[168]*input[168] + w[169]*input[169] +
			  w[170]*input[170] + w[171]*input[171] + w[172]*input[172] + w[173]*input[173] + w[174]*input[174] +
			  w[175]*input[175] + w[176]*input[176] + w[177]*input[177] + w[178]*input[178] + w[179]*input[179] +
			  w[180]*input[180] + w[181]*input[181] + w[182]*input[182] + w[183]*input[183] + w[184]*input[184] +
			  w[185]*input[185] + w[186]*input[186] + w[187]*input[187] + w[188]*input[188] + w[189]*input[189] +
			  w[190]*input[190] + w[191]*input[191] + w[192]*input[192] + w[193]*input[193] + w[194]*input[194] +
			  w[195]*input[195] + w[196]*input[196] + w[197]*input[197] + w[198]*input[198] + w[199]*input[199] +
			  w[200]*input[200] + w[201]*input[201] + w[202]*input[202] + w[203]*input[203] + w[204]*input[204] +
			  w[205]*input[205] + w[206]*input[206] + w[207]*input[207] + w[208]*input[208] + w[209]*input[209] +
			  w[210]*input[210] + w[211]*input[211] + w[212]*input[212] + w[213]*input[213] + w[214]*input[214] +
			  w[215]*input[215] + w[216]*input[216] + w[217]*input[217] + w[218]*input[218] + w[219]*input[219] +
			  w[220]*input[220] + w[221]*input[221] + w[222]*input[222] + w[223]*input[223] + w[224]*input[224] +
			  w[225]*input[225] + w[226]*input[226] + w[227]*input[227] + w[228]*input[228] + w[229]*input[229] +
			  w[230]*input[230] + w[231]*input[231] + w[232]*input[232] + w[233]*input[233] + w[234]*input[234] +
			  w[235]*input[235] + w[236]*input[236] + w[237]*input[237] + w[238]*input[238] + w[239]*input[239] +
			  w[240]*input[240] + w[241]*input[241] + w[242]*input[242] + w[243]*input[243] + w[244]*input[244] +
			  w[245]*input[245] + w[246]*input[246] + w[247]*input[247] + w[248]*input[248] + w[249]*input[249] +
			  w[250]*input[250] + w[251]*input[251] + w[252]*input[252] + w[253]*input[253] + w[254]*input[254] +
			  w[255]*input[255] + w[256]*input[256] + w[257]*input[257] + w[258]*input[258] + w[259]*input[259] +
			  w[260]*input[260] + w[261]*input[261] + w[262]*input[262] + w[263]*input[263] + w[264]*input[264] +
			  w[265]*input[265] + w[266]*input[266] + w[267]*input[267] + w[268]*input[268] + w[269]*input[269] +
			  w[270]*input[270] + w[271]*input[271] + w[272]*input[272] + w[273]*input[273] + w[274]*input[274] +
			  w[275]*input[275] + w[276]*input[276] + w[277]*input[277] + w[278]*input[278] + w[279]*input[279] +
			  w[280]*input[280] + w[281]*input[281] + w[282]*input[282] + w[283]*input[283] + w[284]*input[284] +
			  w[285]*input[285] + w[286]*input[286] + w[287]*input[287] + w[288]*input[288] + w[289]*input[289] +
			  w[290]*input[290] + w[291]*input[291] + w[292]*input[292] + w[293]*input[293] + w[294]*input[294] +
			  w[295]*input[295] + w[296]*input[296] + w[297]*input[297] + w[298]*input[298] + w[299]*input[299];
		t1=1.0/(1.0+exp(-w[300]*sum));//bug 2012-05-24 w[10] --> w[100]
		error=desired_output-t1;
		break;
	case 100:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  w[40]*input[40] + w[41]*input[41] + w[42]*input[42] + w[43]*input[43] + w[44]*input[44] +
			  w[45]*input[45] + w[46]*input[46] + w[47]*input[47] + w[48]*input[48] + w[49]*input[49] +
			  w[50]*input[50] + w[51]*input[51] + w[52]*input[52] + w[53]*input[53] + w[54]*input[54] +
			  w[55]*input[55] + w[56]*input[56] + w[57]*input[57] + w[58]*input[58] + w[59]*input[59] +
			  w[60]*input[60] + w[61]*input[61] + w[62]*input[62] + w[63]*input[63] + w[64]*input[64] +
			  w[65]*input[65] + w[66]*input[66] + w[67]*input[67] + w[68]*input[68] + w[69]*input[69] +
			  w[70]*input[70] + w[71]*input[71] + w[72]*input[72] + w[73]*input[73] + w[74]*input[74] +
			  w[75]*input[75] + w[76]*input[76] + w[77]*input[77] + w[78]*input[78] + w[79]*input[79] +
			  w[80]*input[80] + w[81]*input[81] + w[82]*input[82] + w[83]*input[83] + w[84]*input[84] +
			  w[85]*input[85] + w[86]*input[86] + w[87]*input[87] + w[88]*input[88] + w[89]*input[89] +
			  w[90]*input[90] + w[91]*input[91] + w[92]*input[92] + w[93]*input[93] + w[94]*input[94] +
			  w[95]*input[95] + w[96]*input[96] + w[97]*input[97] + w[98]*input[98] + w[99]*input[99];
		t1=1.0/(1.0+exp(-w[100]*sum));//bug 2012-05-24 w[10] --> w[100]
		error=desired_output-t1;
		break;
	case 50:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39] +
			  w[40]*input[40] + w[41]*input[41] + w[42]*input[42] + w[43]*input[43] + w[44]*input[44] +
			  w[45]*input[45] + w[46]*input[46] + w[47]*input[47] + w[48]*input[48] + w[49]*input[49];
		t1=1.0/(1.0+exp(-w[50]*sum));
		error=desired_output-t1;
		break;
	case 40:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19] +
			  w[20]*input[20] + w[21]*input[21] + w[22]*input[22] + w[23]*input[23] + w[24]*input[24] +
			  w[25]*input[25] + w[26]*input[26] + w[27]*input[27] + w[28]*input[28] + w[29]*input[29] +
			  w[30]*input[30] + w[31]*input[31] + w[32]*input[32] + w[33]*input[33] + w[34]*input[34] +
			  w[35]*input[35] + w[36]*input[36] + w[37]*input[37] + w[38]*input[38] + w[39]*input[39];
		t1=1.0/(1.0+exp(-w[40]*sum));
		error=desired_output-t1;
		break;
	case 20:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9] +
			  w[10]*input[10] + w[11]*input[11] + w[12]*input[12] + w[13]*input[13] + w[14]*input[14] +
			  w[15]*input[15] + w[16]*input[16] + w[17]*input[17] + w[18]*input[18] + w[19]*input[19];
		t1=1.0/(1.0+exp(-w[20]*sum));
		error=desired_output-t1;
		break;
	case 10:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7] + w[8]*input[8] + w[9]*input[9];
		t1=1.0/(1.0+exp(-w[10]*sum));
		error=desired_output-t1;
		break;
	case 8:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4] +
			  w[5]*input[5] + w[6]*input[6] + w[7]*input[7];
		t1=1.0/(1.0+exp(-w[8]*sum));
		error=desired_output-t1;
		break;
	case 5:
		sum = w[0]*input[0] + w[1]*input[1] + w[2]*input[2] + w[3]*input[3] + w[4]*input[4];
		t1=1.0/(1.0+exp(-w[5]*sum));
		error=desired_output-t1;
		break;
	default:
		sum=0.0;
		for (i=0; i<n_nodes; i++){
			sum+=w[i]*input[i];
		}
		t1=1.0/(1.0+exp(-w[n_nodes]*sum));//bug 2012-05-24 w[100] --> w[n_nodes]
		error=desired_output-t1;
		break;
	}

	return (error*error);
}
#endif
