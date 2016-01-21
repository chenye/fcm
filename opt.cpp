
#include "opt.h"
#include "fcm.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void COpt::write_result(char filename[])
{
	FILE *fp=fopen(filename, "a+");

	fprintf(fp, "%% NEW RESULTS %%\n");
	
#ifdef WIN32
	char str_date[10];
	char str_time[10];
	fprintf(fp, "Current_Time = '%s, %s'\n\n", _strdate(str_date), _strtime(str_time));
#else
	time_t rawtime;
  	struct tm * timeinfo;
  	char buffer [80];
  	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
  	strftime (buffer,80,"%c",timeinfo);
	fprintf(fp, "Current_Time = '%s'\n\n", buffer);
#endif
	
	fprintf(fp, "seed = %ld;\n", seed);
	fprintf(fp, "fcm_num_node = %d;\n", fcm->n_nodes);
	fprintf(fp, "map_density = %f;\n", fcm->map_density);
	fprintf(fp, "len_data = %d;\n", fcm->len_data);
	fprintf(fp, "num_seq = %d;\n", fcm->num_seq);
	fprintf(fp, "lambda = %f;\n", fcm->lambda[0]);

	fprintf(fp, "dimension = %d;\n", dimension);
	fprintf(fp, "max_iter = %d;\n", max_iter);
	fprintf(fp, "lambda_lowerbound = %f;\n", solution_lowerbound[fcm->n_nodes]);
	fprintf(fp, "lambda_upperbound = %f;\n", solution_upperbound[fcm->n_nodes]);
	fprintf(fp, "sparse_penalty = %f;\n", sparse_penalty);
	// algorithm related parameters:
	write_parameters(fp);
	fprintf(fp, "pop_size = %d;\n", pop_size);
	// fprintf(fp, "omega = %f;\n", pso->omega);
	// fprintf(fp, "phi1 = %f;\n", pso->phi1);
	// fprintf(fp, "phi2 = %f;\n", pso->phi2);
	
	fprintf(fp, "result_data_error = %f;\n", composed_bestfit);
	fprintf(fp, "result_model_error= %f;\n", fcm->get_model_error(composed_best_solution));
	
	fcm->get_undirected_error(composed_best_solution, 0.05);
	
	///////////////////////////////////////////////////////////////////////////
	fprintf(fp, "true_pos=%d;\n", fcm->true_pos);
	fprintf(fp, "true_neg=%d;\n", fcm->true_neg);
	fprintf(fp, "false_pos=%d;\n", fcm->false_pos);
	fprintf(fp, "false_neg=%d;\n", fcm->false_neg);
	fprintf(fp, "true_pos_s=%d;\n", fcm->true_pos_s);
	fprintf(fp, "true_neg_s=%d;\n", fcm->true_neg_s);
	fprintf(fp, "false_pos_s=%d;\n", fcm->false_pos_s);
	fprintf(fp, "false_neg_s=%d;\n", fcm->false_neg_s);
	//float spec=(float)acs->true_pos_s/(acs->true_pos_s+acs->false_neg_s);
	//float sens=(float)acs->true_neg_s/(acs->true_neg_s+acs->false_pos_s);
	fprintf(fp, "specificity=%f;\n", fcm->spec);
	fprintf(fp, "sensitivity=%f;\n", fcm->sens);
	fprintf(fp, "SS_mean=%f;\n", fcm->ss_mean);

	fprintf(fp, "accuracy = %f;\n", (float)(fcm->true_neg+fcm->true_pos)/(fcm->true_neg+fcm->true_pos+fcm->false_neg+fcm->false_pos));

	float auroc, aupr;

	fcm->get_auc(composed_best_solution, &auroc, &aupr);
	fprintf(fp, "auroc = %f;\n", auroc);
	fprintf(fp, "aupr = %f;\n", aupr);

	fprintf(fp, "%%%lf\t%lf\t%f\t%f\t%f\t%f\t%f\n", composed_bestfit,
											fcm->get_model_error(composed_best_solution),
											fcm->spec, fcm->sens, fcm->ss_mean,
											auroc, aupr );

	fprintf(fp, "weight=[\n");
	int ind=0;
	for (int i=0; i<fcm->n_nodes; i++){
		for (int j=0; j<fcm->n_nodes; j++){
			fprintf(fp, "%f\t", composed_best_solution[ind++]);
		}
		fprintf(fp, ";\n");
	}
	fprintf(fp, "];\n");

	//write lambda
	fprintf(fp, "lambda=[\n");
	for (int i=0; i<fcm->n_nodes; i++){
		fprintf(fp, "%f\t", composed_best_solution[ind++]);
	}
	fprintf(fp, "];\n");

	//target_weight
	fprintf(fp, "target_weight=[\n");
	for (int i=0; i<fcm->n_nodes; i++){
		for (int j=0; j<fcm->n_nodes; j++){
			fprintf(fp, "%f\t", fcm->target_weight[i][j]);
		}
		fprintf(fp, ";\n");
	}
	fprintf(fp, "];\n");
	
	fclose(fp);
}

void COpt::write_simple_result(char filename[], bool write_new_line_only)
{
	float auroc, aupr;

	FILE *fp=fopen(filename, "a+");
	
	if (write_new_line_only==true)
		fprintf(fp, "\n");
	else
	{
		//acor->fcm->get_undirected_error(acor->bestkant[acor->bestkantrank[0]], 0.05);
		fcm->get_undirected_error(composed_best_solution, 0.05);
		fcm->get_auc(composed_best_solution, &auroc, &aupr);

		fprintf(fp, "%lf\t%lf\t%f\t%f\t%f\t%f\t%f\n", composed_bestfit,
											fcm->get_model_error(composed_best_solution),
											fcm->spec, fcm->sens, fcm->ss_mean,
											auroc,aupr );
	}

	fclose(fp);
}

void COpt::init_intermediate_result()
{
	intermediate_fit = new float[max_iter];
	intermediate_model_error = new float[max_iter];

	int i;
	for (i=0; i<max_iter; i++){
		intermediate_fit[i]=0.0;
		intermediate_model_error[i]=0.0;
	}
}

void COpt::free_intermediate_result()
{
	delete[] intermediate_fit;
	delete[] intermediate_model_error;
}

void COpt::save_iter_best_solution(int cur_iter, int cur_node, float *solu)
{
	float temp_fit, temp_model_error;
	temp_fit = fcm->fcm_fitness_partial(solu, cur_node);

	int i;
	if (use_sparse){
		for (i=0; i<dimension-1; i++){
			temp_fit += sparse_penalty * fabs(solu[i]/(dimension-1));
		}
	}

	temp_model_error=0.0;
	for (i=0; i<fcm->n_nodes; i++){
		temp_model_error += fabs(fcm->target_weight[i][cur_node]-solu[i]);
	}
	temp_model_error = temp_model_error/(fcm->n_nodes);

	intermediate_fit[cur_iter] += temp_fit;
	intermediate_model_error[cur_iter] += temp_model_error;
}

void COpt::write_intermediate_result_to_file(char filename[])
{
	int i;

	FILE *fp=fopen(filename, "a+");
	for (i=0; i<max_iter; i++){
		fprintf(fp, "%f\t%f\n", intermediate_fit[i], intermediate_model_error[i]/(fcm->n_nodes));
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}
