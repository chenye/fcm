#ifndef _MAIN_OPT
#define _MAIN_OPT

void try_one_parameter(int method, char seq_name[], float sparse_penalty = 0.0,float lambda_min=5.0, float lambda_max=5.0, int normalization_method = 0, float lambda_in_normalization = 2.0, bool save_intermediate_results = false);
void try_diff_parameters(int method, char seq_name[]);
void try_diff_parameters(int method, char seq_name[], float sparse_penalty, float lambda_min, float lambda_max, int normalization_method);

#endif
