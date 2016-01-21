#include "acor.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include "rand_gen.h"

using namespace std;

CAcoR::CAcoR(long seed, bool use_sparse, float sparse_penalty)
{
	this->seed=seed;
	srand(seed);
	init_random(seed);
	//solution_upperbound=1.0;
	//solution_lowerbound=-1.0;
	this->use_sparse=use_sparse;
	this->sparse_penalty=sparse_penalty;
	return;
}

CAcoR::CAcoR()
{
	n_parameters = 4;
	strcpy(method_name, "ACOR");
}

CAcoR::~CAcoR()
{
	return;
}

void CAcoR::allocate_memory()
{
	int i;
	int n_nodes=fcm->n_nodes;
	
	dimension=n_nodes+1;

	bestkant=new float*[k_best];
	for (i=0; i<k_best; i++){
		bestkant[i]=new float[dimension];
	}
	bestkantfit=new float[k_best];
	bestkantrank=new int[k_best];

	ant=new float*[pop_size];
	for (i=0; i<pop_size; i++){
		ant[i]=new float[dimension];
	}
	antfit=new float[pop_size];
	antrank=new int[pop_size];

	best_ant_weight=new float[k_best];

	solution_lowerbound=new float[dimension];
	solution_upperbound=new float[dimension];
	solution_range = new float[dimension];
	for (i=0; i<dimension; i++){
		solution_lowerbound[i]=-1.0;
		solution_upperbound[i]=1.0;
		solution_range[i] = 2.0;
	}

	composed_best_solution=new float[dimension*n_nodes];
}

void CAcoR::release_memory()
{
	int i;

	for (i=0; i<k_best; i++){
		delete[] bestkant[i];
	}
	delete[] bestkant;
	delete[] bestkantfit;
	delete[] bestkantrank;

	for (i=0; i<pop_size; i++){
		delete[] ant[i];
	}
	delete[] ant;
	delete[] antfit;
	delete[] antrank;

	delete[] best_ant_weight;
	
	delete[] solution_lowerbound;
	delete[] solution_upperbound;
	delete[] solution_range;

	delete[] composed_best_solution;
}

void CAcoR:: write_parameters(FILE *fp)
{
	fprintf(fp, "k_best = %d\n", k_best);
	fprintf(fp, "q = %f\n", q);
	fprintf(fp, "xi = %f\n", xi);
	fprintf(fp, "p0 = %f\n", p0);
}

void CAcoR::set_parameters(float para[])
{
	k_best = (int) (para[0]); //*pop_size);
	q = para[1];
	xi = para[2];
	p0 = para[3];
}

void CAcoR::init_acor()
{
	init_random(seed);
}

int CAcoR::get_iter_best_id()
{
	int i, id;
	float iter_best_fit = antfit[0];
	id = 0;
	for (i=1; i<pop_size; i++){
		if (antfit[i]<iter_best_fit){
			iter_best_fit = antfit[i];
			id = i;
		}
	}

	return id;
}

void CAcoR::init_ant_rand(int n)
{
	int i, j;

	for (i=0; i<k_best; i++){
		for (j=0; j<dimension; j++){
			// bestkant[i][j]=((double)rand()/RAND_MAX)*(solution_upperbound[j]-solution_lowerbound[j]) + solution_lowerbound[j];
			bestkant[i][j]=rand0to1()*(solution_upperbound[j]-solution_lowerbound[j]) + solution_lowerbound[j];
		}
		bestkantfit[i]=evaluate(bestkant[i], n);
	}

	sort_solutions(bestkantfit, bestkantrank, k_best);
}

void CAcoR::calculate_best_ant_weight()
{
	int i;
	float sqrt2pi=sqrt(2*PI);

	sum_best_ant_weight=0.0;
	for (i=0; i<k_best; i++){
		best_ant_weight[i]=1.0/(q*k_best*sqrt2pi)*exp(-i*i/(2.0*q*q*k_best*k_best));
		sum_best_ant_weight=sum_best_ant_weight+best_ant_weight[i];
	}
}

float CAcoR::evaluate(float *solu, int cur_node)
{
	/*
	int i;
	float sum;

	sum=0.0;
	for (i=0; i<dimension; i++){
		sum=sum+solu[i]*solu[i];
	}

	//return 1.0/(sum+0.001);
	return sum;
	*/
	
	//return fcm->fcm_fitness(solu);
	float fit=fcm->fcm_fitness_partial(solu, cur_node);
	int i;

	if (use_sparse){
		for (i=0; i<dimension-1; i++) {
			fit+= sparse_penalty * fabs(solu[i])/(dimension-1);
		}
	}

	return fit;
}

void CAcoR::sort_solutions(float *fit, int *rank, int n)
{
	int i, j;
	//float max_val;
	//int max_index, temp;
	float min_val;
	int min_index, temp;

	for (i=0; i<n; i++){
		rank[i]=i;
	}

	/*
	for (i=0; i<n; i++){
		max_val=fit[rank[i]];
		max_index=i;
		for (j=i+1; j<n; j++){
			if (max_val<fit[rank[j]]){
				max_val=fit[rank[j]];
				max_index=j;
			}
		}
		temp=rank[i];
		rank[i]=rank[max_index];
		rank[max_index]=temp;
	}
	*/
	for (i=0; i<n; i++){
		min_val=fit[rank[i]];
		min_index=i;
		for (j=i+1; j<n; j++){
			if (min_val>fit[rank[j]]){
				min_val=fit[rank[j]];
				min_index=j;
			}
		}
		temp=rank[i];
		rank[i]=rank[min_index];
		rank[min_index]=temp;
	}
}

void CAcoR::update_archive(/*float *bestkantfit, float *antfit, int *bestkantrank, int *antrank, int n1, int n2*/)
{
	int i, j;
	int reverse_index;
	int ant_index;

	//k_best=n1;
	//num_ants=n2;

	ant_index=0;
	reverse_index=k_best-1;
	for (i=0; i<=reverse_index && ant_index<pop_size; ){//2012-05-05 bug fixed
		//if (bestkantfit[bestkantrank[i]]<antfit[antrank[ant_index]]){
		if (bestkantfit[bestkantrank[i]]>antfit[antrank[ant_index]]){
			bestkantfit[bestkantrank[reverse_index]]=antfit[antrank[ant_index]];
			for (j=0; j<dimension; j++){
				bestkant[bestkantrank[reverse_index]][j]=ant[antrank[ant_index]][j];
			}
			reverse_index--;
			ant_index++;
		}
		else{
			i++;
		}
	}

	sort_solutions(bestkantfit, bestkantrank, k_best);
}

void CAcoR::build_tour(int cur_ant)
{
	int i, j, selected;
	float temp, sum, sigma;
	float temp_solu;

	double debug_t1;
	double debug_t2;

	p0=0.0; /////////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	for (i=0; i<dimension; i++){
		
		//select an ant (gaussian distribution)
		// if ((double)rand()/RAND_MAX<p0){
		if (rand0to1()<p0){
			//select the best one

		}
		else {
			//select randomly according to the weights
			// temp=(double)rand()/RAND_MAX*sum_best_ant_weight;
			temp=rand0to1()*sum_best_ant_weight;
			sum=0.0;
			for (j=0; j<k_best-1; j++){	// if the (k_best-1)-th one did not work, don't need to try the k_best-th choice. In this way, j will always less than k_best
				sum=sum+best_ant_weight[j];
				if (sum>temp)
					break;
			}
			selected=bestkantrank[j];

			//calculate sigma for the selected ant
			sigma=0.0;
			for (j=0; j<k_best; j++){
				sigma=sigma+fabs(bestkant[j][i]-bestkant[selected][i]);
			}
			sigma=xi*sigma/(float)(k_best-1);

			//sampling the selected gaussian
			//sqrt(-2.0*log((double)rand()))*cos(2*PI*rand());
			do{
				do{
					// debug_t1=(double)rand()/RAND_MAX;
					debug_t1=rand0to1();
				}while (debug_t1==0.0);
				// debug_t2=(double)rand()/RAND_MAX;
				debug_t2=rand0to1();
				temp_solu=sigma*sqrt(-2.0*log(debug_t1))*cos(2*PI*debug_t2)+bestkant[selected][i];
				//temp_solu=sigma*sqrt(-2.0*log((double)rand()/RAND_MAX))*cos(2*PI*(double)rand()/RAND_MAX)+bestkant[selected][i];
			}while(temp_solu < solution_lowerbound[i] || temp_solu > solution_upperbound[i]);
			
			ant[cur_ant][i]=temp_solu;
		}
	}
}

void CAcoR::compose_bestant(int n)
{
	int i;
	int n_nodes=fcm->n_nodes;
	for (i=0; i<n_nodes; i++) {
		composed_best_solution[i*n_nodes+n]=bestkant[bestkantrank[0]][i];
	}
	composed_best_solution[n_nodes*n_nodes+n]=bestkant[bestkantrank[0]][n_nodes];

	composed_bestfit+=bestkantfit[bestkantrank[0]];
}

void CAcoR::run()
{
	int cur_iter, cur_ant;

	calculate_best_ant_weight();
	init_ant_rand(0);

	for (cur_iter=0; cur_iter<max_iter; cur_iter++){
		for (cur_ant=0; cur_ant<pop_size; cur_ant++){
			build_tour(cur_ant);
			antfit[cur_ant]=evaluate(ant[cur_ant], 0);
		}
		sort_solutions(antfit, antrank, pop_size);
		update_archive();

		if (cur_iter%500==0)
			cout<<"cur_iter="<<cur_iter<<" best fit="<<bestkantfit[bestkantrank[0]]<<endl;
	}

}

void CAcoR::run_decomposed()
{
	int cur_iter, cur_ant, cur_node;

	init_acor();
	calculate_best_ant_weight();
	composed_bestfit=0.0;

	for (cur_node=0; cur_node<fcm->n_nodes; cur_node++) {
		init_ant_rand(cur_node);
		
		for (cur_iter=0; cur_iter<max_iter; cur_iter++) {
			for (cur_ant=0; cur_ant<pop_size; cur_ant++) {
				build_tour(cur_ant);
				antfit[cur_ant]=evaluate(ant[cur_ant], cur_node);
			}
			sort_solutions(antfit, antrank, pop_size);
			update_archive();

			if (save_intermediate_results==true){
				save_iter_best_solution(cur_iter, cur_node, ant[get_iter_best_id()]);
			}

			if (cur_iter%500==0)
				cout<<"cur_node="<<cur_node<<" cur_iter="<<cur_iter<<" best fit="<<bestkantfit[bestkantrank[0]]<<endl;
		}
		compose_bestant(cur_node);
	}

}
