#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <tuple>
#include <random>
#include <algorithm>


double evaluate_sphere(float* individual);
float* generate_random(float* individual);
void swap(float ***ptr_a, float ***ptr_b);

//int partition (float** population , int low, int high,double* fitness_scores); 

double uniform_random(double range_from, double range_to);
double random_random();
double sphere_function(float* sub_individual, float* old_individual, float func_value);
void mySort(float** population, double* fitness_scores);
float** generate_pop();
double* rank_pop(float** population, double* fitness_scores, float* old_individual, double func_value);
float** crossover(float** individual, float* old_individual, float func_value);
float* mutate(float* individual, int iteration_number);
float** breed(float** individual, float* old_individual, double func_value, int iteration_number);
int roulette(double* inclusive_scan_result);
double* compute_inclusive_scan(double* fitness_scores);
float** selectFittest(double* fitness_scores, float** population, float** individual);
float** iterate_population(float** population, double* fitness_scores, float* individual, float func_value, int iteration_number);
float** internal_ga(float** individual, double func_value, float **sub_individual);





/***********************************************************/
#define npop 100
#define nVar 100000 // Number of Variables
#define target 0.000001		// Target value for convergence
#define max_value 20.0
#define min_value -10.0
#define niterations 10000	// Number of iterations
#define ncircle 5	// Number of Circles
#define max_iter 4
#define crossover_rate 0.7
#define mutation_rate 0.05

extern int var_num;
extern int start_var;
extern int end_var;
/***********************************************************/
using namespace std;

// Function to compare by the Mth element
template<int M, template<typename> class F = std::less>
struct TupleCompare
{
    template<typename T>
    bool operator()(T const &t1, T const &t2)
    {
        return F<typename tuple_element<M, T>::type>()(get<M>(t1), get<M>(t2));
    }
};

// To generate uniformly distributed variable value, in range (minimum value to maximum value)
double uniform_random(double range_from, double range_to)
{
	random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(range_from, range_to);
	return distribution(generator);
}


// To generate uniformly distributed variable value, in range (0 to 1)
double random_random()
{
	return uniform_random(0, 1);
}

// Defining custom sphere function
double sphere_function(float* sub_individual, float* old_individual, float func_value)
{
	double fitness_scores1 = 0.0;
	double fitness_scores2 = 0.0;
	for (int i = 0; i < var_num; ++i)
		fitness_scores1 += pow(sub_individual[i], 2);
	for (int i = start_var; i < end_var; ++i)
		fitness_scores2 += pow(old_individual[i], 2);

	return fitness_scores1 - fitness_scores2 + func_value;
}

// Function to rank the population according to the fitness scores
// Input
// 		- Population (vector of vector)
// Output
// 		- Sorted population according to the fitness scores, with each individual concatenated with its fitness score
// 		- vector of Tuples <individual_vector, rank> in sorted according to fitness scores
int partition (float **population,double *fitness_scores, int low, int high) 
{ 	

    double pivot = fitness_scores[high];    // pivot 
    
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (fitness_scores[j] <= pivot) 
        { 
            i++;    // increment index of smaller element
    		swap(fitness_scores[i],fitness_scores[j]);
    		swap(population[i],population[j]);
        
        } 
    } 
    swap(fitness_scores[i+1],fitness_scores[high]);
	swap(population[i+1],population[high]);
    
    return (i + 1); 
}  
/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
void quickSort(float **population,double *fitness_scores, int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(population,fitness_scores, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(population,fitness_scores, low, pi - 1); 
        quickSort(population,fitness_scores, pi + 1, high); 
    } 
}

// Function to rank the population according to the fitness scores
// Input
// 		- Population (vector of vector)
// Output
// 		- Sorted population according to the fitness scores, with each individual concatenated with its fitness score
// 		- vector of Tuples <individual_vector, rank> in sorted according to fitness scores
void mySort(float** population, double* fitness_scores)
{
	
		quickSort(population,fitness_scores,0,npop-1);

	
}



// Generate Population
// Input 
// 		- Number of variables for a particular individual in population
// 		- Number of population (To be generated)
// Returns
// 		- Population in format - vector of vector (i.e. Matrix, with each row is a variable and each column is an individual)
float** generate_pop()
{
	float** population = new float* [npop];

  	random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<float> uniform_distribution(min_value, max_value);

	// Each i generate an individual
	for(int i=0; i<npop; i++)
	{
		population[i] = new float[var_num];
		// Initialize the individual, for all variables - (0 to var_num-1)
		for(int k=0; k< var_num; k++)
			population[i][k] = uniform_distribution(generator);
		// cout << i+1 << " individual generated\n";
	}
	return population;
}

double* rank_pop(float** population, double* fitness_scores, float* old_individual, double func_value)
{ 
	
	// Initialize ranked population and concatenate with fitness score of individual
	for(int i = 0; i < npop; i++)
	// Calculate fitness score and store it in fitness_scores list
		fitness_scores[i] = sphere_function(population[i], old_individual, func_value);
	// Sort the concatenated individual vector and fitness score
	mySort(population, fitness_scores);
		
	return fitness_scores;
}

// Crossover Functions
// Input
// 		- Two individuals
// Output
// 		- Returns 2 new individuals that are crossover of given two individuals
float** crossover(float** individual, float* old_individual, float func_value)
{
    // coefficient of crossover of two individuals
	float coefficient[] = {0.5, 0.5, -1.5, 1.5};
    
    // Result of crossover of two individuals
	float** crossover_result = new float*[5];
	if(crossover_result == NULL)
		cout << "Memory OverFlow\n";
	crossover_result[0] = individual[0];
	crossover_result[1] = individual[1];

	for (int i = 2; i < 5; ++i)
	{
		crossover_result[i] = new float[var_num];
		if(crossover_result[i] == NULL)
			cout << "Memory OverFlow\n";
	}

    vector<tuple<double, int>> fitness_scores(5);

    fitness_scores[0] = make_tuple(sphere_function(individual[0], old_individual, func_value), 0);
    fitness_scores[1] = make_tuple(sphere_function(individual[1], old_individual, func_value), 1);

	// Crossover the individuals to generate new individuals
	for (int j = 0; j < 3; ++j)
	{
	    for(int i = 0; i < var_num; i++)
	    	crossover_result[j+2][i] = coefficient[j]*individual[0][i] + coefficient[j+1]*individual[1][i];

	    // Calculate the fitness scores of all the individuals in candidate
        fitness_scores[j+2] = make_tuple(sphere_function(crossover_result[j+2], old_individual, func_value), j+2);
	}

	sort(fitness_scores.begin(), fitness_scores.end(), TupleCompare<0>());

	int index1 = get<1>(fitness_scores[0]);
	int index2 = get<1>(fitness_scores[1]);

	individual[0] = crossover_result[index1];
	individual[1] = crossover_result[index2];

	for (int i = 0; i < 5; ++i)
		if(index1 != i && index2 != i)
			delete [] crossover_result[i];
		
	delete [] crossover_result;
    // Return top two crossover generated individuals
    return individual;
}

// Mutate function
// If random number between 0 and 1 is less than mutation rate then mutate the individual.
// Input
// 		- Individual (vector)
// Output : Returns
// 		- Mutated individual if random number is less than mutation rate
// 		- Same individual, otherwise
float* mutate(float* individual, int iteration_number)
{
	// If random number between 0 and 1 is less than mutation rate
    if (random_random() < mutation_rate)
    {
    	// Uniform random number between -1 and 1
        double tau = uniform_random(-1.0, 1.0);
        double exponent = 1 - pow(random_random(), 1 - iteration_number/max_iter);
        double diff = max_value - min_value;

        // calculate mutated individual
        for(int i = 0; i < var_num; i++)
            individual[i] = individual[i] + tau*diff*exponent;

        return individual;
    }
    // Return mutated individual or same individual
    return individual;
}

// Breed function
// Input
// 		- Two individuals
// Output
// 		- Returns the breed of the two given individuals
float** breed(float** individual, float* old_individual, double func_value, int iteration_number)
{
	if (random_random() < crossover_rate)
		individual = crossover(individual, old_individual, func_value);

	individual[0] = mutate(individual[0], iteration_number);
	individual[1] = mutate(individual[1], iteration_number);

	return individual;
}

// Roulette Function
// Input
// 		- Sorted list of Fitness Scores for all individuals
// Output
// 		- Index which divide fitness scores in two halves, separated by a uniform random number between 0 and total sum of fitness scores
int roulette(double* inclusive_scan_result)
{
	double random_number = 0.0;
	int index = 0;
	// Calculate uniform random between 0 and total sum of fitness scores
	random_number = uniform_random(0, inclusive_scan_result[npop-1]);
	index = lower_bound(inclusive_scan_result, inclusive_scan_result+npop, random_number)-inclusive_scan_result;

	// Return index of the separation of the list of fitness scores
	return index;
}

// Need to optimize here
// Calculate inclusive sum scan of given input array and returns a vector
double* compute_inclusive_scan(double* fitness_scores)
{
	double *inclusive_scan_result = new double [npop];
	inclusive_scan_result[0] = fitness_scores[0];
	for (int i = 1; i < npop; ++i)
		inclusive_scan_result[i] = inclusive_scan_result[i-1] + fitness_scores[i];

	return inclusive_scan_result;
}

// SelectFittest Function
// Input
// 	 	 - Population and their fitness scores
// Output
// 		 - Return two individuals whose index is selected by Roulette in the population (tuple of vector)
float** selectFittest(double* fitness_scores, float** population, float** individual)
{ 
	int i=0;

	int index1, index2;
	if(individual == NULL)
		cout << "Memory OverFlow\n";
	// Calculate two indices to select the individuals

	double* inclusive_scan_result = compute_inclusive_scan(fitness_scores);
	while(1)
	{

		index1 = roulette(inclusive_scan_result);
		index2 = roulette(inclusive_scan_result);
		


		
		if(index1 != index2)
			break;

	}
	

	if(index1 >= npop || index2 >= npop)
		cout << "OverFlow\n";
	else
	{
		individual[0] = new float[var_num];
		individual[1] = new float[var_num];
		for (int i = 0; i < var_num; ++i)
		{
			// 1st selected individual
			individual[0][i] = population[index1][i];
			//cout<<individual[0][i]<<" ";
			// 2nd selected individual
			individual[1][i] = population[index2][i];
		}
	}
	return individual;
}


// Iterate population
// Generate new generation's population
// Input
// 		- Ranked population
// Output
// 		- Returns new generation Population that has top 1/15th of the previous generation population
float** iterate_population(float** population, double* fitness_scores, float* individual, float func_value, int iteration_number)
{
	;
	float** new_population = new float*[npop];
	if(new_population == NULL)
		cout << "Memory OverFlow\n";
	// New population retains top 1/15th of the previous generation population
	for (int i = 0; i < npop/15; ++i)
		new_population[i] = population[i];
	int newpop_size = npop/15;

	// New individuals are added till new population size becomes equal to npop
	float** tuple_individual = new float*[2];
	while(newpop_size < npop)
	{
		// Take the fittest among population
		// Allocate new memory and return its address
		tuple_individual = selectFittest(fitness_scores, population, tuple_individual);
         
		// Allow breeding of the individuals
		tuple_individual = breed(tuple_individual, individual, func_value, iteration_number);
       /*  for(int i=0;i<newpop_size;i++)
         {
         	for(int j=start_var;j<end_var;j++)
         	cout<<"tuple_individual[i][j]";
         }*/
		// Push the individuals into the new generation population
		new_population[newpop_size++] = tuple_individual[0];
		if(newpop_size < npop)
			new_population[newpop_size++] = tuple_individual[1];
	}
	delete[] tuple_individual;
	for (int i = npop/15 ; i < npop; ++i)
		delete [] population[i];
	delete [] population;

	return new_population;
}

float** internal_ga(float** individual, double func_value, float **sub_individual)
{
	bool solution_found = false;
	int iteration_number = 1;
	
	// Generating npop-2 number of individuals as 2 individuals are coming from previous population
	float **population = generate_pop();

	// Extract the required chromosome from the individual1 and individual2
	// Chromosomes which are processing right now
	// Store the extracted individual in the population
	for (int i = 0; i < var_num; ++i)
	{
		population[npop-2][i] = individual[0][start_var+i];
		population[npop-1][i] = individual[1][start_var+i];
	}

	double* fitness_scores = new double [npop];
	fitness_scores = rank_pop(population, fitness_scores, individual[0], func_value);

	// Loop till solution is not found or the loop iteration exceeds the maximum iterations
	while(solution_found == false && iteration_number <= max_iter)
	{
		// Generate new generation population that retains top 1/15th of old generation population
		population = iterate_population(population, fitness_scores, individual[0], func_value, iteration_number);
		// Rank the new population
		fitness_scores = rank_pop(population, fitness_scores, individual[0], func_value);

		// calculate the score of the first individual of population
		double value = fitness_scores[0];
		if(value <= target)
			solution_found = true;
		// Increment the Loop Counter
		cout << "\tInner iteration number: " << iteration_number << " Value = "<< value << endl;
		iteration_number++;
	}

	sub_individual[0] = population[0];
	sub_individual[1] = population[1];

	for (int i = 2; i < npop; ++i)
		delete [] population[i];
	delete[] population;
	delete[] fitness_scores;

	// return two individual's chromosome value
	return sub_individual;
}

using namespace std;

int var_num = 0;
int start_var = 0;
int end_var = 0;

double evaluate_sphere(float* individual)
{
	double result = 0.0;
	for (int i = 0; i < nVar; ++i)
		result += pow(individual[i], 2);
	return result;
}

// Generate an individual, for all variables - (0 to nVar-1)
float* generate_random(float* individual)
{
	// To generate uniform random number between range_from to range_to
	random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<float> uniform_distribution(min_value, max_value);

	// Each i generate a variable value in individual
	for (int i = 0; i < nVar; ++i){
		individual[i] = uniform_distribution(generator);
	    }
	return individual;
}

void swap(float ***ptr_a, float ***ptr_b)
{
	float **temp = *ptr_a;
	*ptr_a = *ptr_b;
	*ptr_b = temp;
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	clock_t start = clock();
	
	// Number of variables per circle
	int var_circle = nVar / ncircle;
	bool solution_found = false;
	// Outer Iteration Initialization
	int outer_iteration = 1;
		
	// Initializing candidate
	float **individual = new float*[2];
	individual[0] = new float[nVar];
	individual[1] = new float[nVar];
	float** new_individual = new float*[2];
	new_individual[0] = new float[nVar];
	new_individual[1] = new float[nVar];
	float ** sub_individual = new float* [2];

	individual[0] = generate_random(individual[0]);
	for (int i = 0; i < nVar; ++i)
		individual[1][i] = individual[0][i];
	
	// Function Value
	double func_value = evaluate_sphere(individual[0]);
     cout<<func_value;
	// Beginning of chakravyuh evaluation
	while(solution_found == false && outer_iteration < niterations)
	{
		cout << "\n\n\nOuter Iteration: " << outer_iteration << endl;

		// Each circle is numbered from 0 onwards
		// Beginning of optimization within loop
		for (int circle_num = 0; circle_num < ncircle; ++circle_num)
		{
			cout << "Circle number: " << circle_num+1 << endl;
			// var_num stores locally, number of true variables that need to be solved first
			if (circle_num*var_circle + var_circle > nVar)
				var_num = circle_num*var_circle + var_circle;
			else
				var_num = var_circle;
			// Starting value of number of chromosome
			start_var = circle_num * var_circle;
			// End value of number of chromosome (excluded)
			end_var = start_var + var_num;

			// Calling internal GA function
			// Allocates memory and return the address
			sub_individual = internal_ga(individual, func_value, sub_individual);

			for (int i = start_var; i < end_var; ++i)
			{
				new_individual[0][i] = sub_individual[0][i - start_var];
				new_individual[1][i] = sub_individual[1][i - start_var];
			}
			delete[] sub_individual[0];
			delete[] sub_individual[1];
		}
		swap(&individual, &new_individual);
		func_value = evaluate_sphere(individual[0]);

		outer_iteration ++;
		cout << "Cost: " << func_value << endl;
		if(func_value <= target)
			solution_found = true;
	}

	delete[] sub_individual;
	delete[] new_individual[0];
	delete[] new_individual[1];
	delete[] new_individual;
	delete[] individual[0];
	delete[] individual[1];
	delete[] individual;

	clock_t stop = clock();
	cout << "\nTime taken: " << (double)(stop - start)/CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
