//  main.c
//  HMSEM_Simulation


/***** Libraries *****/
#include    <stdio.h>
#include    <math.h>
#include    <string.h>
#include    <gsl/gsl_math.h>
#include    <gsl/gsl_randist.h>
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_vector.h>
#include    <gsl/gsl_matrix.h>
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_linalg.h>
#include    <gsl/gsl_sf_exp.h>
#include    <gsl/gsl_sort_vector.h>

/***** Definitions *****/
#define     CLUSTER 1

/***** Function Prototype *****/
int distance(const gsl_matrix *x, int n, gsl_matrix *result);


int main(int argc, char *argv[]) {

    int N         = atoi(argv[1]);
    int ALPHA     = atoi(argv[2]);
    int LAMBDA    = atoi(argv[3]);
    double BETA_1 = atof(argv[4]);
    double BETA_2 = atof(argv[5]);
    int RAND_SEED = atoi(argv[6]);
    int P         = atoi(argv[7]);
    char ID[100];
    strcpy(ID, argv[8]);
    /****** RNG Environment ******/
    int seed = RAND_SEED;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r ,seed);



    int j;                                              //initialise counter
    gsl_matrix * x         = gsl_matrix_alloc(N, 2);    //Matrix of points
    if(CLUSTER == 1){
        //Uniform points on grid
        for(j = 0; j < N; j++){
            gsl_matrix_set(x, j, 0, gsl_rng_uniform(r));    //Set Random Coordinates of points
            gsl_matrix_set(x, j, 1, gsl_rng_uniform(r));
        }
    }
    if(CLUSTER == 0){
        //two cluster
        for(j = 0; j < N/2; j++){
            gsl_matrix_set(x, j, 0, gsl_ran_gaussian_ziggurat(r, 0.1) + 0.5);    //Set Random Coordinates of points
            gsl_matrix_set(x, j, 1, gsl_ran_gaussian_ziggurat(r, 0.1) + 0.5);
        }
         for(j = N/2; j < N; j++){
            gsl_matrix_set(x, j, 0, gsl_ran_gaussian_ziggurat(r, 0.25) + 1.5);    //Set Random Coordinates of points
            gsl_matrix_set(x, j, 1, gsl_ran_gaussian_ziggurat(r, 0.3) + 1.5);
        }
    }


    /***** Initialise Vectors and Matrices *****/
    gsl_vector * s                          = gsl_vector_alloc(2*N);    //Vector with number of susceptibles
    gsl_vector * i                          = gsl_vector_alloc(2*N);    //Vector with number of infecteds
    gsl_vector * time                       = gsl_vector_alloc(2*N);    //Vector of event times
    gsl_vector * infecteds                  = gsl_vector_alloc(N);      //Vector with 1 if individual has ever been infected
    gsl_vector * current_infecteds          = gsl_vector_alloc(N);      //Vector with 1 if individual is currently infected
    gsl_matrix * times                      = gsl_matrix_alloc(2, N);   //Matrix with infection and recovery times for each individual
    gsl_matrix * dist_mat                   = gsl_matrix_alloc(N, N);   //Distance Matrix
    gsl_vector * R                          = gsl_vector_alloc(N);      //Vector with upcoming recovery times
    gsl_vector * ones                       = gsl_vector_alloc(N);      //Vector where every element is 1.
    gsl_vector * type                       = gsl_vector_alloc(N);      //Vector of types
    gsl_vector * beta_values                = gsl_vector_alloc(P);

    FILE * f_beta = fopen("beta_values.txt", "r");
    gsl_vector_fscanf(f_beta, beta_values);
    fclose(f_beta);


    gsl_vector_set_all(R, 100000);                                      //Set upcoming recovery times to be impossibly large
    gsl_vector_set_all(ones, 1);                                        //Initialise vector of ones
    gsl_vector_set(current_infecteds, 0, 0);

    /***** Initial Conditions ****/
    gsl_vector_set(i, 0, 1);                                            //Initially 1 infected
    gsl_vector_set(s, 0, N-1);                                          //Initially N-1 susceptibles
    long i0;                                                            //Uniformly choose first infected
    i0 = gsl_rng_uniform_int(r, N);                                     //1st infected
    //printf("i0 = %ld\n", i0);
    gsl_vector_set(infecteds, i0, 1);                                   //Set values in infected vectors for i0
    gsl_vector_set(current_infecteds, i0, 1);
    gsl_matrix_set(times, 0, i0, 0.0);
    gsl_matrix_set(times, 1, i0, gsl_ran_gamma(r, ALPHA,  1.0/LAMBDA)); //Generate Recovery Time for i0
    gsl_vector_set(R, i0, gsl_matrix_get(times, 1, i0));
    distance(x, N, dist_mat);                                           //Calculate Distance between each point

    for(j = 0; j < N; j++)
      gsl_vector_set(type, j, gsl_rng_uniform_int(r, P));  //Assign Types

    for(j = 0; j < N/2.0; j++)
      gsl_vector_set(type, j, 0);  //Assign Types
   for(j = N/2.0; j < N; j++)
        gsl_vector_set(type, j, 1);  //Assign Types

    //printf("infection time = %f \t recovery time = %f\n", gsl_matrix_get(times, 0, i0), gsl_matrix_get(times, 1, i0));

    int k, m, q;
    gsl_matrix * beta_all = gsl_matrix_alloc(N, N);
    gsl_matrix * beta_all_type = gsl_matrix_alloc(N, N);
    for(k = 0; k < N; k++){
      for(m = 0; m < N; m++){
        for(q = 0; q < P; q++)
         if(gsl_vector_get(type, k) == q){
           gsl_matrix_set(beta_all_type, k, m, gsl_vector_get(beta_values, q));            //Calculate beta value
         }
       }
     }
    for(k = 0; k < N; ++k)
        for(m = 0; m < N; ++m)
          if(gsl_vector_get(type, m) == 0){
            gsl_matrix_set(beta_all, m, k, exp(-BETA_1*gsl_matrix_get(dist_mat, m, k)));            //Calculate beta value
          } else {
            gsl_matrix_set(beta_all, m, k, exp(-BETA_2*gsl_matrix_get(dist_mat, m, k)));            //Calculate beta value
          }
    gsl_matrix_mul_elements(beta_all, beta_all_type);

    /***** Generate Times *****/
    j = 0;
    int r_count1 = 0, r_count2 = 0, i_count = 1;                        //Initialise Counters
    while(gsl_vector_get(i, j) > 0){
        //If there are still infecteds....
        if(gsl_vector_get(s, j) > 0){
            //If there are still susceptibles....



            //Generate beta matrix
            gsl_matrix * beta_work = gsl_matrix_alloc(N, N);
            gsl_matrix * beta = gsl_matrix_alloc(N, N);
            gsl_matrix_memcpy(beta, beta_all);
            for(k = 0; k < N; k++){
                    for(m = 0; m < N; m++){
                        if(m == k){
                            gsl_matrix_set(beta_work, m, k, 0);
                        }else if(gsl_vector_get(infecteds, m) == 1){
                            gsl_matrix_set(beta_work, m, k, 0);
                        }else if (gsl_vector_get(current_infecteds, m) == 0 && gsl_vector_get(current_infecteds, k) == 1){
                            gsl_matrix_set(beta_work, m, k, 1);
                        } else {
                            gsl_matrix_set(beta_work, m, k, 0);
                        }
                    }
            }
            gsl_matrix_mul_elements(beta, beta_work);
            gsl_matrix_free(beta_work);

            gsl_vector * col_sums = gsl_vector_alloc(N);                                                       //Vector for column sums
            gsl_vector * col_work = gsl_vector_alloc(N);
            for(k = 0; k < N; k++){
                //Calculate sum of each column in beta matrix
                double  sum_result;
                gsl_matrix_get_col(col_work, beta, k);
                gsl_blas_ddot(ones, col_work, &sum_result);
                gsl_vector_set(col_sums, k, sum_result);

            }
            gsl_vector * prop_times                 = gsl_vector_alloc(N);                                    //Vector for proposed time for each column
            gsl_vector * infection_probabailities   = gsl_vector_alloc(N);                                    //Probability each individual is infected
            gsl_vector * cumulative_probabilites    = gsl_vector_alloc(N);
            long infector;
            int next_infected;
            for(k = 0; k < N; k++){
                if(gsl_vector_get(col_sums, k) == 0){
                    gsl_vector_set(prop_times, k, 10000);
                } else{
                    gsl_vector_set(prop_times, k, gsl_ran_exponential(r, 1.0/gsl_vector_get(col_sums, k)));
                }
            }
            //gsl_vector_fprintf(stdout, prop_times, "%f");
            infector = gsl_vector_min_index (prop_times);
            gsl_matrix_get_col(infection_probabailities, beta, infector);
            gsl_vector_set(cumulative_probabilites, 0, gsl_vector_get(infection_probabailities, 0));
            for(k = 1; k < N; k++){
                gsl_vector_set(cumulative_probabilites, k, gsl_vector_get(cumulative_probabilites, k-1) + gsl_vector_get(infection_probabailities, k));
            }
            gsl_vector_scale(cumulative_probabilites, 1.0/gsl_vector_get(cumulative_probabilites, N-1));
            //printf("Standardised Probabilites\n");
            //gsl_vector_fprintf(stdout, cumulative_probabilites, "%f");
            double u = gsl_rng_uniform(r);
            int ga_counter = 0;
            while(u > gsl_vector_get(cumulative_probabilites, ga_counter)){
                ga_counter = ga_counter + 1;

            }
            next_infected = ga_counter;
            double time_new;
            time_new = gsl_vector_min(prop_times) + gsl_vector_get(time, j);



            /***** Check if infection or recovery occurs *****/
            if(time_new < gsl_vector_min(R)){
                //An Infection Occurs...
                //printf("next infected =  %d\n", next_infected);
                gsl_vector_set(s, j+1, gsl_vector_get(s, j) - 1);
                gsl_vector_set(i, j+1, gsl_vector_get(i, j) + 1);
                gsl_matrix_set(times, 0, next_infected, time_new);
                gsl_matrix_set(times, 1, next_infected, time_new + gsl_ran_gamma(r, ALPHA, 1.0/LAMBDA));
                gsl_vector_set(R, next_infected, gsl_matrix_get(times, 1, next_infected));
                gsl_vector_set(current_infecteds, next_infected, 1);
                gsl_vector_set(infecteds, next_infected, 1);
                gsl_vector_set(time, j+1, time_new);
                i_count = i_count + 1;
                //printf("Individual = %d \t Next Infection = %f \t Recovers at = %f\n", next_infected, gsl_matrix_get(times, 0, next_infected), gsl_matrix_get(times, 1, next_infected));
            } else{
                //A Recovery Occurs...
                gsl_vector_set(s, j+1, gsl_vector_get(s, j));
                gsl_vector_set(i, j+1, gsl_vector_get(i, j) - 1);
                long next_recovery = gsl_vector_min_index(R);
                gsl_vector_set(time, j+1, gsl_vector_get(R, next_recovery));
                gsl_vector_set(current_infecteds, next_recovery, 0);
                gsl_vector_set(R, next_recovery, 10000);
                r_count1 = r_count1 + 1;

            }
        } else {
            //If there are no susceptibles...
            gsl_vector_set(s, j+1, gsl_vector_get(s, j));
            gsl_vector_set(i, j+1, gsl_vector_get(i, j) - 1);
            long next_recovery = gsl_vector_min_index(R);
            gsl_vector_set(time, j+1, gsl_vector_get(R, next_recovery));
            gsl_vector_set(current_infecteds, next_recovery, 0);
            gsl_vector_set(R, next_recovery, 1000);
            r_count2 = r_count2 + 1;

        }
        j = j + 1;
    }
    //gsl_vector_fprintf(stdout, infecteds, "%f");
    for(k = 0; k < N; k++){
        if(gsl_vector_get(infecteds, k) == 0){
            gsl_matrix_set(times, 0, k, GSL_POSINF);
            gsl_matrix_set(times, 1, k, GSL_POSINF);

        }
    }

    /***** Translate Times *****/
    gsl_vector * i_times = gsl_vector_alloc(N);
    gsl_vector * r_times = gsl_vector_alloc(N);
    gsl_matrix_get_row(i_times, times, 0);                        // Split Times Matrix into two vectors
    gsl_matrix_get_row(r_times, times, 1);
    gsl_vector_add_constant(i_times, -1*gsl_vector_min(r_times)); // Set r[0] = 0
    gsl_vector_add_constant(r_times, -1*gsl_vector_min(r_times));



    /***** Output Data *****/
    char str_in1[100]  = "";
    char str_in2[100]  = "";
    char str_in3[100]  = "";
    char str_in4[100]  = "";
    char str_in5[100]  = "";



    strcat(str_in1, "i_times");
    strcat(str_in1, ID);
    strcat(str_in1, ".dat");
    strcat(str_in2, "r_times");
    strcat(str_in2, ID);
    strcat(str_in2, ".dat");
    strcat(str_in3, "infecteds");
    strcat(str_in3, ID);
    strcat(str_in3, ".dat");
    strcat(str_in4, "positions");
    strcat(str_in4, ID);
    strcat(str_in4, ".dat");
    strcat(str_in5, "types");
    strcat(str_in5, ID);
    strcat(str_in5, ".dat");



    printf("infecteds = %d \t recovereds = %d\n", i_count, r_count1 + r_count2);
    FILE * f1 = fopen(str_in1, "w");
    FILE * f2 = fopen(str_in2, "w");
    FILE * f3 = fopen(str_in3, "w");
    FILE * f4 = fopen(str_in4, "w");
    FILE * f5 = fopen(str_in5, "w");
    gsl_vector_fprintf(f1, i_times, "%f");
    gsl_vector_fprintf(f2, r_times, "%f");
    gsl_vector_fprintf(f3, infecteds, "%f");
    gsl_matrix_fprintf(f4, x, "%f");
    gsl_vector_fprintf(f5, type, "%f");
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    /***** Free Matrices *****/
    gsl_matrix_free(x);
    gsl_vector_free(s);
    gsl_vector_free(i);
    gsl_vector_free(time);
    gsl_vector_free(R);
    gsl_vector_free(current_infecteds);
    gsl_vector_free(infecteds);
    gsl_matrix_free(times);
    gsl_matrix_free(dist_mat);



    return 0;
}

int distance(const gsl_matrix *x, int n, gsl_matrix *result){
    /* Distance Matrix where M_{j, k} is the distance between x_j and x_k */
    /*
     *  x       matrix of positions (n*2)
     *  n       length of position matrix
     *  result  distance matrix (n*n)
     */
    int j, k;
    double diff1, diff2;
    gsl_matrix * temp = gsl_matrix_alloc(n, n);
    gsl_matrix_set_all(result, 0.0);
    gsl_matrix_set_all(temp, 0.0);
    for(j = 0; j < n; j++){
        for(k = 0; k < (j+1); k++){
            diff1 = pow(gsl_matrix_get(x, j, 0) - gsl_matrix_get(x, k, 0), 2); /*Calculate difference between x coords*/
            diff2 = pow(gsl_matrix_get(x, j, 1) - gsl_matrix_get(x, k, 1), 2); /*Calculate difference between y coords*/
            gsl_matrix_set(temp, j, k, sqrt(diff1 + diff2)); /* Euclidean Distance*/
        }
    }
    gsl_matrix_memcpy(result, temp); /*Copy matrix for transposition */
    gsl_matrix_transpose(result); /*Transpose Matrix to upper triangular*/
    gsl_matrix_add(result, temp); /*Sum matrices for full distance matrix*/
    gsl_matrix_free(temp);
    return 0;

}
