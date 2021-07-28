/******************************************************
 *  GP Inference for Heterogeneously Mixing Epidemic  *
 *         with Projected Process Aproximation        *
 *                                                    *
 ******************************************************/




/***** Libraries *****/
#include    <stdio.h>
#include    <stdlib.h>
#include    <string.h>
#include    <math.h>
#include    <time.h>

/***** GSL Libraries *****/
#include    <gsl/gsl_math.h>
#include    <gsl/gsl_randist.h>
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_vector.h>
#include    <gsl/gsl_matrix.h>
#include    <gsl/gsl_blas.h>
#include    <gsl/gsl_linalg.h>
#include    <gsl/gsl_sf_exp.h>
#include    <gsl/gsl_sort_vector.h>

#define     BURN 1000
#define     n_update 70


/***** Function Prototypes *****/
double  double_sum(const gsl_vector * i, const gsl_vector * r, const gsl_vector * types, const gsl_vector * infecteds, const gsl_matrix * beta, gsl_matrix * times, int n, int N1, int N2, int t);
double  likelihood(const gsl_vector * i, const gsl_vector * r, const gsl_vector * types, const gsl_vector * infecteds, const gsl_matrix * beta, double d, double sum, double log_sum, double alpha, double gamma, int N, int n, int t);
double  chol_det(const gsl_matrix * A, int N);
double  log_normal_pdf(const gsl_vector * x, const gsl_matrix * A_inverse, const double log_A_det, int N);
void    beta_vector_to_matrix(const gsl_vector * v, gsl_matrix * m, const gsl_matrix * index, const gsl_vector * types, const int N1, const int N2, const int t);
void    distance(const gsl_matrix *x, int n, gsl_matrix *result);
void    sq_exp(const gsl_vector *x, const gsl_vector *y, double alph, double ell, const int n1, const int n2, gsl_matrix * result);
void    rmvnorm_chol(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);

int         main        (int argc, char *argv[]){

    printf("Reading Inputs....\n");
    int N           = atoi(argv[1]);
    int ALPHA       = atoi(argv[2]);
    int M           = atoi(argv[3]);
    int RAND_SEED   = atoi(argv[4]);
    double INF      = atof(argv[5]);
    double VAR      = atof(argv[6]);
    double ELL      = atof(argv[7]);
    double DELTA    = atof(argv[8]);
    int S           = atoi(argv[9]);
    int P           = atoi(argv[10]);
    char ID[100];
    strcpy(ID, argv[11]);

    /****** RNG Environment ******/
    int seed = RAND_SEED;
    const gsl_rng_type * T;
    gsl_rng * random;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    random = gsl_rng_alloc (T);
    gsl_rng_set (random, seed);

    int     j, k;                            //set for loop counters
    int     beta_count[2],  i_count = 0;    //set MH counters
    beta_count[0] = 0; beta_count[1] = 0;
    double  diff_t;                          //start timer
    time_t  start, stop;
    time(&start);

 /*****Import Simulation *****/
    printf("Importing Data...\n");
    gsl_matrix * x                  = gsl_matrix_alloc(N, 2);   //Matrix of points
    gsl_vector * infecteds          = gsl_vector_alloc(N);      //Matrix with 1 if individual was infected, 0 otherwise
    gsl_vector * i                  = gsl_vector_alloc(N);      //Vector of infection times
    gsl_vector * r                  = gsl_vector_alloc(N);      //Vector of recovery times
    gsl_vector * types_full         = gsl_vector_alloc(N);      //Vector of recovery times


    char str_in2[100]  = "";
    char str_in3[100]  = "";
    char str_in4[100]  = "";
    char str_in5[100]  = "";

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

    FILE * f1 = fopen(str_in4, "r");
    FILE * f3 = fopen(str_in2, "r");
    FILE * f4 = fopen(str_in3, "r");
    FILE * f5 = fopen(str_in5, "r");                                  //matrix of point coordinates
    gsl_matrix_fscanf(f1, x);         //matrix of point coordinates
    gsl_vector_fscanf(f3, r);         //vector of recovery times
    gsl_vector_fscanf(f4, infecteds); //vector of binary values, 0 => not infected, 1 => Infected
    gsl_vector_fscanf(f5, types_full);
    fclose(f1);
    fclose(f3);
    fclose(f4);
    fclose(f5);
    gsl_vector_set_all(types_full, 1.0);
    for(j = 0; j < (N/2.0); j++)
      gsl_vector_set(types_full, j, 0.0);


    if(INF == 0){
        for(j = 0; j < N; j++)
            if(gsl_vector_get(infecteds, j) == 1){
                gsl_vector_set(i, j, gsl_vector_get(r, j) - 4.0 + gsl_ran_ugaussian(random));   //Initial infection times are r - 3
            } else {
                gsl_vector_set(i, j, GSL_POSINF);
           }
    }

    if(INF == 1){
        char str_in1[100]  = "";
        strcat(str_in1, "i_times");
        strcat(str_in1, ID);
        strcat(str_in1, ".dat");
        FILE * f2 = fopen(str_in1, "r");
        gsl_vector_fscanf(f2, i);
        fclose(f2);
    }

    int n = 0;
    for(j = 0; j < N; j++){
        n = gsl_vector_get(infecteds, j) + n;
    }

/***** Sort Types ******/
    printf("Sorting Types...\n");
    gsl_vector * types = gsl_vector_alloc(P);
    for(j = 0; j < P; j++)
      gsl_vector_set(types, j, j);
    gsl_matrix * type_mat_index_work = gsl_matrix_alloc(N, N);
    gsl_matrix * type_mat_index = gsl_matrix_alloc(N, N);
    gsl_matrix_set_all(type_mat_index, 0.0);
    k = 1;
    double m, index;
    for(j = 0; j < N; j++){
        for(k = 0; k < N; k ++){
            double w = gsl_vector_get(types_full, j);
            for(m = 0 ; m < P; m++){
                if(gsl_vector_get(types, m) == w){
                    index = m;
                }
            }
            gsl_matrix_set(type_mat_index, j, k, index);
        }
    }

    int n_type[P];
    for(k = 0; k < P; k++)
      n_type[k] = 0;

    for(k = 0; k < P; k++){
      for(j = 0; j < N; j++){
        if(gsl_vector_get(types_full, j) == k & gsl_finite(gsl_vector_get(i, j)) == 1)
          n_type[k]++;
      }
    }


    int N_type[P];
    for(k = 0; k < P; k++)
      N_type[k] = 0;

    for(k = 0; k < P; k++){
      for(j = 0; j < N; j++){
        if(gsl_vector_get(types_full, j) == k)
          N_type[k]++;
      }
    }


/***** Calculate Distance ******/
    printf("Sorting Distances...\n");
    gsl_matrix * dist_mat           = gsl_matrix_alloc(N, N);       //Matrix of distances
    gsl_vector * dist_vec           = gsl_vector_alloc(N*N);        //Vector of distances
    gsl_vector * dist_vec_work      = gsl_vector_alloc(N*N);        //Vector of distances
    gsl_vector * dist_vec_sorted    = gsl_vector_alloc(N*(N-1)/2);  //Vector of distances
    gsl_matrix * dist_mat_index     = gsl_matrix_alloc(N, N);       //Matrix of distance indices
    gsl_matrix * dist_mat_index_work= gsl_matrix_alloc(N, N);       //Matrix of distance indices
    distance(x, N, dist_mat);

    for(j = 0; j < N; j++){                                         //Turn Distance Matrix into Vector
        for(k = 0; k < j; k++){
            gsl_vector_set(dist_vec, j*N + k, gsl_matrix_get(dist_mat, j, k));
        }
    }
    gsl_vector_memcpy(dist_vec_work, dist_vec);
    gsl_sort_vector(dist_vec_work);                                  //sort distances

    for(j = 0; j < N*(N-1)/2; j++)
        gsl_vector_set(dist_vec_sorted, j, gsl_vector_get(dist_vec_work, N*N - (N*(N-1))/2 + j)); //remove 0s and duplicates


/***** Import Distances Per Type *****/
    char str_dist0[100]  = "";
    strcat(str_dist0, "distances0_");
    strcat(str_dist0, ID);
    strcat(str_dist0, ".txt");
    char str_dist1[100]  = "";
    strcat(str_dist1, "distances1_");
    strcat(str_dist1, ID);
    strcat(str_dist1, ".txt");
    char str_dist3[100]  = "";
    strcat(str_dist3, "distances0_sparse");
    strcat(str_dist3, ID);
    strcat(str_dist3, ".txt");
    char str_dist4[100]  = "";
    strcat(str_dist4, "distances1_sparse");
    strcat(str_dist4, ID);
    strcat(str_dist4, ".txt");
    char str_dist_lengths[100]  = "";
    strcat(str_dist_lengths, "lengths_");
    strcat(str_dist_lengths, ID);
    strcat(str_dist_lengths, ".txt");
    char str_index_mat_0[100]  = "";
    strcat(str_index_mat_0, "index_mat0_");
    strcat(str_index_mat_0, ID);
    strcat(str_index_mat_0, ".txt");
    char str_index_mat_1[100]  = "";
    strcat(str_index_mat_1, "index_mat1_");
    strcat(str_index_mat_1, ID);
    strcat(str_index_mat_1, ".txt");


    int unique[2];
    FILE * unique_lengths = fopen(str_dist_lengths, "r");
    fscanf(unique_lengths, "%d \n%d", &unique[0], &unique[1]);
    gsl_vector * dist_0_sorted = gsl_vector_alloc(unique[0]);
    gsl_vector * dist_1_sorted = gsl_vector_alloc(unique[1]);
    gsl_vector * dist_sparse_0 = gsl_vector_alloc(S);
    gsl_vector * dist_sparse_1 = gsl_vector_alloc(S);
    gsl_matrix * dist_mat_index_0 = gsl_matrix_alloc(N, N_type[0]);
    gsl_matrix * dist_mat_index_1 = gsl_matrix_alloc(N, N_type[1]);
    FILE * f_d0 = fopen(str_dist0, "r");
    FILE * f_d1 = fopen(str_dist1, "r");
    FILE * f_d0_sparse = fopen(str_dist3, "r");
    FILE * f_d1_sparse = fopen(str_dist4, "r");
    FILE * f_d_mat0 = fopen(str_index_mat_0, "r");
    FILE * f_d_mat1 = fopen(str_index_mat_1, "r");

    gsl_vector_fscanf(f_d0, dist_0_sorted);
    gsl_vector_fscanf(f_d1, dist_1_sorted);
    gsl_vector_fscanf(f_d0_sparse, dist_sparse_0);
    gsl_vector_fscanf(f_d1_sparse, dist_sparse_1);
    gsl_matrix_fscanf(f_d_mat0, dist_mat_index_0);
    gsl_matrix_fscanf(f_d_mat1, dist_mat_index_1);

    fclose(f_d0);
    fclose(f_d1);
    fclose(f_d0_sparse);
    fclose(f_d1_sparse);
    fclose(unique_lengths);
    fclose(f_d_mat0);
    fclose(f_d_mat1);
    //gsl_matrix_transpose(dist_mat_index_0);
    //gsl_matrix_transpose(dist_mat_index_1);


/***** Covariance Matrix *****/
    printf("Constructing Covariance Matrix...\n");
    gsl_matrix * sigma_0            = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1            = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_0_prop       = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1_prop       = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_0_chol       = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1_chol       = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_0_chol_prop  = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1_chol_prop  = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_0_inverse    = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1_inverse    = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_0_inverse_prop    = gsl_matrix_alloc(S, S);
    gsl_matrix * sigma_1_inverse_prop    = gsl_matrix_alloc(S, S);
    gsl_matrix * id                      = gsl_matrix_alloc(S, S);
    gsl_matrix * k_s_d_0                 = gsl_matrix_alloc (unique[0], S);
    gsl_matrix * k_s_d_1                 = gsl_matrix_alloc (unique[1], S);


    sq_exp(dist_sparse_0, dist_sparse_0, VAR, ELL, S, S, sigma_0);    //Compute Square Exponential Covariance Function
    sq_exp(dist_0_sorted, dist_sparse_0, VAR, ELL, unique[0], S, k_s_d_0);
    sq_exp(dist_sparse_1, dist_sparse_1, VAR, ELL, S, S, sigma_1);    //Compute Square Exponential Covariance Function
    sq_exp(dist_1_sorted, dist_sparse_1, VAR, ELL, unique[1], S, k_s_d_1);




    gsl_matrix_set_identity(id);
    gsl_matrix_scale(id, 0.000001);
    gsl_matrix_add(sigma_0, id);                                                          //Add correction term
    gsl_matrix_add(sigma_1, id);
    gsl_matrix_memcpy(sigma_0_chol, sigma_0);
    gsl_linalg_cholesky_decomp(sigma_0_chol);                                            //Compyte Cholesky Decomposition
    gsl_matrix_memcpy(sigma_0_inverse, sigma_0_chol);
    gsl_linalg_cholesky_invert(sigma_0_inverse);
    gsl_matrix_memcpy(sigma_1_chol, sigma_1);
    gsl_linalg_cholesky_decomp(sigma_1_chol);                                            //Compyte Cholesky Decomposition
    gsl_matrix_memcpy(sigma_1_inverse, sigma_1_chol);
    gsl_linalg_cholesky_invert(sigma_1_inverse);
/***** Initialise Matrices and Vectors *****/
    printf("Preparing for MCMC...\n");
    gsl_vector * beta_vector        = gsl_vector_alloc(N*(N-1)/2);
    gsl_vector * beta_vector_0      = gsl_vector_alloc(unique[0]);
    gsl_vector * beta_vector_1      = gsl_vector_alloc(unique[1]);
    gsl_vector * beta_vector_0_prop = gsl_vector_alloc(unique[0]);
    gsl_vector * beta_vector_1_prop = gsl_vector_alloc(unique[1]);
    gsl_matrix * beta_matrix_0      = gsl_matrix_alloc(N, N_type[0]);
    gsl_matrix * beta_matrix_0_prop = gsl_matrix_alloc(N, N_type[0]);
    gsl_matrix * beta_matrix_1      = gsl_matrix_alloc(N, N_type[1]);
    gsl_matrix * beta_matrix_1_prop = gsl_matrix_alloc(N, N_type[1]);
    gsl_vector * f                  = gsl_vector_alloc(N*(N-1)/2);                //Store current values of PPA GP
    gsl_vector * f_0                = gsl_vector_alloc(unique[0]);                //Store current values of PPA GP
    gsl_vector * f_1                = gsl_vector_alloc(unique[1]);                //Store current values of PPA GP
    gsl_vector * f_prop             = gsl_vector_alloc(N*(N-1)/2);                //Store current values of PPA GP
    gsl_vector * f_0_prop           = gsl_vector_alloc(unique[0]);                //Store current values of PPA GP
    gsl_vector * f_1_prop           = gsl_vector_alloc(unique[1]);                //Store current values of PPA GP
    gsl_vector * f_bar              = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_work             = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_0_bar            = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_1_bar            = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_bar_prop         = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_0_bar_prop       = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_1_bar_prop       = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_add              = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_add_prop         = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * f_temp             = gsl_vector_alloc(S);                //Store current values of PPA GP
    gsl_vector * mu1                = gsl_vector_calloc(S);
    gsl_vector * diff_vec           = gsl_vector_alloc(S);
    gsl_vector * diff_vec_prop      = gsl_vector_alloc(S);
    gsl_vector * temp_vec           = gsl_vector_alloc(S);
    gsl_vector * mu2                = gsl_vector_alloc(N*(N-1)/2);
    gsl_vector * mu3                = gsl_vector_alloc(N*(N-1)/2);
    gsl_vector * i_props            = gsl_vector_alloc(N);
    gsl_vector * beta_vector_prop   = gsl_vector_alloc(N*(N-1)/2);
    gsl_vector * prior_draw  	      = gsl_vector_alloc(S);
    gsl_matrix * d_sum_times_0      = gsl_matrix_calloc(N, N_type[0]);
    gsl_matrix * d_sum_times_1      = gsl_matrix_calloc(N, N_type[1]);
    gsl_matrix * d_sum_times_0_prop = gsl_matrix_calloc(N, N_type[0]);
    gsl_matrix * d_sum_times_1_prop = gsl_matrix_calloc(N, N_type[1]);
    gsl_vector * mean_0             = gsl_vector_calloc(unique[0]);
    gsl_vector * mean_1             = gsl_vector_calloc(unique[1]);
    gsl_vector * mean_0_work        = gsl_vector_calloc(unique[0]);
    gsl_vector * mean_1_work        = gsl_vector_calloc(unique[1]);
    gsl_vector * var_0              = gsl_vector_calloc(unique[0]);
    gsl_vector * var_1              = gsl_vector_calloc(unique[1]);
    gsl_vector * var_0_work         = gsl_vector_calloc(unique[0]);
    gsl_vector * var_1_work         = gsl_vector_calloc(unique[1]);

    // lengthscale  set up
    double ell_current[2], ell_prop[2], det[2], det_prop[2];
    ell_current[0] = ELL;
    ell_current[1] = ELL;
    det[0] = chol_det(sigma_0_chol, S);
    det[1] = chol_det(sigma_1_chol, S);


    double log_prior_ratio_ell, log_p_acc, log_prior_ratio_f, loglike_ratio, log_prior_ratio_v, v_dot, v_prop_dot, ell_count;
    gsl_vector_set_all(mu3, 0.0);
    ell_count = 0;

    double gamma_current, i_sum, d_sum[2], d_sum_prop[2], s_sum[2], s_sum_prop[2], l_sum[2], l_sum_prop[2], q_prob, loglike_prop[2], loglike[2], i_prop, p_acc, bet_sum;
    double delta1 = sqrt(1-(DELTA*DELTA));
    int i_type;

    gamma_current = 1.0;
    s_sum[0] = 0.0;
    s_sum[1] = 0.0;
    l_sum[0] = 0.0;
    l_sum[1] = 0.0;

    for(k = 0; k < N; k++){
       if(gsl_vector_get(infecteds, k)==1){
         if(gsl_vector_get(types_full, k) == 0){
            s_sum[0] = s_sum[0] + gsl_vector_get(r, k)-gsl_vector_get(i, k);             //compute Sum(r-i)
            l_sum[0] = l_sum[0] + log(gsl_vector_get(r, k)-gsl_vector_get(i, k));             //compute Sum(r-i)
          } else {
             s_sum[1] = s_sum[1] + gsl_vector_get(r, k)-gsl_vector_get(i, k);             //compute Sum(r-i)
             l_sum[1] = l_sum[1] + log(gsl_vector_get(r, k)-gsl_vector_get(i, k));             //compute Sum(r-i)
           }
        }
      }
    gsl_vector_set_all(f_0_bar, log(0.001));
    gsl_vector_set_all(f_add, log(0.001));
    gsl_vector_memcpy(f_1_bar, f_add);
    gsl_vector_add(f_1_bar, f_0_bar);

    gsl_blas_dgemv (CblasNoTrans, 1.0, sigma_0_inverse, f_0_bar, 0.0, f_work); //Project f_0_bar onto full datatset
    gsl_blas_dgemv (CblasNoTrans, 1.0, k_s_d_0, f_work, 0.0, f_0);
    for(j = 0; j < unique[0]; j++)
        gsl_vector_set(beta_vector_0, j, exp(gsl_vector_get(f_0, j)));

    gsl_blas_dgemv (CblasNoTrans, 1.0, sigma_1_inverse, f_1_bar, 0.0, f_work); //Project f_1_bar onto full datatset
    gsl_blas_dgemv (CblasNoTrans, 1.0, k_s_d_1, f_work, 0.0, f_1);
    for(j = 0; j < unique[1]; j++)
        gsl_vector_set(beta_vector_1, j, exp(gsl_vector_get(f_1, j)));

    beta_vector_to_matrix(beta_vector_0, beta_matrix_0, dist_mat_index_0, types_full, N_type[0], N, 0);
    beta_vector_to_matrix(beta_vector_1, beta_matrix_1, dist_mat_index_1, types_full, N_type[1], N, 1);

    d_sum[0]    = double_sum(i, r, types_full, infecteds, beta_matrix_0, d_sum_times_0, n_type[0], N_type[0], N, 0);                        //Compute initial Double Sum
    d_sum[1]    = double_sum(i, r, types_full, infecteds, beta_matrix_1, d_sum_times_1, n_type[1], N_type[1], N, 1);                        //Compute initial Double Sum
    loglike[0]  = likelihood(i, r, types_full, infecteds, beta_matrix_0, d_sum[0], s_sum[0], l_sum[0], ALPHA, 1.0, N, n_type[0], 0);//Compute intial likelihood
    loglike[1]  = likelihood(i, r, types_full, infecteds, beta_matrix_1, d_sum[1], s_sum[1], l_sum[1], ALPHA, 1.0, N, n_type[1], 1);//Compute intial likelihood


    //Track i times
    double i_array[N];
    if(INF == 0){
      for(k = 0; k < N; k++){
         i_array[k] = gsl_vector_get(i, k);
      }
    }

    //Open files to save MCMC output
    char str_local1[100]  = "";
    char str_local2[100]  = "";
    char str_local3[100]  = "";
    strcat(str_local1, "i_sum");
    strcat(str_local1, ID);
    strcat(str_local1, ".txt");
    strcat(str_local2, "gamma");
    strcat(str_local2, ID);
    strcat(str_local2, ".txt");
    strcat(str_local3, "ell");
    strcat(str_local3, ID);
    strcat(str_local3, ".txt");
    FILE * f_i_sum = fopen(str_local1, "w");
    FILE * f_gamma = fopen(str_local2, "w");
    FILE * f_ell   = fopen(str_local3, "w");
    int k1;
    double log_f_prior, counter;
    double a[n_update], b[N];                //for updating a smaller number of infection times
    for (k = 0; k < N; k++){
	       b[k] = (double) k;
    }

/***** MCMC Loop *****/
    printf("MCMC Loop...\n");
    for(j = 0; j < M; j++){

        gamma_current = gsl_ran_gamma(random, (n_type[0] + n_type[1])*ALPHA + 1, 1.0/(s_sum[0] + s_sum[1] + 0.01));          //Gibbs step for gamma from infectious period distribution
        loglike[0]  = likelihood(i, r, types_full, infecteds, beta_matrix_0, d_sum[0], s_sum[0], l_sum[0], ALPHA, gamma_current, N, n_type[0], 0);//Compute intial likelihood
        loglike[1]  = likelihood(i, r, types_full, infecteds, beta_matrix_1, d_sum[1], s_sum[1], l_sum[1], ALPHA, gamma_current, N, n_type[1], 1);//Compute intial likelihood


        //Sample length scale
        ell_prop[0] = gsl_ran_gaussian(random, 0.05) + ell_current[0];
        ell_prop[1] = gsl_ran_gaussian(random, 0.05) + ell_current[1];
        //printf("ell' = %f\n", ell_prop);
        if(ell_prop[0] > 0){
          sq_exp(dist_sparse_0, dist_sparse_0, VAR, ell_prop[0], S, S, sigma_0_prop);    //Compute Square Exponential Covariance Function
          gsl_matrix_add(sigma_0_prop, id);                                                          //Add correction term
          gsl_matrix_memcpy(sigma_0_chol_prop, sigma_0_prop);
          gsl_linalg_cholesky_decomp(sigma_0_chol_prop);                                            //Compyte Cholesky Decomposition
          gsl_matrix_memcpy(sigma_0_inverse_prop, sigma_0_chol_prop);
          gsl_linalg_cholesky_invert(sigma_0_inverse_prop);
          det_prop[0]          = chol_det(sigma_0_chol_prop, S);
          log_prior_ratio_ell  = 0.01*(ell_current[0] - ell_prop[0]);
          log_prior_ratio_f    = log_normal_pdf(f_0_bar, sigma_0_inverse_prop, det_prop[0], S)  - log_normal_pdf(f_0_bar, sigma_0_inverse, det[0], S);
          if(log(gsl_ran_flat(random, 0.0, 1.0)) < (log_prior_ratio_f + log_prior_ratio_ell)){
            ell_current[0]        = ell_prop[0];
            det[0]                = det_prop[0];
            gsl_matrix_memcpy (sigma_0, sigma_0_prop);
            gsl_matrix_memcpy (sigma_0_chol, sigma_0_chol_prop);
            gsl_matrix_memcpy (sigma_0_inverse, sigma_0_inverse_prop);
            sq_exp(dist_0_sorted, dist_sparse_0, VAR, ell_current[0], unique[0], S, k_s_d_0);



            if(j > BURN)
              ell_count++;
            //printf("ell = %f\n", ell_current);
          }
        }

        if(ell_prop[1] > 0){
          sq_exp(dist_sparse_1, dist_sparse_1, VAR, ell_prop[1], S, S, sigma_1_prop);    //Compute Square Exponential Covariance Function
          gsl_matrix_add(sigma_1_prop, id);                                                          //Add correction term
          gsl_matrix_memcpy(sigma_1_chol_prop, sigma_1_prop);
          gsl_linalg_cholesky_decomp(sigma_1_chol_prop);                                            //Compyte Cholesky Decomposition
          gsl_matrix_memcpy(sigma_1_inverse_prop, sigma_1_chol_prop);
          gsl_linalg_cholesky_invert(sigma_1_inverse_prop);
          det_prop[1]             = chol_det(sigma_1_chol_prop, S);
          log_prior_ratio_ell  = 0.01*(ell_current[1] - ell_prop[1]);
          log_prior_ratio_f    = log_normal_pdf(f_add, sigma_1_inverse_prop, det_prop[1], S)  - log_normal_pdf(f_add, sigma_1_inverse, det[1], S);
          if(log(gsl_ran_flat(random, 0.0, 1.0)) < (log_prior_ratio_f + log_prior_ratio_ell)){
            ell_current[1]        = ell_prop[1];
            det[1]                = det_prop[1];
            gsl_matrix_memcpy (sigma_1, sigma_1_prop);
            gsl_matrix_memcpy (sigma_1_chol, sigma_1_chol_prop);
            gsl_matrix_memcpy (sigma_1_inverse, sigma_1_inverse_prop);
            sq_exp(dist_1_sorted, dist_sparse_1, VAR, ell_current[1], unique[1], S, k_s_d_1);



            //if(j > BURN)
              //ell_count++;
            //printf("ell = %f\n", ell_current);
          }
        }

        // // //Sample beta
        for(counter = 0; counter < 5; counter++){
        rmvnorm_chol(random, S, mu1, sigma_0_chol, prior_draw);                    //sample from prior
        gsl_vector_memcpy(f_0_bar_prop, f_0_bar);
        gsl_vector_scale(prior_draw, DELTA);
        gsl_vector_scale(f_0_bar_prop, delta1);
        gsl_vector_add(f_0_bar_prop, prior_draw);

        rmvnorm_chol(random, S, mu1, sigma_1_chol, prior_draw);                    //sample from prior
        gsl_vector_memcpy(f_add_prop, f_add);
        gsl_vector_scale(prior_draw, DELTA);
        gsl_vector_scale(f_add_prop, delta1);
        gsl_vector_add(f_add_prop, prior_draw);
        gsl_vector_memcpy(f_1_bar_prop, f_add_prop);
        gsl_vector_add(f_1_bar_prop, f_0_bar_prop);


        gsl_blas_dgemv (CblasNoTrans, 1.0, sigma_0_inverse, f_0_bar_prop, 0.0, f_work);
        gsl_blas_dgemv (CblasNoTrans, 1.0, k_s_d_0, f_work, 0.0, f_0_prop);
        gsl_blas_dgemv (CblasNoTrans, 1.0, sigma_1_inverse, f_1_bar_prop, 0.0, f_work);
        gsl_blas_dgemv (CblasNoTrans, 1.0, k_s_d_1, f_work, 0.0, f_1_prop);


        //turn into beta
        for(k = 0; k < unique[0]; k++)
            gsl_vector_set(beta_vector_0_prop, k, exp(gsl_vector_get(f_0_prop, k)));

        for(k = 0; k < unique[1]; k++)
            gsl_vector_set(beta_vector_1_prop, k, exp(gsl_vector_get(f_1_prop, k)));

        beta_vector_to_matrix(beta_vector_0_prop, beta_matrix_0_prop, dist_mat_index_0, types_full, N_type[0], N, 0);   //turn proposal from vector to matrix
        beta_vector_to_matrix(beta_vector_1_prop, beta_matrix_1_prop, dist_mat_index_1, types_full, N_type[1], N, 1);   //turn proposal from vector to matrix
        d_sum_prop[0]    = double_sum(i, r, types_full, infecteds, beta_matrix_0_prop, d_sum_times_0, n_type[0], N_type[0], N, 0); //double_sum_beta_update(beta_matrix_prop, d_sum_times, N);               //compute proposal double sum and likelihood
        d_sum_prop[1]    = double_sum(i, r, types_full, infecteds, beta_matrix_1_prop, d_sum_times_1, n_type[1], N_type[1], N, 1); //double_sum_beta_update(beta_matrix_prop, d_sum_times, N);               //compute proposal double sum and likelihood
        loglike_prop[0]  = likelihood(i, r, types_full, infecteds, beta_matrix_0_prop, d_sum_prop[0], s_sum[0], l_sum[0], ALPHA, gamma_current, N, n_type[0], 0);
        loglike_prop[1]  = likelihood(i, r, types_full, infecteds, beta_matrix_1_prop, d_sum_prop[1], s_sum[1], l_sum[1], ALPHA, gamma_current, N, n_type[1], 1);
        //printf("log p_acc = %f\n", loglike_prop[1] -loglike[1] + loglike_prop[0] - loglike[0]);
        //printf("loglike[0] = %f\t loglike[1] = %f\n", loglike_prop[0], loglike_prop[1]);

        if(log(gsl_ran_flat(random, 0.0, 1.0)) < (loglike_prop[0] - loglike[0] + loglike_prop[1] - loglike[1])){                   //MH Step
            //If accepted then update values...
            gsl_vector_memcpy(beta_vector_0, beta_vector_0_prop);
            gsl_matrix_memcpy(beta_matrix_0, beta_matrix_0_prop);
            gsl_vector_memcpy(f_0, f_0_prop);
            gsl_vector_memcpy(f_0_bar, f_0_bar_prop);
            gsl_vector_memcpy(beta_vector_1, beta_vector_1_prop);
            gsl_matrix_memcpy(beta_matrix_1, beta_matrix_1_prop);
            gsl_vector_memcpy(f_1, f_1_prop);
            gsl_vector_memcpy(f_1_bar, f_1_bar_prop);
            gsl_vector_memcpy(f_add, f_add_prop);
            memcpy(loglike, loglike_prop, sizeof(loglike));
            memcpy(d_sum, d_sum_prop, sizeof(d_sum));
            beta_count[0]++;
        }
        }

       /* INFECTION TIMES MH STEP */
        /* This section updates each infection time using a Metropolis Hastings algorithm. Foe each infection time
        it proposes a new infection time by i_{new} = r - Gamma(alpha, gamma). It then computes the acceptance probability. */
        if(INF == 0){
        gsl_ran_choose (random, a, n_update, b, N, sizeof (double));           //choose infection times to update
            for(k = 0; k < n_update; k++){
                if(gsl_vector_get(infecteds, a[k]) == 1){ //If individual k is infected...
		                i_prop = gsl_vector_get(r, a[k]) - gsl_ran_gamma(random, ALPHA, 1.0/gamma_current);     //propose new infection time
                    gsl_vector_memcpy(i_props, i);
                    gsl_vector_set(i_props, a[k], i_prop);
                    if(a[k] < (N/2.0)){
                      s_sum_prop[0] = s_sum[0] + gsl_vector_get(i, a[k]) - i_prop;                               //compute proposed sums and likelihoods
                      l_sum_prop[0] = l_sum[0] - log(gsl_vector_get(r, a[k]) - gsl_vector_get(i, a[k])) + log(gsl_vector_get(r, a[k]) - i_prop);
                      s_sum_prop[1] = s_sum[1];
                      l_sum_prop[1] = l_sum[1];
                    } else {
                      s_sum_prop[1] = s_sum[1] + gsl_vector_get(i, a[k]) - i_prop;                               //compute proposed sums and likelihoods
                      l_sum_prop[1] = l_sum[1] - log(gsl_vector_get(r, a[k]) - gsl_vector_get(i, a[k])) + log(gsl_vector_get(r, a[k]) - i_prop);
                      s_sum_prop[0] = s_sum[0];
                      l_sum_prop[0] = l_sum[0];
                    }                                             //replace current time by proposed time
                    d_sum_prop[0]    = double_sum(i_props, r, types_full, infecteds, beta_matrix_0, d_sum_times_0_prop, n_type[0], N_type[0], N, 0); //double_sum_beta_update(beta_matrix_prop, d_sum_times, N);               //compute proposal double sum and likelihood
                    d_sum_prop[1]    = double_sum(i_props, r, types_full, infecteds, beta_matrix_1, d_sum_times_1_prop, n_type[1], N_type[1], N, 1); //double_sum_beta_update(beta_matrix_prop, d_sum_times, N);               //compute proposal double sum and likelihood
                    loglike_prop[0]  = likelihood(i_props, r, types_full, infecteds, beta_matrix_0, d_sum_prop[0], s_sum_prop[0], l_sum_prop[0], ALPHA, gamma_current, N, n_type[0], 0);
                    loglike_prop[1]  = likelihood(i_props, r, types_full, infecteds, beta_matrix_1, d_sum_prop[1], s_sum_prop[1], l_sum_prop[1], ALPHA, gamma_current, N, n_type[1], 1);
                    q_prob = gsl_ran_gamma_pdf(gsl_vector_get(r, a[k]) - gsl_vector_get(i, a[k]), ALPHA, 1.0/gamma_current)/gsl_ran_gamma_pdf(gsl_vector_get(r, a[k]) - i_prop, ALPHA, 1.0/gamma_current);
                    p_acc = loglike_prop[0] + loglike_prop[1] - loglike[0] - loglike[1] + log(q_prob);                                           //compute new acceptance probabilities
                    //printf("log(p_acc) = %f\n", p_acc);
                    if(log(gsl_ran_flat(random, 0.0, 1.0)) < p_acc){
                        //If accepted update values...
                        gsl_vector_memcpy(i, i_props);
                        gsl_matrix_memcpy(d_sum_times_0, d_sum_times_0_prop);
			                  gsl_matrix_memcpy(d_sum_times_1, d_sum_times_1_prop);
                        d_sum[0] = d_sum_prop[0];
                        d_sum[1] = d_sum_prop[1];
                        loglike[0] = loglike_prop[0];
                        loglike[1] = loglike_prop[1];
                        s_sum[0] = s_sum_prop[0];
                        s_sum[1] = s_sum_prop[1];
                        l_sum[0] = l_sum_prop[0];
                        l_sum[1] = l_sum_prop[1];
                        i_count = i_count + 1;
                    }
                }
            }
            //Compute sum of infection times and save to file
            i_sum = 0.0;
            for(k = 0; k < N; k++){
                if(gsl_vector_get(infecteds, k)==1)
                    i_sum = i_sum + gsl_vector_get(i, k);
            }

            fprintf(f_i_sum, "%f\n", i_sum);
        }
        //printf("loglike[0] = %f \t loglike[1] = %f\n", loglike[0], loglike[1]);
        //Save MCMC output
        if(j % 1000 == 0)
           printf("%.2f %%...\n", (double) j/(M)*100);
        if(j == (double) M / 2.0 || j == (double) M / 4.0)
            printf("Current Beta Acceptance Rate = %.2f %% \nCurrent i Acceptance Rate  = %.2f %%\n", ((double)beta_count[0]/(j*20))*100, ((double)i_count/((j-BURN)*n_update))*100);

        if(j > BURN){
          gsl_vector_add(mean_0, beta_vector_0);
          gsl_vector_add(mean_1, beta_vector_1);

          gsl_vector_memcpy(var_0_work, beta_vector_0);
          gsl_vector_mul(var_0_work, beta_vector_0);
          gsl_vector_add(var_0, var_0_work);
          gsl_vector_memcpy(var_1_work, beta_vector_1);
          gsl_vector_mul(var_1_work, beta_vector_1);
          gsl_vector_add(var_1, var_1_work);
          fprintf(f_gamma, "%f\n", gamma_current);
          fprintf(f_ell, "%f \t %f \n", ell_current[0], ell_current[1]);

        }


    }
    fclose(f_gamma);
    char str_mean0[100]  = "";
    char str_mean1[100]  = "";
    char str_var_0[100]  = "";
    char str_var_1[100]  = "";
    strcat(str_mean0, "mean0_");
    strcat(str_mean1, "mean1_");
    strcat(str_mean0, ID);
    strcat(str_mean1, ID);
    strcat(str_mean0, ".txt");
    strcat(str_mean1, ".txt");
    strcat(str_var_0, "var0_");
    strcat(str_var_0, ID);
    strcat(str_var_0, ".txt");
    strcat(str_var_1, "var1_");
    strcat(str_var_1, ID);
    strcat(str_var_1, ".txt");
    gsl_vector_scale(mean_0, 1.0/(M-BURN));
    gsl_vector_scale(mean_1, 1.0/(M-BURN));
    gsl_vector_scale(var_0, 1.0/(M-BURN));
    gsl_vector_scale(var_1, 1.0/(M-BURN));
    FILE * f_mean0 = fopen(str_mean0, "w");
    FILE * f_mean1 = fopen(str_mean1, "w");
    gsl_vector_fprintf(f_mean0, mean_0, "%f");
    gsl_vector_fprintf(f_mean1, mean_1, "%f");
    fclose(f_mean0);
    fclose(f_mean1);
    gsl_vector_memcpy(mean_0_work, mean_0);
    gsl_vector_mul(mean_0, mean_0_work);
    gsl_vector_sub(var_0, mean_0);
    gsl_vector_memcpy(mean_1_work, mean_1);
    gsl_vector_mul(mean_1, mean_1_work);
    gsl_vector_sub(var_1, mean_1);
    FILE * f_var_0 = fopen(str_var_0, "w");
    FILE * f_var_1 = fopen(str_var_1, "w");
    gsl_vector_fprintf(f_var_0, var_0, "%f");
    gsl_vector_fprintf(f_var_1, var_1, "%f");
    fclose(f_var_0);
    fclose(f_var_1);
    printf("Beta Acceptance Rate = %f Infection Times Acceptance Rate = %f %% \n", ((double)beta_count[0]/(M*20))*100, ((double)i_count/(M*N))*100);
    time(&stop);
    printf("Elapsed Time %f \n", difftime(stop, start));

 return 0;
}


double double_sum(const gsl_vector * i, const gsl_vector * r, const gsl_vector * types, const gsl_vector * infecteds, const gsl_matrix * beta, gsl_matrix * times, int n, int N1, int N2, int t){
    /* Computes Double Sum in Likelihood Function */
    /* i, r         -- vectors with infection and removal times
     * infecteds    -- vector with 0 if susceptible 1 if infected
     * beta         -- matrix of beta_{i,j} values
     * times        -- matrix of min differences in times
     * n,           -- number of infecteds
     * N            -- population size
     */

    double j, k, result;
    gsl_matrix_set_all(times, 0.0);
    gsl_matrix * work = gsl_matrix_alloc(N2, N1);
    for(j = 0; j < N2; j++){
        if(gsl_vector_get(infecteds, j) == 1){
            for(k = 0; k < N1; k++){
                gsl_matrix_set(times, j, k , GSL_MIN(gsl_vector_get(r, j), gsl_vector_get(i, k+(t*N2/2.0))) - GSL_MIN(gsl_vector_get(i, j), gsl_vector_get(i, k+(t*N2/2.0))));
            }
        }
    }
    gsl_matrix_memcpy(work, times);
    gsl_matrix_mul_elements(work, beta);
    result = 0.0;
    for(j = 0; j < N2; j++){
        for(k = 0; k < N1; k++){
            result = result + gsl_matrix_get(work, j, k);
        }
    }
    gsl_matrix_free(work);

    return result;

}



double  likelihood(const gsl_vector * i, const gsl_vector * r, const gsl_vector * types, const gsl_vector * infecteds, const gsl_matrix * beta, double d, double sum, double log_sum, double alpha, double gamma, int N, int n, int t){
    /* i, r         -- vectors with infection and removal times
     * infecteds    -- vector with 0 if susceptible 1 if infected
     * beta         -- matrix of beta_{i,j} values
     * d            -- double sum value
     * sum          -- sum(r-i)
     * log_sum      -- sum(log(r-i))
     * alpha        -- rate parameter of infectious period distribution
     * gamma        -- shape parameter of infectious period distribution
     * n,           -- number of infecteds
     * N            -- population size
     * t            -- type (0 or 1)
     */


     gsl_vector * y = gsl_vector_alloc(N/2.0);
     gsl_vector_set_all(y, 1.0);
     int j, k;
     double bet_sum = 0.0;
     double loglike;


     //Compute pressure exerted by current infecteds on each susceptible
     for(j = 0; j < N/2.0; j++){
         if(gsl_vector_get(infecteds, j + t*(N/2.0)) == 1){
           gsl_vector_set(y, j, 0.0);
             for(k = 0; k < N; k++){
                 if(gsl_vector_get(i, k) < gsl_vector_get(i, j + t*(N/2.0)) && gsl_vector_get(r, k) > gsl_vector_get(i, j + t*(N/2.0))){
                     gsl_vector_set(y, j, gsl_vector_get(y, j) + gsl_matrix_get(beta, k, j));
                 }
             }
         }
         else
             gsl_vector_set(y, j, 1.0);                      //if not infected set to 1 to preserve sum of logs
     }
     for(j = 0; j < N/2.0; j ++){
          gsl_vector_set(y, j, log(gsl_vector_get(y, j)));    //log pressure values
     }
     int min_index = gsl_vector_min_index(i);
     if(min_index >= (N/2.0)*t && min_index < ((N/2.0) + t*(N/2)))
        gsl_vector_set(y, min_index-(t*N/2.0), 0.0);
     //compute sum of log values
     for(j = 0; j < N/2.0; j++)
         bet_sum = bet_sum + gsl_vector_get(y, j);

     loglike =  bet_sum - d + (alpha - 1)*log_sum  - gamma*sum;
     //printf("%f \t %f \t %f\t %f \n", bet_sum,  d, (alpha - 1)*log_sum, gamma*sum);

     gsl_vector_free(y);
     return loglike;
}




void distance(const gsl_matrix *x, int n, gsl_matrix *result){
    /* Distance Matrix where M_{j, k} is the distance between x_j and x_k */
    /*
     *  x       matrix of positions (n*2)
     *  n       length of position matrix
     *  result  distance matrix (n*n)
     */
    int j, k;
    double diff1, diff2;
    gsl_matrix * temp = gsl_matrix_alloc(n, n);
    for(j = 0; j < n; j++){
        for(k = 0; k < (j+1); k++){
            diff1 = pow(gsl_matrix_get(x, j, 0) - gsl_matrix_get(x, k, 0), 2); /*Calculate difference between x coords*/
            diff2 = pow(gsl_matrix_get(x, j, 1) - gsl_matrix_get(x, k, 1), 2); /*Calculate difference between y coords*/
            gsl_matrix_set(temp, j, k, sqrt(diff1 + diff2));                   /* Euclidean Distance*/
        }
    }
    gsl_matrix_memcpy(result, temp); /*Copy matrix for transposition */
    gsl_matrix_transpose(result); /*Transpose Matrix to upper triangular*/
    gsl_matrix_add(result, temp); /*Sum matrices for full distance matrix*/
    gsl_matrix_free(temp);

}

void beta_vector_to_matrix(const gsl_vector * v, gsl_matrix * m, const gsl_matrix * index, const gsl_vector * types, const int N1, const int N2, const int t){
    /* Turns beta vecotr from GP into N*N matrix
        * v           -- vector of beta_{i,j}
        * m           -- beta_{i,j} matrix to be outputed
        * index       -- matrix showing whcih element of vector corresponds to which matrix element
        * N           -- Population Size
    */

    int j, k;
    for(j = 0; j < N2; j++){
        for(k = 0; k < N1; k++){
            if(j == k){
                gsl_matrix_set(m, j, k, 0.0);
            } else{
                gsl_matrix_set(m, j, k, gsl_vector_get(v, gsl_matrix_get(index, j, k)));
            }
        }
    }
}


void sq_exp(const gsl_vector *x, const gsl_vector *y, double alph, double ell, const int n1, const int n2, gsl_matrix * result){
    /* Squared Exponential Covariance Function */
    /*
     *  x, y        -- input vectors
     *  alph, ell   -- GP Parameters
     *  n1, n2      -- length of x and y
     *  results     -- output matrix of dimension n1*n2 with the covariance matrix
     */

    int count1, count2;
    double covar;
    for(count1 = 0; count1 < n1; count1++){
        for(count2 = 0; count2 < n2; count2++){
            covar = pow(alph, 2)*exp(-1*pow(gsl_vector_get(x, count1) - gsl_vector_get(y, count2), 2)/pow(ell, 2));
            gsl_matrix_set(result, count1, count2, covar);
        }
    }
}

void rmvnorm_chol(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
     /* Multivariate Normal Distribution Sample */
     /*
      *  r          -- GSL random number initialiseresult*result
      *  n          -- length of sample
      *  var        -- variance matrix (n x n)
      *  result     -- vector where sample will be stored (n x 1)
      */
     int k;
     gsl_matrix * work = gsl_matrix_alloc(n,n);
     gsl_matrix_memcpy          (work,var);                     //Preserve variance matrix
     for(k=0; k<n; k++)
         gsl_vector_set(result, k, gsl_ran_ugaussian(r));      //Generate standard Gaussian variates

     gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, work, result);
     gsl_vector_add (result, mean);
     gsl_matrix_free(work);
}

double chol_det(const gsl_matrix * A, int N){

  int j;

  double result;
  result = 0.0;

  for(j = 0; j < N; j++){
    result = result + log(gsl_matrix_get(A, j, j));
  }

  result = 2*result;
  return result;
}


double log_normal_pdf(const gsl_vector * x, const gsl_matrix * A_inverse, const double log_A_det, int N){

  double result, temp;
  gsl_vector * work = gsl_vector_alloc(N);

  gsl_vector_memcpy(work, x);
  gsl_blas_dgemv(CblasNoTrans, 1.0, A_inverse, x, 0.0, work);
  gsl_blas_ddot(x, work, &temp);

  result = -0.5*N*log(2*M_PI) - 0.5*log_A_det - 0.5*temp;
  return result;

}
