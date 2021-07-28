######################################
# Simulation, Inference and Plotting #
#           for GP MCMC              #
######################################

#This bash script simulates a disease with a heterogenously mixing
#infection rate. It then uses GPs and MCMC to infer the infection
#rate and data augmentation to infer the unobserved infection times.
#It then plots the results and attaches these plots to an email.


SECONDS=0
## Set Variables
N=500    #Population Size
beta_1=3
beta_2=2
beta=(0.01 0.007)    #array of beta values
alpha=4         #Infectious Period Distribution Rate Parameter
gamma=3         #Infectious Period Distribution Shape Parameter
M=10000   	    #Number of MCMC Iterations
seed=156         #RNG Seed
var=6
ell=4
inf=0           #0 => Infer infection times, 1 => Infection times known
ID=5            #Run ID
delta=0.15
S=30
rho=0.9
echo ${beta[*]} > /local/pmxrs4/Type_dist/input/beta_values.txt
P=$(echo ${beta[*]} | wc -w)

## Compile c Code and Move
echo 'Compiling Code...'
gcc Simulation.c -lm -lgsl -lgslcblas -o Multi_Output_Sim.o
gcc Multi_Output_Inference.c -g -lm -fopenmp -lgsl -lgslcblas -o Multi_Output_Inf.o

## Run Code
echo 'Simulating Outbreak...'
./Multi_Output_Sim.o $N $alpha $gamma $beta_1 $beta_2 $seed $P $ID
echo 'Processing data in R'
Rscript Multi_Output_Distance_Preperation.R $N $ID $S
echo 'MCMC Inference...'
./Multi_Output_Inf.o $N $alpha $M $seed $inf $var $ell $delta $rho $S $P $ID
rm Multi_Output_Sim.o
rm Multi_Output_Inf.o
Rscript Multi_Output_plot.R $N $beta_1 $beta_2 $inf $gamma $ID


