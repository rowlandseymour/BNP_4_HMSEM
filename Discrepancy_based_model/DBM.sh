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
N=1000    #Population Size
beta_1=3
beta_2=2
beta=(0.005 0.001)    #array of beta values
alpha=6         #Infectious Period Distribution Rate Parameter
gamma=3         #Infectious Period Distribution Shape Parameter
M=10000   	    #Number of MCMC Iterations
seed=15666         #RNG Seed
var=6
ell=4
inf=1           #0 => Infer infection times, 1 => Infection times known
ID=6            #Run ID
delta=0.01
S=30
echo ${beta[*]} >  beta_values.txt
P=$(echo ${beta[*]} | wc -w)

## Compile c Code and Move
echo 'Compiling Code...'
cd ~/HMSEM/Multi_Type
gcc Simulation.c -lm -lgsl -lgslcblas -o Type_dist_Sim.o
gcc DBM.c -g -lm -fopenmp -lgsl -lgslcblas -o Type_dist_Inf.o

## Run Code
echo 'Simulating Outbreak...'
./Type_dist_Sim.o $N $alpha $gamma $beta_1 $beta_2 $seed $P $ID
echo 'Processing data in R'
Rscript Multi_Output_Distance_Preperation.R $N $ID $S
echo 'MCMC Inference...'
./Type_dist_Inf.o $N $alpha $M $seed $inf $var $ell $delta $S $P $ID
rm Type_dist_Sim.o
rm Type_dist_Inf.o
Rscript DBM_plot.R $N $beta_1 $beta_2 $inf $gamma $ID