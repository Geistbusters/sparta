#PBS -S /bin/csh
#PBS -N SPANT
# This example uses the Sandy Bridge nodes
# User job can access ~31 GB of memory per Sandy Bridge node.
# A memory intensive job that needs more than ~1.9 GB
# per process should use less than 16 cores per node
# to allow more memory per MPI process. This example
# asks for 32 nodes and 8 MPI processes per node.
# This request implies 32x8 = 256 MPI processes for the job.
## David uses: broadwell, 28 nodes
#PBS -q devel
#PBS -l select=1:ncpus=1:model=bro
#PBS -l min_walltime=1:00:00,max_walltime=2:00:00
#PBS -W group_list=e1305
#PBS -m e
#PBS -j oe

module load comp-intel/2020.4.304  mpi-hpe/mpt.2.28_25Apr23_rhel87

./spa_rome_2 < in.test 

