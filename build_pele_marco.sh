module load gcc
module load cuda
module load cmake
module load python

export BASE=/ccs/home/cyrush/WORKSCRATCH/2021_11_cyrush_pelelm/
export AMREX_HOME=${BASE}/amrex/
export PELELM_HOME=${BASE}/PeleLM
export PELE_PHYSICS_HOME=${BASE}/PelePhysics
export IAMR_HOME=${BASE}/IAMR
export PELEMP_HOME=${BASE}/PeleMP
export AMREX_HYDRO_HOME=${BASE}/AMReX-Hydro

echo $BASE


make USE_MPI=TRUE USE_CUDA=TRUE  AMREX_HOME=../amrex/


