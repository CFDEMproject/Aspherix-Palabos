#!/bin/bash

#===================================================================#
# allrun script for testcase
#===================================================================#
# make sure the path to Palabos is set correctly in the Makefile
# compile with: make clean / make
# run this script

trap cleanup 1 2 3 6 15

# start DEM separately
bash DEMrun.sh &
BASHPID=$!

#source DEMrun.sh&

function cleanup ()
{
    pkill -P $BASHPID
    kill -15 $BASHPID
}

# start CFD ./sedimentingSphere N uMax rho_f mu_f v_inf maxT outDir
mpirun -np 2 ./sedimentingSphere 5 0.1 970 0.373 2.5 1 ./log 2>&1 | tee log_CFD
#./sedimentingSphere 10 1 970 0.373 2.5 1 ./log 2>&1 | tee log_CFD

cleanup

echo "clean up? (otherwise ctrl-C)"
read
rm logtmp/*
rm post/*
rm logpost/*
touch logtmp/.gitignore
