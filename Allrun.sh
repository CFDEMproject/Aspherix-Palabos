#!/bin/bash

#===================================================================#
# allrun script for testcase
#===================================================================#
# make sure the path to Palabos is set correctly in the Makefile
# compile with: make clean / make
# run this script

trap cleanup 1 2 3 6 15

function cleanup ()
{
    pkill -P $DEMPID
    pkill -P $CFDPID
#    kill -15 $DEMPID
#    kill -15 $CFDPID
}

# start DEM separately
bash DEMrun.sh &
DEMPID=$!

#source DEMrun.sh&

# start CFD ./sedimentingSphere N uMax rho_f mu_f v_inf maxT outDir
#bash -c "mpirun -np 2 ./sedimentingSphere 5 0.1 970 0.373 2.5 1 ./log " & # 2>&1 | tee log_CFD" &
./sedimentingSphere 10 1 970 0.373 2.5 0.1 ./log #2>&1 | tee log_CFD
CFDPID=$!

while [ $(kill -0 $DEMPID 2> /dev/null; echo $?) == "0" ] && [ $(kill -0 $CFDPID 2> /dev/null; echo $?) == "0" ]
do
    sleep 1s
done

trap - 1 2 3 6 15

echo "clean up? (otherwise ctrl-C)"
read
rm logtmp/*
rm post/*
rm logpost/*
touch logtmp/.gitignore
