#!/bin/bash

#===================================================================#
# allrun script for testcase
#===================================================================#

#mpirun -np 2 $CFDEM_ASX_EXEC -in in.lbdem 2>&1 | tee log_DEM
$CFDEM_ASX_EXEC -in in.lbdem # 2>&1 | tee log_DEM
