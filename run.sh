#!/bin/bash

g=$1
networkSize=$2
ensembleSize=$3
machine=$4
coreNum=$5

spg run ${machine} ./mBFW.sh ${g} ${networkSize} ${ensembleSize} ${coreNum}