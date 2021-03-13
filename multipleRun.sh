#!/bin/bash

machine=$1
fromCore=$2
toCore=$((fromCore+4))

g=0.5
networkSize=160000
ensembleSize=80000

coreNum=$((fromCore+1))
while [ $coreNum -le $toCore ]
do
    # echo $coreNum
    spg run $machine ./generate.sh $g $networkSize $ensembleSize $coreNum
    coreNum=$((coreNum+1))
done
