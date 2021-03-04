#!/bin/bash

networkSizeList=(10000 20000 40000 80000 160000 320000 640000 1280000 2560000 5120000 10240000)
g=0.2

for networkSize in ${networkSizeList[@]}
do
    ./bin/data.out ${g} ${networkSize}
    # echo ${networkSize} ${g}
done
