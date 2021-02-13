#!/bin/bash

g=$1
networkSize=$2
ensembleSize=$3
coreNum=$4

name=N${networkSize}G${g}E${ensembleSize}C${coreNum}

g++ -O3 -march=native -std=c++17 -o bin/${name}.out main-generate.cpp

./bin/${name}.out ${networkSize} ${g} ${ensembleSize} ${coreNum} >> log.txt
rm bin/${name}.out