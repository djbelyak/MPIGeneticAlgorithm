#!/bin/sh
for i in `seq 1 32`;
do
	mpirun -n $i MPIGeneticAlgorithm |tee $i.txt
done 