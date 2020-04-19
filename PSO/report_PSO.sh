#! /bin/zsh

g++ -o report PSO_simulation100.cpp

for dim in 2 5 20
do
    ./report $dim
done

rm report
