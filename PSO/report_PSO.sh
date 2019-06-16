#!/bin/zsh
## sh エミュレーションモード
emulate -R sh
## do
eval "g++ PSO_simulation100.cpp"
for var in 2 5 20
do eval "./a.out $var"
done