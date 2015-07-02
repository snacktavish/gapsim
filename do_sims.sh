#!/bin/bash
set -x
n=$1
i=1
while true
do
        dirn="rep$i"
        mkdir "${dirn}"
        cd "${dirn}"
        python ../gap_sim.py -l 100000 > sim.nex 2>sim-err.txt || exit
        paup -n ../master.nex
        paired_invariants_cull.py --p-inv=0.5 nontransposed.nex > culled.nex || exit
        paup -n ../masterculled.nex
        
        cd ..
        i=$(expr 1 + $i)
        if test $i -gt $n
        then
                exit
        fi
done