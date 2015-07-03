#!/bin/bash
i=1
totincc=0
totcc=0
while true
do
        dirn="rep$i"
        incdistfromtrue=$(grep -E '^[2]   [0-2]+ ' "${dirn}/sim.log" | awk '{print $2}' | head -n1)
        cdistfromtrue=$(grep -E '^[2]   [0-2]+ ' "${dirn}/sim.log" | awk '{print $2}' | tail -n1)
        if test $incdistfromtrue -eq 0
        then
                totincc=$(expr 1 + $totincc)
        fi
        if test $cdistfromtrue -eq 0
        then
                totcc=$(expr 1 + $totcc)
        fi
        i=$(expr 1 + $i)
        ndirn="rep$i"
        if ! test -d $ndirn
        then
                break
        fi
done
echo $(expr $i - 1) " runs"
echo "$totincc correct from raw analysis"
echo "$totcc correct from culled analysis"