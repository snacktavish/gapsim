#!/usr/bin/env python
import random
import itertools
# 4 taxon case
# 2 states
# 81 possible patterns
states = ['-',0, 1]

poss_states = set()

import itertools
s=[ states, states, states, states ]
poss=list(itertools.product(*s))


\[Lambda] = 0.5;
\[Mu] = 1;
ta = .01
tb = 1
ti = ta/100
tc = tb
td = ta


# Tree:
#   A         C
#    \       /
#     \     /
#      \___/
#     /    \
#    B      D

#Expected site pattern probabilities.

#calculate distance A-B
#A-C
#A-D


fmt = '''
#NEXUS
begin taxa;
 dimensions ntax = 4;
 taxlabels A B C D ;
end;

begin distances;
    Format triangle=upper;
    matrix
        A 0.0 {dab} {dac} {dad} 
        B     0.0 {dbc} {dbd} 
        C         0.0 {dcd}
        D             0.0;

end;
'''
.format(dab=dab, dac=dac, dad=dad, dbc=dbc, dbd=dbd, dcd=dcd)