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

# Tree:
#   A           C
#    \         /
#     \       /
#      \i___j/
#      /     \
#     B      D


# Sequence length is length 
seqlen = 10 # seqlen is length of observed sequence at tip a
insrate = 0.5;
delrate = 1;
ta = 1 #long branch
tb = .01 #short branch
ti = ta/100
tc = ta
td = tb
pinvar = 0.5


#generator will 

#for a base in 'a' simulate other tips - walk through the tree branch by branch

#simulate realization of columns.


#Expected site pattern probabilities.

#calculate distance A-B
#A-C
#A-D


#generator that 


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