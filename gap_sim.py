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

def gen_column(n):
    for index in xrange(n):
        column_set = simulate_columns_from_A()
        for col in column_set:
            yield col
            
ins+n(ins+del))
def simulate_columns_from_A():
  if random.random( )

#simulate changes on all branches from root 'a'
def branch_sim(im,base,bl):



def sim_mat(out, seqlen):
    for column in gen_column(seqlen):
        out.write(column)
        out.write('\n')

sim_mat(sys.stdout, seqlen)

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