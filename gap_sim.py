#!/usr/bin/env python
import random
import itertools
import math
import sys
import time
# 4 taxon case
# 2 states
# 81 possible patterns
states = ['-',0, 1]

RNG=random.Random()
import os
if 'SIMSEED' in os.environ:
    s = int(os.environ['SIMSEED'])
else:
    s = RNG.randint(2, 293485729298)
   
sys.stderr.write("Seed is {}\n".format(s))
RNG.seed(s)


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
seqlen = int(sys.argv[1]) # seqlen is length of observed sequence at tip a
insrate = .4 #E-10
delrate = .5 
ta = 1.0 #long branch
tb = .01 #short branch
ti = tb
tc = ta
td = tb
pinvar = 0.5
one_minus_pinv = 1.0 - pinvar

epsilon = 1.0E-10

subst = {
   'A' : 'GCT',
   'G' : 'ACT',
   'C' : 'AGT',
   'T' : 'AGC'
}
#generator will 

#for a base in 'a' simulate other tips - walk through the tree branch by branch

#simulate realization of columns.


#Expected site pattern probabilities.

#calculate distance A-B
#A-C
#A-D


## Notes make debugging optional,
## Write some tests..


def gen_column(n):
    for index in xrange(n):
        column_set = simulate_columns_from_A()
        for col in column_set:
            yield col
            
def simulate_columns_from_A():
  if RNG.random() < pinvar:
     return [RNG.choice('AGCT') * 4]
  subcol = [RNG.choice('AGCT') ]+['-']*5  #initial 
  subseq = [subcol]
  branch_sim(subseq, 0, 4, ta/one_minus_pinv)
  branch_sim(subseq, 4, 1, tb/one_minus_pinv)
  branch_sim(subseq, 4, 5, ti/one_minus_pinv)
  branch_sim(subseq, 5, 2, tc/one_minus_pinv)
  branch_sim(subseq, 5, 3, td/one_minus_pinv)
  c = [''.join(lis[:4]) for lis in subseq]
  return c

VERBOSE = False
#simulate changes on all branches from root 'a'
def branch_sim(subseq, start, end, blen):
    nondel=[]
    for ite, subcol in enumerate(subseq):
       if subcol[start]!='-':
           nondel.append(ite)
       subcol[end] = subcol[start]
    n = len(nondel)
    indelrate = insrate+delrate
    eventrate = insrate+(n*indelrate)
    currpoint = 0.0
    if VERBOSE:
        print("start", nondel)
    while True:
        wt = RNG.expovariate(eventrate)
        currpoint+=wt
        if currpoint >= blen:
              break
        pimmlink = (epsilon+insrate)/eventrate        
        pdel = (n*delrate)/eventrate
        u = RNG.random()
        if u < pimmlink:
            newcol=['-']*6
            newcol[end]=RNG.choice('AGCT')
            subseq.insert(0,newcol)
            nondel = [0] + [ite + 1 for ite in nondel]
            eventrate += indelrate
            if VERBOSE:
                print("immortal {}".format(eventrate),nondel)  
        else:
            ndi = RNG.randrange(0, len(nondel))
            u -= pimmlink   
            if u < pdel:
               ci = nondel.pop(ndi)
               subseq[ci][end]='-'
               eventrate -= indelrate
               if VERBOSE:
                    print("deletion {}".format(eventrate),nondel)  
            else:
                newcol=['-']*6
                newcol[end]=RNG.choice('AGCT')
                ci = nondel[ndi]
                subseq.insert(ci+1,newcol)
                newnondel = [ci+1]
                for ite in  nondel:
                   if ite > ci:
                      newnondel.append(ite+1)
                   else: 
                      newnondel.append(ite)
                nondel = newnondel
                eventrate += indelrate
                if VERBOSE:
                    print("insertion {}".format(eventrate),nondel)  
    subprob = 0.75 - 0.75*math.exp((-4.0*blen)/3.0)
    for ci in nondel:
        if RNG.random() < subprob:
            anc = subseq[ci][end]
            desc = RNG.choice(subst[anc])
            subseq[ci][end] = desc
    return subseq
        
def sim_mat(out, seqlen):
    n = 0
    for column in gen_column(seqlen):
        if n > 0:
            out.write('\n')
        out.write(column)
        n += 1
    return n
from cStringIO import StringIO
buffer = StringIO()
nc = sim_mat(buffer, seqlen)
numbered = ['c_{n} {i}'.format(n=n, i=i) for n, i in enumerate(buffer.getvalue().split('\n'))]
out = sys.stdout
out.write('''#NEXUS
begin data ;
    dimensions ntax = 4 nchar = {c} ;
    taxlabels A B C D ;
    format transpose datatype=dna gap = '-';
matrix 
{m}
;
end;
begin trees;
tree true = [&U] (A,B,(C,D));
end;
begin paup;
    set crite = like;
    lset nst = 1 basefreq = eq pinv = est;
    lscore 1;
    describe / brl ;
    alltrees ;
    showtree ;
end;
'''.format(c=nc, m='\n'.join(numbered)))

#generator that 

"""
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
'''.format(dab=dab, dac=dac, dad=dad, dbc=dbc, dbd=dbd, dcd=dcd)"""