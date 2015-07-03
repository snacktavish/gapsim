#!/usr/bin/env python
import random
import itertools
import math
import sys
import time
import argparse

# 4 taxon case
# 2 states

RNG=random.Random()
import os
if 'SIMSEED' in os.environ:
    s = int(os.environ['SIMSEED'])
else:
    s = RNG.randint(2, 293485729298)
   
sys.stderr.write("Seed is {}\n".format(s))
RNG.seed(s)

sys.stderr.write('''
# Tree:
#   A  A1    C  C1
#    \/I1     \/I4
#     \       /
#      \i___j/
#      /     \ 
#     B       D
#
#
#
#
# Tree:
#   0  1     3  4
#    \/       \/
#    6\       /9
#      \ 7 _8/
#      /     \ 
#     2       5
#
''')

parser = argparse.ArgumentParser(description='Process some arguments.')
parser.add_argument('-l','--seqlen', help='average sequence length (if avg column indel = 1)', default = 100, type=int)
parser.add_argument('-ins','--insrate', help='insertion rate', default = 0.5, type=float)
parser.add_argument('-del','--delrate', help='deletion rate',  default = 1.0, type=float)
parser.add_argument('-invar','--pinvar', help='proportion of invariant block',  default = 0.5, type=float)
parser.add_argument('-lb','--lb', help='length of path from A to MRCA of A,A1|B,C,C1,D',  default = 0.6, type=float)
parser.add_argument('-sb','--sb', help='branch length at b',  default = 0.01, type=float)

args = vars(parser.parse_args())


# Sequence length is length 
seqlen = args['seqlen']
insrate = args['insrate']
delrate = args['delrate'] 
pinvar = args['pinvar']
sb = args['sb']
lb = args['lb']
one_minus_pinv = 1.0 - pinvar
probzlen = (1 - (insrate/delrate))
geomprobins = (insrate/delrate)
epsilon = 1.0E-10


sys.stderr.write("Avg seqlen should be {} (if avg indel rate is 1)\n".format(seqlen))
sys.stderr.write("Insertion rate is {}\n".format(insrate))
sys.stderr.write("Deletion rate is {}\n".format(delrate))
sys.stderr.write("Pinvar is {}\n".format(pinvar))
sys.stderr.write("Short branch length is {} \n".format(sb))
sys.stderr.write("Long branchlength is {}\n".format(lb))



subst = {
   'A' : 'GCT',
   'G' : 'ACT',
   'C' : 'AGT',
   'T' : 'AGC'
}

## Notes make debugging optional,
## Write some tests..


def gen_column(n):
    for index in xrange(n):
        column_set = simulate_columns_from_A()
        for col in column_set:
            yield col

evolve_along_tree = True
def simulate_columns_from_A():
  if RNG.random() < pinvar:
     return [RNG.choice('AGCT') * 6]
    #initial 
  subseq = []
  u = RNG.random()
  u -= probzlen
  currprob = probzlen
  while u > 0:
    subcol = [RNG.choice('AGCT') ]+['-']*9
    subseq.append(subcol)
    currprob *= geomprobins
    u -= currprob 

  if evolve_along_tree:
    #  sys.stderr.write("A starting length is {}\n".format(len(subseq)))
    rem_lb = lb - sb
    # A to internal I1
    subseq = branch_sim(subseq, 0, 6, sb/one_minus_pinv)
    # I1 to A1
    subseq = branch_sim(subseq, 6, 1, sb/one_minus_pinv)
    # I1 to i
    subseq = branch_sim(subseq, 6, 7, rem_lb/one_minus_pinv)
    # i to B
    subseq = branch_sim(subseq, 7, 2, sb/one_minus_pinv)
    # i to j
    subseq = branch_sim(subseq, 7, 8, sb/one_minus_pinv)
    # j to D
    subseq = branch_sim(subseq, 8, 5, sb/one_minus_pinv)
    # j to I2
    subseq = branch_sim(subseq, 8, 9, rem_lb/one_minus_pinv)
    # i4 to C
    subseq = branch_sim(subseq, 9, 3, sb/one_minus_pinv)
    # i4 to C
    subseq = branch_sim(subseq, 9, 4, sb/one_minus_pinv)
  c = [''.join(lis[:6]) for lis in subseq]
  return c

VERBOSE = False
#simulate changes on all branches from root 'a'
def branch_sim(subseq, start, end, blen):
    tracking = False
    if VERBOSE:
      if (len(subseq)==1):
        tracking = True
        sys.stderr.write("tracking one base\n")
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
        pdel = (len(nondel)*delrate)/eventrate
        u = RNG.random()
        if u < pimmlink:
            newcol=['-']*10
            newcol[end]=RNG.choice('AGCT')
            subseq.insert(0,newcol)
            nondel = [0] + [ite + 1 for ite in nondel]
            eventrate += indelrate
            if tracking:
                sys.stderr.write("u is {}, insertion at immortal link, {} bases\n".format(u, len(nondel)))
        else:
            ndi = RNG.randrange(0, len(nondel))
            u -= pimmlink
            if tracking:
                sys.stderr.write("pdel is now {}\n".format(pdel))
            if u < pdel:
               ci = nondel.pop(ndi)
               subseq[ci][end]='-'
               eventrate -= indelrate
               if tracking:
                  sys.stderr.write("u is {}, deletion, {} bases left\n".format(u, len(nondel)))
            else:
                newcol=['-']*10
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
                if tracking:
                  sys.stderr.write("u is {}, insertion! now {} bases\n".format(u, len(nondel)))
    subprob = 0.75 - 0.75*math.exp((-4.0*blen)/3.0)
    for ci in nondel:
        if RNG.random() < subprob:
            anc = subseq[ci][end]
            desc = RNG.choice(subst[anc])
            subseq[ci][end] = desc
    #remove columns with all gaps
    allgap = ['-']*10
    subseq = [ i for i in subseq if i != allgap ]
    return subseq
    

states=('-','N')
s = [states]*6
counts = {''.join(stat):0 for stat in itertools.product(*s)}


def gap_or_no(a):
   if a in 'ATGC':
      return('N')
   elif a == '-':
      return('-')
   else:
       return(a)

def states_trans(col):
   counts["".join([gap_or_no(i) for i in col])] += 1
   

def sim_mat(out, seqlen):
    n = 0
    for column in gen_column(seqlen):
        states_trans(column)
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
    dimensions ntax = 6 nchar = {c} ;
    taxlabels A A1 B C C1 D ;
    format transpose datatype=dna gap = '-';
matrix 
{m}
;
end;
'''.format(c=nc, m='\n'.join(numbered)))


sortke = list(counts.keys())
sortke.sort()
for item in sortke:
    sys.stderr.write("Count of {it} was {count} in sequence of final length {c}\n".format(it=item,count=counts[item],c=nc))

for ii, tip in enumerate(['A','B','C','D']):
    tip_seq_len = sum([counts[key] for key in counts.keys() if key[ii]=='N'])# removing invarioable sites
    sys.stderr.write("lenth of non-gap seq at tip {} is {}\n".format(tip, tip_seq_len))


