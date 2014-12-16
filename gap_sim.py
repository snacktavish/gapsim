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
#   A           C
#    \         /
#     \       /
#      \i___j/
#      /     \ 
#     B       D

''')

parser = argparse.ArgumentParser(description='Process some arguments.')
parser.add_argument('-l','--seqlen', help='average sequence length (if avg column indel = 1)', default = 100, type=int)
parser.add_argument('-ins','--insrate', help='insertion rate', default = 0.5, type=float)
parser.add_argument('-del','--delrate', help='deletion rate',  default = 1.0, type=float)
parser.add_argument('-invar','--pinvar', help='proportion of invariant sites',  default = 0.5, type=float)
parser.add_argument('-ta','--ta', help='branch length at a',  default = 0.6, type=float)
parser.add_argument('-tb','--tb', help='branch length at b',  default = 0.01, type=float)
parser.add_argument('-tc','--tc', help='branch length at a',  default = 0.6, type=float)
parser.add_argument('-td','--td', help='branch length at b',  default = 0.01, type=float)
parser.add_argument('-ti','--ti', help='internal branch length',  default = 0.1, type=float)

args = vars(parser.parse_args())


# Sequence length is length 
seqlen = args['seqlen']
insrate = args['insrate']
delrate = args['delrate'] 
pinvar = args['pinvar']
ta = args['ta']
tb = args['tb']
ti = args['ti']
tc = args['tc']
td = args['td']
one_minus_pinv = 1.0 - pinvar
probzlen = (1 - (insrate/delrate))
geomprobins = (insrate/delrate)
epsilon = 1.0E-10


sys.stderr.write("Avg seqlen should be {} (if avg indel rate is 1)\n".format(seqlen))
sys.stderr.write("Insertion rate is {}\n".format(insrate))
sys.stderr.write("Deletion rate is {}\n".format(delrate))
sys.stderr.write("Pinvar is {}\n".format(pinvar))
sys.stderr.write("Branch A length is {} \n".format(ta))
sys.stderr.write("Branch B length is {}\n".format(tb))
sys.stderr.write("Branch C length is {} \n".format(tc))
sys.stderr.write("Branch D length is {}\n".format(td))
sys.stderr.write("Internal branch length is {}\n".format(ti))



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


def simulate_columns_from_A():
  if RNG.random() < pinvar:
     return [RNG.choice('AGCT') * 4]
    #initial 
  subseq = []
  u = RNG.random()
  u -= probzlen
  currprob = probzlen
  while u > 0:
    subcol = [RNG.choice('AGCT') ]+['-']*5
    subseq.append(subcol)
    currprob *= geomprobins
    u -= currprob 
#  sys.stderr.write("A starting length is {}\n".format(len(subseq)))
  subseq = branch_sim(subseq, 0, 4, ta/one_minus_pinv)
  subseq = branch_sim(subseq, 4, 1, tb/one_minus_pinv)
  subseq = branch_sim(subseq, 4, 5, ti/one_minus_pinv)
  subseq = branch_sim(subseq, 5, 2, tc/one_minus_pinv)
  subseq = branch_sim(subseq, 5, 3, td/one_minus_pinv)
  c = [''.join(lis[:4]) for lis in subseq]
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
            newcol=['-']*6
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
                if tracking:
                  sys.stderr.write("u is {}, insertion! now {} bases\n".format(u, len(nondel)))
    subprob = 0.75 - 0.75*math.exp((-4.0*blen)/3.0)
    for ci in nondel:
        if RNG.random() < subprob:
            anc = subseq[ci][end]
            desc = RNG.choice(subst[anc])
            subseq[ci][end] = desc
    #remove columns with all gaps
    allgap = ['-']*6
    subseq = [ i for i in subseq if i != allgap ]
    return subseq
    

states=('-','N')
s = [states]*4
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
begin _paup;
    set crite = like;
    lset nst = 1 basefreq = eq pinv = est;
    lscore 1;
    describe / brl ;
    alltrees ;
    showtree ;
    describetrees;
    matrixrep brlens=yes file=tmp.nex replace;
end;
'''.format(c=nc, m='\n'.join(numbered)))


sortke = list(counts.keys())
sortke.sort()
for item in sortke:
    sys.stderr.write("Count of {it} was {count} in sequence of final length {c}\n".format(it=item,count=counts[item],c=nc))

for ii, tip in enumerate(['A','B','C','D']):
    tip_seq_len = sum([counts[key] for key in counts.keys() if key[ii]=='N'])# removing invarioable sites
    sys.stderr.write("lenth of non-gap seq at tip {} is {}\n".format(tip, tip_seq_len))


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