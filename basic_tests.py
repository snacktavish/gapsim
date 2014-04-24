#! /usr/bin/env python
import unittest
import subprocess

#Things I need to test
#Am I geting the expected sequence length in the absence of invariant sites
args=["python","gap_sim.py", "-invar", "0"]
testruns = {}
for i in range(15):
    output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
    testruns[i] = {'output':output.split('\n'), 'error':error.split('\n')}

assert(testruns[i]['error'][35].startswith("lenth of non-gap seq at tip A"))
lenA = [int(testruns[i]['error'][35].split()[-1]) for i in testruns]
assert( 90 < sum(lenA)/len(lenA) < 110)

site_patterns={}
for i in testruns:
   for lin in testruns[i]['error'][19:35]:
       ki = lin.split()[2]
       val = int(lin.split()[4])
       if ki not in site_patterns:
            site_patterns[ki] = [val]
       else:
            site_patterns[ki].append(val)
            
assert(sum(site_patterns['--N-'])*0.9 < sum(site_patterns['N---']) < sum(site_patterns['--N-'])*1.1)
assert(sum(site_patterns['-N--'])*0.9 < sum(site_patterns['---N']) < sum(site_patterns['-N--'])*1.1)



#Are the branch lengths appropriate? i.e. is the simulation working correctly

#Under what conditions do we get the wrong tree?

#

#class TestGapsim(unittest.TestCase):
#    def testSimple(self):
        
#if __name__ == "__main__":
#    unittest.main()
    
