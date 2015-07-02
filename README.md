gapsim
======
Simulation of the paired-invariants model using JC+TFK91+Invar

prereqs for analysis:

  1. clone https://github.com/mtholder/DendroBites and put the dendrobites subdir on your PATH

  2. clone https://github.com/jeetsukumaran/DendroPy checkout and install the charset-export branch

  3. get paup on your PATH

## running

    bash do_sims.sh 10

should run 10 replicates creating `rep1`, `rep2`, `rep3`... `rep10` subdirectories.

Parsing the `sim.log` in each should show you:

  1. whether or not JC+pinv model on the simulated data gets the correct tree (or the LBA tree), and

  2. whether or not JC model on the culled data gets the correct tree (or the LBA tree)

