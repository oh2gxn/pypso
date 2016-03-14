# PyPSO
Python program for running Particle Swarm Optimization (PSO) of 
Yagi-Uda antennas simulated with NEC2

The student ham radio station OH2TI had an old and damaged 8-element yagi
antenna for the 6 meter (50 MHz) frequency band. Like so many tech student
and ham projects, also this one got carried away...

My earlier studies in machine learning and artificial intelligence
intrigued me into trying out PSO, which is essentially an intuitive search
heuristic for optimizing anything, where: 
- a design can be parametrized as a D-dimensional real vector, and 
- the cost / goodness of a design can be evaluated on a univariate 
continuous scale.

So why not take the element positions and lengths of a yagi antenna,
evaluate it using NEC2, and see if tweaking any of the lengths would
produce better antenna designs. Additional challenges were added by
considering a dual band yagi for both 6m band and the new 4m band.

## Parts

`yagi.py`: isolates all aspects specific to the antenna designs. It can
operate as a standalone CLI program for turning a simple construction plan
(element dimensions and positions) into a NEC simulation. Just write the 
list of element lengths, positions, and diameters into a CSV text file and 
convert into a NEC file with yagi.py.

Example:
```
$ yagi.py -o 6-el.nec 6-el.csv
$ nec2c -i 6-el.nec -o 6-el.out
$ xnecview 6-el.nec 6-el.out
```

`pypso`: has a graphical user interface for
- setting PSO parameters and optimization criterion
- loading a Yagi design from a file, and evaluating it
- drawing a randomized swarm of slightly different Yagis ("particles")
- running a particle swarm optimization step, searching for better designs
- plotting the evaluations for each of the particles in the process

Note, that just clicking "Randomize" over and over again will effectively
perform a random walk. Clicking "Next" will perform a PSO iteration.

See the file TODO.org for more details on what is still missing and why...

Copyright: Janne Toivola, OH2GXN, 2015

License: [GNU General Public License v3.0](http://www.gnu.org/licenses/gpl-3.0.en.html)
