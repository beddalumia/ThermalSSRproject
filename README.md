## `Scripts&Data` for our project on thermal entanglement under SSR
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/beddalumia/ThermalSSRproject/main?urlpath=%2Fdoc%2Ftree%2FGenerate_plots.ipynb)

- The `U=0/` directory contains MATLAB scripts computing the 2-site RDMs
  for noninteracting lattice systems, leveraging on the formalism by 
  Chung, Peschel and Eisler [1,2,3], and then evaluate measures of the 
  full and the superselected fermionic entanglement. The calculations
  are carried in the large-size limit.

- The `U>0/` directory contains ED data for the same systems, in presence
  of a Hubbard interaction. Due to the exponential numerical effort, the
  systems are solved in the small-size limit. The ED calculations are 
  performed exploiting the recently introduced _trie-ranking_, for performance [4]. 
  We include also a [jupyter notebook](./Generate_plots.ipynb) showcasing the data analysis 
  and reproduce the published plots. You can also [run it interactively on binder](https://mybinder.org/v2/gh/beddalumia/ThermalSSRproject/main?urlpath=%2Fdoc%2Ftree%2FGenerate_plots.ipynb).


[1] Ming-Chiang Chung and Ingo Peschel, [Phys. Rev. B 64, 064412 (2001)](https://doi.org/10.1103/PhysRevB.64.064412)    
[2] Ingo Peschel, [J. Phys. A: Math. Gen. 36 L205 (2003)](https://doi.org/10.1088/0305-4470/36/14/101)    
[3] Ingo Peschel and Viktor Eisler, [J. Phys. A: Math. Theor. 42 504003 (2009)](https://doi.org/10.1088/1751-8113/42/50/504003)    
[4] Markus Wallerberger and Karsten Held, [Phys. Rev. Research 4, 033238 (2022)](https://doi.org/10.1103/PhysRevResearch.4.033238)    
