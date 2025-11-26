## `Scripts&Data` for our project on thermal entanglement under SSR

- The `U=0/` directory contains MATLAB scripts computing the 2-site RDMs
  for noninteracting lattice systems, leveraging on the formalism by 
  Chung, Peschel and Eisler [1,2,3], and then evaluate measures of the 
  full and the superselected fermionic entanglement. The calculations
  are carried in the large-size limit.

- The `U>0/` directory contains ED data for the same systems, in presence
  of a Hubbard interaction. Due to the exponential numerical effort, the
  systems are solved in the small-size limit. The ED calculations are 
  performed with the ??? code [4,5]. We include also custom Julia scripts
  to showcase the data analysis (numerical computation of reduced density
  matrices and their entanglement) and reproduce the published plots.


[1] Ming-Chiang Chung and Ingo Peschel, [Phys. Rev. B 64, 064412 (2001)](https://doi.org/10.1103/PhysRevB.64.064412)    
[2] Ingo Peschel, [J. Phys. A: Math. Gen. 36 L205 (2003)](https://doi.org/10.1088/0305-4470/36/14/101)    
[3] Ingo Peschel and Viktor Eisler, [J. Phys. A: Math. Theor. 42 504003 (2009)](https://doi.org/10.1088/1751-8113/42/50/504003)    
[4] Markus Wallerberger and Karsten Held, [Phys. Rev. Research 4, 033238 (2022)](https://doi.org/10.1103/PhysRevResearch.4.033238)    
[5] ^I imagine it is related(?). Is the Julia code "published" somewhere?