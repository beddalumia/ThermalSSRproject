## Two-site entanglement 

Here we include ED data for the Hubbard model at half filling (and $U>0$), for different finite-size geometries.    
This includes rings of 2-8 sites and 2x3 and 3x3 clusters.    
Each geometry has its own folder containing a subfolder `U%d_n_%.1f_mu%.1f_beta%g_tp0.0`, for each set of parameter values another. `U` is the Hubbard interaction strenght, in units of the
nearest neighbor hopping $t$, `n` is the electron density, `mu` is the chemical potential $\mu/t$, `beta` is the inverse temperature $t\beta$ and `tp` is the next-nearest-neighbor hopping,
here set to $t'=t$.     
Each of these directories contains a `.csv` file with the corresponding 2-site reduced density matrix, for nearest neighbors `rho_(1, 0).csv`.     
The matrix basis follows [arXiv:2312.14275](https://arxiv.org/abs/2312.14275)
