# PowerSpectrum21cm
Code for computing the power spectrum of the 21-cm line during the dark ages. As in 1506.04152 (J.B. Munoz, Y.Ali-Haimoud, and M. Kamionkowski), improved to include full-sky calculation at low ell and flat-sky at high ell for efficiency. 

The file Full21cmCells.cpp is a C++ code to compute C_\ell, and C_{\ell,\ell+j}, given a redshift z, and a bandwidth Deltanu.
It includes the possibility of integrating with a function f(k) (given by fkswitch), to add noise as in 1603.01206 (M. Shiraishi, J.B. Munoz, M. Kamionkowski, and A. Raccanelli) (given by noiseswitch), and to remove the low ells to speed up the code (given by lowellswitch).

Output is ell, C_\ell, and C_{\ell, \ell+n}, for n<=4, at some redshift and with some bandwidth.
