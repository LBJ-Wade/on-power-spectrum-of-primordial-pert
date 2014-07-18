on-power-spectrum-of-primordial-pert
===========================================================================

For any form of potential of inflaton, numerically solve the Mukharnov equation, i.e.,  the e.o.m. of the primordial scalar perturbation; and compute the power spectrum of it.


---

If you use the original definition of the z in Mukharnov equation, then 10 times of the time will be spent to obtain the same result, comparing with the one using slow-rolling parameters in the master branch. And there exist the risk of inaccuracy of the nunerical result. If you insist, you may have to change the working precision of the "NDSolve[]" used in the file "modules.m".

For using the one in master branch, the slow-rolling condition shall first be checked.
