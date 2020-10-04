README
-------------------------------------------------------------------------------------------------------------------------------------------
This code is for C. Y. Wu and J. J. Ding, ¡§Nonconvex Approach for Sparse and Low-Rank Constrained models with Dual Momentum,"

including 

main.m : As main program, the parameters are explained in the comment, we can get the results in Section VI.C.1.

The folder nonconvex_funs containing different types of nonconvex shrinkage functions

linear_sg.m : Proposed linear piecewise shrinkage.

rpca: Origial RPCA from Lib-ADMM Canyi Lu, Zhouchen Lin, Shuicheng Yan "A Unified Alternating Direction Method of Multipliers by Majorization Minimization," accepted by TPAMI

rpca2: convex RPCA + dual momentum adapted from Lib-ADMM

rpca3: nonconvex RPCA adapted from Lib-ADMM

rpca4: nonconvex RPCA + dual momentum adapted from Lib-ADMM
 