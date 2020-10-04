README
-------------------------------------------------------------------------------------------------------------------------------------------
This code is for C. Y. Wu and J. J. Ding, ¡§Nonconvex Approach for Sparse and Low-Rank Constrained models with Dual Momentum,"

including 

main.m : As main program, the parameters are explained in the comment, we can get the results in Section VI.B.1.

The folder nonconvex_funs containing different types of nonconvex shrinkage functions

linear_sg.m : Proposed linear piecewise shrinkage.

lrr: Origial LRR from Lib-ADMM Canyi Lu, Zhouchen Lin, Shuicheng Yan "A Unified Alternating Direction Method of Multipliers by Majorization Minimization," accepted by TPAMI

lrr2: convex LRR + dual momentum adapted from Lib-ADMM

lrr3: nonconvex LRR adapted from Lib-ADMM

lrr4: nonconvex LRR + dual momentum adapted from Lib-ADMM

showW2.m: the spectral clustering adapted from Lib-ADMM

AccMeasure.m: calculating the clustering accuracy from Praisan Padungweang