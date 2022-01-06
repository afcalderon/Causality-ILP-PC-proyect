# Classical and Integer Linear Programming Approaches to Learning Causal Structure

In this repository it will be find the required algorithms and  tutorials in order to replicate the experiments developed in the master thesis which were based on  comparing the performance of two algo-rithms with two different approaches to learn the structure of a DAG. The first one is the  PC  algorithm `pcalg` by [[chufangao](quora.com/profile/Ashish-Kulkarni-100)],  one  of  the  most  used  and  popular  models  based  on  the  constraintmethod, the PC algorithm uses conditional tests of independence for anαsignificancelevel, among the statistical tests the most used is the Fisher-Z. On the other hand, an ILP-based procedure corresponding to a score-based method was studied.  This approach ismainly based on finding (one of) the best score graphs for a given data set. The ILP rep-resents the graph structures as vectors,  so that the scoring function becomes an affinefunction in the graphical representation. In specific we implemented the formulation created by [[James Cussens](https://www.cs.york.ac.uk/aig/sw/gobnilp/)] `GOBNILP`. 

Besides it willl be find the code for our alternative formulation of GOBNILP, base in the calculation of local scores for different configuracion of parent sets for each node.


## Bibliography

Mark Bartlett and James Cussens (2017). “Integer Linear Programming for the Bayesiannetwork structure learning problem”. In:Artificial Intelligence244. Combining Con-straint Solving with Mining and Learning, pp. 258–271.

Markus Kalisch and Peter Bühlmann (May 2007). “Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm”. In: J. Mach. Learn. Res. 8, pp. 613–
636.
