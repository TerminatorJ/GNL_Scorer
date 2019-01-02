# GNL_Scorer
Introduction
-----------------
this is a software for sgRNA activity prediction with great generalization  
The efficiency of sgRNA is the genome editing capabilities to target site. all the sgRNA we designed to target the interesting site in the genome of specific cells,organism,species should cut the upstream sequence of Protospacer Adjacent Motif sequence. However, not all the sgRNA can be activity in such a process. So, the prediction of the efficiency of sgRNA is urgent when applying the knockout(KO) or knockin(KI) experiments in practice. although there are many prediction algorithms developed by many groups, and they can perform well in specific datasets. However, most of the models can not generalized well in the non-seen datasets, which means the models have bad generalization. Here, we advocate using the dataset using the "powerful" detected method to measure the activity of sgRNA, We trained it using the Bayesian Ridge Regression. this model can has a great performance in practice in spite of the cell types,organism and species. Maybe it can be used for the new species without developing the algorithm for the sgRNA efficiency, and also suit for the human, mouse, zebra fish, drosophila, elegans with a general performance.   


Publications:




Usage:

