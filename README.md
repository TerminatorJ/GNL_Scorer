# GNL_Scorer
Introduction
-----------------
This software is built for sgRNA activities prediction.
The efficiency of sgRNAs is the genome editing capabilities to target sites. All the sgRNAs we design to target the potential sites in the genome of specific cells, organisms, species should cut the upstream sequences of Protospacer Adjacent Motif (PAM). However, not all of the sgRNAs can be imported in such a process. Therefore, the prediction of the efficiency of sgRNAs is urgent when applying the knockout (KO) or knockin (KI) experiments in practice. Although, there are quite a few predictive algorithms developed by some groups, but they can just perform well in specific datasets or specific species. Most of the models don't have considerable generalization in the non-seen datasets and even different species data. Here, we suggest using the model trained by "powerful" dataset, which is measured by sequencing. We trained the model by Bayesian Ridge Regression. which could have a great performance in practice in spite of the cell types, organisms and species. Maybe it can be used for those species without developed sgRNA predictive algorithms, and suit for species used most frequently, like human, mouse, zebra fish, drosophila, elegans.
 


Publications:

Wang J, Xiang X, Bolund L, et al. GNL-Scorer: A generalized model for predicting CRISPR on-target activity by machine learning and featurization[J]. Journal of Molecular Cell Biology, 2020.

# Note:
This tool has been refreshed into a new version, please get accesse to https://github.com/TerminatorJ/CRISPR-TRAP-seq to get a solid result to facilitate your experiment.

Usage:

Adding the "oligoarrayaux-3.8" in the file of cal_deltaG to you linux environment.
```Bash
export PATH=your_path/cal_deltaG/oligoarrayaux-3.8/bin :$PATH
```
Then you can use the programme of hybrid-ss-min directly
Make sure the file structure under work direction is the same as below:
 
```Bash
your_script.py  
cal_deltaG  
        oligoarrayaux-3.8  
                bin  
                ...  
                ...  
        deltaG_result.txt.dG  
        for_deltaG.txt  
example  
        example_30mer.fasta  
        example_30mer.txt  
featurization  
        __pycache__  
        featurization.log  
        featurization.py  
model  
        hela_BRR.pickle  
        hela_and_hct116_BRR.pickle  
prediction  
        prediction.py  
readme.md  
```

Make sure these packages were installed in your environment:
```Bash
python3
python3 package: os sys pandas=0.23.4 numpy=1.15.3 time sklearn=0.19.2 Bio=1.72 pickle itertools  
```
An example script can be seen in the test.py and you can follow it step by step!
There are two modes you can use in this software:
1.If you want to input several sequences, you can run the following codes: In your own python script, you can write ".py" file as:

 ```python
from featurization import featurization
import sys
import os
from prediction import prediction
#dir processing
root_dir=os.path.abspath(".")
cal_deltaG_dir=os.path.join(root_dir,"cal_deltaG")
list_input=[]
for argv in sys.argv[1:]:
   assert len(argv)==30, "the sequence input must be 30mer"
   list_input.append(argv)
deltaG_result_dir=featurization.cal_deltaG(list_input,cal_deltaG_dir,input_seq=True)
this_input=featurization.get_input(list_input,deltaG_result_dir,input_seq=True)
##if your sequence is from human run the script below
#model_dir=os.path.join(root_dir,"model","hela_and_hct116_BRR.pickle")

##if your sequence is from other species run the script below
#model_dir=os.path.join(root_dir,"model","hela_BRR.pickle")

predict_result=prediction.predict("no_epi",model_dir,this_input)
for item,result in enumerate(predict_result):
    print("the %dth result is : %s\n" % (item,result))
```
In your shell/bash, you can run your_script.py:
```Bash
your_script.py AACCATGTGACTGTGCATGCTGTACGGCTC ACTCGTGACTGACTAGCTAGGGACTGGCTA
```
Note that: the sequences you input following as the ".py" have to be as long as 30nt, and the 26th-27th bases should be "GG".

2) If your input sequences are set in a fasta file:
```python
from featurization import featurization
import sys
import os
from prediction import prediction
#dir processing
root_dir=os.path.abspath(".")
cal_deltaG_dir=os.path.join(root_dir,"cal_deltaG")
fasta_dir=sys.argv[1]
with open(fasta_dir,"r") as f1:
    fasta_file=f1.read()
deltaG_result_dir=featurization.cal_deltaG(fasta_file,cal_deltaG_dir,input_seq=False)
print(deltaG_result_dir)
this_input=featurization.get_input(fasta_file,deltaG_result_dir)
##if your sequence is from human run the script below
#model_dir=os.path.join(root_dir,"model","hela_and_hct116_BRR.pickle")

##if your sequence is from other species run the script below
#model_dir=os.path.join(root_dir,"model","hela_BRR.pickle")

predict_result=prediction.predict("no_epi",model_dir,this_input)
for item,result in enumerate(predict_result):
    print("The %dth result is : %s\n" % (item+1,result))
```

In your shell/bash env, you can run your_script.py:
```Bash
your_script.py your_fastafile.fasta
```


Output:
The 1th result is : 2.4192648021191783
...
...









