# GNL_Scorer
Introduction
-----------------
This is a software for sgRNA activity prediction with great generalization  
The efficiency of sgRNA is the genome editing capabilities to target site. All the sgRNA we designed to target the interesting site in the genome of specific cells,organism,species should cut the upstream sequence of Protospacer Adjacent Motif(PAM). However, not all the sgRNA can be activity in such a process. So, prediction of the efficiency of sgRNA is urgent when applying the knockout(KO) or knockin(KI) experiments in practice. Although there are many prediction algorithms developed by many groups, and they can perform well in specific datasets. However, most of the models can not generalized well in the non-seen datasets, which means the models have bad generalization. Here, we advocate using the dataset using the "powerful" detected method to measure the activity of sgRNA, We trained it used the Bayesian Ridge Regression. This model could have a great performance in practice in spite of the cell types,organisms and species. Maybe it can be used for the new species without developing the algorithm for the sgRNA efficiency, and also suit for the human, mouse, zebra fish, drosophila, elegans with a general performance.   


Publications:




Usage:

Add the "oligoarrayaux-3.8" in the cal_deltaG file to you linux environment.
```Bash
export PATH=your_path/cal_deltaG/oligoarrayaux-3.8/bin :$PATH
```
 then run the test.py
 
 if you input file is sevearl sequences run like below:
 in your python script write your_script.py:
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
model_dir=os.path.join(root_dir,"model","hela_and_hct116_1_datasets_BR_model_80_train.pickle")
predict_result=prediction.predict("no_epi",model_dir,this_input)
for item,result in enumerate(predict_result):
    print("the %dth result is : %s\n" % (item,result))
```
in your shell/bash run your_script.py:
```Bash
your_script.py AACCATGTGACTGTGCATGCTGTACGGCTC ACTCGTGACTGACTAGCTAGGGACTGGCTA
```
note that: the sequence you input following the ".py" show have the length of 30, and the 26th-27th should be the "GG"

if your input is a file of fasta:
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
model_dir=os.path.join(root_dir,"model","hela_and_hct116_1_datasets_BR_model_80_train.pickle")
predict_result=prediction.predict("no_epi",model_dir,this_input)
for item,result in enumerate(predict_result):
    print("The %dth result is : %s\n" % (item,result))
```






