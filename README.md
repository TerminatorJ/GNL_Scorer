# GNL_Scorer
Introduction
-----------------
This is a software for sgRNA activity prediction with great generalization  
The efficiency of sgRNA is the genome editing capabilities to target site. All the sgRNA we designed to target the interesting site in the genome of specific cells,organism,species should cut the upstream sequence of Protospacer Adjacent Motif(PAM). However, not all the sgRNA can be activity in such a process. So, prediction of the efficiency of sgRNA is urgent when applying the knockout(KO) or knockin(KI) experiments in practice. Although there are many prediction algorithms developed by many groups, and they can perform well in specific datasets. However, most of the models can not generalized well in the non-seen datasets, which means the models have bad generalization. Here, we advocate using the dataset using the "powerful" detected method to measure the activity of sgRNA, We trained it used the Bayesian Ridge Regression. This model could have a great performance in practice in spite of the cell types,organisms and species. Maybe it can be used for the new species without developing the algorithm for the sgRNA efficiency, and also suit for the human, mouse, zebra fish, drosophila, elegans with a general performance.   


Publications:

（Add the paper）


Usage:

Add the "oligoarrayaux-3.8" in the cal_deltaG file to you linux environment.
```Bash
export PATH=your_path/cal_deltaG/oligoarrayaux-3.8/bin :$PATH
```
Then you can use the hybrid-ss-min directly

Make sure the structure of your direction is like below:  
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

Make sure these packages were installed in your python3:
```Bash
python package: os sys pandas=0.23.4 numpy=1.15.3 time sklearn=0.19.2 Bio=1.72 pickle itertools  
```




An examlple script was shown in the test.py and you can copy it to your own file of your_script.py!

There are two modes for you to use this software:  
1) If your want to input sevearl sequences run like below:
In your own python script write ".py" as :
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
In your shell/bash run your_script.py:
```Bash
your_script.py AACCATGTGACTGTGCATGCTGTACGGCTC ACTCGTGACTGACTAGCTAGGGACTGGCTA
```
Note that: the sequence you input following the ".py" show have the length of 30, and the 26th-27th should be the "GG"

2) If your input is a file of fasta:
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

In your shell/bash run your_script.py:
```Bash
your_script.py your_fastafile.fasta
```


Output:
The 1th result is : 2.4192648021191783
...
...









