from featurization import featurization
import sys
import os 
from prediction import prediction
#dir processing
root_dir=os.path.abspath(".")
cal_deltaG_dir=os.path.join(root_dir,"cal_deltaG")

#print(cal_deltaG_dir)




##if you input the fasta file(if your input is sequence, please note the lines liake we do here )

#fasta_dir=sys.argv[1]
#with open(fasta_dir,"r") as f1:
#    fasta_file=f1.read()
#deltaG_result_dir=featurization.cal_deltaG(fasta_file,cal_deltaG_dir,input_seq=False)
#print(deltaG_result_dir)
#this_input=featurization.get_input(fasta_file,deltaG_result_dir)





##if you just input the sequence to the 
list_input=[]
for argv in sys.argv[1:]:
    assert len(argv)==30, "the sequence input must be 30mer"
    argv=argv.upper()
    list_input.append(argv)
deltaG_result_dir=featurization.cal_deltaG(list_input,cal_deltaG_dir,input_seq=True)
this_input=featurization.get_input(list_input,deltaG_result_dir,input_seq=True)




#load the model
model_dir=os.path.join(root_dir,"model","hela_and_hct116_1_datasets_BR_model_80_train.pickle")


#predict the sequence
predict_result=prediction.predict("no_epi",model_dir,this_input)


#show the result
for item,result in enumerate(predict_result):
   
    print("the %dth result is : %s\n" % (item,result))
