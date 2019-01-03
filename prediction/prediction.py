import pickle
import os


def predict(model_type,model_dir,this_input):
    #root_dir=os.path.abspath("..")
    #model_dir=os.path.join(root_dir,"hela_and_hct116_1_datasets_BR_model_80_train.pickle")
    if model_type=="no_epi":
        with open(model_dir,"rb") as f1:
            model=pickle.load(f1)
            prediction=model.predict(this_input)
            return prediction
