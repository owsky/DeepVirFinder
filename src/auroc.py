import numpy as np
from keras.models import load_model
from sklearn.metrics import roc_auc_score
import optparse

options = optparse.OptionParser()
options.add_option("-m", "--model", action="store", type="string", dest="model_path")

model_path = options.model_path

model = load_model(model_path)
test_data = []
y_pred = model.predict(test_data)
y_true = [] # true labels for test data
auroc = roc_auc_score(y_true, y_pred)
print("AUROC:", auroc)
