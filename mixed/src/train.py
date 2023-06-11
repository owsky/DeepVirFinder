import os
import sys
from sklearn.metrics import roc_auc_score
import numpy as np
from keras.models import load_model, Model
from keras.layers import Dense, Dropout, Input, Conv1D, GlobalMaxPooling1D, BatchNormalization
from keras.layers.merge import Average
from keras.optimizers import Adam
from keras.callbacks import ModelCheckpoint, EarlyStopping
from custom_layers import RowNormalization, ColumnNormalization

base_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),"..")
out_dir = os.path.join(base_path, "models")
os.makedirs(out_dir, exist_ok=True)
training_encodings = os.path.join(base_path, "data", "tr", "encoded")
validation_encodings = os.path.join(base_path, "data", "val", "encoded")

######## loading data for training, validation, testing ##########
print("...loading data...")

## phage RefSeq
print("...loading virus data...")
# training
filename_codetrfw = [ x for x in os.listdir(training_encodings) if 'fw.npy' in x and 'virus' in x  ][0]
print("data for training " + filename_codetrfw)
phageRef_codetrfw = np.load(os.path.join(training_encodings, filename_codetrfw), allow_pickle=True)
phageRef_codetrbw = np.load(os.path.join(training_encodings, filename_codetrfw.replace('fw', 'bw')), allow_pickle=True)
# validation
filename_codevalfw = [ x for x in os.listdir(validation_encodings) if 'fw.npy' in x and 'virus' in x  ][0]
print("data for validation " + filename_codevalfw)
phageRef_codevalfw = np.load(os.path.join(validation_encodings, filename_codevalfw), allow_pickle=True)
phageRef_codevalbw = np.load(os.path.join(validation_encodings, filename_codevalfw.replace('fw', 'bw')), allow_pickle=True)


## host RefSeq
print("...loading host data...")
# training
filename_codetrfw = [ x for x in os.listdir(training_encodings) if 'fw.npy' in x and 'host' in x  ][0]
print("data for training " + filename_codetrfw)
hostRef_codetrfw = np.load(os.path.join(training_encodings, filename_codetrfw), allow_pickle=True)
hostRef_codetrbw = np.load(os.path.join(training_encodings, filename_codetrfw.replace('fw', 'bw')), allow_pickle=True)
# validation
filename_codevalfw = [ x for x in os.listdir(validation_encodings) if 'fw.npy' in x and 'host' in x  ][0]
print("data for validation " + filename_codevalfw)
hostRef_codevalfw = np.load(os.path.join(validation_encodings, filename_codevalfw), allow_pickle=True)
hostRef_codevalbw = np.load(os.path.join(validation_encodings, filename_codevalfw.replace('fw', 'bw')), allow_pickle=True)

######## combine V and H, shuf training data ##########
print("...combining V and H...")
### training V+B
Y_tr = np.concatenate((np.repeat(0, hostRef_codetrfw.shape[0]), np.repeat(1, phageRef_codetrfw.shape[0])))
X_trfw = np.concatenate((hostRef_codetrfw, phageRef_codetrfw), axis=0)
del hostRef_codetrfw, phageRef_codetrfw
X_trbw = np.concatenate((hostRef_codetrbw, phageRef_codetrbw), axis=0)
del hostRef_codetrbw, phageRef_codetrbw
print("...shuffling training data...")
#size, seq_len, channel_num = X_tr.shape
index_trfw = list(range(0, X_trfw.shape[0]))
np.random.shuffle(index_trfw)
X_trfw_shuf = X_trfw[np.ix_(index_trfw, range(X_trfw.shape[1]), range(X_trfw.shape[2]))]
del X_trfw
X_trbw_shuf = X_trbw[np.ix_(index_trfw, range(X_trbw.shape[1]), range(X_trbw.shape[2]))]
del X_trbw
Y_tr_shuf = Y_tr[index_trfw]

### validation V+B
Y_val = np.concatenate((np.repeat(0, hostRef_codevalfw.shape[0]), np.repeat(1, phageRef_codevalfw.shape[0])))
X_valfw = np.concatenate((hostRef_codevalfw, phageRef_codevalfw), axis=0)
del hostRef_codevalfw, phageRef_codevalfw
X_valbw = np.concatenate((hostRef_codevalbw, phageRef_codevalbw), axis=0)
del hostRef_codevalbw, phageRef_codevalbw

######### training model #############
# parameters
POOL_FACTOR = 1
dropout_cnn = 0.1
dropout_pool = 0.1
dropout_dense = 0.1
learningrate = 0.001
batch_size=int(X_trfw_shuf.shape[0]/(1000*1000/1000)) ## smaller batch size can reduce memory
pool_len1 = int((1000-500+1)/POOL_FACTOR)

channel_num = 4
epochs = 10
filter_len1 = 10
nb_filter1 = 500
nb_dense = 500

modPattern = 'model_siamese_varlen_'+'k_fl'+str(filter_len1)+'_fn'+str(nb_filter1)+'_dn'+str(nb_dense)
modName = os.path.join( out_dir, modPattern + '.h5')
checkpointer = ModelCheckpoint(filepath=modName, verbose=1,save_best_only=True)
earlystopper = EarlyStopping(monitor='val_acc', min_delta=0.0001, patience=5, verbose=1)


##### build model #####

def get_output(input_layer, hidden_layers):
    output = input_layer
    for hidden_layer in hidden_layers:
        output = hidden_layer(output)
    return output


print("...building model...")
## if model exists
if os.path.isfile(modName):
  model = load_model(modName) # , {norm_layer.__name__: norm_layer})
  print("...model exists...")
else :
  ## siamese
  forward_input = Input(shape=(None, channel_num))
  reverse_input = Input(shape=(None, channel_num))
  hidden_layers = [
    Conv1D(filters = nb_filter1, kernel_size = filter_len1, activation='relu'),
    GlobalMaxPooling1D(),
    Dropout(dropout_pool),
    Dense(nb_dense, activation='relu'),
    Dropout(dropout_dense),
    Dense(1, activation='sigmoid')
  ]
  forward_output = get_output(forward_input, hidden_layers)     
  reverse_output = get_output(reverse_input, hidden_layers)
  output = Average()([forward_output, reverse_output])
  model = Model(inputs=[forward_input, reverse_input], outputs=output)
  model.compile(Adam(lr=learningrate), 'binary_crossentropy', metrics=['accuracy'])

print("...fitting model...")
print(str(filter_len1)+'_fn'+str(nb_filter1)+'_dn'+str(nb_dense)+'_ep'+str(epochs))
model.fit(
   x = [X_trfw_shuf, X_trbw_shuf], y = Y_tr_shuf, 
    batch_size=batch_size, epochs=epochs, verbose=2, 
    validation_data=([X_valfw, X_valbw], Y_val), 
    callbacks=[checkpointer, earlystopper]
)

## Final evaluation AUC ###

## train data
type = 'tr'
X_fw = X_trfw_shuf
X_bw = X_trbw_shuf
Y = Y_tr_shuf
print("...predicting "+type+"...\n")
Y_pred = model.predict([X_fw, X_bw], batch_size=1)
auc_tr = roc_auc_score(Y, Y_pred)
del Y, X_fw, X_bw


# val data
type = 'val'
X_fw = X_valfw
X_bw = X_valbw
Y = Y_val
print("...predicting "+type+"...\n")
Y_pred = model.predict([X_fw, X_bw], batch_size=1)
auc_val = roc_auc_score(Y, Y_pred)

with open(os.path.join(out_dir, modPattern + '_' + type + 'auroc.txt'), "w") as f:
    f.write(f"Training AUROC:\t{str(auc_tr)}\nValidation AUROC:\t{str(auc_val)}")
    f.close()

del Y, X_fw, X_bw