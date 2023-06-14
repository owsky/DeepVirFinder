import numpy as np
from keras.models import load_model
from sklearn.metrics import roc_auc_score
import optparse
import os

def main():
    parser = optparse.OptionParser()
    parser.add_option("-m", "--model", action="store", type="string", dest="model_path")
    parser.add_option("-t", "--test", action="store", type="string", dest="test_path")

    (options, _) = parser.parse_args()

    model_path = options.model_path
    model_dir = os.path.dirname(model_path)
    model_name = os.path.basename(model_path).split(".")[0]
    model = load_model(model_path)

    test_path = options.test_path

    # Loading test data
    lengths = ["0.15k", "0.3k", "0.5k", "1.0k"]

    for l in lengths:
        host = np.load(os.path.join(test_path, f"host_{l}_codefw.npy"))
        hostbw = np.load(os.path.join(test_path, f"host_{l}_codebw.npy"))
        virus = np.load(os.path.join(test_path, f"virus_{l}_codefw.npy"))
        virusbw = np.load(os.path.join(test_path, f"virus_{l}_codebw.npy"))

        ### validation V+B
        y_true = np.concatenate((np.repeat(0, host.shape[0]), np.repeat(1, virus.shape[0])))

        X_valfw = np.concatenate((host, virus), axis=0)
        del host, virus

        X_valbw = np.concatenate((hostbw, virusbw), axis=0)
        del hostbw, virusbw

        y_pred = model.predict([X_valfw, X_valbw], batch_size=1)
        auroc = roc_auc_score(y_true, y_pred)
        
        res_path = os.path.join(model_dir, model_name + f"_{l}_auroc.txt")
        with open(res_path, "w") as file:
            file.write(f"{auroc}\n")
        file.close()

if __name__ == "__main__":
    main()