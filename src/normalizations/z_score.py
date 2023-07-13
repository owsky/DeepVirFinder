import numpy as np

def encodeSeq(seq: str) -> list: 
    seq_code = list()
    for pos in range(len(seq)) :
        letter = seq[pos]
        if letter in ['A', 'a'] :
            code = [1,0,0,0]
        elif letter in ['C', 'c'] :
            code = [0,1,0,0]
        elif letter in ['G', 'g'] :
            code = [0,0,1,0]
        elif letter in ['T', 't'] :
            code = [0,0,0,1]
        else :
            code = [1/4, 1/4, 1/4, 1/4]
        seq_code.append(code)
    return seq_code 

def normalize_sequence(sequence: str) -> list:
    encoded = encodeSeq(sequence)
    data_array = np.array(encoded)
    
    # Calculate the mean and standard deviation along each column
    mean = np.mean(data_array, axis=0)
    std_dev = np.std(data_array, axis=0)
    
    # Perform z-score normalization
    normalized_data = (data_array - mean) / std_dev
    
    return normalized_data.tolist()