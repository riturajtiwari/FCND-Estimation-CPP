import numpy as np


def load_file(fname):
    result = np.genfromtxt(fname, skip_header=1, delimiter=",")
    return result


base_path = "/Users/rituraj_tiwari/Documents/Udacity/FCND/FCND-Estimation-CPP/config/log/"
files = ["Graph1.txt", "Graph2.txt"]
for file in files:
    file_path = base_path + file
    data = load_file(file_path)
    print("File: %s, std dev: %f" % (file_path, np.std(data[:, 1])))
