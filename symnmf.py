import math
import sys
import pandas as pd
import numpy as np
from numpy import random
import mysymnmfsp
np.random.seed(0)


# Function to average the matrix entries
def avg_matrix(M):
    cols = len(M[0])
    rows = len(M)
    sum = 0
    for i in range (rows):
        for j in range (cols):
            sum += M[i][j]
    print (sum/(cols*rows))
    return (sum/(cols*rows))

# Function to create the output by printing the final centroids
def create_output(vectors_array):
    for vector in vectors_array:
        for i in range(len(vector)):
            vector[i] = '{:.4f}'.format(vector[i])
        print(','.join(vector))


def main():
    dataframe = None
    data = None
    state = 0
    lst_of_arguments = sys.argv
    # Check the number of command-line arguments
    if (len(lst_of_arguments) != 4 ):
        print("An Error Has Occurred")
        return 1
    first_argument = lst_of_arguments[1]
    # Check if the first argument is a valid number of clusters
    if first_argument.isdigit():
        first_argument = int (first_argument)
    else:
        print("An Error Has Occurred")
        return 1
    second_argument = lst_of_arguments[2]
    # Check if the third argument exists and validate it
    if (second_argument == "symnmf"):
        state = 1
    elif (second_argument == "sym"):
        state = 2
    elif (second_argument == "ddg"):
        state = 3
    elif (second_argument == "norm"):
        state = 4
    else:
        print("An Error Has Occurred")
        return 1
    file = lst_of_arguments[3]
    if (file[-4:] == ".txt"):
        try:
            data = np.genfromtxt(file, dtype=float, encoding=None, delimiter=",")
        except:
            print("An Error Has Occurred")
            return 1
    else:
        print("An Error Has Occurred")
        return 1
    list = data.tolist()
    n = len(list)
    k = first_argument
    H = [[0.0 for i in range(k)] for j in range(n)]
    #   Initializing H   #
    ###########################################
    if (state == 1):
        W = mysymnmfsp.fit(list,H,4,k)
        m = avg_matrix(W)
        max_val_interval = math.sqrt(m / k) * 2
        H_min = [0 for i in range(k)]
        H_max = [max_val_interval for i in range(k)]
        H = np.random.uniform(low=H_min, high=H_max, size=(n, k))
        H = H.tolist()
        X = mysymnmfsp.fit(list,H,1,k)
        create_output(X)
        return 0
    else:
        X = mysymnmfsp.fit(list, H, state, k)
        create_output(X)
        return 0
    

    ###########################################







if __name__ == "__main__":
    main()
