import math
import sys
import pandas as pd
import numpy as np
from numpy import random
np.random.seed(0)

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
        print("An error has occurred")
        return 1
    first_argument = lst_of_arguments[1]
    # Check if the first argument is a valid number of clusters
    if first_argument.isdigit():
        first_argument = int (first_argument)
    else:
        print("Invalid number of clusters!")
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
        print("Invalid goal name!")
        return 1
    file = lst_of_arguments[3]
    if (file[-4:] == ".txt"):
        try:
            data = np.genfromtxt(file, dtype=float, encoding=None, delimiter=",")
            #len_of_text = len(data[0])
            #cols = [i for i in range(len_of_text)]
            #dataframe = pd.DataFrame(data, columns=cols)
        except:
            print("Invalid file's name!")
            return 1
    else:
        print("Invalid file's name!")
        return 1
    list = data.tolist()







if __name__ == "__main__":
    main()
