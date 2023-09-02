import math
import sys
from sklearn.metrics import silhouette_score
import numpy as np
from numpy import random
import mysymnmfsp

np.random.seed(0)


def kmeans(K, iter, input_data):
    epsilon = 0.0001
    num_of_clusters = K
    file = input_data
    # Initialize array_of_centroids with the first K vectors in the input file
    array_of_centroids = build_table(num_of_clusters, input_data)[
        0]  # initialize a table with the first k vectors in the input file
    N = build_table(num_of_clusters, input_data)[1]
    # Check if the number of clusters is valid
    if num_of_clusters >= N:
        print("Invalid number of clusters!")
        return
    iter_cnt = 0
    delta = epsilon + 1
    # Perform k-means iterations until convergence or maximum iterations
    while (iter_cnt < iter and delta >= epsilon):
        clusters = [[] for i in range(num_of_clusters)]
        iter_cnt += 1
        # Read the input file and assign vectors to the nearest centroid
        interactive_file = open(file, 'r')
        lines = interactive_file.readlines()
        line_index = 0
        for line in lines:
            line_index += 1
            # For the first iteration, assign the first K vectors as centroids
            if (iter_cnt == 1) and (line_index <= num_of_clusters):
                vector = list(line.split(","))
                vector1 = create_vector(vector)
                clusters[line_index - 1].append(vector1)
                continue
            vector = list(line.split(","))
            vector1 = create_vector(vector)
            d = len(vector1)
            # Find the nearest centroid for the vector
            min_distance = math.inf
            min_distance_index = -1
            for i in range(num_of_clusters):
                curr_distance = calc_distance(vector1, array_of_centroids[i])
                if curr_distance < min_distance:
                    min_distance = curr_distance
                    min_distance_index = i
            clusters[min_distance_index].append(vector1)
        # Update the centroids based on the assigned vectors
        array_of_centroids = update_centroid(clusters, num_of_clusters, d)
        # Calculate the change in centroids and check for convergence
        if iter_cnt != 1:
            delta = calc_delta(array_of_centroids, prev_array_of_centroids)
        prev_array_of_centroids = array_of_centroids.copy()
        line_index = 0
        interactive_file.close()
    # Create output by printing the final centroids
    cluster_assignments = []
    for line in lines:
        vector = list(line.split(","))
        vector1 = create_vector(vector)
        min_distance = math.inf
        min_distance_index = -1
        for i in range(num_of_clusters):
            curr_distance = calc_distance(vector1, array_of_centroids[i])
            if curr_distance < min_distance:
                min_distance = curr_distance
                min_distance_index = i
        cluster_assignments.append(min_distance_index)

    return array_of_centroids, cluster_assignments



# Function to calculate the change in centroids
def calc_delta(new_centroids, prev_centroids):
    delta = 0
    for i in range(len(new_centroids)):
        distance = calc_distance(new_centroids[i], prev_centroids[i])
        if distance > delta:
            delta = distance
    return delta


# Function to build the initial table of centroids
def build_table(num_of_clusters, input_data):
    array_of_centroids = [[] for i in range(num_of_clusters)]
    file = input_data
    interactive_file = open(file, 'r')
    lines = interactive_file.readlines()
    line_index = 0
    for line in lines:
        vector = list(line.split(","))
        vector1 = create_vector(vector)
        if line_index < num_of_clusters:
            array_of_centroids[line_index] = vector1
            line_index += 1
        else:
            break
    N = len(lines)
    interactive_file.close()
    return [array_of_centroids, N]


# Function to create a vector from a list of values
def create_vector(vector):
    vector[-1] = vector[-1][:len(vector[-1]) - 1]
    for j in range(len(vector)):
        vector[j] = float(vector[j])
    return vector


# Function to calculate the Euclidean distance between two vectors
def calc_distance(vector1, vector2):
    i = 0
    summ = 0
    for i in range(len(vector1)):
        summ += pow((vector1[i] - vector2[i]), 2)
    return math.sqrt(summ)


# Function to update the centroids based on the assigned vectors
def update_centroid(clusters, num_of_clusters, d):
    new_centroids = [[] for i in range(num_of_clusters)]
    for i in range(num_of_clusters):
        new_centroids[i] = [0 for i in range(d)]
        num_of_points = len(clusters[i])
        for point in clusters[i]:
            for coord in range(d):
                new_centroids[i][coord] += point[coord] / num_of_points
    return new_centroids


# Function to create the output by printing the final centroids

def avg_matrix(M):
    cols = len(M[0])
    rows = len(M)
    sum = 0
    for i in range(rows):
        for j in range(cols):
            sum += M[i][j]
    return (sum / (cols * rows))


def create_output(vectors_array):
    for vector in vectors_array:
        for i in range(len(vector)):
            vector[i] = '{:.4f}'.format(vector[i])
        print(','.join(vector))


def main():
    lst_of_arguments = sys.argv
    # Check the number of command-line arguments
    if (len(lst_of_arguments) != 3):
        print("An error has occurred")
        return 0
    first_argument = lst_of_arguments[1]
    # Check if the first argument is a valid number of clusters
    if first_argument.isdigit():
        first_argument = int(first_argument)
    else:
        print("An error has occurred")
        return 0
    file = lst_of_arguments[2]
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
    # Check if the third argument exists and validate it
    n = len(list)
    k = first_argument
    H = [[0.0 for i in range(k)] for j in range(n)]
    #   Initializing H   #
    ###########################################
    W = mysymnmfsp.fit(list, H, 4, k)
    m = avg_matrix(W)
    max_val_interval = math.sqrt(m / k) * 2
    H_min = [0 for i in range(k)]
    H_max = [max_val_interval for i in range(k)]
    H = np.random.uniform(low=H_min, high=H_max, size=(n, k))
    H = H.tolist()
    finalH = mysymnmfsp.fit(list, H, 1, k)

    finalH_np = np.array(finalH)
    cluster_array_symnmf = np.argmax(finalH, axis=1)
    kmeans_table, cluster_array_kmeans = kmeans(k, 300, file)
    symnmf_score = silhouette_score(list, cluster_array_symnmf)
    kmeans_score = silhouette_score(list, cluster_array_kmeans)
    formatted_string_symnmf = "{:.4f}".format(symnmf_score)
    formatted_string_kmeans = "{:.4f}".format(kmeans_score)
    print("nmf:", formatted_string_symnmf)
    print("kmeans:", formatted_string_kmeans)


if __name__ == "__main__":
    main()
