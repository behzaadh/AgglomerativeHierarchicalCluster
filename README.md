# Agglomerative Hierarchical Cluster (AHC)
This repository contains an implementation of the Agglomerative Hierarchical Clustering algorithm in C++. The algorithm merges clusters iteratively until a stopping criterion is met, resulting in a set of clusters. 
Got inspired by [this project](https://github.com/OlaPietka/Agglomerative-Hierarchical-Clustering-from-scratch/tree/main).

# Algorithm
During the clustering process, we iteratively aggregate the most similar two clusters, until there are 
 clusters left. For initialization, each data point forms its own cluster.

 # Requirements
 [Eigen library](https://eigen.tuxfamily.org/) is used for matrix operations.

 # Running the Code
1. Make sure to install Eigen library.
2. Compile and run the main.cpp file to perform agglomerative hierarchical clustering on your data.


Feel free to explore and modify the code as needed for your projects or research. If you have any questions or suggestions, please feel free to reach out (bdhosseinzadeh@gmail.com)!