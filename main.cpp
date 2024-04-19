#include <fstream>
#include "AgglomerativeHierarchicalCluster.h"

using namespace Eigen;
using namespace std;

/**
 * @brief The main function of the program.
 * 
 * This function reads data from a file, performs agglomerative hierarchical clustering on the data,
 * and prints the resulting clusters.
 * 
 * @return 0 if the program executed successfully, 1 otherwise.
 */
int main() {
    // Define the dimensions of the matrix
    int rows = 2;
    int cols = 100;

    // Create an Eigen::MatrixXd to hold the data
    MatrixXd data(rows, cols);

    // Open the input file
    ifstream inputFile("../../sample_input.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open input file!" << endl;
        return 1;
    }

    // Read data from the file into the Eigen::MatrixXd
    for (int col = 0; col < cols; ++col) {
        for (int row = 0; row < rows; ++row) {
            if (!(inputFile >> data(row, col))) {
                cerr << "Error: Failed to read data from input file!" << endl;
                return 1;
            }
        }
    }

    // Close the input file
    inputFile.close();


    int K = 4; // Number of clusters desired
    int M = 0; // Use single-linkage distance measure

    // Create an instance of AgglomerativeHierarchicalClustering
    AgglomerativeHierarchicalClustering clustering(data, K, M);

    // Run the clustering algorithm
    clustering.run_algorithm();

    // Print the resulting clusters
    clustering.print();
    clustering.print_label();

    //print cluster index for each data point
    auto cluster_idx = clustering.cluster_idx();

    for (int i = 0; i < cluster_idx.size(); ++i) {
        cout << "Data point [" << data(0,i) << "," << data(1,i) << "] is in cluster " << cluster_idx[i] << endl;
    }

    return 0;
}
