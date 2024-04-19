#include <iostream>
#include <map>
#include <vector>
#include "Measures.h"

/**
 * @class AgglomerativeHierarchicalClustering
 * @brief Represents the Agglomerative Hierarchical Clustering algorithm.
 *
 * The AgglomerativeHierarchicalClustering class performs agglomerative hierarchical clustering on a given input data matrix.
 * It merges clusters iteratively until a stopping criterion is met, resulting in a set of clusters.
 * The class provides methods to initialize clusters, find the closest clusters, merge clusters, and run the clustering algorithm.
 * It also provides methods to print the resulting clusters and the labels of the data points.
 */
class AgglomerativeHierarchicalClustering {
public:
    /**
     * @brief Constructs an instance of the AgglomerativeHierarchicalClustering class.
     *
     * @param data The input data matrix.
     * @param K The number of clusters to be formed.
     * @param M The linkage method to be used for clustering.
     */
    AgglomerativeHierarchicalClustering(const Eigen::MatrixXd& data, int K, int M)
        : _data(data), _N(data.cols()), _K(K) {
        _measure = get_distance_measure(M);
        init_clusters();
    }

    void init_clusters() {
        Eigen::VectorXi index(1);
        for (int i = 0; i < _data.cols(); ++i) {
            _clusters[i] = _data.col(i);
            index[0] = i;
            _label[i] = index;
        }
    }

    /**
     * Finds the closest clusters in the Agglomerative Hierarchical Cluster algorithm.
     *
     * @return A pair of integers representing the indices of the closest clusters.
     */
    std::pair<int, int> find_closest_clusters() {
        double min_dist = std::numeric_limits<double>::infinity();
        std::pair<int, int> closest_clusters;

        auto cluster_ids = get_cluster_ids();

        for (size_t i = 0; i < cluster_ids.size() - 1; ++i) {
            for (size_t j = i + 1; j < cluster_ids.size(); ++j) {
                int ci_id = cluster_ids[i];
                int cj_id = cluster_ids[j];
                double dist = _measure(_clusters[ci_id], _clusters[cj_id]);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_clusters = {ci_id, cj_id};
                }
            }
        }
        return closest_clusters;
    }

    /**
     * Merges two clusters identified by their IDs and forms a new cluster.
     *
     * @param ci_id The ID of the first cluster to merge.
     * @param cj_id The ID of the second cluster to merge.
     */
    void merge_and_form_new_clusters(int ci_id, int cj_id) {
        Eigen::MatrixXd merged_cluster(_data.rows(), _clusters[ci_id].cols() + _clusters[cj_id].cols());
        merged_cluster << _clusters[ci_id], _clusters[cj_id];
        Eigen::VectorXi merged_label(_label[ci_id].size() + _label[cj_id].size());
        merged_label << _label[ci_id], _label[cj_id];

        _clusters.erase(cj_id);
        _clusters[ci_id] = merged_cluster;
        _label.erase(cj_id);
        _label[ci_id] = merged_label;
    }

    /**
     * Runs the agglomerative hierarchical clustering algorithm.
     * This function performs the steps of the agglomerative hierarchical clustering algorithm.
     * It merges the clusters iteratively until a stopping criterion is met.
     * The resulting clusters are stored in the data structure of the class.
     */
    void run_algorithm() {
        int loop_count = 0;
        while (_clusters.size() > static_cast<size_t>(_K)) {
            auto closest_clusters = find_closest_clusters();
            merge_and_form_new_clusters(closest_clusters.first, closest_clusters.second);
            loop_count++;
        }
    }

    Eigen::VectorXi cluster_idx() const {
        Eigen::VectorXi res(_N);
        int cluster_id = 0;
        for (const auto& [_, idx] : _label) {
            for (int i = 0; i < idx.rows(); ++i) {
                res[idx(i)] = cluster_id;
            }
            cluster_id++;
        }

        return res;
    }

    /**
     * Prints the Cluster object.
     */
    void print() const {
        int cluster_id = 0;
        for (const auto& [id, points] : _clusters) {
            std::cout << "Cluster: " << cluster_id++ << ", Cluster size: " << points.cols() <<std::endl;
            for (int i = 0; i < points.cols(); ++i) {
                std::cout << "    ";
                for (int j = 0; j < points.rows(); ++j) {
                    std::cout << points(j, i) << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    /**
     * Prints the label of the object.
     */
    void print_label() const {
        int label_id = 0;
        for (const auto& [_, idx] : _label) {
            std::cout << "Cluster: " << label_id++ << ", Cluster size: " << idx.rows() <<std::endl;
            std::cout << "    ";
            for (int i = 0; i < idx.rows(); ++i) {
                std::cout << idx(i) << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    /**
     * Returns the cluster IDs of the agglomerative hierarchical cluster.
     *
     * @return A vector of integers representing the cluster IDs.
     */
    std::vector<int> get_cluster_ids() const {
        std::vector<int> ids;
        for (const auto& [id, _] : _clusters) {
            ids.push_back(id);
        }
        return ids;
    }

private:
    Eigen::MatrixXd _data;
    int _N;
    int _K;
    DistanceFunction _measure;
    std::map<int, Eigen::MatrixXd> _clusters;
    std::map<int, Eigen::VectorXi> _label;
};
