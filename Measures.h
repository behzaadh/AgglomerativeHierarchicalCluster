#include <Eigen/Dense>
#include <limits>

// Define distance function using Eigen::Array
/**
 * Calculates the distance between two arrays.
 *
 * @param p The first array.
 * @param q The second array.
 * @return The distance between the two arrays.
 */
double distance(const Eigen::ArrayXd& p, const Eigen::ArrayXd& q) {
    return (p - q).matrix().norm();
}

// Define distance measures using Eigen::MatrixXd
/**
 * Calculates the single-linkage distance between two clusters.
 *
 * @param ci The first cluster represented as a matrix.
 * @param cj The second cluster represented as a matrix.
 * @return The single-linkage distance between the two clusters.
 */
double single_link(const Eigen::MatrixXd& ci, const Eigen::MatrixXd& cj) {
    double min_dist = std::numeric_limits<double>::infinity();
    for (int i = 0; i < ci.cols(); ++i) {
        for (int j = 0; j < cj.cols(); ++j) {
            double dist = distance(ci.col(i), cj.col(j));
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
    }
    return min_dist;
}

// Get distance measure based on M using function pointers
typedef double (*DistanceFunction)(const Eigen::MatrixXd&, const Eigen::MatrixXd&);

/**
 * Retrieves the distance measure based on the given parameter M.
 *
 * @param M The parameter used to determine the distance measure.
 * @return The distance function corresponding to the given parameter.
 */
DistanceFunction get_distance_measure(int M) {
    if (M == 0)
        return single_link;
    else
        return nullptr; // For other measures like complete_link or average_link
}
