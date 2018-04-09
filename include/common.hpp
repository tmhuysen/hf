#ifndef HF_COMMON_HPP
#define HF_COMMON_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

namespace hf{


/*
 * TYPEDEFS
 */

typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> DensityFunction;
typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &, const Eigen::Tensor<double, 4> &)> TwoElectronMatrixFunction;

}


#endif //HF_COMMON_HPP
