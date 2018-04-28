// This file is part of GQCG-hf.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-hf is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-hf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-hf.  If not, see <http://www.gnu.org/licenses/>.
#ifndef HF_COMMON_HPP
#define HF_COMMON_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>



namespace hf {


/*
 * TYPEDEFS
 */

typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> DensityFunction;
typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &, const Eigen::Tensor<double, 4> &)> TwoElectronMatrixFunction;


}  // namespace hf



#endif  // HF_COMMON_HPP
