//
// Utils.hpp
//
// Created by Wei Chen on 9/1/21
//

#ifndef Utils_hpp
#define Utils_hpp

#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

namespace fem3d{

// a static class implementing some assitant functions
class Utils{
public:
    static bool readMesh(const std::string& filePath, int& nvet, int& nele,
                         Eigen::MatrixXd& vet, Eigen::MatrixXi& ele);
    static bool writeU(const std::string& filePath, const Eigen::VectorXd& U);
    static void subMatrix(const Eigen::MatrixXd& H, const Eigen::VectorXi& rowInd,
                          const Eigen::VectorXi& colInd, Eigen::MatrixXd& sub);
};

} // namespace

#endif // Utils_hpp
