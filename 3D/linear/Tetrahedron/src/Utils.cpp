//
// Utils.cpp
//
// Created by Wei Chen on 9/1/21
//

#include "Utils.hpp"

namespace fem3d{

bool Utils::readMesh(const std::string& filePath, int& nvet, int& nele,
              Eigen::MatrixXd& vet, Eigen::MatrixXi& ele){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        return false;
    }

    fscanf(in, "%d%d", &nvet, &nele);
    vet.resize(nvet, 3);
    ele.resize(nele, 4);

    for(int vI=0; vI<nvet; ++vI){
        fscanf(in, "%lf%lf%lf", &vet(vI, 0), &vet(vI, 1), &vet(vI, 2));
    }
    for(int eleI=0; eleI<nele; ++eleI){
        fscanf(in, "%d%d%d%d", &ele(eleI, 0), &ele(eleI, 1), &ele(eleI, 2), &ele(eleI, 3));
    }
    ele.array() -= 1;

    fclose(in);
    return true;
}

bool Utils::writeU(const std::string& filePath, const Eigen::VectorXd& U){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
        return false;
    }

    for(int i=0; i<U.size(); ++i){
        fprintf(out, "%.20f\n", U(i));
    }

    fclose(out);
    return true;
}

// return sub-matrix(3*3) of H(4*4)
void Utils::subMatrix(const Eigen::MatrixXd& H, const Eigen::VectorXi& rowInd,
                      const Eigen::VectorXi& colInd, Eigen::MatrixXd& sub){
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            sub(i, j) = H(rowInd(i), colInd(j));
        }
    }
}

} // namespace
