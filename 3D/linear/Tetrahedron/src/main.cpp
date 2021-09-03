//
// main.cpp
//
// Created by Wei Chen on 8/30/21
//

#include "Mesh.hpp"

#include <Eigen/Eigen>
#include <iostream>

int main(){
    fem3d::Mesh mesh = fem3d::Mesh();
    mesh.setLoad();
    std::cout << "Logging: [fem3d] set load" << std::endl;
    mesh.setDBC();
    std::cout << "Logging: [fem3d] set DBC" << std::endl;

    mesh.computeK();
    std::cout << "Logging: [fem3d] assembly stiffness matrix" << std::endl;
    mesh.boundaryK();

    mesh.solve();
    std::cout << "Logging: [fem3d] solve" << std::endl;
    mesh.write();
    std::cout << "Logging: [fem3d] write displacement to output/U.txt" << std::endl;
    return 0;
}
