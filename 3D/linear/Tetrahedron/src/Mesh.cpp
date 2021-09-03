//
// Mesh.cpp
//
// Created by Wei Chen on 9/1/21
//

#include "Mesh.hpp"
#include "Utils.hpp"

#include <iostream>
#include <cmath>
#include <cstdlib>

namespace fem3d{

Mesh::Mesh(void){
    E = 1.0;
    nu = 0.3;
    D.resize(6, 6);
    D << 1.0-nu,     nu,     nu,              0.0,              0.0,              0.0,
             nu, 1.0-nu,     nu,              0.0,              0.0,              0.0,
             nu,     nu, 1.0-nu,              0.0,              0.0,              0.0,
            0.0,    0.0,    0.0, (1.0-2.0*nu)/2.0,              0.0,              0.0,
            0.0,    0.0,    0.0,              0.0, (1.0-2.0*nu)/2.0,              0.0,
            0.0,    0.0,    0.0,              0.0,              0.0, (1.0-2.0*nu)/2.0;
    double tmp = E / (1.0 + nu) / (1.0 - 2.0 * nu);
    D = D * tmp;

    if(!Utils::readMesh("../input/cube.txt", nvet, nele, vet, ele)){
        std::cout << "Logging: [fem3d] error on reading vertices" << std::endl;
        exit(-1);
    }
    ndof = nvet * 3;

    F.setZero(ndof);
    U.setZero(ndof);

    computeFeatures();
    linSysSolver = new EigenLibSolver();
    linSysSolver->set_pattern(vNeighbor);
    linSysSolver->analyze_pattern();
    std::cout << "Logging: [fem3d] mesh constructed" << std::endl;
}

Mesh::~Mesh(void){
    delete linSysSolver;
}

void Mesh::computeFeatures(void){ // compute vNeighbor
    vNeighbor.resize(0);
    vNeighbor.resize(nvet);

    for(int vI=0; vI<nvet; ++vI){
        vNeighbor[vI].insert(vI);
    }

    for(int eleI=0; eleI<nele; ++eleI){
        const Eigen::Matrix<int, 1, 4>& eleVInd = ele.row(eleI);
        for(int i=0; i<4; ++i){
            for(int j=i+1; j<4; ++j){
                vNeighbor[eleVInd(i)].insert(eleVInd(j));
                vNeighbor[eleVInd(j)].insert(eleVInd(i));
            }
        }
    }
}

void Mesh::computeK(void){ // assembly stiffness matrix
    linSysSolver->setZero();
    
    for(int eleI=0; eleI<nele; ++eleI){
        const Eigen::Matrix<int, 1, 4>& eleVInd = ele.row(eleI);
        Eigen::Matrix<double, 4, 3> X;
        X.row(0) = vet.row(eleVInd(0));
        X.row(1) = vet.row(eleVInd(1));
        X.row(2) = vet.row(eleVInd(2));
        X.row(3) = vet.row(eleVInd(3));

        int num = 12;
        Eigen::MatrixXd Ke(num, num);
        computeKe(X, Ke);

        Eigen::VectorXd edof(num);
        edof << eleVInd(0)*3, eleVInd(0)*3+1, eleVInd(0)*3+2,
                eleVInd(1)*3, eleVInd(1)*3+1, eleVInd(1)*3+2,
                eleVInd(2)*3, eleVInd(2)*3+1, eleVInd(2)*3+2,
                eleVInd(3)*3, eleVInd(3)*3+1, eleVInd(3)*3+2;
        for(int i=0; i<num; ++i){
            for(int j=0; j<num; ++j){
                linSysSolver->addCoeff(edof(i), edof(j), Ke(i, j));
            }
        }
    }
}

// compute element stiffness matrix
void Mesh::computeKe(const Eigen::Matrix<double, 4, 3>& X, Eigen::MatrixXd& Ke){
    Eigen::MatrixXd H(4, 4);
    H << 1.0, X(0, 0), X(0, 1), X(0, 2),
         1.0, X(1, 0), X(1, 1), X(1, 2),
         1.0, X(2, 0), X(2, 1), X(2, 2),
         1.0, X(3, 0), X(3, 1), X(3, 2);
    double V6 = H.determinant();

    Eigen::VectorXi rowInd(3);
    Eigen::VectorXi colInd(3);
    Eigen::MatrixXd sub(3, 3);

    colInd << 0, 2, 3;
    rowInd << 1, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double b1 = -sub.determinant();
    rowInd << 0, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double b2 = sub.determinant();
    rowInd << 0, 1, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double b3 = -sub.determinant();
    rowInd << 0, 1, 2;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double b4 = sub.determinant();

    colInd << 0, 1, 3;
    rowInd << 1, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double c1 = sub.determinant();
    rowInd << 0, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double c2 = -sub.determinant();
    rowInd << 0, 1, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double c3 = sub.determinant();
    rowInd << 0, 1, 2;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double c4 = -sub.determinant();


    colInd << 0, 1, 2;
    rowInd << 1, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double d1 = -sub.determinant();
    rowInd << 0, 2, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double d2 = sub.determinant();
    rowInd << 0, 1, 3;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double d3 = -sub.determinant();
    rowInd << 0, 1, 2;
    Utils::subMatrix(H, rowInd, colInd, sub);
    double d4 = sub.determinant();

    Eigen::MatrixXd B(6, 12);
    B << b1,  0,  0, b2,  0,  0, b3,  0,  0, b4,  0,  0,
          0, c1,  0,  0, c2,  0,  0, c3,  0,  0, c4,  0,
          0,  0, d1,  0,  0, d2,  0,  0, d3,  0,  0, d4,
         c1, b1,  0, c2, b2,  0, c3, b3,  0, c4, b4,  0,
          0, d1, c1,  0, d2, c2,  0, d3, c3,  0, d4, c4,
         d1,  0, b1, d2,  0, b2, d3,  0, b3, d4,  0, b4;
    B = B / V6;

    Ke = V6 / 6.0 * (B.transpose() * D * B);
}

void Mesh::boundaryK(void){
    for(int i=0; i<DBCV.size(); ++i){
        int DBCI = DBCV(i);
        linSysSolver->setZeroCol(DBCI);
        linSysSolver->setUnitRow(DBCI);
    }
}

void Mesh::solve(void){
    linSysSolver->factorize();
    linSysSolver->solve(F, U);
}

void Mesh::write(void){
    if(!Utils::writeU("../output/U.txt", U)){
        std::cout << "Logging: [fem3d] error on writing displacement U" << std::endl;
        exit(-1);
    }
}

// <<< interfaces
void Mesh::setLoad(void){ // user defined load
    for(int vI=0; vI<nvet; ++vI){
        if(abs(vet(vI, 2)-1.0)<eps && abs(vet(vI, 0))<eps){
            F(3*vI) = -0.01;
        }
    }
}

void Mesh::setDBC(void){ // user defined fixed dofs
    DBCV.resize(0);
    for(int vI=0; vI<nvet; ++vI){
        if(abs(vet(vI, 2))<eps){
            DBCV.conservativeResize(DBCV.size()+3);
            DBCV(DBCV.size()-3) = vI * 3;
            DBCV(DBCV.size()-2) = vI * 3 + 1;
            DBCV(DBCV.size()-1) = vI * 3 + 2;
        }
    }
}
// >>> interfaces

} // namespace
