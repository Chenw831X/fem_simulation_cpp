//
// Mesh.hpp
//
// Created by Wei Chen on 9/1/21
//

#ifndef Mesh_hpp
#define Mesh_hpp

#include "EigenLibSolver.hpp"

#include <Eigen/Eigen>
#include <vector>
#include <set>

namespace fem3d{

class Mesh{
public: // owned data
    int nvet, nele, ndof;
    Eigen::MatrixXd vet; // vertices coordinates
    Eigen::MatrixXi ele; // vertice index of each tetrahedron
    double E, nu;
    Eigen::MatrixXd D; // constitutive matrix

public: // owned features
    double eps = 1e-8;
    Eigen::VectorXd F; // load of each dof
    Eigen::VectorXd U; // dofs' displacement to be computed
    Eigen::VectorXi DBCV; // dofs in DBC
    std::vector<std::set<int>> vNeighbor; // records all vertices' indices adjacent to each vertice
    LinSysSolver* linSysSolver;

public: // constructor
    Mesh(void);
    ~Mesh(void);

public:
    void computeFeatures(void);
    void computeK(void);
    void computeKe(const Eigen::Matrix<double, 4, 3>& X, Eigen::MatrixXd& Ke);
    void boundaryK(void);
    void solve(void);
    void write(void);

public: // interface: set load and DBC
    void setLoad(void);
    void setDBC(void);
};

} // namespace

#endif // Mesh_hpp
