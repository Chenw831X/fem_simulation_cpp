//
// EigenLibSolver.hpp
//
// Created by Wei Chen on 9/2/21
//

#ifndef EigenLibSolver_hpp
#define EigenLibSolver_hpp

#include "LinSysSolver.hpp"

namespace fem3d{

class EigenLibSolver : public LinSysSolver{
    typedef LinSysSolver Base;

protected:
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    virtual void set_pattern(const std::vector<std::set<int>>& vNeighbor);

    virtual void analyze_pattern(void);
    virtual void factorize(void);
    virtual void solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& res);

    virtual void setCoeff(int rowI, int colI, double val);
    virtual void addCoeff(int rowI, int colI, double val);
    virtual void setZero(void);
    virtual void setUnitRow(int rowI);
    virtual void setZeroCol(int colI);
};

} // namespace

#endif // EigenLibSolver_hpp