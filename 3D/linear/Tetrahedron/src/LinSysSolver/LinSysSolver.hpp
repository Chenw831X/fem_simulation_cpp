//
// LinSysSolver.hpp
//
// Created by Wei Chen on 9/2/21
//

#ifndef LinSysSolver_hpp
#define LinSysSolver_hpp

#include <Eigen/Eigen>
#include <iostream>
#include <set>
#include <map>
#include <vector>

namespace fem3d{

class LinSysSolver{
protected:
    int numRows, nnz;
    Eigen::VectorXi ia, ja;
    std::vector<std::map<int, int>> IJ2aI;
    Eigen::VectorXd a;

public:
    virtual ~LinSysSolver(void){};

    virtual void set_pattern(const std::vector<std::set<int>>& vNeighbor){
        int nvet = static_cast<int>(vNeighbor.size());
        numRows = nvet * 3;
        nnz = 0;
        for(int vI=0; vI<nvet; ++vI){
            nnz += static_cast<int>(vNeighbor[vI].size()) * 9;
        }

        ia.setZero(numRows+1);
        ja.setZero(nnz);
        IJ2aI.resize(0);
        IJ2aI.resize(numRows);
        a.setZero(nnz);

        for(int vI=0; vI<nvet; ++vI){
            int num = static_cast<int>(vNeighbor[vI].size()) * 3;
            for(int dof=0; dof<3; ++dof){
                int row = vI * 3 + dof;
                ia(row+1) = ia(row) + num;
                int cnt = 0;
                for(const auto& nbVI : vNeighbor[vI]){
                    ja(ia(row)+cnt) = nbVI * 3;
                    ja(ia(row)+cnt+1) = nbVI * 3 + 1;
                    ja(ia(row)+cnt+2) = nbVI * 3 + 2;
                    IJ2aI[row][nbVI*3] = ia(row) + cnt;
                    IJ2aI[row][nbVI*3+1] = ia(row) + cnt + 1;
                    IJ2aI[row][nbVI*3+2] = ia(row) + cnt + 2;
                    cnt += 3;
                }
            }
        }
    }

    virtual void analyze_pattern(void) = 0;
    virtual void factorize(void) = 0;
    virtual void solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& res) = 0;

    virtual void setCoeff(int rowI, int colI, double val){
        assert(rowI<numRows);
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder!=IJ2aI[rowI].end());
        a[finder->second] = val;
    }

    virtual void addCoeff(int rowI, int colI, double val){
        assert(rowI<numRows);
        const auto finder = IJ2aI[rowI].find(colI);
        assert(finder!=IJ2aI[rowI].end());
        a[finder->second] += val;
    }

    virtual void setZero(void){
        a.setZero();
    }

    virtual void multiply(const Eigen::VectorXd& x, Eigen::VectorXd& Ax){
        assert(x.size()==numRows);
        Ax.setZero(numRows);
        for(int rowI=0; rowI<numRows; ++rowI){
            for(const auto& colI : IJ2aI[rowI]){
                Ax(rowI) += a[colI.second] * x(colI.first);
            }
        }
    }

    virtual void setUnitRow(int rowI){
        assert(rowI<numRows);
        for(const auto& colI : IJ2aI[rowI]){
            a[colI.second] = (rowI==colI.first);
        }
    }

    virtual void setZeroCol(int colI){
        assert(colI<numRows);
        for(int rowI=0; rowI<numRows; ++rowI){
            const auto finder = IJ2aI[rowI].find(colI);
            if(finder != IJ2aI[rowI].end()){
                a[finder->second] = 0.0;
            }
        }
    }
};

} // namespace

#endif // LinSysSolver_hpp
