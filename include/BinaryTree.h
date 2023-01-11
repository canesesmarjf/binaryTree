#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <armadillo>

// Header file for binary search tree

class node
{
    public:
    
    // Data stored in node:
    int x_count; // data counts in this node
    std::vector<uint> ix; // indices of data in this node
    
    // Node attributes:
    double x_center;
    double x_left;
    double x_right;
    
    // Constructor:
    node();
    node(double x_left, double x_right, int depth, int depth_max);
    
    // Methods:
    void Insert(uint i, arma::vec * r);
    node * Find(double xq);
    
    private:
        
    // Tree depth level attributes:
    int depth;
    int depth_max;
    
    // Subnodes:
    node * node_left;
    node * node_right;
        
    // Methods:
    bool IsPointInsideBoundary(double p);
    bool HasNodeReachMaxDepth();
    int WhichSubNodeDoesItBelongTo(double p);
    void CreateSubNodeIfItDoesNotExist(double p);
    void InsertPointIntoSubNode(uint i, arma::vec * r);
};

#endif
