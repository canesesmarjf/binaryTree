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
    node * GetSubNode(int index);

    private:
        
    // Tree depth level attributes:
    int depth;
    int depth_max;
    
    // Subnodes within this node:
    // subnode[0] : right_node
    // subnode[1] : left_node
    //   +------------------+------------------+
    //   |  left_node = 1   |  right_node = 0  |
    std::vector<node *> subnode;

        
    // Methods:
    bool IsPointInsideBoundary(double p);
    bool HasNodeReachMaxDepth();
    int WhichSubNodeDoesItBelongTo(double p);
    bool DoesSubNodeExist(int subNode);
    void CreateSubNode(int subNode);
};

#endif
