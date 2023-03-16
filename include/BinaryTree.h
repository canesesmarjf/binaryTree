#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <armadillo>

// This header file contains two main classes:
// BinaryTree_TYP and node
// The BinaryTree_TYP makes use of node type but all binary tree operations are performed by node

class node
{
public:
  // Variables:
  int x_count;            // Data counts in this node
  std::vector<uint> ix;   // Indices of data in this node
  double x_center;        // Center position of the node
  double x_left;          // Left boundary of node
  double x_right;

  // Constructor:
  node();
  node(double x_left, double x_right, int depth, int depth_max);

  // Methods:
  void insert(uint i, arma::vec * r); // Insert the ith element of vector *r
  void insert_all(arma::vec * r);     // Insert all elements of vector *r
  node * find(double xq);             // Find and return pointer of node corresponding to position xq
  node * get_subnode(int index);       // ?

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

class BinaryTree_TYP
{
public:
  // Constructor:
  BinaryTree_TYP();
  BinaryTree_TYP(double x_left, double x_right, int depth_max);

  // Variables:
  std::vector<node *> node_list;  // List of pointers to nodes at maxmimum depth

  // Methods:
  void insert_all(arma::vec * r);
  int get_num_nodes();
  arma::vec get_node_centers();
  int get_max_depth();

private:
  // Variables:
  int depth_max; // Maximum depth of tree
  int num_nodes; // Number of nodes in tree
  double dx; // Distance between adjacent node centers
  arma::vec node_centers; // vector containing list of node centers at maximum depth
  node * root; // Root node of tree

  // Methods:
  void assemble_node_list();
};

#endif
