#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <armadillo>
using namespace std;

// This header file contains three main classes:
// - tree_params_TYP
// - binaryTree_TYP
// - node_TYP

// Is important to node that a tree is composed of nodes; thus a tree as a whole will have parameters such as:
// - Maximum depth
// - Left and right boundaries
// - Number of nodes in the final layer
// Moreover, each node will also have its own parameters:
// - current depth
// - Left and right boundaries
// - Data stored in node

// =====================================================================================
class tree_params_TYP
{
public:
  uint dimensionality;
  arma::vec min;
  arma::vec max;
  arma::uvec max_depth;
  arma::uvec dim_levels; // Determines when the tree needs to change to the next dimension
  uint num_nodes; // Number of nodes based on 2^depth_max

  // Interface variables:
  // double x_left;  // Left boundary of tree
  // double x_right; // Right boundary of tree
  // int depth_max; // Maximum depth of tree

  // Derived variables:
  // double length; // Length of entire domain
  // double dx;
  // int num_elem; // Total number of elements to hold in tree

  // int elem_per_node; // number of elements per node based on a uniform distribution
  // arma::vec node_centers; // Vector containing the center locations of the nodes

  // Constructor:
  tree_params_TYP(){};
  tree_params_TYP(double x_left, double x_right, int depth_max, int num_elem);
};

// =====================================================================================
class node_TYP
{
public:
  // Data variables:
  int x_count;            // Number of points indexed in the current node
  std::vector<uint> ix;   // Indices of data appended to this node

  // Node parameters:
  uint depth;
  uint dim;
  arma::vec min;
  arma::vec max;
  arma::vec center;

  // Constructor:
  node_TYP(){};
  node_TYP(arma::vec min, arma::vec max, uint depth, tree_params_TYP * tree_params);

  // double x_center;        // Center position of the node
  // double x_left;          // Left boundary of node
  // double x_right;         // Right boundary of node

  // Constructor:
  // node_TYP(double x_left, double x_right, int depth, tree_params_TYP * tree_params);

  // Methods:
  void insert(uint i, vector<arma::vec *> data, bool write_data); // Insert the ith element of data
  void insert_all(vector<arma::vec *> data); // Insert all elements of the data
  node_TYP * find(uint i, vector<arma::vec *> data);
  // node_TYP * find(double xq); // Find and return pointer of node corresponding to position xq
  // node_TYP * get_subnode(int index); // ?

private:
  // Variables:
  tree_params_TYP * tree_params; // Pointer to tree parameters

  // Subnodes within this node:
  // subnode[0] : right_node
  // subnode[1] : left_node
  //   +------------------+------------------+
  //   |  left_node = 1   |  right_node = 0  |
  std::vector<node_TYP *> subnode;

  // Methods:
  bool IsPointInsideBoundary(double p);
  bool HasNodeReachMaxDepth();
  int WhichSubNodeDoesItBelongTo(double p);
  bool DoesSubNodeExist(int subNode);
  void CreateSubNode(int subNode);
};

// =====================================================================================
class binaryTree_TYP
{
public:
  // Constructor:
  binaryTree_TYP();
  binaryTree_TYP(tree_params_TYP * tree_params);
  // binaryTree_TYP(double x_left, double x_right, int depth_max, int num_elem);

  // Variables:
  node_TYP * root; // Root node of tree
  // std::vector<node_TYP *> node_list; // List of pointers to nodes at maxmimum depth
  tree_params_TYP * tree_params; // Pointer to tree attributes

  // Methods:
  void insert_all(vector<arma::vec *> data);
  // int get_num_nodes();
  // arma::vec get_node_centers();
  // int get_max_depth();
  // void clear_all();
  // void print_info(int ii);
  // void save_data_all(string prefix);

private:
  // Variables:
  // node_TYP * root; // Root node of tree

  // Methods:
  // void assemble_node_list();
  // void save_data(int ii, string prefix);
};

#endif
