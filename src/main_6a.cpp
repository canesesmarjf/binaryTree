#include <iostream>
#include <armadillo>
#include "BinaryTree.h"
#include "Vranic.h"
#include <chrono>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <filesystem>

using namespace std::chrono;
using namespace std;
using namespace arma;
namespace fs = filesystem;

/* OBJECTIVE:
Test basic operation of multi=dim binary tree
*/

int main()
{
  // Load input data:
  // ===================================================================
  arma::vec x = {0.1,0.2,0.3};
  arma::vec y = {0.3,0.4,0.5};
  arma::vec z = {0.5,0.6,0.7};

  // Start with a 3D binary tree:
  // ======================================================================
  tree_params_TYP tree_params;

  tree_params.dimensionality = 3;
  tree_params.min       = {-1 , -1, -1};
  tree_params.max       = {+1 , +1, +1};
  tree_params.max_depth = {1  , 1 , 1 };

  // Assemble data into format needed for tree:
  // ===================================================================
  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&x, &y,&z};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 2D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Count how many leaf nodes were populated:
  // ======================================================================
  int k = tree.count_leaf_nodes();
  cout << "total number of leaf nodes populated: " << k << endl;

  // Count how many leaf points were inserted:
  // ======================================================================
  k = tree.count_leaf_points();
  cout << "total number of leaf points inserted: " << k << endl;

  return 0;
}
