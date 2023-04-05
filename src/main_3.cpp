#include <iostream>
#include <armadillo>
#include "BinaryTree.h"
#include <chrono>
#include <string>
#include <sstream>

using namespace std::chrono;
using namespace std;
using namespace arma;

int main()
{
  // Load input data:
  // ===================================================================
  arma::vec x_p;
  arma::mat v_p;
  arma::vec a_p;

  x_p.load("input_files/Step_1_x_p.csv",csv_ascii);
  v_p.load("input_files/Step_1_v_p.csv",csv_ascii);
  a_p.load("input_files/Step_1_a_p.csv",csv_ascii);

  cout << "columns of v_p = " << v_p.n_cols << endl;
  cout << "rows of v_p = " << v_p.n_rows << endl;

  // Start with a 3D binary tree:
  // ======================================================================
  tree_params_TYP tree_params;

  double v_max = 1.60E6;
  double x_max = 1;

  tree_params.dimensionality = 1;
  tree_params.min       = {-x_max, -v_max, 0     };
  tree_params.max       = {+x_max, +v_max, +v_max};
  tree_params.max_depth = {4     ,5     , 4     };

  // Assemble data into format needed for tree:
  // ===================================================================

  // Copy the data to arma::vec types:
  auto start = high_resolution_clock::now();
  arma::vec v_par = v_p.col(0);
  arma::vec v_per = v_p.col(1);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "t2 = " << duration.count() << endl;

  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&x_p, &v_par,&v_per};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 2D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Test the find option:
  // ======================================================================
  // int i = 10318;
  int i = 9500-1;
  node_TYP * result = tree.root->find(i,data);

  // Save data to csv:
  // ======================================================================
  arma::uvec ix = conv_to<arma::uvec>::from(result->ix);
  ix.save("ix_main_3.csv", arma::csv_ascii);
  
  // NEED TO DEVELOP A WAY TO CLEAN TREE!


  // // Test clearing the tree from data:
  // // ===================================================================
  // tree.print_info(43);
  //
  // tree.clear_all();
  //
  // tree.print_info(43);

  return 0;
}
