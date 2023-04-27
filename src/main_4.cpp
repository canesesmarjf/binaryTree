#include <iostream>
#include <armadillo>
#include "BinaryTree.h"
#include <chrono>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>

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

  tree_params.dimensionality = 3;
  tree_params.min       = {-x_max, -v_max, 0     };
  tree_params.max       = {+x_max, +v_max, +v_max};
  tree_params.max_depth = {4     ,5      , 4     };

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

  // Create variables:
  int i;
  int search_dimensionality;
  node_TYP * leaf = NULL;
  arma::uvec ip;

  // Folder where data is to be stored
  string file_root = "./output_files/main_4/";

  // Get data:
  // ======================================================================
  i = 1811-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4a.csv", arma::csv_ascii);

  // Get data:
  // ======================================================================
  i = 18857-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4b.csv", arma::csv_ascii);

  // Get data:
  // ======================================================================
  i = 16491-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4c.csv", arma::csv_ascii);

  // Get data:
  // ======================================================================
  i = 9071-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4d.csv", arma::csv_ascii);

  // Get data:
  // ======================================================================
  i = 1520-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4e.csv", arma::csv_ascii);

  // Get data:
  // ======================================================================
  i = 7526-1;
  search_dimensionality = 3;
  leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_4f.csv", arma::csv_ascii);

  return 0;
}
