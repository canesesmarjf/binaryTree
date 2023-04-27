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

  // Test the first find option:
  // ======================================================================
  // int i = 10318;
  // int i = 9500-1;
  int i = 9165-1;
  int search_dimensionality = 3;
  node_TYP * leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  arma::uvec ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save("ip_main_3.csv", arma::csv_ascii);

  // Test the second find option:
  // ======================================================================
  // This one is intended to be used with the node centers so that we can assemble a node list
  // node_x below will contain all the data stored undernearth that node
  double x_m = 0.2;
  node_TYP * node_x = tree.find(x_m);
  cout << "node_x = " << node_x << endl;
  cout << "node_x coordinates: " << endl;
  node_x->center.print();
  cout << "node_x->ip.size() = " << node_x->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(node_x->ip);
  ip.save("ip_main_3b.csv", arma::csv_ascii);

  // Lets take a random index from that node_x and get all the data stored in the 3D node:
  int p_count = node_x->p_count;
  std::srand(std::time(nullptr));
  int random_index = rand()%p_count;
  int iq = node_x->ip[random_index];

  node_TYP * node_v = tree.root->find(iq,data,3);
  ip = conv_to<arma::uvec>::from(node_v->ip);
  ip.save("ip_main_3c.csv", arma::csv_ascii);

  // Test clearing the tree from data:
  // ===================================================================
  cout << "Testing clearing the tree: " << endl;
  x_m = 0.3;
  node_x = tree.find(x_m);
  cout << "Prior to clearing, node_x->ip.size() = " << node_x->ip.size() << endl;

  // i = 1;
  // leaf = tree.root->find(i,data,search_dimensionality);
  // cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  tree.clear_all();

  x_m = 0.3;
  node_x = tree.find(x_m);
  cout << "After clearing, node_x->ip.size() = " << node_x->ip.size() << endl;

  // leaf = tree.root->find(i,data,search_dimensionality);
  // cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  return 0;
}
