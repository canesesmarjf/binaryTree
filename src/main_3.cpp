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

/* OBJECTIVE:
The Objetive of this script is to demonstrate the use the binary tree for 3D and show how to use the find(xq) and find(i,data,dimensionality) options
*/

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
  auto start = high_resolution_clock::now();

  // Copy the data to arma::vec types:
  arma::vec v_par = v_p.col(0);
  arma::vec v_per = v_p.col(1);

  auto stop = high_resolution_clock::now();

  // Check the time it took to copy data over:
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "t2 = " << duration.count() << endl;

  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&x_p, &v_par,&v_per};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 3D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Folder where data is to be stored
  string file_root = "./output_files/main_3/";

  // Test the first find option:
  // ======================================================================
  int i = 9165-1;
  int search_dimensionality = 3;
  node_TYP * leaf = tree.find(i,data,search_dimensionality);

  cout << "leaf = " << leaf << endl;
  cout << "node coordinates: " << endl;
  leaf->center.print();
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  arma::uvec ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_3a.csv", arma::csv_ascii);

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
  ip.save(file_root + "ip_main_3b.csv", arma::csv_ascii);

  // Lets take a random index from that node_x and get all the data stored in the 3D node:
  int p_count = node_x->p_count;
  // std::srand(std::time(nullptr));
  int random_index = rand()%p_count;
  int iq = node_x->ip[random_index];

  cout << "random integer from x_m find is: " << iq << endl;

  node_TYP * node_v = tree.root->find(iq,data,3);
  ip = conv_to<arma::uvec>::from(node_v->ip);
  ip.save(file_root + "ip_main_3c.csv", arma::csv_ascii);

  // Test clearing the tree from data:
  // ===================================================================
  cout << "Testing clearing the tree: " << endl;
  x_m = 0.3;
  node_x = tree.find(x_m);
  cout << "Prior to clearing, node_x->ip.size() = " << node_x->ip.size() << endl;

  tree.clear_all();

  x_m = 0.3;
  node_x = tree.find(x_m);
  cout << "After clearing, node_x->ip.size() = " << node_x->ip.size() << endl;

  // Test deleteing node:
  // ===================================================================
  cout << "Testing deleting all nodes" << endl;
  cout << "Prior to deleting ..." << endl;
  cout << "Address stored in node_x = "  << node_x << endl;
  cout << "Capacity of node_x->ip = " << node_x->ip.capacity() << endl;

  tree.delete_nodes();

  cout << "After deleting nodes ..." << endl;
  cout << "Address strored in node_x = "  << node_x << endl;
  cout << "Capacity of node_x->ip = " << node_x->ip.capacity() << endl;

  return 0;
}
