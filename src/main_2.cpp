#include <iostream>
#include <armadillo>
#include "BinaryTree.h"
#include <chrono>
#include <string>
#include <sstream>

using namespace std::chrono;
using namespace std;
using namespace arma;

/* OBJECTIVE:
The Objetive of this script is to demonstrate the use the binary tree for 2D and also to demonstrate the use of the dimensionality option
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

  // Start with a 2D binary tree for velocity space:
  // ======================================================================
  tree_params_TYP tree_params;

  double v_max = 1.60E6;
  double v_par_min = round(arma::min(v_p.col(0))/1e4)*1e4;
  double v_par_max = round(arma::max(v_p.col(0))/1e4)*1e4;
  double v_per_min = 0;
  double v_per_max = round(arma::max(v_p.col(1))/1e4)*1e4;

  // Normalizing variable:
  double v_norm = v_per_max;

  tree_params.dimensionality = 2;
  tree_params.min       = {-v_max, 0     };
  tree_params.max       = {+v_max, +v_max};
  tree_params.max_depth = {5     , 4     };

  // Assemble data into format needed for tree:
  // ===================================================================
  // The data in the current form of PICOS++ is inside a arma::mat. In order to supply
  // only arma::vec to the binary tree, we need to divide this mat. The simple way is to
  // use arma::vec vpar_p = IONS->at(ss).v_p.col(0); however, this method relies on copying data
  // This means that are no longer accessing the originarl data but a copy of it.
  // there seems to be at least two solutions to this right now. the first one is to modify PICOS++
  // so that we have vpar_p and vper_p as arma::vec and not together in arma::mat.
  // The other solution that we found out by using chatGPT is to use arma::subview_col<double>.
  // This will produce reference to the original data; however, it is not a arma::vec and thus this
  // affects how we define the type of the argument.
  // Based on this discussion it seems that the best pathway to have a general case is to
  // describe all data using arma::vec throughout PICOS++ so that we can use the binary tree
  // with arma::vec x_p and arma::vec vpar_p and arma::vec vper_p

  // Copy the data to arma::vec types:
  auto start = high_resolution_clock::now();
  arma::vec v_par = v_p.col(0);
  arma::vec v_per = v_p.col(1);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "t2 = " << duration.count() << endl;

  // The following demonstrates that the copying process takes considerable time compared
  // to just using pointers:
  start = high_resolution_clock::now();
  arma::subview_col<double> v_par_sub = v_p.col(0);
  arma::subview_col<double> v_per_sub = v_p.col(1);
  stop = high_resolution_clock::now();
  duration = duration_cast<microseconds>(stop - start);
  cout << "t2 = " << duration.count() << endl;

  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&v_par,&v_per};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 2D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Folder where data is to be stored
  string file_root = "./output_files/main_2/";

  // Search tree by index:
  // ======================================================================
  // Initialize variables:
  int i;
  int search_dimensionality;
  node_TYP * leaf = NULL;
  arma::uvec ip;

  i = 10318-1;

  /* How is search dimensionality used?
  Consider how multi dimensional tree is organized and how it searches. It first searches along the first dimension which has an associated depth. Once it reaches the final depth of the 1st dimension it can move on to the next dimension and progress in depth and so on.
  At each final depth of a corresponding dimension, the local leaf node stores all the data that is to be used in subsequent depths. So the search dimensionality can be used to determine the depth of the search in terms of the dimensions.

  Thus consider a 3D search, a find with a dimensionlity of 3 returns a set of indices in a space space cube. A dimensionality of 2 returns a slice along dim1 and dim 2.
  */

  search_dimensionality = 1;
  leaf = tree.root->find(i,data,search_dimensionality);
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip.resize(leaf->ip.size());
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_2a.csv", arma::csv_ascii);

  search_dimensionality = 2;
  leaf = tree.root->find(i,data,search_dimensionality);
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip.resize(leaf->ip.size());
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_2b.csv", arma::csv_ascii);

  // Search tree by position of first node:
  // ======================================================================
  double x_m = 0.25e6;
  leaf = tree.root->find(x_m);
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  // Save data to csv:
  ip.resize(leaf->ip.size());
  ip = conv_to<arma::uvec>::from(leaf->ip);
  ip.save(file_root + "ip_main_2c.csv", arma::csv_ascii);

  // Test clearing the tree from data:
  // ===================================================================
  i = 1;
  leaf = tree.root->find(i,data,search_dimensionality);
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  tree.clear_all();

  leaf = tree.root->find(i,data,search_dimensionality);
  cout << "leaf->ip.size() = " << leaf->ip.size() << endl;

  return 0;
}
