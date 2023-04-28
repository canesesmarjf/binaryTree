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
  // Generate random sample:
  // ===================================================================
  arma_rng::set_seed_random();
  int N_CP = 1e+4;
  double mean_r = 0;
  double x_min = -2;
  double x_max = +2;
  double std_r = (x_max - x_min)/15;
  vec r = randn(N_CP)*std_r + mean_r;

  // Save to test data:
  r.save("output_files/main_1/data.txt",raw_ascii);

  // 1D binary tree parameters:
  // ======================================================================
  tree_params_TYP tree_params;

  tree_params.dimensionality = 1;
  tree_params.min       = {+x_min};
  tree_params.max       = {+x_max};
  tree_params.max_depth = {6     };

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Assemble data into format needed for tree:
  // ===================================================================
  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&r};

  // Insert 2D points into tree:
  // ======================================================================
  auto start = high_resolution_clock::now();
  tree.insert_all(data);
  auto stop = high_resolution_clock::now();

  // Create a vector with points inside each cell:
  int Nx = pow(2,tree_params.max_depth[0]);
  double dx_m = (x_max-x_min)/(Nx);
  arma::vec x_m(Nx);
  for (int ii = 0; ii < Nx; ii++ )
  {
    x_m[ii] = x_min + dx_m/2 + ii*dx_m;
  }

  // Find data in each node:
  node_TYP * leaf = NULL;
  arma::uvec ip;
  string fileName;
  for (int ii = 0; ii < Nx; ii++ )
  {
    leaf = tree.find(x_m[ii]);

    if (leaf != NULL)
    {
      // Allocate memory:
      ip.resize(leaf->ip.size());

      // Put data into armadillo containers:
      ip = arma::conv_to<arma::uvec>::from(leaf->ip);

      // Save file:
      fileName = "output_files/main_1/result_" + to_string(ii) + ".txt";
      cout << "Saving data with fileName = " << fileName << endl;
      ip.save(fileName,arma::raw_ascii);
    }
    else
    {
        cout << "no points in node i =  " << ii << endl;
    }

  }

  // Create a binary tree empty of data:
  // ===================================================================
  // Tree parameters:
  // double x_left  = -2;
  // double x_right = +2;
  // int depth_max = 6;
  // int num_elem  = N_CP;
  // binaryTree_TYP tree(x_left,x_right,depth_max,N_CP);

  // tree.print_info(20);

  // Populate binary tree with data by inserting entire dataset "r" into binary tree structure:
  // ===================================================================
  // auto start = high_resolution_clock::now();
  // tree.insert_all(&r);
  // auto stop = high_resolution_clock::now();

  // Calculate time it took to assemble:
  // ===================================================================
  // auto duration = duration_cast<microseconds>(stop - start);
  // cout << duration.count() << endl;

  // Save data to file:
  // ===================================================================
  // tree.save_data_all("output_files/main_1/" + "result_");

  /*
  if (0)
  {
    // Create string stream object:
    stringstream sso;

    for (int i = 0; i < tree.get_num_nodes(); i++)
    {
        // Print data to CLI
        if (tree.node_list[i]->ix.empty() == false)
        {
            cout << "x_left  = " << tree.node_list[i]->x_left  << endl;
            cout << "x_right = " << tree.node_list[i]->x_right << endl;
            cout << "x_count = " << tree.node_list[i]->x_count << endl;

            // Put data into armadillo containers:
            arma::uvec ix = conv_to<arma::uvec>::from(tree.node_list[i]->ix);
            arma::vec rix = r.elem(ix);

            // Print data to CLI:
            rix.print("r[ix] = ");
        }
        else
        {
            cout << "no points in node i =  " << i << endl;
        }
    }
  }
  */

  // Test clearing the tree from data:
  // ===================================================================
  // tree.print_info(43);

  tree.clear_all();

  // tree.print_info(43);

  return 0;
}
