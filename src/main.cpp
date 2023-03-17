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
  r.save("output_files/data.txt",raw_ascii);

  // Create a binary tree empty of data:
  // ===================================================================
  // Tree parameters:
  double x_left  = -2;
  double x_right = +2;
  int depth_max = 6;
  int num_elem  = N_CP;
  binaryTree_TYP tree(x_left,x_right,depth_max,N_CP);

  tree.print_info(20);

  // Populate binary tree with data by inserting entire dataset "r" into binary tree structure:
  // ===================================================================
  auto start = high_resolution_clock::now();
  tree.insert_all(&r);
  auto stop = high_resolution_clock::now();

  // Calculate time it took to assemble:
  // ===================================================================
  auto duration = duration_cast<microseconds>(stop - start);
  cout << duration.count() << endl;

  // Save data to file:
  // ===================================================================
  tree.save_data_all("result_");

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
  tree.print_info(43);

  tree.clear_all();

  tree.print_info(43);

  return 0;
}
