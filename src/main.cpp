#include <iostream>
#include <armadillo>
#include "H5Cpp.h"
#include "BinaryTree.h"
#include <chrono>
#include <string>
#include <sstream>

using namespace std::chrono;

//using namespace H5;
using namespace std;
using namespace arma;

int main()
{
  // Setup binary tree:
  // ===================================================================
  // Tree attributes:
  double x_left  = -2;
  double x_right = +2;
  int depth_max = 6;

  // Create instance of binary tree:
  BinaryTree_TYP tree(x_left,x_right,depth_max);

  // Generate random sample:
  // ===================================================================
  arma_rng::set_seed_random();
  int N_CP = 1e+4;
  double mean_r = 0.5;
  double std_r = (x_right - x_left)/7;
  vec r = randn(N_CP)*std_r + mean_r;

  // Save to test data:
  // ===================================================================
  r.save("output_files/data.txt",raw_ascii);

  // Assemble binary tree using the input data:
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
  node * result = NULL;

  // Create string stream object:
  stringstream sso;

  for (int i = 0; i < tree.get_num_nodes(); i++)
  {
      // Get data from binary tree:
      result = tree.node_list[i];

      // Print data to CLI
      if (NULL != result)
      {
          cout << "x_left  = " << result->x_left  << endl;
          cout << "x_right = " << result->x_right << endl;
          cout << "x_count = " << result->x_count << endl;

          // Print data using std::vector:
          if (0)
          {
              for (int j = 0; j < result->x_count ; j++ )
              {
                  cout << "r[ix] = " << r[result->ix[j]] << " , ix = " << result->ix[j] << endl;
              }
              arma::mat set;
          }

          // Put data into armadillo containers:
          arma::uvec ix = conv_to<arma::uvec>::from(result->ix);
          arma::vec rix = r.elem(ix);

          // Print data to CLI:
          rix.print("r[ix] = ");

          // Save data to file:
          sso << i;
          string fileName = "result_";
          fileName = "output_files/" + fileName + sso.str() + ".txt";
          cout << "Saving data with fileName = " << fileName << endl;
          ix.save(fileName,raw_ascii);
          sso.str("");
      }
      else
      {
          cout << "no points in node: " << endl;
      }

  }

  return 0;
}
