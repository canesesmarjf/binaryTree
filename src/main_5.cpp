#include <iostream>
#include <armadillo>
#include "BinaryTree.h"
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
The Objetive of this script is to develop a C++ 3D vranic method that acts on the data produced by a 3D binary tree.
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

  // Start with a 3D binary tree:
  // ======================================================================
  tree_params_TYP tree_params;

  double v_max = 1.60E6;
  double x_max = 1;

  tree_params.dimensionality = 3;
  tree_params.min       = {-x_max, -v_max, 0     };
  tree_params.max       = {+x_max, +v_max, +v_max};
  tree_params.max_depth = {4     ,6      , 4     };

  // Assemble data into format needed for tree:
  // ===================================================================
  // Copy the data to arma::vec types:
  auto start = high_resolution_clock::now();
  arma::vec v_par = v_p.col(0);
  arma::vec v_per = v_p.col(1);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "time required for copy = " << duration.count() << endl;

  // Create a vector with pointers to the data:
  vector<arma::vec *> data = {&x_p, &v_par,&v_per};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 2D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Count how many nodes where populated:
  // ======================================================================
  int k = tree.count_nodes();
  cout << "total number of nodes populated: " << k << endl;

  // Folder where data is to be stored:
  // ======================================================================
  string file_root = "./output_files/main_5/";
  if (fs::is_directory(file_root) == false)
  {
    fs::create_directory(file_root);
  }

  // Create variables for find operations:
  // ======================================================================
  double xq;
  int i;
  int search_dimensionality;
  node_TYP * leaf = NULL;
  node_TYP * leaf_cube = NULL;
  arma::uvec ip;

  // Find data based on find(xq):
  // ======================================================================
  xq = 0.3;
  leaf = tree.find(xq);

  // Get data from random index from leaf:
  // ======================================================================
  i = leaf->ip[0];
  cout << "Random index from leaf->ip: " << i << endl;

  search_dimensionality = 3;
  leaf_cube = tree.find(i,data,search_dimensionality);

  cout << "leaf_cube = " << leaf_cube << endl;
  cout << "node coordinates: " << endl;
  leaf_cube->center.print();
  cout << "leaf_cube->ip.size() = " << leaf_cube->ip.size() << endl;

  // Save data to csv:
  ip = conv_to<arma::uvec>::from(leaf_cube->ip);
  ip.save(file_root + "ip_main_5a.csv", arma::csv_ascii);

  // Down sample data using vranic method:
  // ======================================================================
  // Select subset N:
  int N = round(ip.n_elem/1);
  cout << "N = " << N << endl;

  // Create aliases to data:
  arma::uvec rng = ip.subvec(0,N-1);
  arma::vec wi = a_p.elem(rng);
  arma::vec xi = x_p.elem(rng);
  arma::vec yi = v_par.elem(rng);
  arma::vec zi = v_per.elem(rng);

  /* Alternatively, we can use the following two options:

  arma::vec& wi = a_p;
  arma::vec& xi = x_p;
  arma::vec& yi = v_par;
  arma::vec& zi = v_per;

  arma::vec& wi = a_p;
  arma::vec& xi = *(data[0]);
  arma::vec& yi = *(data[1]);
  arma::vec& zi = *(data[2]);
  */

  // Calculate merge-cell statistics:
  double wt_N = sum(wi);
  arma::vec ri = wi/wt_N;

  // Expectation values:
  double E_x = dot(ri,xi);
  double E_y = dot(ri,yi);
  double E_z = dot(ri,zi);
  // double E_r = sqrt(pow(E_x,2) + pow(E_y,2) + pow(E_z,2));

  // Calculate deltas:
  arma::vec dx = xi - E_x;
  arma::vec dy = yi - E_y;
  arma::vec dz = zi - E_z;
  arma::vec dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

  // Calculate skewness vector for deltas:
  int n = 2;
  double mu3_dx = dot(ri,dx%pow(dr,n));
  double mu3_dy = dot(ri,dy%pow(dr,n));
  double mu3_dz = dot(ri,dz%pow(dr,n));
  double mu3_dr = sqrt(pow(mu3_dx,2) + pow(mu3_dy,2) + pow(mu3_dz,2));

  cout << "mu3_dx = " << mu3_dx << endl;
  cout << "mu3_dy = " << mu3_dy << endl;
  cout << "mu3_dz = " << mu3_dz << endl;

  // Calculate the merge-cell coordinate system:
  arma::vec x_hat = {1,0,0};
  arma::vec x_hat_prime = {mu3_dx/mu3_dr,mu3_dy/mu3_dr,mu3_dz/mu3_dr};
  arma::vec z_hat_prime = cross(x_hat,x_hat_prime);
  arma::vec y_hat_prime = cross(z_hat_prime,x_hat_prime);

  cout << "coordinate system test, dot(x_hat_prime,y_hat_prime) = " << dot(x_hat_prime,y_hat_prime) << endl;
  x_hat_prime.print("x_hat_prime");
  y_hat_prime.print("y_hat_prime");

  // Merge-cell coordinate system matrix:
  arma::mat e_prime(3,3);
  e_prime.col(0) = x_hat_prime;
  e_prime.col(1) = y_hat_prime;
  e_prime.col(2) = z_hat_prime;

  cout << "e_prime.n_rows() = " << e_prime.n_rows << endl;
  cout << "e_prime.n_cols() = " << e_prime.n_cols << endl;

  // Rotation matrix:
  arma::mat R = e_prime.t();

  // Convert deltas from standard coordinate to merge-cell coordinate system:
  arma::vec xi_prime = R(0,0)*xi + R(0,1)*yi + R(0,2)*zi;
  arma::vec yi_prime = R(1,0)*xi + R(1,1)*yi + R(1,2)*zi;
  arma::vec zi_prime = R(2,0)*xi + R(2,1)*yi + R(2,2)*zi;

  // Calculate merge-cell statistics in new coordinate system:
  // Expectation values:
  double E_x_prime = dot(ri,xi_prime);
  double E_y_prime = dot(ri,yi_prime);
  double E_z_prime = dot(ri,zi_prime);

  // Deltas:
  arma::vec dx_prime = xi_prime - E_x_prime;
  arma::vec dy_prime = yi_prime - E_y_prime;
  arma::vec dz_prime = zi_prime - E_z_prime;
  arma::vec dr_prime = sqrt(pow(dx_prime,2) + pow(dy_prime,2) + pow(dz_prime,2));

  // Standard deviation:
  double sigma_x_prime = sqrt(dot(ri,pow(dx_prime,2)));
  double sigma_y_prime = sqrt(dot(ri,pow(dy_prime,2)));
  double sigma_z_prime = sqrt(dot(ri,pow(dz_prime,2)));
  double sigma_r_prime = sqrt(dot(ri,pow(dr_prime,2)));

  cout << "sigma_r_prime = " << sigma_r_prime << endl;

  // Calculate deltas of new set M:
  int M = 6;
  arma::vec dx_prime_M(M);
  arma::vec dy_prime_M(M);
  arma::vec dz_prime_M(M);

  cout << "sigma_x_prime = " << sigma_x_prime << endl;
  cout << "sigma_y_prime = " << sigma_y_prime << endl;
  cout << "sigma_z_prime = " << sigma_z_prime << endl;

  dx_prime_M(0) = + sqrt(M/2)*sigma_x_prime;
  dx_prime_M(1) = - sqrt(M/2)*sigma_x_prime;
  dx_prime_M(2) = 0;
  dx_prime_M(3) = 0;
  dx_prime_M(4) = 0;
  dx_prime_M(5) = 0;

  dy_prime_M(0) = 0;
  dy_prime_M(1) = 0;
  dy_prime_M(2) = + sqrt(M/2)*sigma_y_prime;
  dy_prime_M(3) = - sqrt(M/2)*sigma_y_prime;
  dy_prime_M(4) = 0;
  dy_prime_M(5) = 0;

  dz_prime_M(0) = 0;
  dz_prime_M(1) = 0;
  dz_prime_M(2) = 0;
  dz_prime_M(3) = 0;
  dz_prime_M(4) = + sqrt(M/2)*sigma_z_prime;
  dz_prime_M(5) = - sqrt(M/2)*sigma_z_prime;

  dx_prime_M.print("dx_prime_M = ");
  dy_prime_M.print("dy_prime_M = ");
  dz_prime_M.print("dz_prime_M = ");

  // Convert deltas of M into standard frame:
  R = e_prime;
  arma::vec dx_M = R(0,0)*dx_prime_M + R(0,1)*dy_prime_M + R(0,2)*dz_prime_M;
  arma::vec dy_M = R(1,0)*dx_prime_M + R(1,1)*dy_prime_M + R(1,2)*dz_prime_M;
  arma::vec dz_M = R(2,0)*dx_prime_M + R(2,1)*dy_prime_M + R(2,2)*dz_prime_M;

  dx_M.print("dx_M = ");
  dy_M.print("dy_M = ");
  dz_M.print("dz_M = ");

  // Produce vectors for M set:
  arma::vec wi_M = arma::ones<vec>(M)*wt_N/M;
  arma::vec xi_M = E_x + dx_M;
  arma::vec yi_M = E_y + dy_M;
  arma::vec zi_M = E_z + dz_M;

  // Checks on M set vectors:
  double wt_M = sum(wi_M);
  ri = wi_M/wt_N;

  // Expectation values:
  double E_x_M = dot(ri,xi_M);
  double E_y_M = dot(ri,yi_M);
  double E_z_M = dot(ri,zi_M);

  // Calculate deltas:
  dx_M = E_x_M - xi_M;
  dy_M = E_y_M - yi_M;
  dz_M = E_z_M - zi_M;
  arma::vec dr_M = sqrt(pow(dx_M,2) + pow(dy_M,2) + pow(dz_M,2));
  double sigma_r_M = sqrt(dot(ri,pow(dr_M,2)));

  xi_M.print("xi_M = ");
  yi_M.print("yi_M = ");
  zi_M.print("zi_M = ");

  cout << "wt_M = " << wt_M << endl;
  cout << "wt_N = " << wt_N << endl;

  cout << "E_x_M = " << E_x_M << endl;
  cout << "E_x   = " << E_x   << endl;

  cout << "E_y_M = " << E_y_M << endl;
  cout << "E_y   = " << E_y   << endl;

  cout << "E_z_M = " << E_z_M << endl;
  cout << "E_z   = " << E_z   << endl;

  cout << "sigma_r_prime = " << sigma_r_prime << endl;
  cout << "sigma_r_M     = " << sigma_r_M     << endl;

  // Save data to postprocess in MATLAB:
  wi_M.save(file_root + "wi_M.csv", arma::csv_ascii);
  xi_M.save(file_root + "xi_M.csv", arma::csv_ascii);
  yi_M.save(file_root + "yi_M.csv", arma::csv_ascii);
  zi_M.save(file_root + "zi_M.csv", arma::csv_ascii);


  /* notes:
  arma::vec x1 = x_p.elem(in) creates a copy of x_p subview and x1 is effectively another vector.
  Consider that on each cube produced by the tree, we are not copying the entire data set but rather just a small fraction so perhaps it is ok to do a copy operation. What we want to avoid is to do a copy of the entire N_CP_MPI.

  The other option, is to use the aliases to the entire data set and just use the indices vector to select subsets as follows:

  arma::vec& x1 = x_p;
  double p = dot(x1.elem(in),x_p.elem(in));

  */

  return 0;
}
