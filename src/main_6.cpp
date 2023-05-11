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
The Objetive of this script is to carry out an entire distribution resampleing using 3D vranic method that acts on the data produced by a 3D binary tree.
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

  // Normalize data:
  double v_max = 1.85E6;
  double x_max = 1;

  double x_norm = x_max;
  double y_norm = v_max;
  double z_norm = v_max;

  x_p = x_p/x_max;
  v_p = v_p/v_max;

  // Copy data to be used in down-sampling method:
  arma::vec a_p_copy = a_p;
  arma::vec x_p_copy = x_p;
  arma::mat v_p_copy = v_p;

  // Start with a 3D binary tree:
  // ======================================================================
  tree_params_TYP tree_params;

  tree_params.dimensionality = 3;
  tree_params.min       = {-1 , -1, 0 };
  tree_params.max       = {+1 , +1, +1};
  tree_params.max_depth = {5  , 6 , 5 };

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
  vector<arma::vec *> data = {&x_p, &v_par, &v_per, &a_p};

  // Create tree based on parameters:
  // ======================================================================
  binaryTree_TYP tree(&tree_params);

  // Insert 2D points into tree:
  // ======================================================================
  tree.insert_all(data);

  // Count how many nodes where populated:
  // ======================================================================
  int k = tree.count_leaf_nodes();
  cout << "total number of leaf nodes populated: " << k << endl;

  // Count how many points where inserted:
  // ======================================================================
  k = tree.count_leaf_points();
  cout << "total number of leaf points inserted: " << k << endl;

  // Folder where data is to be stored:
  // ======================================================================
  string file_root = "./output_files/main_6/";
  if (fs::is_directory(file_root) == false)
  {
    fs::create_directory(file_root);
  }

  // Create grid for each dimension in binary tree:
  // ======================================================================
  // This section probably needs to be part of the tree_params_TYP class:
  // tree_params.calculate_grids();

  double Lx = tree_params.max[0] - tree_params.min[0];
  double Ly = tree_params.max[1] - tree_params.min[1];
  double Lz = tree_params.max[2] - tree_params.min[2];

  int Nx = pow(2,tree_params.max_depth[0]);
  int Ny = pow(2,tree_params.max_depth[1]);
  int Nz = pow(2,tree_params.max_depth[2]);

  double dx = Lx/Nx;
  double dy = Ly/Ny;
  double dz = Lz/Nz;

  arma::vec xq = tree_params.min[0] + dx/2 + regspace(0,Nx-1)*dx;
  arma::vec yq = tree_params.min[1] + dy/2 + regspace(0,Ny-1)*dy;
  arma::vec zq = tree_params.min[2] + dz/2 + regspace(0,Nz-1)*dz;

  // Create pointers for leaf nodes:
  // ======================================================================
  vector<node_TYP *> leaf_x(Nx);
  node_TYP * leaf_y = NULL;
  node_TYP * leaf_z = NULL;

  // Calculate leaf_x list:
  // ======================================================================
  arma::ivec p_count(Nx);
  for (int xx = 0; xx < Nx ; xx++)
  {
   cout << "xq(xx) = " << xq(xx) << endl;
   leaf_x[xx] = tree.find(xq(xx));
   p_count[xx] = leaf_x[xx]->p_count;
  }

  p_count.save(file_root + "p_count_" + ".csv", arma::csv_ascii);

  // Save contents of distribution contained in all the tree:
  // ======================================================================
  arma::imat v_count(Ny,Nz,fill::zeros);

  // Loop over leaf_x:
  cout << "leaf_x.size() = " << leaf_x.size() << endl;
  cout << "leaf_x.size() = " << Nx << endl;

  for (int xx = 0; xx < Nx ; xx++)
  {
    // Check if current leaf_x is in surplus:
    if (1)
    {
      // Loop over y:
      for (int yy = 0; yy < Ny; yy++)
      {
       leaf_y = leaf_x[xx]->find(yq[yy],1);
       if (leaf_y == NULL)
       {
         v_count(yy,span::all).zeros();
       }
       else
       {
         // Loop over z:
         for (int zz = 0; zz < Nz; zz++)
         {
           leaf_z = leaf_y->find(zq[zz],2);
           if (leaf_z == NULL)
           {
             v_count(yy,zz) = 0;
           }
           else
           {
             v_count(yy,zz) = leaf_z->p_count;
           }
         }
       }

      }
    }
    // Save velocity space density for every leaf_x node:
    cout << "xq[xx] = " << xq[xx] << " , xx = " << to_string(xx) << endl;
    v_count.save(file_root + "v_count_" + to_string(xx) + ".csv", arma::csv_ascii);
    v_count.zeros();
  }

  // Save contents of distribution contained ONLY in surplus region of tree:
  // ======================================================================
  std::vector<double> x_surplus;
  int mean_pcount = mean(p_count);

  // Loop over leaf_x:
  for (int xx = 0; xx < Nx ; xx++)
  {
    // Clear v_count:
    v_count.fill(0);

    // Check if current leaf_x is in surplus:
    int surplus = p_count[xx] - mean_pcount; //
    if (surplus <= 0)
      continue;

    // Record surplus x location:
    x_surplus.push_back(xq[xx]);

    // Loop over y:
    for (int yy = 0; yy < Ny; yy++)
    {
      // Get current leaf_y:
      leaf_y = leaf_x[xx]->find(yq[yy],1);

       if (leaf_y == NULL)
         continue;

       // Loop over z:
       for (int zz = 0; zz < Nz; zz++)
       {
         leaf_z = leaf_y->find(zq[zz],2);

         if (leaf_z == NULL)
          continue;

         // Record total number of particles in leaf_z cube:
         v_count(yy,zz) = leaf_z->p_count;

       } // Loop over z

    } // Loop over y

    // Save velocity space density for every leaf_x node:
    v_count.save(file_root + "surplus_v_count_" + to_string(xx) + ".csv", arma::csv_ascii);

  } // Loop over x

  // Save x_surplus to file:
  arma::vec x_out(x_surplus.data(),x_surplus.size());
  x_out.save(file_root + "x_surplus.csv", arma::csv_ascii);

  // Apply Vranic resampling method::
  // ======================================================================
  // The following loop does the following:
  // - Loop over all "x" x_nodes
  // Select only surplus "x" nodes
  // Loop over y-z slice and assemble v_count
  // One v_count is assembled, sort cubes by descending order
  // Loop over cubes and down sample until total number of particles is equal to desired value

  vranic_TYP vranic;
  std::vector<uint> ip_free;
  int N_min = 10;
  int N = 0;
  int M = 6;
  int particle_surplus = 0;
  int exit_flag_z;

  // Input:
  // N_min, N, M
  // ip_free
  // data
  // leaf_x from tree: leaf_x[xx]->p_count
  // xq, yq, zq

  // Output:
  // data
  // ip_free

  // Loop over leaf_x:
  for (int xx = 0; xx < Nx ; xx++)
  {
    // Clear flags:
    exit_flag_z = 0;

    // Calculate particle surplus
    particle_surplus = p_count[xx] - mean_pcount;

    // Check if current leaf_x is in surplus:
    if (particle_surplus > 0)
    {
      // Record surplus x location:
      x_surplus.push_back(xq[xx]);

      // Loop over y:
      for (int yy = 0; yy < Ny; yy++)
      {
        if (exit_flag_z == 1)
          break;

        // Get current leaf_y:
        leaf_y = leaf_x[xx]->find(yq[yy],1);

         if (leaf_y != NULL)
         {
           // Loop over z:
           for (int zz = 0; zz < Nz; zz++)
           {
             leaf_z = leaf_y->find(zq[zz],2);
             if (leaf_z != NULL)
             {
               // Total number of particles in leaf_z cube:
               N = leaf_z->p_count;

               if (N >= N_min)
               {
                 // Adjust N so that surplus count does not become negative:
                 // Negative means that leaf_x will be in deficit.
                 if (particle_surplus - N < 0)
                {
                  N = particle_surplus;

                  // Check that new set N is not smaller than M:
                  if (N < (M+1))
                  {
                    exit_flag_z = 1;
                    break;
                  }
                }

                 // Create objects for down-sampling:
                 merge_cell_TYP set_N(N);
                 merge_cell_TYP set_M(M);

                 // First N indices of particles in current leaf_z:
                 arma::uvec ip = conv_to<uvec>::from(leaf_z->ip);

                 // Create set_N:
                 set_N.xi = data[0]->elem(ip.head(N));
                 set_N.yi = data[1]->elem(ip.head(N));
                 set_N.zi = data[2]->elem(ip.head(N));
                 set_N.wi = data[3]->elem(ip.head(N));

                 // Calculate set M:
                 // Pass by reference
                 vranic.down_sample(&set_N, &set_M);

                 // Print statistics:
                 cout << "Set N: " << endl;
                 vranic.print_stats(&set_N);

                 // Print statistics:
                 cout << "Set M: " << endl;
                 vranic.print_stats(&set_M);

                 for (int ii = 0; ii < N; ii++)
                 {
                   // Get global index:
                   int jj = ip(ii);

                   if (ii < set_M.n_elem)
                   {
                     // Down-sample global distribution:
                     a_p_copy(jj)   = set_M.wi(ii);
                     x_p_copy(jj)   = set_M.xi(ii);
                     v_p_copy(jj,0) = set_M.yi(ii);
                     v_p_copy(jj,1) = set_M.zi(ii);
                   }
                   else
                   {
                     // Record global indices that correspond to memory locations that are to be repurposed:
                     ip_free.push_back(jj);

                     // Clear values:
                     a_p_copy(jj)   = -1;
                     x_p_copy(jj)   = -1;
                     v_p_copy(jj,0) = -1;
                     v_p_copy(jj,1) = -1;
                   }

                   // Decrement particle surplus:
                   particle_surplus--;
                 }
               }
             }
           } // Loop over z
         }
      } // Loop over y
    } // Surplus check
  } // Loop over x

  // Save data:
  arma::file_type format = arma::csv_ascii;
  a_p_copy.save(file_root + "a_p_copy" + ".csv", format);
  x_p_copy.save(file_root + "x_p_copy" + ".csv", format);
  v_p_copy.save(file_root + "v_p_copy" + ".csv", format);

  return 0;


  // Find data based on find(xq):
  // ======================================================================
  // xq = 0.39/x_norm;
  // leaf = tree.find(xq);

  /*
  // Get data from random index from leaf:
  // ======================================================================
  for (int j = 0; j < leaf->ip.size(); j++)
  {
    i = leaf->ip[j];
    if (data[2]->at(i) > 0.333)
      break;
  }
  cout << "Random index from leaf->ip: " << i << endl;
  cout << "vper(i) = " << data[2]->at(i) << endl;

  search_dimensionality = 3;
  leaf_cube = tree.find(i,data,search_dimensionality);

  cout << "leaf_cube = " << leaf_cube << endl;
  cout << "node coordinates: " << endl;
  leaf_cube->center.print();
  cout << "leaf_cube->ip.size() = " << leaf_cube->ip.size() << endl;

  // Save data to csv:
  ip_cube = conv_to<arma::uvec>::from(leaf_cube->ip);
  ip_cube.save(file_root + "ip_main_5a.csv", arma::csv_ascii);

  // Calculate down-sampled set using vranic method:
  // ======================================================================
  // Select subset N:
  int N = round(ip_cube.n_elem);
  int M = 6;
  arma::uvec rng = ip_cube.subvec(0,N-1);

  // Number of particles in set_M
  cout << "N = " << N << endl;

  // Create objects for down-sampling:
  merge_cell_TYP set_N(N);
  merge_cell_TYP set_M(M);
  vranic_TYP vranic;

  // Populate set N:
  // This process requires a copy operation
  set_N.wi = a_p.elem(rng);
  set_N.xi = x_p.elem(rng);
  set_N.yi = v_par.elem(rng);
  set_N.zi = v_per.elem(rng);

  // Calculate set M:
  // Pass by reference
  vranic.down_sample(&set_N, &set_M);

  // Print statistics:
  cout << "Set N: " << endl;
  vranic.print_stats(&set_N);

  // Print statistics:
  cout << "Set M: " << endl;
  vranic.print_stats(&set_M);

    // Save data to postprocess in MATLAB:
  set_M.wi.save(file_root + "wi_M.csv", arma::csv_ascii);
  set_M.xi.save(file_root + "xi_M.csv", arma::csv_ascii);
  set_M.yi.save(file_root + "yi_M.csv", arma::csv_ascii);
  set_M.zi.save(file_root + "zi_M.csv", arma::csv_ascii);

  for (int ii = 0; ii < ip_cube.n_elem; ii++)
  {
    // Get global index:
    int jj = ip_cube(ii);

    if (ii < set_M.n_elem)
    {
      // Down-sample global distribution:
      a_p_copy(jj)   = set_M.wi(ii);
      x_p_copy(jj)   = set_M.xi(ii);
      v_p_copy(jj,0) = set_M.yi(ii);
      v_p_copy(jj,1) = set_M.zi(ii);
    }
    else
    {
      // Record global indices that correspond to memory locations that are to be repurposed:
      ip_free.push_back(jj);

      // Clear values:
      a_p_copy(jj)   = -1;
      x_p_copy(jj)   = -1;
      v_p_copy(jj,0) = -1;
      v_p_copy(jj,1) = -1;
    }
  }

  // Save data to csv:
  arma::uvec free_indices = conv_to<arma::uvec>::from(ip_free);
  free_indices.save(file_root + "ip__free_main_5a.csv", arma::csv_ascii);

  // Save copy of distribution:
  a_p_copy.save(file_root + "a_p_copy_main_5.csv", arma::csv_ascii);
  x_p_copy.save(file_root + "x_p_copy_main_5.csv", arma::csv_ascii);
  v_p_copy.save(file_root + "v_p_copy_main_5.csv", arma::csv_ascii);

  */
}
