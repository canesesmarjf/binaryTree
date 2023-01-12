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
    
    // Generate computational domain:
    double x_left = -4;
    double x_right = +4;
    
    // Generate random sample:
//    arma_rng::set_seed_random();
    int N_CP = 1e+4;
    double mean_r = 0;
    double std_r = (x_right - x_left)/10;
    vec r = randn(N_CP)*std_r + mean_r;
    
    // Save to test data:
    r.save("output_files/data.txt",raw_ascii);
    
    // Define maximum depth of tree:
    int depth_max = 7;
    
    // Root node depth:
    int depth_root = 0;

    // Create instance of binary tree:
    node BinaryTree(x_left,x_right,depth_root,depth_max);
            
    //Build binary tree:
    auto start = high_resolution_clock::now();
    for (int i = 0; i < r.size(); i++)
    {
        BinaryTree.Insert(i,&r);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;

    // Find data in tree:
    int queries = 3;
    node * result = NULL;

    // Create string stream object:
    stringstream sso;

    for (int i = 0; i < queries; i++)
    {
        // Prompt and gather data from user:
        cout << "Enter a number: " << endl;
        double xq;
        cin >> xq;

        // Get data from binary tree:
        result = BinaryTree.Find(xq);

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
