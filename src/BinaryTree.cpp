#include "BinaryTree.h"
#include <iostream>

using namespace std;

/*
// Constructors:
// ================================================================================================================
tree_params_TYP::tree_params_TYP(double x_left, double x_right, int depth_max, int num_elem)
{
  this->x_left = x_left;
  this->x_right = x_right;
  this->depth_max = depth_max;
  this->num_elem = num_elem;

  this->num_nodes = pow(2,depth_max);
  this->length = x_right - x_left;
  this->dx = (x_right - x_left)/this->num_nodes;
  this->elem_per_node = round(num_elem/this->num_nodes);
  this->node_centers = arma::linspace(x_left,(x_right-dx),this->num_nodes) + dx/2;
}
*/

// ================================================================================================================
binaryTree_TYP::binaryTree_TYP(tree_params_TYP * tree_params)
{
  // Update tree_params:
  tree_params->num_nodes = pow(2,sum(tree_params->max_depth));
  tree_params->dim_levels = cumsum(tree_params->max_depth);

  // Store pointer to tree parameters
  this->tree_params = tree_params;

  // Create root node:
  uint depth_root = 0;
  arma::vec min = tree_params->min;
  arma::vec max = tree_params->max;
  root = new node_TYP(min,max,depth_root,tree_params);
}

/*
// ================================================================================================================
binaryTree_TYP::binaryTree_TYP(double x_left, double x_right, int depth_max, int num_elem)
{
  // Interface variables:
  tree_params = new tree_params_TYP(x_left, x_right, depth_max, num_elem);

  // Allocate memory to node list, which is a vector of pointers:
  node_list.resize(tree_params->num_nodes);

  // Create root node:
  int depth_root = 0;
  root = new node_TYP(tree_params->x_left,tree_params->x_right,depth_root,tree_params);

  // Traverse entire tree:
  // This is done by inserting node_centers so as to create all the node
  // at the final depth while not writting any data to the nodes.
  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    bool write_data = false;
    root->insert(nn,&tree_params->node_centers,write_data);
  }

  // Assemble list of nodes into a vector of pointers:
  assemble_node_list();
}
*/

// ================================================================================================================
void binaryTree_TYP::insert_all(vector<arma::vec *> data)
{
  // Insert points into nodes:
  this->root->insert_all(data);
}

/*
// ================================================================================================================
int binaryTree_TYP::get_num_nodes()
{
  return tree_params->num_nodes;
}

// ================================================================================================================
arma::vec binaryTree_TYP::get_node_centers()
{
  return tree_params->node_centers;
}

// ================================================================================================================
int binaryTree_TYP::get_max_depth()
{
  return tree_params->depth_max;
}

// ================================================================================================================
void binaryTree_TYP::assemble_node_list()
{
  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    // Get the location of the nth node center:
    double xq = tree_params->node_centers.at(nn);

    // Copy the pointer of the node containing xq into the node list:
    node_list.at(nn) = root->find(xq);
  }
}
*/

// ================================================================================================================
void binaryTree_TYP::clear_all()
{
  if (NULL == root)
  {
    cout << "Tree has NULL root node" << endl;
    return;
  }

  root->clear_node();
}

/*
// ================================================================================================================
void binaryTree_TYP::print_info(int ii)
{
  cout << " " << endl;
  cout << "tree.node_list[ii]->x_left        = " << node_list[ii]->x_left << endl;
  cout << "tree.node_list[ii]->x_center      = " << node_list[ii]->x_center << endl;
  cout << "tree.node_list[ii]->x_right       = " << node_list[ii]->x_right << endl;
  cout << "tree.node_list[ii]->depth         = " << node_list[ii]->depth << endl;
  cout << "tree.node_list[ii]->ip.capacity() = " << node_list[ii]->ip.capacity() << endl;
  cout << "tree.node_list[ii]->p_count       = " << node_list[ii]->p_count << endl;
  cout << "tree.node_list[ii]->ip.size()     = " << node_list[ii]->ip.size() << endl;
  cout << " " << endl;
}

// ================================================================================================================
void binaryTree_TYP::save_data(int ii, string prefix)
{
  // Put data into armadillo containers:
  arma::uvec ip = arma::conv_to<arma::uvec>::from(node_list[ii]->ip);

  // Save data to file:
  string fileName = "output_files/" + prefix + to_string(ii) + ".txt";
  cout << "Saving data with fileName = " << fileName << endl;
  ip.save(fileName,arma::raw_ascii);
}

// ================================================================================================================
void binaryTree_TYP::save_data_all(string prefix)
{
  for (int ii = 0; ii < get_num_nodes(); ii++)
  {
      if (node_list[ii]->ip.empty() == false)
      {
          // Save data for node ii:
          save_data(ii, prefix);
      }
      else
      {
          cout << "no points in node i =  " << ii << endl;
      }
  } // for
}
*/

// ================================================================================================================
node_TYP::node_TYP(arma::vec min, arma::vec max, uint depth, tree_params_TYP * tree_params)
{
  // Node attributes:
  this->center  = (min + max)/2;
  this->min     = min;
  this->max     = max;
  this->depth       = depth;
  this->tree_params = tree_params;

  // Allocate memory for subnodes:
  this->subnode.reserve(2);
  this->subnode[0] = NULL;
  this->subnode[1] = NULL;
  this->p_count    = 0;
}

// insert method:
// ================================================================================================================
void node_TYP::insert(uint i, vector<arma::vec *> data, bool write_data)
{
    // Objective:
    // insert point into a subnode of current node. When maximum depth is reached, append point to node.

    // Determine dimension to operate on:
    int dim = sum(depth > tree_params->dim_levels);

    // Current data point:
    double p = arma::as_scalar(data[dim]->at(i)); // make it a function that takes in i, data,dim

    // Check if data is within node's boundaries:
    // ===============================================
    if (!IsPointInsideBoundary(p,dim))
    {
        // Warning message:
        cout << "point " << p <<" is outside domain" << endl;
        return;
    }

    // Check if maximum tree depth has been reached:
    // ===============================================
    // If yes, append point's index to leaf and RETURN up the stack
    if (HasNodeReachMaxDepth(dim))
    {
        // Append point:
        if (write_data)
        {
          this->ip.push_back(i);
          this->p_count++;
        }

        // Return control to calling stack IF max_depth for all dimensions has been reached:
        uint total_dims = tree_params->dimensionality;
        if (depth >= tree_params->dim_levels(total_dims-1))
        {
          return;
        }
        else
        {
          // Evaluate data at next dimension:
          dim++;
          p = arma::as_scalar(data[dim]->at(i));
        }
    }

    // Determine which subnode to insert point:
    // ========================================
    int node_index = WhichSubNodeDoesItBelongTo(p,dim);

    // Check if subnode needs to be created:
    // ====================================
    if ( !DoesSubNodeExist(node_index) )
    {
        CreateSubNode(node_index,dim);
    }

    // Insert point into subnode:
    // ==========================
    this->subnode[node_index]->insert(i,data,write_data);

} // node_TYP::insert

// insert_all method:
// ================================================================================================================
void node_TYP::insert_all(vector<arma::vec *> data)
{
  for (int i = 0; i < data[0]->size(); i++)
  {
      // if ( i == 9500 - 1 )
      // {
      //   cout << i << endl;
      // }
      this->insert(i,data,true);
  }
}


// =================================================================================================================
bool node_TYP::IsPointInsideBoundary(double p, int dim)
{
    // Objective:
    // if p is inside the boundaries of the node, return true, otherwise false

    // Define boundaries of node:
    double x_left  = this->min[dim];
    double x_right = this->max[dim];

    // Create boolean result:
    bool flag = ((p >= x_left) && (p <= x_right));

    return flag;
}

// ================================================================================================================
bool node_TYP::HasNodeReachMaxDepth(int dim)
{
    int depth = this->depth;

    if (depth >= tree_params->dim_levels[dim]) // Greater only? equal only?
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// ================================================================================================================
int node_TYP::WhichSubNodeDoesItBelongTo(double p, int dim)
{
    //   +------------------+------------------+
    //   |  node_left = 1  |  node_right = 0  |

    // Center location of present node:
    double x_center = this->center[dim];

    // Number associated with subnode:
    int node_index;

    // Determine which location point belongs to:
    if ( p < x_center)
    {
        node_index = 1; // left
    }
    else
    {
        node_index = 0; // right
    }

    return node_index;
}

// find method:
// =================================================================================================================
node_TYP * binaryTree_TYP::find(int i,vector<arma::vec *> data,int search_dimensionality)
{
  return this->root->find(i,data,search_dimensionality);
}

// find method:
// =================================================================================================================
node_TYP * binaryTree_TYP::find(double xq)
{
  return this->root->find(xq);
}

// find method:
// =================================================================================================================
node_TYP * node_TYP::find(double xq)
{
    // 1D data input inplies the following:
    int search_dimensionality = 1;

    // Calculate the relevant dimension based on the depth of the present node:
    int dim = sum(depth > tree_params->dim_levels);

    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq,dim))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth(dim))
    {
        return this;
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq,dim);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(xq);

}

// find method:
// =================================================================================================================
node_TYP * node_TYP::find(uint i, vector<arma::vec *> data, int search_dimensionality)
{
    // Check the value of search _dimensionality:
    assert(search_dimensionality <= tree_params->dimensionality && "search_dimensionality needs to be less or equal to dimensionality of tree");

    // Calculate the relevant dimension based on the depth of the present node:
    int dim = sum(depth > tree_params->dim_levels);

    // Extract the data based on the dimension being considered:
    double xq = arma::as_scalar(data[dim]->at(i));

    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq,dim))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth(dim))
    {
      // Return control to calling stack IF maximum search dimension has been reached:
      if (depth >= tree_params->dim_levels(search_dimensionality-1))
      {
        return this;
      }
      else
      {
        // Evaluate data at next dimension:
        dim++;
        xq = arma::as_scalar(data[dim]->at(i));
      }
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq,dim);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(i,data,search_dimensionality);
}

// ================================================================================================================
bool node_TYP::DoesSubNodeExist(int node_index)
{
    if (NULL == subnode[node_index])
    {
        // Does not exist:
        return 0;
    }
    else
    {
        // It already exists:
        return 1;
    }
}

// ================================================================================================================
void node_TYP::CreateSubNode(int node_index, int dim)
{
    // Attributes for new subnode:
    uint depth = this->depth + 1;
    arma::vec min = this->min; //x_left;
    arma::vec max = this->max; //x_right;

    switch (node_index)
    {
    case 0: // Right subnode:
        {
            // Attributes for new subnode:
            min[dim] = this->center[dim];
            max[dim] = this->max[dim];

            // Exit:
            break;
        }
    case 1: // Left subnode:
        {
            // Attributes for new subnode:
            min[dim] = this->min[dim];
            max[dim] = this->center[dim];

            // Exit:
            break;
        }
    }

    // Create new subnode:
    this->subnode[node_index] = new node_TYP(min, max, depth, tree_params);
}

// ================================================================================================================
void node_TYP::clear_node()
{
  this->p_count = 0;
  this->ip.clear();

  if (NULL != this->subnode[0])
  {
    this->subnode[0]->clear_node();
  }

  if (NULL != this->subnode[1])
  {
    this->subnode[1]->clear_node();
  }
}
