#include "BinaryTree.h"
#include <iostream>

// Improvements:
// [DONE] 1 - Remove the else statement in the ::Insert() method
// [DONE] 2 - Use vector<node*> subnodes insteads of node_left and node_right
// [DONE] 3 - Pass i and r insteads of sometimes p or i,p and r, otherwise, somehow save the value of the data's pointer &r
// 4 - Make getter and setter functions to hide all variables as private
// 5- Improve code execution by optimizing dynamic memory allocation via (a) insert for counting and then (b) insert for allocating + appending


using namespace std;

// Constructors:
// ==================================================================================================================
node::node()
{
    cout << "default constructor" << endl;
}

node::node(double x_left, double x_right, int depth, int depth_max)
{
    // Node attributes:
    this->x_center   = (x_left + x_right)/2;
    this->x_left     = x_left;
    this->x_right    = x_right;
    this->depth      = depth;
    this->depth_max  = depth_max;
    
    // Allocate memory for subnodes:
    this->subnode.reserve(2);
    this->subnode[0] = NULL;
    this->subnode[1] = NULL;
    this->x_count    = 0;
    
    // To improve performance, we need to preallocate memory so that std::vector does not need to resize at every insertion:
    // This may be acheived in a two step fashion, where we first count how many points are inserted in each end node, and
    // then we loop again over all particles again and allocate memory at the very first insertion occurance.
    // this->ix.reserve(500); // Reserve memory
}

// Insert method:
// ==================================================================================================================
void node::Insert(uint i, arma::vec * r)
{
    // Objective:
    // Insert point into a subnode of current node. When maximum depth is reached, append point to node.
    
    // Current data point:
    double p = arma::as_scalar(r->at(i));
    
    // Check if data is within node's boundaries:
    // ===============================================
    if (!IsPointInsideBoundary(p))
    {
        // Warning message:
        cout << "point " << p <<" is outside domain" << endl;
        return;
    }
    
    // Check if maximum tree depth has been reached:
    // ===============================================
    // If yes, append point's index to leaf and RETURN up the stack
    if (HasNodeReachMaxDepth())
    {
        // Append point:
        this->ix.push_back(i);
        this->x_count++;

        // Diagnostics:
        if (0)
        {
            cout << "i = " << i << endl;
            cout << "size of ix = " << this->ix.size() << endl;
            cout << "this->x_count = " << this->x_count << endl;
            cout << "x_left = " << this->x_left << endl;
            cout << "x_right = " << this->x_right << endl;
            for (int j = 0; j < this->ix.size(); j++)
            {
                cout << "j = " << j << ", ix[j] = " << this->ix[j] <<", value = " << r->at(this->ix[j]) <<endl;
            }
        }
        
        // Return control to calling stack:
        return;
    }
    
    // Determine which subnode to insert point:
    // ========================================
    int node_index = WhichSubNodeDoesItBelongTo(p);
    
    // Check if subnode needs to be created:
    // ====================================
    if ( !DoesSubNodeExist(node_index) )
    {
        CreateSubNode(node_index);
    }
    
    // Insert point into subnode:
    // ==========================
    this->subnode[node_index]->Insert(i,r);
        
} // node::Insert


// ==================================================================================================================
bool node::IsPointInsideBoundary(double p)
{
    // Objective:
    // if p is inside the boundaries of the node, return true, otherwise false
    
    // Define boundaries of node:
    double x_left  = this->x_left;
    double x_right = this->x_right;
        
    // Create boolean result:
    bool flag = ((p >= x_left) && (p <= x_right));
    
    
    return flag;
}

// ==================================================================================================================
bool node::HasNodeReachMaxDepth()
{
    int depth     = this->depth;
    int depth_max = this->depth_max;
        
    if (depth >= depth_max)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// ==================================================================================================================
int node::WhichSubNodeDoesItBelongTo(double p)
{
    //   +------------------+------------------+
    //   |  node_left = 1  |  node_right = 0  |
    
    // Center location of present node:
    double x_center = this->x_center;
    
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

// Find method:
// ==================================================================================================================
node * node::Find(double xq)
{
    
    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }
    
    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth())
    {
        return this;
    }
    
    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq);
        
    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }
       
    // Drill further into subnode:
    return this->subnode[node_index]->Find(xq);
    
}

// ==================================================================================================================
bool node::DoesSubNodeExist(int node_index)
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

// ==================================================================================================================
void node::CreateSubNode(int node_index)
{
    // Attributes for new subnode:
    int depth     = this->depth + 1;
    int depth_max = this->depth_max;
    double x_left;
    double x_right;
    
    switch (node_index)
    {
    case 0: // Right subnode:
        {
            // Attributes for new subnode:
            x_left  = this->x_center;
            x_right = this->x_right;
                
            // Exit:
            break;
        }
    case 1: // Left subnode:
        {
            // Attributes for new subnode:
            x_left  = this->x_left;
            x_right = this->x_center;

            // Exit:
            break;
        }
    }
    
    // Create new subnode:
    this->subnode[node_index] = new node(x_left, x_right, depth, depth_max);
    
}

// ==================================================================================================================
node * node::GetSubNode(int index)
{
    return this->subnode[index];
}
