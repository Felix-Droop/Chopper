#pragma once

#include <chopper/union/hyperloglog.hpp>

struct clustering_node 
{
    // children in the tree
    size_t left;
    size_t right;
    // hll sketch of the union if the node is still a root
    hyperloglog hll;
};