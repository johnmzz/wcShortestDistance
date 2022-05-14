#pragma once
#ifndef _WGraph_H
#define _WGraph_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
//#include <assert.h>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "Utils.h"

using namespace std;

class WGraph {
public:
    AdjList graph;  // adj list of edges, eg. [[v1,v2,v3], [v2,v3,v5], ...]
    AdjList wlist;  // weight of edges, eg. [[1,2,2], [4,2,2], ...]
    uint32_t nsize;
    uint32_t msize;
    uint32_t wsize;
    uint32_t max_w;
    vector<uint32_t> ndegree;  // degree of each node?

    ////for building index
    vector<vector<uint32_t>> labeling_v;
    vector<vector<uint32_t>> offset;
    vector<vector<uint32_t>> labeling_d;
    vector<vector<uint32_t>> labeling_w;

    // tree decomposition
    AdjList hgraph_;
    vector<set<uint32_t>> hgraph;
    vector<uint32_t> decomp_order;

    WGraph();
    WGraph(string, string, bool);

    void read_graph(ifstream&);
    void read_edgelist(ifstream&);
    void save_graph_to_bin(string, string);
    void get_graph_statistics();
    void print_graph();
    void preprocess(int);
    void sort_by_degree();
    void check_degree_order();

    //********************online solutions************************//
    uint32_t min_distance(vector<uint32_t>&, boost::dynamic_bitset<>&);
    uint32_t dijkstra(AdjList&, uint32_t, uint32_t);
    uint32_t constrained_shortest_distance_naive(uint32_t, uint32_t, uint32_t, double&);
    uint32_t constrained_shortest_distance_dijkstra(uint32_t, uint32_t, uint32_t, double&);
    uint32_t constrained_shortest_distance_plus(uint32_t, uint32_t, uint32_t, double&);

    //****************2-hop labeling based solutions*******************//
    // query with index
    uint32_t query_while_indexing_vertex(uint32_t, uint32_t, uint32_t);
    uint32_t query_while_indexing_vertex_V2(uint32_t, uint32_t, uint32_t, vector<vector<pair<uint32_t, uint32_t>>>&);
    uint32_t query_while_indexing_vertex_V3(uint32_t, uint32_t, uint32_t, vector<vector<uint32_t>>&, vector<vector<uint32_t>>&);
    uint32_t query_while_indexing_vertex_V4(uint32_t, uint32_t, uint32_t, vector<vector<pair<uint32_t, uint32_t>>>&, vector<uint32_t>& , vector<uint32_t>&);
    uint32_t query_while_indexing_vertex_V5(uint32_t, uint32_t, uint32_t, vector<vector<pair<uint32_t, uint32_t>>>&, vector<uint32_t>& , vector<uint32_t>&, vector<uint32_t>&);
    uint32_t query_with_index_vertex(uint32_t, uint32_t, uint32_t, double&);
    uint32_t query_with_index_vertex_skip(uint32_t, uint32_t, uint32_t, uint32_t);
    uint32_t query_for_pll(uint32_t, uint32_t);

    // index construction algorithms
    void build_index(string, int);
    void vertex_prioritized_indexing(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&);
    void vertex_prioritized_indexing_plus(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&, boost::dynamic_bitset<>&);
    void vertex_prioritized_indexing_V2(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&, boost::dynamic_bitset<>&);
    void vertex_prioritized_indexing_V3(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&, boost::dynamic_bitset<>&);
    void vertex_prioritized_indexing_V4(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&, boost::dynamic_bitset<>&);
    void vertex_prioritized_indexing_V5(uint32_t, boost::dynamic_bitset<>&, vector<uint32_t>&, boost::dynamic_bitset<>&, vector<uint32_t>&, vector<uint32_t>&);

    void pruned_landmark_labeling(uint32_t, vector<uint32_t>&);

    void print_index();
    void print_index_pll();
    void get_index_size();
    void check_minimality();

    // tree decomposition
    void tree_decomp();
    void print_hgraph();
    void print_vec(vector<uint32_t> &vec);
    void tree_decomp_partial(float);
};
#endif
#pragma once