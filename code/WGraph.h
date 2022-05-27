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
    // graph storage
    AdjList graph;  // adj list of edges, eg. [[v1,v2,v3], [v2,v3,v5], ...]
    AdjList wlist;  // weight of edges, eg. [[1,2,2], [4,2,2], ...]
    uint32_t nsize;
    uint32_t msize;
    uint32_t wsize;
    uint32_t max_w;
    vector<uint32_t> ndegree;  // degree of each node

    // index storage
    vector<vector<uint32_t>> labeling_v;
    vector<vector<uint32_t>> offset;
    vector<vector<uint16_t>> labeling_d;
    vector<vector<uint8_t>> labeling_w;
    vector<vector<pair<uint16_t, uint8_t>>> labeling_dw;

    // tree decomposition
    vector<set<uint32_t>> hgraph;

    WGraph();
    WGraph(string, string, bool);

    //******************** graph processing ************************//
    void read_graph(ifstream&);
    void read_edgelist(ifstream&);
    void save_graph_to_bin(string, string);
    void get_graph_statistics();
    void print_graph();
    void sort_by_degree();
    void check_degree_order();

    //******************** index processing ************************//
    void print_index();
    void get_index_size();
    void check_minimality();

    //******************** tree decomposition ************************//
    void tree_decomp();
    void tree_decomp_partial(float);
    void print_hgraph();

    //********************online solutions************************//
    uint32_t min_distance(vector<uint32_t>&, boost::dynamic_bitset<>&);
    uint32_t dijkstra(AdjList&, uint32_t, uint32_t);
    uint32_t constrained_shortest_distance_naive(uint32_t, uint32_t, uint32_t, double&);
    uint32_t constrained_shortest_distance_dijkstra(uint32_t, uint32_t, uint32_t, double&);
    uint32_t constrained_shortest_distance_plus(uint32_t, uint32_t, uint32_t, double&);

    //************************** query ****************************//
    uint32_t query(uint32_t, uint32_t, uint32_t, double&);
    uint32_t query_skip(uint32_t, uint32_t, uint32_t, uint32_t);

    bool query_while_indexing_vertex_V1(uint32_t, uint32_t, uint16_t, uint8_t);
    bool query_while_indexing_vertex_V8(uint32_t, uint32_t, uint16_t, uint8_t, vector<vector<pair<uint16_t, uint8_t>>>&);

    //******************** index construction ************************//
    void build_index(string, int);
    void vertex_prioritized_indexing_V1(uint32_t, boost::dynamic_bitset<>&, vector<uint8_t>&, boost::dynamic_bitset<>&);
    void vertex_prioritized_indexing_V8(uint32_t, boost::dynamic_bitset<>&, vector<uint8_t>&, boost::dynamic_bitset<>&);
};
#endif
#pragma once