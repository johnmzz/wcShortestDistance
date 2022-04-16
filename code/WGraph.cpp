#include "WGraph.h"

WGraph::WGraph() {
    graph = AdjList();
    wlist = AdjList();
    nsize = 0;
    msize = 0;
    wsize = 0;
    ndegree = vector<uint32_t>();

    labeling_v = vector<vector<uint32_t>>();
    offset = vector<vector<uint32_t>>();
    labeling_d = vector<vector<uint32_t>>();
    labeling_w = vector<vector<uint32_t>>();
}

WGraph::WGraph(string path, string graph, bool isBin) {
    if (isBin) {
        ifstream in(path + graph + string("/") + string("graph.bin"), ios_base::binary);
        if (!in) {
            cout << "Cannot open graph edgelist!\n"
                << endl;
            cout << path + graph + string("/") + string("graph.bin") << endl;
            exit(0);
        }
        read_graph(in);
    }
    else {
        ifstream in(path + graph + string("/") + string("graph.txt"));
        if (!in) {
            cout << "Cannot open graph edgelist!\n"
                << endl;
            cout << path + graph + string("/") + string("graph.txt") << endl;
            exit(0);
        }
        read_edgelist(in);
    }
}

// Read graph from bin format
void WGraph::read_graph(ifstream& input) {
    input.read(reinterpret_cast<char*>(&nsize), 1 * sizeof(uint32_t));
    input.read(reinterpret_cast<char*>(&msize), 1 * sizeof(uint32_t));
    input.read(reinterpret_cast<char*>(&wsize), 1 * sizeof(uint32_t));
    input.read(reinterpret_cast<char*>(&max_w), 1 * sizeof(uint32_t));

    ndegree = vector<uint32_t>(nsize);
    input.read(reinterpret_cast<char*>(&ndegree[0]), nsize * sizeof(uint32_t));

    graph = AdjList(nsize);
    for (uint32_t u = 0; u < nsize; u++) {
        if (ndegree[u] == 0) continue;
        vector<uint32_t> vec(ndegree[u]);
        input.read(reinterpret_cast<char*>(&vec[0]), ndegree[u] * sizeof(uint32_t));
        graph[u].assign(vec.begin(), vec.end());
    }
    wlist = AdjList(nsize);
    for (uint32_t u = 0; u < nsize; u++) {
        if (ndegree[u] == 0) continue;
        vector<uint32_t> vec(ndegree[u]);
        input.read(reinterpret_cast<char*>(&vec[0]), ndegree[u] * sizeof(uint32_t));
        wlist[u].assign(vec.begin(), vec.end());
    }
}

// read graph by edge list (txt)
void WGraph::read_edgelist(ifstream& in) {
    nsize = 0;
    msize = 0;
    wsize = 0;
    max_w = 0;
    string buf, from, to, sig;
    uint32_t from_id, to_id, sig_id;
    uint32_t idx;
    uint32_t cnt = 0;

    while (getline(in, buf)) {
        if (buf.empty()) {
            continue;
        }
        if (cnt == 0) {
            idx = buf.find(" ");
            from = buf.substr(0, idx);
            istringstream(from) >> nsize;
            buf.erase(0, idx + 1);
            idx = buf.find(" ");
            to = buf.substr(0, idx);
            istringstream(to) >> msize;
            buf.erase(0, idx + 1);
            idx = buf.find(" ");
            sig = buf.substr(0, idx);
            istringstream(sig) >> wsize;

            cnt++;
            graph = AdjList(nsize);
            wlist = AdjList(nsize);
            ndegree = vector<uint32_t>(nsize, 0);
        }
        else {
            idx = buf.find(" ");
            from = buf.substr(0, idx);
            istringstream(from) >> from_id;  /// read starting node of the edge
            buf.erase(0, idx + 1);
            idx = buf.find(" ");
            to = buf.substr(0, idx);
            istringstream(to) >> to_id;  /// read starting node of the edge
            buf.erase(0, idx + 1);
            idx = buf.find(" ");
            sig = buf.substr(0, idx);
            istringstream(sig) >> sig_id;  /// read weight of the edge
            graph[from_id].push_back(to_id);
            graph[to_id].push_back(from_id);
            wlist[from_id].push_back(sig_id);
            wlist[to_id].push_back(sig_id);
            ndegree[from_id]++;
            ndegree[to_id]++;
            // msize++;
            if (sig_id > max_w) {
                max_w = sig_id;
            }
        }
    }
}

//////save graph to bin, bin
void WGraph::save_graph_to_bin(string path, string file) {
    ofstream outfile(path + file + "/" + string("graph.bin"), ofstream::binary);
    if (!outfile) {
        cout << "Cannot open file for saving!" << endl;
        cout << path + file + "/" + string("graph.bin") << endl;
        exit(0);
    }

    outfile.write(reinterpret_cast<char*>(&nsize), 1 * sizeof(uint32_t));
    outfile.write(reinterpret_cast<char*>(&msize), 1 * sizeof(uint32_t));
    outfile.write(reinterpret_cast<char*>(&wsize), 1 * sizeof(uint32_t));
    outfile.write(reinterpret_cast<char*>(&max_w), 1 * sizeof(uint32_t));
    outfile.write(reinterpret_cast<char*>(&ndegree[0]), nsize * sizeof(uint32_t));
    //// save graph into bin
    for (uint32_t u = 0; u < nsize; u++) {
        outfile.write(reinterpret_cast<char*>(&graph[u][0]), ndegree[u] * sizeof(uint32_t));
    }

    //// save wlist into bin
    for (uint32_t u = 0; u < nsize; u++) {
        if (ndegree[u] == 0) {
            continue;
        }
        outfile.write(reinterpret_cast<char*>(&wlist[u][0]), ndegree[u] * sizeof(uint32_t));
    }
    cout << "Finish saving binary graph!" << endl;
}

void WGraph::get_graph_statistics() {
    cout << "the number of vertices: " << nsize << endl;
    cout << "the number of edges: " << msize << endl;
    cout << "the number of weights: " << wsize << endl;
    cout << "the maximum weight value: " << max_w << endl;
}

void WGraph::print_graph() {
    // for (uint32_t u = 0; u < nsize; u++) {
    //     cout << u << ": ";
    //     for (uint32_t idx = 0; idx < ndegree[u]; idx++) {
    //         cout << "<" << graph[u][idx] << ", " << wlist[u][idx] << ">, ";
    //     }
    //     cout << endl;
    // }

    for (uint32_t u = 0; u < nsize; u++) {
        cout << u << ": ";
        for (uint32_t idx = 0; idx < ndegree[u]; idx++) {
            cout << "<" << graph[u][idx] << ", " << wlist[u][idx] << ">, ";
        }
        cout << endl;
    }
}

uint32_t WGraph::min_distance(vector<uint32_t>& dist, boost::dynamic_bitset<>& visited) {
    uint32_t min_idx = INVALID_VALUE;
    uint32_t min_dist = INVALID_VALUE;
    for (uint32_t i = 0; i < dist.size(); i++) {
        if (!visited[i] && dist[i] < min_dist) {
            min_idx = i;
            min_dist = dist[i];
        }
    }
    return min_idx;
}

///// O(n^)
uint32_t WGraph::dijkstra(AdjList& newg, uint32_t s, uint32_t t) {
    vector<uint32_t> dist = vector<uint32_t>(newg.size(), INVALID_VALUE);
    boost::dynamic_bitset<> visited(newg.size());
    dist[s] = 0;
    // Find shortest path for all vertices
    for (uint32_t count = 0; count < nsize - 1; count++) {
        // Pick the minimum distance vertex from the set of vertices not yet processed. u is always equal to src in the first iteration.
        uint32_t u = min_distance(dist, visited);
        if (u == INVALID_VALUE) {
            continue;
        }
        // Mark the picked vertex as processed
        visited[u] = 1;
        if (u == t) {
            return dist[t];
        }
        // Update dist value of the adjacent vertices of the picked vertex.
        for (auto w : newg[u])
            if (!visited[w] && dist[u] + 1 < dist[w]) {  // dist[u] + weight for weighted edge.
                dist[w] = dist[u] + 1;
            }
    }
    return dist[t];
}

uint32_t WGraph::constrained_shortest_distance_naive(uint32_t s, uint32_t t, uint32_t r, double& qtime) {
    ////// compute shortest distance under the constraint r
    auto t1 = std::chrono::high_resolution_clock::now();
    /// step 1: filter graph and leave all the weights that are larger than the constraint r
    vector<vector<uint32_t>> newg = vector<vector<uint32_t>>(nsize);
    for (uint32_t u = 0; u < nsize; u++) {
        for (uint32_t nidx = 0; nidx < ndegree[u]; nidx++) {
            if (wlist[u][nidx] >= r) {
                newg[u].push_back(graph[u][nidx]);
            }
        }
    }

    // cout << "newg is: " << endl;
    // for (uint32_t u = 0; u < nsize; u++) {
    //     cout << u << ":";
    //     for (auto w: newg[u]) {
    //         cout << " " << w;
    //     }
    //     cout << endl;
    // }

    /// step 2: compute the shortest distance based on the filtered graph
    queue<uint32_t> myqueue;
    myqueue.emplace(s);
    vector<uint32_t> visited = vector<uint32_t>(nsize, INVALID_VALUE);
    visited[s] = 0;
    while (!myqueue.empty()) {
        uint32_t u = myqueue.front();
        myqueue.pop();
        for (auto w : newg[u]) {
            if (visited[w] == INVALID_VALUE) {  // use visited for both visited info and store disntance
                visited[w] = visited[u] + 1;
                myqueue.emplace(w);
                // cout << "check 1, " << u << ", " << w << endl;
                // cout << "check 2, " << visited[u] << ", " << visited[w] << endl;
                if (w == t) {
                    auto t2 = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
                    qtime += fp_ms.count();
                    // cout << "naive bfs query time: " << fp_ms.count() << " ms." << endl;
                    return visited[w];
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    qtime += fp_ms.count();
    // cout << "naive bfs query time: " << fp_ms.count() << " ms." << endl;
    return INVALID_VALUE;
}

uint32_t WGraph::constrained_shortest_distance_dijkstra(uint32_t s, uint32_t t, uint32_t r, double& qtime) {
    ////// compute shortest distance under the constraint r
    auto t1 = std::chrono::high_resolution_clock::now();
    /// step 1: filter graph and leave all the weights that are larger than the constraint r
    vector<vector<uint32_t>> newg = vector<vector<uint32_t>>(nsize);
    for (uint32_t u = 0; u < nsize; u++) {
        for (uint32_t nidx = 0; nidx < ndegree[u]; nidx++) {
            if (wlist[u][nidx] >= r) {
                newg[u].push_back(graph[u][nidx]);
            }
        }
    }
    // cout << "newg is: " << endl;
    // for (uint32_t u = 0; u < nsize; u++) {
    //     cout << u << ":";
    //     for (auto w: newg[u]) {
    //         cout << " " << w;
    //     }
    //     cout << endl;
    // }

    /// step 2: compute the shortest distance based on the filtered graph
    uint32_t d = dijkstra(newg, s, t);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    qtime += fp_ms.count();
    // cout << "dijkstra query time: " << fp_ms.count() << " ms." << endl;
    return d;
}

uint32_t WGraph::constrained_shortest_distance_plus(uint32_t s, uint32_t t, uint32_t r, double& qtime) {
    ////// compute shortest distance through constrained bfs
    auto t1 = std::chrono::high_resolution_clock::now();
    queue<uint32_t> myqueue;
    myqueue.emplace(s);
    vector<uint32_t> visited = vector<uint32_t>(nsize, INVALID_VALUE);
    visited[s] = 0;
    while (!myqueue.empty()) {
        uint32_t u = myqueue.front();
        myqueue.pop();
        for (uint32_t idx = 0; idx < ndegree[u]; idx++) {
            uint32_t w = graph[u][idx];
            if (visited[w] == INVALID_VALUE && wlist[u][idx] >= r) {
                visited[w] = visited[u] + 1;
                myqueue.emplace(w);
                if (w == t) {
                    auto t2 = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
                    qtime += fp_ms.count();
                    // cout << "constrained bfs query time: " << fp_ms.count() << " ms." << endl;
                    return visited[w];
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    qtime += fp_ms.count();
    // cout << "constrained bfs query time: " << fp_ms.count() << " ms." << endl;
    return INVALID_VALUE;
}

//// query algorithms based on the index
uint32_t WGraph::query_with_index_vertex(uint32_t s, uint32_t t, uint32_t r, double& qtime) {
    auto t1 = std::chrono::high_resolution_clock::now();
    uint32_t min_dist = INVALID_VALUE;
    uint32_t i = 0, j = 0;

    while (i < labeling_v[s].size() && j < labeling_v[t].size()) {
        uint32_t curr_u = labeling_v[s][i];
        uint32_t curr_w = labeling_v[t][j];
        if (curr_u < curr_w) {
            i++;
        }
        else if (curr_u > curr_w) {
            j++;
        }
        else {
            uint32_t d1 = INVALID_VALUE, d2 = INVALID_VALUE;
            for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
                if (labeling_w[s][p] >= r) {
                    d1 = labeling_d[s][p];
                    break;
                }
            }
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_w[t][q] >= r) {
                    d2 = labeling_d[t][q];
                    break;
                }
            }

            if (d1 != INVALID_VALUE && d2 != INVALID_VALUE) {
                if (d1 + d2 < min_dist) {
                    min_dist = d1 + d2;
                }
            }
            i++;
            j++;
        }
    }

    while (i < labeling_v[s].size()) {
        uint32_t curr_u = labeling_v[s][i];
        if (curr_u == t) {
            uint32_t d1 = INVALID_VALUE;
            for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
                if (labeling_w[s][p] >= r) {
                    d1 = labeling_d[s][p];
                    break;
                }
            }
            if (d1 != INVALID_VALUE && d1 < min_dist) {
                min_dist = d1;
            }
        }
        i++;
    }

    while (j < labeling_v[t].size()) {
        uint32_t curr_w = labeling_v[t][j];
        if (curr_w == s) {
            uint32_t d2 = INVALID_VALUE;
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_w[t][q] >= r) {
                    d2 = labeling_d[t][q];
                    break;
                }
            }
            if (d2 != INVALID_VALUE && d2 < min_dist) {
                min_dist = d2;
            }
        }
        j++;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    qtime += fp_ms.count();
    return min_dist;
}

uint32_t WGraph::query_with_index_vertex_skip(uint32_t s, uint32_t t, uint32_t r, uint32_t dpos) {
    uint32_t min_dist = INVALID_VALUE;
    uint32_t i = 0, j = 0;

    while (i < labeling_v[s].size() && j < labeling_v[t].size()) {
        uint32_t curr_u = labeling_v[s][i];
        uint32_t curr_w = labeling_v[t][j];
        if (curr_u < curr_w) {
            i++;
        }
        else if (curr_u > curr_w) {
            j++;
        }
        else {
            uint32_t d1 = INVALID_VALUE, d2 = INVALID_VALUE;
            for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
                if (p == dpos) {
                    continue;
                }
                if (labeling_w[s][p] >= r) {
                    d1 = labeling_d[s][p];
                    break;
                }
            }
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_w[t][q] >= r) {
                    d2 = labeling_d[t][q];
                    break;
                }
            }

            if (d1 != INVALID_VALUE && d2 != INVALID_VALUE) {
                if (d1 + d2 < min_dist) {
                    min_dist = d1 + d2;
                }
            }
            i++;
            j++;
        }
    }

    while (i < labeling_v[s].size()) {
        uint32_t curr_u = labeling_v[s][i];
        if (curr_u == t) {
            uint32_t d1 = INVALID_VALUE;
            for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
                if (p == dpos) {
                    continue;
                }
                if (labeling_w[s][p] >= r) {
                    d1 = labeling_d[s][p];
                    break;
                }
            }
            if (d1 != INVALID_VALUE && d1 < min_dist) {
                min_dist = d1;
            }
        }
        i++;
    }

    while (j < labeling_v[t].size()) {
        uint32_t curr_w = labeling_v[t][j];
        if (curr_w == s) {
            uint32_t d2 = INVALID_VALUE;
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_w[t][q] >= r) {
                    d2 = labeling_d[t][q];
                    break;
                }
            }
            if (d2 != INVALID_VALUE && d2 < min_dist) {
                min_dist = d2;
            }
        }
        j++;
    }
    return min_dist;
}

uint32_t WGraph::query_while_indexing_vertex(uint32_t s, uint32_t t, uint32_t r) {
    uint32_t min_dist = INVALID_VALUE;
    uint32_t i = 0, j = 0;

    uint32_t isize = labeling_v[s].size();
    uint32_t jsize = labeling_v[t].size();
    while (i < isize && j < jsize) {
        if (labeling_v[s][i] < labeling_v[t][j]) {
            i++;
        }
        else if (labeling_v[s][i] > labeling_v[t][j]) {
            j++;
        }
        else {
            uint32_t l1 = lower_bound(labeling_w[s].begin() + offset[s][i], labeling_w[s].begin() + offset[s][i + 1], r) - labeling_w[s].begin();
            if (l1 != offset[s][i + 1]) {
                uint32_t l2 = lower_bound(labeling_w[t].begin() + offset[t][j], labeling_w[t].begin() + offset[t][j + 1], r) - labeling_w[t].begin();
                if (l2 != offset[t][j + 1]) {
                    if (labeling_d[s][l1] + labeling_d[t][l2] < min_dist) {
                        min_dist = labeling_d[s][l1] + labeling_d[t][l2];
                    }
                }
            }
            i++;
            j++;
        }
    }

    // while (i < isize) {
    //     if (labeling_v[s][i] == t) {
    //         uint32_t l1 = lower_bound(labeling_w[s].begin() + offset[s][i], labeling_w[s].begin() + offset[s][i + 1], r) - labeling_w[s].begin();
    //         if (l1 != offset[s][i + 1] && labeling_d[s][l1] < min_dist) {
    //             min_dist = labeling_d[s][l1];
    //         }
    //     }
    //     i++;
    // }

    // while (j < jsize) {
    //     if (labeling_v[t][j] == s) {
    //         uint32_t l2 = lower_bound(labeling_w[t].begin() + offset[t][j], labeling_w[t].begin() + offset[t][j + 1], r) - labeling_w[t].begin();
    //         if (l2 != offset[t][j + 1] && labeling_d[t][l2] < min_dist) {
    //             min_dist = labeling_d[t][l2];
    //         }
    //     }
    //     j++;
    // }
    return min_dist;
}

////
uint32_t WGraph::query_while_indexing_vertex_V2(uint32_t s, uint32_t t, uint32_t r, vector<vector<pair<uint32_t, uint32_t>>>& u_label) {
    uint32_t min_dist = INVALID_VALUE;

    for (uint32_t j = 0; j < labeling_v[t].size(); j++) {
        uint32_t v = labeling_v[t][j];

        if (v > s) break;

        if (u_label[v].empty()) continue;

        bool updated_flag = false;
        for (uint32_t p = offset[t][j]; p < offset[t][j + 1]; p++) {
            if (updated_flag == true) {
                break;
            }

            uint32_t d = labeling_d[t][p];
            uint32_t w = labeling_w[t][p];
            if (w < r) {
                continue;
            }

            for (uint32_t i = 0; i < u_label[v].size(); i++) {
                if (u_label[v][i].second >= r) {
                    if ((u_label[v][i].first + d) < min_dist) {
                        min_dist = u_label[v][i].first + d;
                        // update turning_point_r[t] =
                        updated_flag = true;
                        break;
                    }
                }
            }
        }
    }
    return min_dist;
}

uint32_t WGraph::query_while_indexing_vertex_V3(uint32_t s, uint32_t t, uint32_t r, vector<vector<uint32_t>>& u_label_d, vector<vector<uint32_t>>& u_label_w) {
    // cout << "Round: s = " << s << ", t = " << t << ", r = " << r << ":" << endl;
    uint32_t min_dist = INVALID_VALUE;

    for (uint32_t j = 0; j < labeling_v[t].size(); j++) {
        uint32_t v = labeling_v[t][j];

        if (v > s) break;   // TODO: break if labeling_v[t][j] > s, because all subsequent v will be larger than current v

        if (u_label_d[v].empty()) continue;

        bool updated_flag = false;
        for (uint32_t p = offset[t][j]; p < offset[t][j + 1]; p++) {
            if (updated_flag == true) {
                break;
            }

            uint32_t d = labeling_d[t][p];
            uint32_t w = labeling_w[t][p];
            if (w < r) {
                continue;
            }
            /////TODO: change it to lowerbound, using binary search instead of linear search
            /////// change it to lowerbound or uperbound function
            for (uint32_t i = 0; i < u_label_w[v].size(); i++) {
                if (u_label_w[v][i] >= r) {
                    // cout << "index = " << i << ", value = " << u_label_w[v][i] << endl;
                    if ((u_label_d[v][i] + d) < min_dist) {
                        min_dist = u_label_d[v][i] + d;
                        // update turning_point_r[t] =
                        updated_flag = true;
                        break;
                    }
                }
            }
            /*
            vector<uint32_t>::iterator low = lower_bound(u_label_w[v].begin(), u_label_w[v].end(), r);
            if (low != u_label_w[v].end()) {
                // cout << "lower bound = " << (low - u_label_w[v].begin()) << ", value = " << u_label_w[v][low-u_label_w[v].begin()] << endl;
                for (uint32_t i = (low - u_label_w[v].begin()); i != u_label_w[v].end() - u_label_w[v].begin(); i++) {
                    if ((u_label_d[v][i] + d) < min_dist) {
                        min_dist = u_label_d[v][i] + d;
                        // update turning_point_r[t] =
                        updated_flag = true;
                        break;
                    }
                }
            }
            */
            /*
            vector<uint32_t>::iterator low = lower_bound(u_label_w[v].begin(), u_label_w[v].end(), r);
            while (low != u_label_w[v].end()) {
                if ((u_label_d[v][low-u_label_w[v].begin()] + d) < min_dist) {
                    min_dist = u_label_d[v][low-u_label_w[v].begin()] + d;
                    updated_flag = true;
                    break;
                }
                else {
                    r++;
                    low = lower_bound(u_label_w[v].begin(), u_label_w[v].end(), r);
                }
            }
            */
        }
    }
    return min_dist;
}

uint32_t WGraph::query_while_indexing_vertex_V4(uint32_t s, uint32_t t, uint32_t r, vector<vector<pair<uint32_t, uint32_t>>>& u_label, vector<uint32_t>& turning_point_d, vector<uint32_t>& turning_point_w) {
    uint32_t min_dist = INVALID_VALUE;

    for (uint32_t j = 0; j < labeling_v[t].size(); j++) {
        uint32_t v = labeling_v[t][j];

        if (v > s) break;

        if (u_label[v].empty()) continue;

        bool updated_flag = false;
        for (uint32_t p = offset[t][j]; p < offset[t][j + 1]; p++) {
            if (updated_flag == true) {
                break;
            }

            uint32_t d = labeling_d[t][p];
            uint32_t w = labeling_w[t][p];
            if (w < r) {
                continue;
            }

            for (uint32_t i = 0; i < u_label[v].size(); i++) {
                if (u_label[v][i].second >= r) {
                    if ((u_label[v][i].first + d) < min_dist) {
                        min_dist = u_label[v][i].first + d;
                        turning_point_w[t] = min(w, u_label[v][i].second);
                        turning_point_d[t] = min_dist;
                        updated_flag = true;
                        break;
                    }
                }
            }
        }
    }
    return min_dist;
}

uint32_t WGraph::query_while_indexing_vertex_V5(uint32_t s, uint32_t t, uint32_t r, vector<vector<pair<uint32_t, uint32_t>>>& u_label, vector<uint32_t>& turning_point_d, vector<uint32_t>& turning_point_w, vector<uint32_t>& visited_tp) {
    uint32_t min_dist = INVALID_VALUE;

    for (uint32_t j = 0; j < labeling_v[t].size(); j++) {
        uint32_t v = labeling_v[t][j];

        if (v > s) break;

        if (u_label[v].empty()) continue;

        bool updated_flag = false;
        for (uint32_t p = offset[t][j]; p < offset[t][j + 1]; p++) {
            if (updated_flag == true) {
                break;
            }

            uint32_t d = labeling_d[t][p];
            uint32_t w = labeling_w[t][p];
            if (w < r) {
                continue;
            }

            for (uint32_t i = 0; i < u_label[v].size(); i++) {
                if (u_label[v][i].second >= r) {
                    if ((u_label[v][i].first + d) < min_dist) {
                        min_dist = u_label[v][i].first + d;
                        turning_point_w[t] = min(w, u_label[v][i].second);
                        turning_point_d[t] = min_dist;
                        visited_tp.push_back(t);
                        updated_flag = true;
                        break;
                    }
                }
            }
        }
    }
    return min_dist;
}

// uint32_t WGraph::query_while_indexing_vertex(uint32_t s, uint32_t t, uint32_t r) {
//     uint32_t min_dist = INVALID_VALUE;
//     uint32_t i = 0, j = 0;
//
//     uint32_t isize = labeling_v[s].size();
//     uint32_t jsize = labeling_v[t].size();
//     while (i < isize && j < jsize) {
//         if (labeling_v[s][i] < labeling_v[t][j]) {
//             i++;
//         } else if (labeling_v[s][i] > labeling_v[t][j]) {
//             j++;
//         } else {
//             uint32_t d1 = INVALID_VALUE, d2 = INVALID_VALUE;
//             for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
//                 if (labeling_w[s][p] >= r) {
//                     d1 = labeling_d[s][p];
//                     break;
//                 }
//             }
//             for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
//                 if (labeling_w[t][q] >= r) {
//                     d2 = labeling_d[t][q];
//                     break;
//                 }
//             }
//
//             if (d1 != INVALID_VALUE && d2 != INVALID_VALUE) {
//                 if (d1 + d2 < min_dist) { min_dist = d1 + d2; }
//             }
//             i++;
//             j++;
//         }
//     }
//
//     while (i < labeling_v[s].size()) {
//         if (labeling_v[s][i] == t) {
//             uint32_t d1 = INVALID_VALUE;
//             for (uint32_t p = offset[s][i]; p < offset[s][i + 1]; p++) {
//                 if (labeling_w[s][p] >= r) {
//                     d1 = labeling_d[s][p];
//                     break;
//                 }
//             }
//             if (d1 != INVALID_VALUE && d1 < min_dist) {
//                 min_dist = d1;
//             }
//         }
//         i++;
//     }
//
//     while (j < labeling_v[t].size()) {
//         if (labeling_v[t][j] == s) {
//             uint32_t d2 = INVALID_VALUE;
//             for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
//                 if (labeling_w[t][q] >= r) {
//                     d2 = labeling_d[t][q];
//                     break;
//                 }
//             }
//             if (d2 != INVALID_VALUE && d2 < min_dist) {
//                 min_dist = d2;
//             }
//         }
//         j++;
//     }
//     return min_dist;
// }

uint32_t WGraph::query_for_pll(uint32_t s, uint32_t t) {
    uint32_t min_dist = INVALID_VALUE;
    uint32_t i = 0, j = 0;

    while (i < labeling_v[s].size() && j < labeling_v[t].size()) {
        uint32_t curr_u = labeling_v[s][i];
        uint32_t curr_w = labeling_v[t][j];
        if (curr_u < curr_w) {
            i++;
        }
        else if (curr_u > curr_w) {
            j++;
        }
        else {
            uint32_t d1 = labeling_d[s][i];
            uint32_t d2 = labeling_d[t][j];
            if (d1 + d2 < min_dist) {
                min_dist = d1 + d2;
            }
            i++;
            j++;
        }
    }

    while (i < labeling_v[s].size()) {
        if (labeling_v[s][i] == t) {
            uint32_t d1 = labeling_d[s][i];
            if (d1 < min_dist) {
                min_dist = d1;
            }
        }
        i++;
    }

    while (j < labeling_v[t].size()) {
        if (labeling_v[t][j] == s) {
            uint32_t d2 = labeling_d[t][j];
            if (d2 < min_dist) {
                min_dist = d2;
            }
        }
        j++;
    }
    return min_dist;
}

void WGraph::build_index(string type, int threshold) {
    labeling_v = vector<vector<uint32_t>>(nsize);
    offset = vector<vector<uint32_t>>(nsize);
    labeling_d = vector<vector<uint32_t>>(nsize);
    labeling_w = vector<vector<uint32_t>>(nsize);

    for (uint32_t u = 0; u < nsize; u++) {
        offset[u].push_back(0);
    }

    boost::dynamic_bitset<> updated = boost::dynamic_bitset<>(nsize);
    vector<uint32_t> visited_r = vector<uint32_t>(nsize, 0);
    boost::dynamic_bitset<> this_visited_flag = boost::dynamic_bitset<>(nsize);
    vector<uint32_t> visited_d = vector<uint32_t>(nsize, INVALID_VALUE);

    vector<uint32_t> turning_point_w = vector<uint32_t>(nsize, 0);
    vector<uint32_t> turning_point_d = vector<uint32_t>(nsize, 0);

    clock_t start, end;
    start = clock();
    if (type == "V0") {
        // clock_t sss,ttt;
        // sss = clock();
        for (uint32_t u = 0; u < nsize; u++) {
            vertex_prioritized_indexing(u, updated, visited_r);
            if (u % 1000 == 0) {
                cout << u << " vertices finished!" << endl;
            }
        }
    }
    else if (type == "V1") {
        for (uint32_t u = 0; u < nsize; u++) {
            vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V2") {
        cout << "use V2 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V2(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V3") {
        cout << "use V3 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V3(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V4") {
        cout << "use V4 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V4(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V5") {
        cout << "use V5 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V5(u, updated, visited_r, this_visited_flag, turning_point_d, turning_point_w);
            }
            else {
                vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V6") {
        tree_decomposition();
        cout << "use V6 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                uint32_t u_ = v_degree[nsize-1-u].first;
                vertex_prioritized_indexing_V4(u, updated, visited_r, this_visited_flag);
            }
            else {
                uint32_t u_ = v_degree[nsize-1-u].first;
                vertex_prioritized_indexing_plus(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "PLL") {
        for (uint32_t u = 0; u < nsize; u++) {
            pruned_landmark_labeling(u, visited_d);
            if (u % 1000 == 0) {
                cout << u << " vertices finished!" << endl;
            }
        }
    }
    end = clock();
    cout << "build index: " << (float)(end - start) / CLOCKS_PER_SEC << " s" << endl;
}

// updated initialized to 000000...
// visited_r initialized to [0,0,0,...]
void WGraph::vertex_prioritized_indexing(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint32_t>& visited_r) {
    vector<uint32_t> visited_vertex;
    priority_queue<pair<uint32_t, uint32_t>, vector<pair<uint32_t, uint32_t>>, descWeightSort> Q1, Q2;  // Asc or Desc??
    Q1.emplace(make_pair(u, INVALID_VALUE));
    // updated[u] = 1;
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);  // use set??
    uint32_t curr_d = 0;
    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.top();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                // visited_r[v] = r;
                // cout << "u, v, r, curr_d, querydist: " << u << ", " << v << ", " << r << ", " << curr_d << ", " << query_while_indexing_dist(u, v, r) << endl;
                if (query_while_indexing_vertex(u, v, r) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                    // cout << "<" << u << "," << curr_d << "," << r << "> to L(" << v << ")" << endl;
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        continue;
                    }
                    Q2.emplace(w, new_r);
                    // cout << curr_d << ": " << w << new_r << endl;
                    visited_r[w] = new_r;
                    visited_vertex.push_back(w);
                }
            }
            curr_d++;
        }
        else {
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.top();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                // cout << "u, v, r, curr_d, querydist: " << u << ", " << v << ", " << r << ", " << curr_d << ", " << query_while_indexing_dist(u, v, r) << endl;
                if (query_while_indexing_vertex(u, v, r) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                    // cout << "<" << u << "," << curr_d << "," << r << "> to L(" << v << ")" << endl;
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        continue;
                    }
                    Q1.emplace(w, new_r);
                    visited_r[w] = new_r;
                    // cout << curr_d << ": " << w << new_r << endl;
                    visited_vertex.push_back(w);
                }
            }
            curr_d++;
        }
    }

    /// reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
}

// change based on vertex_prioritized_indexing but remove the priority queue
void WGraph::vertex_prioritized_indexing_plus(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint32_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint32_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE));
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);
    uint32_t curr_d = 0;
    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            // unordered_set<uint32_t> this_visited;
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex(u, v, r) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    visited_vertex.push_back(w);
                }
            }

            for (auto it : this_visited) {
                Q2.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
        else {
            // unordered_set<uint32_t> this_visited;
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex(u, v, r) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q1.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
    }

    /// reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
}

// V2, same as V1 but changed query method
void WGraph::vertex_prioritized_indexing_V2(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint32_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint32_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE));
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);
    uint32_t curr_d = 0;

    vector<vector<pair<uint32_t, uint32_t>>> u_label = vector<vector<pair<uint32_t, uint32_t>>>(nsize, vector<pair<uint32_t, uint32_t>>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint32_t d = labeling_d[u][p];
            uint32_t r = labeling_w[u][p];
            u_label[v].emplace_back(make_pair(d, r));
        }
    }
    u_label[u].emplace_back(make_pair(0, INVALID_VALUE));

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex_V2(u, v, r, u_label) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    // this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q2.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
        else {
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex_V2(u, v, r, u_label) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    // this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q1.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
    }
    // reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
}

// V3, changed data stucture of mapping (labels of u)
void WGraph::vertex_prioritized_indexing_V3(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint32_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint32_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE));
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);
    uint32_t curr_d = 0;

    vector<vector<uint32_t>> u_label_d = vector<vector<uint32_t>>(nsize, vector<uint32_t>());       // TODO: change size from nsize to u+1?
    vector<vector<uint32_t>> u_label_w = vector<vector<uint32_t>>(nsize, vector<uint32_t>());

    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint32_t d = labeling_d[u][p];
            uint32_t r = labeling_w[u][p];
            u_label_d[v].emplace_back(d);
            u_label_w[v].emplace_back(r);
        }
    }
    u_label_d[u].emplace_back(0);
    u_label_w[u].emplace_back(INVALID_VALUE);

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex_V3(u, v, r, u_label_d, u_label_w) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q2.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
        else {
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                if (query_while_indexing_vertex_V3(u, v, r, u_label_d, u_label_w) <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q1.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
    }
    // reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
}


// V4, optimization that stores turning points of r, initialized each round (each u)
void WGraph::vertex_prioritized_indexing_V4(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint32_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint32_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE));
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);
    uint32_t curr_d = 0;

    vector<uint32_t> turning_point_w = vector<uint32_t>(nsize, 0);
    vector<uint32_t> turning_point_d = vector<uint32_t>(nsize, 0);

    vector<vector<pair<uint32_t, uint32_t>>> u_label = vector<vector<pair<uint32_t, uint32_t>>>(nsize, vector<pair<uint32_t, uint32_t>>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint32_t d = labeling_d[u][p];
            uint32_t r = labeling_w[u][p];
            u_label[v].emplace_back(make_pair(d, r));
        }
    }
    u_label[u].emplace_back(make_pair(0, INVALID_VALUE));

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                uint32_t query_d;

                if (r < turning_point_w[v]) {
                    query_d = turning_point_d[v];
                }
                else {
                    query_d = query_while_indexing_vertex_V4(u, v, r, u_label, turning_point_d, turning_point_w);
                }
                if (query_d <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q2.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
        else {
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                uint32_t query_d;

                if (r < turning_point_w[v]) {
                    query_d = turning_point_d[v];
                }
                else {
                    query_d = query_while_indexing_vertex_V4(u, v, r, u_label, turning_point_d, turning_point_w);
                }
                if (query_d <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q1.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
    }
    // reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
}

// V5, different implementation of V4
void WGraph::vertex_prioritized_indexing_V5(uint32_t u,
                                            boost::dynamic_bitset<>& updated,
                                            vector<uint32_t>& visited_r,
                                            boost::dynamic_bitset<>& this_visited_flag,
                                            vector<uint32_t>& turning_point_w,
                                            vector<uint32_t>& turning_point_d) {
    vector<uint32_t> visited_tp;
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint32_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE));
    visited_r[u] = INVALID_VALUE;
    visited_vertex.push_back(u);
    uint32_t curr_d = 0;

    vector<vector<pair<uint32_t, uint32_t>>> u_label = vector<vector<pair<uint32_t, uint32_t>>>(nsize, vector<pair<uint32_t, uint32_t>>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint32_t d = labeling_d[u][p];
            uint32_t r = labeling_w[u][p];
            u_label[v].emplace_back(make_pair(d, r));
        }
    }
    u_label[u].emplace_back(make_pair(0, INVALID_VALUE));

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint32_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                uint32_t query_d;

                if (r < turning_point_w[v]) {
                    query_d = turning_point_d[v];
                }
                else {
                    query_d = query_while_indexing_vertex_V5(u, v, r, u_label, turning_point_d, turning_point_w, visited_tp);
                }
                if (query_d <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q2.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
        else {
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint32_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint32_t r = p.second;
                uint32_t query_d;

                if (r < turning_point_w[v]) {
                    query_d = turning_point_d[v];
                }
                else {
                    query_d = query_while_indexing_vertex_V5(u, v, r, u_label, turning_point_d, turning_point_w, visited_tp);
                }
                if (query_d <= curr_d) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_d[v].push_back(curr_d);
                        labeling_w[v].push_back(r);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint32_t new_r = min(r, wlist[v][idx]);
                    if (new_r <= visited_r[w]) {
                        if (new_r != 0 || visited_r[w] != 0) {
                            continue;
                        }
                    }
                    visited_r[w] = new_r;
                    if (this_visited_flag[w] & 1) {
                        continue;
                    }
                    this_visited.push_back(w);
                    this_visited_flag[w] = 1;
                    visited_vertex.push_back(w);
                }
            }
            for (auto it : this_visited) {
                Q1.emplace(make_pair(it, visited_r[it]));
                this_visited_flag[it] = 0;
            }
            curr_d++;
        }
    }
    // reset flag vectors
    for (auto w : visited_vertex) {
        updated[w] = 0;
        visited_r[w] = 0;
    }
    for (auto w : visited_tp) {
        turning_point_d[w] = 0;
        turning_point_w[w] = 0;
    }
}

// index without weightings
void WGraph::pruned_landmark_labeling(uint32_t u, vector<uint32_t>& visited_d) {
    vector<uint32_t> visited_vertex;
    queue<uint32_t> Q;
    Q.emplace(u);
    visited_d[u] = 0;
    visited_vertex.push_back(u);
    while (!Q.empty()) {
        uint32_t v = Q.front();
        Q.pop();
        // cout << "u, v, d: " << u << ", " << v << ", " <<  visited_d[v] << ", " << query_for_pll(u,v) << endl;
        if (v != u && query_for_pll(u, v) <= visited_d[v]) {  //// need to implement
            continue;
        }
        else {
            labeling_v[v].push_back(u);
            // offset[v].push_back(offset[v].back() + 1);
            labeling_d[v].push_back(visited_d[v]);
        }
        for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
            uint32_t w = graph[v][idx];
            if (w <= u) {
                continue;
            }  /// only explore the vertex that ranks lower than the starting vertex u
            if (visited_d[w] != INVALID_VALUE) {
                continue;
            }
            visited_d[w] = visited_d[v] + 1;
            visited_vertex.push_back(w);
            Q.emplace(w);
        }
    }
    /// reset flag vectors
    for (auto w : visited_vertex) {
        visited_d[w] = INVALID_VALUE;
    }
}

void WGraph::print_index() {
    for (uint32_t u = 0; u < nsize; u++) {
        cout << u << ": ";
        for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
            for (uint32_t j = offset[u][i]; j < offset[u][i + 1]; j++) {
                cout << "<" << labeling_v[u][i] << "," << labeling_d[u][j] << "," << labeling_w[u][j] << ">, ";
            }
        }
        cout << endl;
    }
}

void WGraph::print_index_pll() {
    for (uint32_t u = 0; u < nsize; u++) {
        cout << u << ": ";
        for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
            cout << "<" << labeling_v[u][i] << "," << labeling_d[u][i] << ">, ";
        }
        cout << endl;
    }
}

void WGraph::get_index_size() {
    cout << "--------------Details of Index Occupation-------------" << endl;
    uint64_t cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
    for (uint32_t u = 0; u < nsize; u++) {
        cnt1 += labeling_v[u].size();
        cnt2 += offset[u].size();
        cnt3 += labeling_d[u].size();
        cnt4 += labeling_w[u].size();
    }

    cout << "labeling v, offset, labeling d, labeling w: " << cnt1 << ", " << cnt2 << ", " << cnt3 << ", " << cnt4 << endl;
    cout << "total number of unsigned ints in the index: " << cnt1 + cnt2 + cnt3 + cnt4 << endl;
}

void WGraph::check_minimality() {
    cout << "check minimality " << endl;
    for (uint32_t u = 0; u < nsize; u++) {
        for (uint32_t wpos = 0; wpos < labeling_v[u].size(); wpos++) {
            for (uint32_t dpos = offset[u][wpos]; dpos < offset[u][wpos + 1]; dpos++) {
                uint32_t d1 = query_with_index_vertex_skip(u, labeling_v[u][wpos], labeling_w[u][dpos], dpos);
                if (d1 == labeling_d[u][dpos]) {
                    cout << "failed at u, v, r, d: " << u << ", " << labeling_v[u][wpos] << ", " << labeling_w[u][dpos] << ", " << labeling_d[u][dpos] << endl;
                }
            }
        }
    }
}

void WGraph::preprocess(int num_intervals) {
    if (num_intervals > wsize) {
        cout << "Number of intervals > size of label set." << endl;
        exit(0);
    }

    std::set<uint32_t> weight_set;
    for (uint32_t i = 0; i < nsize; i++) {
        for (uint32_t j = 0; j < wlist[i].size(); j++) {
            weight_set.insert(wlist[i][j]);
        }
    }
    vector<uint32_t> weight_values;
    std::copy(weight_set.begin(), weight_set.end(), std::back_inserter(weight_values));
    /*
    cout << "Weights" << endl;
    for(auto it = weight_values.begin(); it != weight_values.end(); it++)
    {
        cout << *it << ", ";
    }
    cout << endl;
    */

    vector<vector<uint32_t>> out_vec;
    size_t length = weight_set.size()/num_intervals;
    size_t remain = weight_set.size()%num_intervals;
    size_t begin = 0;
    size_t end = 0;

    for (uint32_t i = 0; i < num_intervals; i++) {
        if (remain > 0) {
            end += (i < remain) ? (length + 1) : length;
        }
        else {
            end += length;
        }
        vector<uint32_t> split(weight_values.begin() + begin, weight_values.begin() + end);
        out_vec.push_back(split);
        begin = end;
    }

    /*
    for (auto i : out_vec) {
        for (auto j : i) {
            cout << j << ",";
        }
        cout << endl;
    }
    */
}

void WGraph::print_hgraph() {
    cout << "h_graph = " << endl;
    for (uint32_t i = 0; i < hgraph.size(); i++) {
        cout << i << ": [";
        for (auto j : hgraph[i]) {
            cout << j << ",";
        }
        cout << "]" << endl;
    }
}

void WGraph::print_vdegree() {
    cout << "v_degree = " << endl;
    for (uint32_t i = 0; i < v_degree.size(); i++) {
        cout << i << ": " << "(" << v_degree[i].first << "," << v_degree[i].second << ")" << endl;
    }
}

void WGraph::print_vindex() {
    cout << "index = [";
    for (auto i : v_index) {
        cout << i << ",";
    }
    cout << "]" << endl;
}

void WGraph::print_beg_end() {
    cout << "beg = [";
    for (auto i : beg) {
        cout << i << ",";
    }
    cout << "]" << endl;

    cout << "end = [";
    for (auto i : end) {
        cout << i << ",";
    }
    cout << "]" << endl;
}

void WGraph::print_decomp_order() {
    cout << "order = [";
    for (auto i : decomp_order) {
        cout << i << ",";
    }
    cout << "]" << endl;
}

void WGraph::tree_decomposition() {
    // *********** copy graph ************
    for (uint32_t i = 0; i < graph.size(); i++) {
        set<uint32_t> s(graph[i].begin(), graph[i].end());
        hgraph.push_back(s);
    }
    // print_hgraph();

    // ************ initialize degree and sort ************
    for (uint32_t i = 0; i < nsize; i++) {
        // cout << i << ": " << ndegree[i] << endl;
        pair<uint32_t, uint32_t> p(i, ndegree[i]);
        v_degree.push_back(p);
    }
    // print_vdegree();
    sort(v_degree.begin(), v_degree.end(), [](auto &left, auto &right) { return left.second < right.second; });
    // print_vdegree();

    // ************ initialize index to degree ************
    v_index = vector<uint32_t>(nsize);
    for (int i = 0; i < nsize; i++) {
        v_index[v_degree[i].first] = i;
    }
    // print_vindex();

    // ************ initialize begin/end ************
    beg = vector<uint32_t>(nsize, INF);
    end = vector<uint32_t>(nsize, INF);
    for (uint32_t i = 0; i < nsize; i++) {
        if (i == 0 || v_degree[i-1].second != v_degree[i].second) {
            beg[v_degree[i].second] = i;
        }
        if (i == (nsize-1) || v_degree[i+1].second != v_degree[i].second) {
            end[v_degree[i].second] = i;
        }
    }
    // print_beg_end();

    // ************ start decomposition ************
    cout << "Start tree decomposition..." << endl;
    for (uint32_t i = 0; i < nsize; i++) {
        uint32_t u = v_degree[i].first;

        for (auto v = hgraph[u].begin(); v != hgraph[u].end(); v++) {
            for (auto w = hgraph[u].begin(); w != hgraph[u].end(); w++) {
                if ((*v != *w) && (*v < *w)) {
                    // cout << *v << "," << *w << endl;
                    hgraph[*v].insert(*w);
                    hgraph[*w].insert(*v);
                }
            }
        }

        for (auto v : hgraph[u]) {
            hgraph[v].erase(u);
            uint32_t idx = v_index[v];
            uint32_t old = v_degree[idx].second;
            uint32_t new_ = hgraph[v].size();

            if (new_ < old) {
                uint32_t pos = max(i+1, beg[old]);

                v_degree[idx] = v_degree[pos];
                v_index[v_degree[pos].first] = idx;

                v_degree[pos] = pair<uint32_t, uint32_t>(v, new_);
                v_index[v] = pos;

                beg[old] = pos + 1;
                if (beg[old] > end[old]) {
                    beg[old] = INF;
                    end[old] = INF;
                }

                if (end[new_] == INF) {
                    beg[new_] = pos;
                    end[new_] = pos;
                } else {
                    end[new_] = pos;
                }
            }
            else if (new_ > old) {
                for (uint32_t j = old; j < new_; j++) {
                    uint32_t pos = end[j];

                    v_degree[idx] = v_degree[pos];
                    v_index[v_degree[pos].first] = idx;

                    v_degree[pos] = pair<uint32_t, uint32_t>(v, j+1);
                    v_index[v] = pos;

                    end[j] = pos - 1;
                    if (end[j] < beg[j]) {
                        beg[j] = INF;
                        end[j] = INF;
                    }

                    if (beg[j+1] == INF) {
                        beg[j+1] = pos;
                        end[j+1] = pos;
                    } else {
                        beg[j+1] = pos;
                    }

                    idx = v_index[v];
                }
            }
        }
    }

    // ************ Finished decomposition ************
    // cout << "Finished decomposition, result: " << endl;
    // print_vdegree();
}

void WGraph::tree_decomposition_2() {
    // *********** copy graph ************
    for (uint32_t i = 0; i < graph.size(); i++) {
        set<uint32_t> s(graph[i].begin(), graph[i].end());
        hgraph.push_back(s);
    }
    //print_hgraph();

    // ************ initialize degree and sort ************
    for (uint32_t i = 0; i < nsize; i++) {
        // cout << i << ": " << ndegree[i] << endl;
        pair<uint32_t, uint32_t> p(i, ndegree[i]);
        v_degree.push_back(p);
    }
    // print_vdegree();
    sort(v_degree.begin(), v_degree.end(), [](auto &left, auto &right) { return left.second < right.second; });
    //print_vdegree();

    // ************ initialize index to degree ************
    v_index = vector<uint32_t>(nsize);
    for (int i = 0; i < nsize; i++) {
        v_index[v_degree[i].first] = i;
    }
    //print_vindex();

    // ************ start decomposition ************
    cout << "Start tree decomposition..." << endl;
    for (uint32_t i = 0; i < nsize; i++) {
        uint32_t u = v_degree[0].first;

        for (auto v = hgraph[u].begin(); v != hgraph[u].end(); v++) {
            for (auto w = hgraph[u].begin(); w != hgraph[u].end(); w++) {
                if ((*v != *w) && (*v < *w)) {
                    // cout << *v << "," << *w << endl;
                    hgraph[*v].insert(*w);
                    hgraph[*w].insert(*v);
                }
            }
        }

        for (auto v : hgraph[u]) {
            hgraph[v].erase(u);
            v_degree[v].second = hgraph[v].size();
        }

        decomp_order.push_back(u);
        v_degree.erase(v_degree.begin());
        sort(v_degree.begin(), v_degree.end(), [](auto &left, auto &right) { return left.second < right.second; });

        for (int j = 0; j < v_degree.size(); j++) {
            v_index[v_degree[j].first] = j;
        }
    }

    // ************ Finished decomposition ************
    print_decomp_order();
}