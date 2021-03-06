#include "WGraph.h"

WGraph::WGraph() {
    graph = AdjList();
    hgraph = vector<set<uint32_t>>();
    wlist = AdjList();
    nsize = 0;
    msize = 0;
    wsize = 0;
    ndegree = vector<uint32_t>();

    labeling_v = vector<vector<uint32_t>>();
    offset = vector<vector<uint32_t>>();
    labeling_d = vector<vector<uint16_t>>();
    labeling_w = vector<vector<uint8_t>>();
    labeling_dw = vector<vector<pair<uint16_t, uint8_t>>>();
    labeling_bit = vector<vector<uint32_t>>();
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

// save graph to bin, bin
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

    uint32_t max_degree = *max_element(ndegree.begin(), ndegree.end());
    cout << "the maximum degree: " << max_degree << endl;
}

void WGraph::print_graph() {
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

// O(n^)
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

// default query algorithm
uint32_t WGraph::query(uint32_t s, uint32_t t, uint32_t r, double& qtime) {
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
                if (labeling_dw[s][p].second >= r) {
                    d1 = labeling_dw[s][p].first;
                    break;
                }
            }
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_dw[t][q].second >= r) {
                    d2 = labeling_dw[t][q].first;
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
                if (labeling_dw[s][p].second >= r) {
                    d1 = labeling_dw[s][p].first;
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
                if (labeling_dw[t][q].second >= r) {
                    d2 = labeling_dw[t][q].first;
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

uint32_t WGraph::query_skip(uint32_t s, uint32_t t, uint32_t r, uint32_t dpos) {
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
                if (labeling_dw[s][p].second >= r) {
                    d1 = labeling_dw[s][p].first;
                    break;
                }
            }
            for (uint32_t q = offset[t][j]; q < offset[t][j + 1]; q++) {
                if (labeling_dw[t][q].second >= r) {
                    d2 = labeling_dw[t][q].first;
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
                if (labeling_dw[s][p].second >= r) {
                    d1 = labeling_dw[s][p].first;
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
                if (labeling_dw[t][q].second >= r) {
                    d2 = labeling_dw[t][q].first;
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

template<typename T>
bool compare_w(const T &a,const T &b){
    return a.second<b.second;
}

bool WGraph::query_while_indexing_vertex_V1(uint32_t s, uint32_t t, uint16_t curr_d, uint8_t r) {
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
            uint32_t l1 = lower_bound(labeling_dw[s].begin() + offset[s][i], labeling_dw[s].begin() + offset[s][i + 1], make_pair(0,r), compare_w<pair<uint16_t, uint8_t>>) - labeling_dw[s].begin();
            if (l1 != offset[s][i + 1]) {
                uint32_t l2 = lower_bound(labeling_dw[t].begin() + offset[t][j], labeling_dw[t].begin() + offset[t][j + 1], make_pair(0,r), compare_w<pair<uint16_t, uint8_t>>) - labeling_dw[t].begin();
                if (l2 != offset[t][j + 1]) {
                    if (labeling_dw[s][l1].first + labeling_dw[t][l2].first <= curr_d) {
                        return true;
                    }
                }
            }
            i++;
            j++;
        }
    }
    return false;
}

bool WGraph::query_while_indexing_vertex_V2(uint32_t s, uint32_t t, uint16_t curr_d, uint8_t r) {
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
            uint32_t l1 = lower_bound(labeling_bit[s].begin() + offset[s][i], labeling_bit[s].begin() + offset[s][i + 1], ((uint32_t) r << 16)) - labeling_bit[s].begin();
            if (l1 != offset[s][i + 1]) {
                uint32_t l2 = lower_bound(labeling_bit[t].begin() + offset[t][j], labeling_bit[t].begin() + offset[t][j + 1], ((uint32_t) r << 16)) - labeling_bit[t].begin();
                if (l2 != offset[t][j + 1]) {
                    if ((uint16_t)labeling_bit[s][l1] + (uint16_t)labeling_bit[t][l2] <= curr_d) {
                        return true;
                    }
                }
            }
            i++;
            j++;
        }
    }
    return false;
}

bool WGraph::query_while_indexing_vertex_V8(uint32_t s, uint32_t t, uint16_t curr_d, uint8_t r, vector<vector<pair<uint16_t, uint8_t>>>& u_label) {
    for (uint32_t i = 0; i < labeling_v[t].size(); i++) {
        uint32_t v = labeling_v[t][i];

        if (v > s) break;

        if (u_label[v].empty()) continue;

        uint32_t l1 = lower_bound(labeling_dw[t].begin() + offset[t][i], labeling_dw[t].begin() + offset[t][i + 1], make_pair(0,r), compare_w<pair<uint16_t, uint8_t>>) - labeling_dw[t].begin();
        if (l1 != offset[t][i + 1]) {
            uint32_t l2 = lower_bound(u_label[v].begin(), u_label[v].end(), make_pair(0,r), compare_w<pair<uint16_t, uint8_t>>) - u_label[v].begin();
            if (l2 != u_label[v].size()) {
                uint16_t d = labeling_dw[t][l1].first;
                uint16_t u_d = u_label[v][l2].first;

                if ((u_d + d) <= curr_d) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool WGraph::query_while_indexing_vertex_V9(uint32_t s, uint32_t t, uint16_t curr_d, uint8_t r, vector<vector<uint32_t>>& u_label) {
    for (uint32_t i = 0; i < labeling_v[t].size(); i++) {
        uint32_t v = labeling_v[t][i];

        if (v > s) break;

        if (u_label[v].empty()) continue;

        uint32_t l1 = lower_bound(labeling_bit[t].begin() + offset[t][i], labeling_bit[t].begin() + offset[t][i + 1], ((uint32_t) r << 16)) - labeling_bit[t].begin();
        if (l1 != offset[t][i + 1]) {
            uint32_t l2 = lower_bound(u_label[v].begin(), u_label[v].end(), ((uint32_t) r << 16)) - u_label[v].begin();
            if (l2 != u_label[v].size()) {
                uint16_t d = (uint16_t) labeling_bit[t][l1];
                uint16_t u_d = (uint16_t) u_label[v][l2];

                if ((u_d + d) <= curr_d) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool WGraph::query_while_indexing_vertex_V10(uint32_t s, uint32_t t, uint16_t curr_d, uint8_t r, vector<vector<uint16_t>>& new_label) {
    for (uint32_t i = 0; i < labeling_v[t].size(); i++) {
        uint32_t v = labeling_v[t][i];

        if (v > s) break;

        if (new_label[v].empty()) continue;

        uint16_t u_d = new_label[v][r];
        if (u_d != INVALID_VALUE_16bit) {
            uint32_t l1 = lower_bound(labeling_dw[t].begin() + offset[t][i], labeling_dw[t].begin() + offset[t][i + 1], make_pair(0,r), compare_w<pair<uint16_t, uint8_t>>) - labeling_dw[t].begin();
            if (l1 != offset[t][i+1]) {
                uint16_t d = labeling_dw[t][l1].first;
                if ((u_d + d) <= curr_d) {
                    return true;
                }
            }
        }
    }
    return false;
}

void WGraph::build_index(string type, int threshold) {
    labeling_v = vector<vector<uint32_t>>(nsize);
    offset = vector<vector<uint32_t>>(nsize);
    labeling_dw = vector<vector<pair<uint16_t, uint8_t>>>(nsize);
    labeling_bit = vector<vector<uint32_t>>(nsize);

    for (uint32_t u = 0; u < nsize; u++) {
        offset[u].push_back(0);
    }

    boost::dynamic_bitset<> updated = boost::dynamic_bitset<>(nsize);
    vector<uint8_t> visited_r = vector<uint8_t>(nsize, 0);
    boost::dynamic_bitset<> this_visited_flag = boost::dynamic_bitset<>(nsize);

    clock_t start, end;
    start = clock();
    if (type == "V8") {
        cout << "use V8 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V8(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_V1(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V9") {
        cout << "use V9 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V9(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_V2(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    else if (type == "V10") {
        cout << "use V10 for first " << threshold << "\% vertices..." << ", which is " << (int)(nsize*threshold/100) << endl;
        for (uint32_t u = 0; u < nsize; u++) {
            if (u < (int)(nsize*threshold/100)) {
                vertex_prioritized_indexing_V10(u, updated, visited_r, this_visited_flag);
            }
            else {
                vertex_prioritized_indexing_V1(u, updated, visited_r, this_visited_flag);
            }
            if (u % 2000 == 0) {
                cout << u << " vertices finished!" << endl;
                clock_t temp = clock();
                cout << "this round took " << (float)(temp - start) / CLOCKS_PER_SEC << " s" << endl;
            }
        }
    }
    end = clock();
    cout << "build index: " << (float)(end - start) / CLOCKS_PER_SEC << " s" << endl;
}

void WGraph::vertex_prioritized_indexing_V1(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint8_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint8_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE_8bit));
    visited_r[u] = INVALID_VALUE_8bit;
    visited_vertex.push_back(u);
    uint16_t curr_d = 0;

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint8_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                if (query_while_indexing_vertex_V1(u, v, curr_d, r)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint8_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                if (query_while_indexing_vertex_V1(u, v, curr_d, r)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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

void WGraph::vertex_prioritized_indexing_V2(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint8_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint8_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE_8bit));
    visited_r[u] = INVALID_VALUE_8bit;
    visited_vertex.push_back(u);
    uint16_t curr_d = 0;

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint8_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                if (query_while_indexing_vertex_V2(u, v, curr_d, r)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_bit[v].push_back((uint32_t) r << 16 | curr_d);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_bit[v].push_back((uint32_t) r << 16 | curr_d);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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
            vector<uint32_t> this_visited;
            while (!Q2.empty()) {
                pair<uint32_t, uint8_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                if (query_while_indexing_vertex_V2(u, v, curr_d, r)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_bit[v].push_back((uint32_t) r << 16 | curr_d);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_bit[v].push_back((uint32_t) r << 16 | curr_d);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }  /// only explore the vertex that ranks lower than the starting vertex u
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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

void WGraph::vertex_prioritized_indexing_V8(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint8_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint8_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE_8bit));
    visited_r[u] = INVALID_VALUE_8bit;
    visited_vertex.push_back(u);
    uint16_t curr_d = 0;

    vector<uint8_t> max_pruned_r = vector<uint8_t>(nsize, 0);

    register vector<vector<pair<uint16_t, uint8_t>>> u_label = vector<vector<pair<uint16_t, uint8_t>>>(nsize, vector<pair<uint16_t, uint8_t>>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint16_t d = labeling_dw[u][p].first;
            uint8_t r = labeling_dw[u][p].second;
            u_label[v].emplace_back(make_pair(d, r));
        }
    }
    u_label[u].emplace_back(make_pair(0, INVALID_VALUE_8bit));

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint8_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                uint16_t query_d;

                // if (r < max_pruned_r[v]) {
                //     continue;
                // } else if
                if (query_while_indexing_vertex_V8(u, v, curr_d, r, u_label)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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
                pair<uint32_t, uint8_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                uint16_t query_d;

                // if (r < max_pruned_r[v]) {
                //     continue;
                // } else if
                if (query_while_indexing_vertex_V8(u, v, curr_d, r, u_label)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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

void WGraph::vertex_prioritized_indexing_V9(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint8_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint8_t>> Q1, Q2;
    Q1.emplace(make_pair(u, INVALID_VALUE_8bit));
    visited_r[u] = INVALID_VALUE_8bit;
    visited_vertex.push_back(u);
    uint16_t curr_d = 0;

    vector<uint8_t> max_pruned_r = vector<uint8_t>(nsize, 0);

    register vector<vector<uint32_t>> u_label = vector<vector<uint32_t>>(nsize, vector<uint32_t>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            u_label[v].emplace_back(labeling_bit[u][p]);
        }
    }
    u_label[u].emplace_back(((uint32_t) 0 << 16) | INVALID_VALUE);

    while (!Q1.empty() || !Q2.empty()) {
        if (!Q1.empty()) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint8_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                uint16_t query_d;

                // if (r < max_pruned_r[v]) {
                //     continue;
                // } else if
                if (query_while_indexing_vertex_V9(u, v, curr_d, r, u_label)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_bit[v].push_back((uint32_t) r << 16 | curr_d);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_bit[v].push_back(((uint32_t) r << 16) | curr_d);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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
                pair<uint32_t, uint8_t> p = Q2.front();
                Q2.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                uint16_t query_d;

                // if (r < max_pruned_r[v]) {
                //     continue;
                // } else if
                if (query_while_indexing_vertex_V9(u, v, curr_d, r, u_label)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_bit[v].push_back(((uint32_t) r << 16) | curr_d);
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_bit[v].push_back(((uint32_t) r << 16) | curr_d);
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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

void WGraph::vertex_prioritized_indexing_V10(uint32_t u, boost::dynamic_bitset<>& updated, vector<uint8_t>& visited_r, boost::dynamic_bitset<>& this_visited_flag) {
    vector<uint32_t> visited_vertex;
    queue<pair<uint32_t, uint8_t>> Q1;
    Q1.emplace(make_pair(u, INVALID_VALUE_8bit));
    visited_r[u] = INVALID_VALUE_8bit;
    visited_vertex.push_back(u);
    uint16_t curr_d = 0;

    vector<uint8_t> max_pruned_r = vector<uint8_t>(nsize, 0);

    register vector<vector<pair<uint16_t, uint8_t>>> u_label = vector<vector<pair<uint16_t, uint8_t>>>(nsize, vector<pair<uint16_t, uint8_t>>());
    for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
        uint32_t v = labeling_v[u][i];
        for (uint32_t p = offset[u][i]; p < offset[u][i + 1]; p++) {
            uint16_t d = labeling_dw[u][p].first;
            uint8_t r = labeling_dw[u][p].second;
            u_label[v].emplace_back(make_pair(d, r));
        }
    }
    u_label[u].emplace_back(make_pair(0, INVALID_VALUE_8bit));

    vector<vector<uint16_t>> new_label(nsize, vector<uint16_t>(max_w+1, INVALID_VALUE_16bit));
    for (uint32_t v = 0; v < u_label.size(); v++) {
        for (int i = max_w; i >= 0; i--) {
            int idx = lower_bound(u_label[v].begin(), u_label[v].end(), make_pair(0,i), compare_w<pair<uint16_t, uint8_t>>) - u_label[v].begin();
            if (idx != u_label[v].size()) {
                new_label[v][i] = u_label[v][idx].first;
            }
        }
    }

    // cout << "u = " << u << ": " << "u_label = " << endl;
    // for (auto u_label_v : u_label) {
    //     cout << "[";
    //     for (auto p : u_label_v) {
    //         cout << "(" << p.first << "," << (int)p.second << "),";
    //     }
    //     cout << "]" << endl;
    // }

    // cout << "new_label = " << endl;
    // for (auto new_label_v : new_label) {
    //     cout << "[";
    //     for (auto new_v : new_label_v) {
    //         cout << new_v << ",";
    //     }
    //     cout << "]" << endl;
    // }
    // cout << endl;

    while (!Q1.empty()) {
        uint32_t q_size = Q1.size();
        for (uint32_t i = 0; i < q_size; i++) {
            vector<uint32_t> this_visited;
            while (!Q1.empty()) {
                pair<uint32_t, uint8_t> p = Q1.front();
                Q1.pop();
                uint32_t v = p.first;
                uint8_t r = p.second;
                uint16_t query_d;

                // if (r < max_pruned_r[v]) {
                //     continue;
                // } else if
                if (query_while_indexing_vertex_V10(u, v, curr_d, r, new_label)) {
                    continue;
                }
                else {
                    if (updated[v] & 1) {
                        offset[v].back() = offset[v].back() + 1;
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                    }
                    else {
                        labeling_v[v].push_back(u);
                        offset[v].push_back(offset[v].back() + 1);
                        labeling_dw[v].push_back(make_pair(curr_d, r));
                        updated[v] = 1;
                    }
                }
                for (uint32_t idx = 0; idx < ndegree[v]; idx++) {
                    uint32_t w = graph[v][idx];
                    if (w <= u) {
                        continue;
                    }
                    uint8_t new_r = min(r, (uint8_t)wlist[v][idx]);
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

void WGraph::print_index() {
    for (uint32_t u = 0; u < nsize; u++) {
        cout << u << ": ";
        for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
            for (uint32_t j = offset[u][i]; j < offset[u][i + 1]; j++) {
                cout << "<" << labeling_v[u][i] << "," << (int)labeling_dw[u][j].first << "," << (int)labeling_dw[u][j].second << ">, ";
            }
        }
        cout << endl;
    }
}

void WGraph::print_index_2() {
    for (uint32_t u = 0; u < nsize; u++) {
        cout << u << ": ";
        for (uint32_t i = 0; i < labeling_v[u].size(); i++) {
            for (uint32_t j = offset[u][i]; j < offset[u][i + 1]; j++) {
                cout << "<" << labeling_v[u][i] << "," << (uint16_t)labeling_bit[u][j] << "," << (labeling_bit[u][j] >> 16) << ">, ";
            }
        }
        cout << endl;
    }
}

void WGraph::get_index_size() {
    cout << "--------------Details of Index Occupation-------------" << endl;
    uint64_t cnt1 = 0, cnt2 = 0, cnt3 = 0;
    for (uint32_t u = 0; u < nsize; u++) {
        cnt1 += labeling_v[u].size();
        cnt2 += offset[u].size();
        cnt3 += labeling_dw[u].size();
    }

    cout << "labeling v, offset, labeling d, labeling w: " << cnt1 << ", " << cnt2 << ", " << cnt3 << endl;
    cout << "total number of unsigned ints 32 in the index: " << cnt1 + cnt2 << ", total number of unsigned int 16+8: " << cnt3 << endl;
}

void WGraph::check_minimality() {
    cout << "check minimality " << endl;
    for (uint32_t u = 0; u < nsize; u++) {
        for (uint32_t wpos = 0; wpos < labeling_v[u].size(); wpos++) {
            for (uint32_t dpos = offset[u][wpos]; dpos < offset[u][wpos + 1]; dpos++) {
                uint32_t d1 = query_skip(u, labeling_v[u][wpos], labeling_dw[u][dpos].second, dpos);
                if (d1 == labeling_dw[u][dpos].first) {
                    cout << "failed at u, v, r, d: " << u << ", " << labeling_v[u][wpos] << ", " << labeling_dw[u][dpos].second << ", " << labeling_dw[u][dpos].first << endl;
                }
            }
        }
    }
}

void WGraph::print_hgraph() {
    cout << "h_graph size = " << hgraph.size() << endl;
    for (uint32_t i = 0; i < hgraph.size(); i++) {
        cout << i << ": [";
        for (auto j : hgraph[i]) {
            cout << j << ",";
        }
        cout << "]" << endl;
    }
}

void print_v_degree_ordered(const vector<pair<uint32_t, uint32_t>> &vec) {
    for (auto v : vec) {
        cout << "(" << v.first << "," << v.second << "),";
    }
    cout << endl;
}

// check if vertices of graph are ordered by degree (decreasing order)
void WGraph::check_degree_order() {
    for (uint32_t i = 0; i < graph.size(); i++) {
        // cout << i << "," << ndegree[i] << endl;

        if (i != 0 && ndegree[i] > ndegree[i-1]) {
            cout << "Error: ndegree[" << i << "] > ndegree[" << i-1 << "]" << endl;
        }
    }
}

void WGraph::sort_by_degree() {
    // sort all vertices by degree
    vector<pair<uint32_t, uint32_t>> v_degree_ordered(nsize);
    for (uint32_t i = 0; i < graph.size(); i++) {
        v_degree_ordered[i] = make_pair(i, ndegree[i]);
    }
    std::sort(v_degree_ordered.begin(), v_degree_ordered.end(), [](auto &left, auto &right) {return left.second > right.second;});

    // stor offset of each v in v_degree_ordered
    vector<uint32_t> offset(nsize);
    for (uint32_t i = 0; i < nsize; i++) {
        offset[v_degree_ordered[i].first] = i;
    }

    // construct new graph
    AdjList new_graph(nsize);
    AdjList new_wlist(nsize);
    for (uint32_t i = 0; i < nsize; i++) {
        uint32_t u = v_degree_ordered[i].first;
        for (uint32_t j = 0; j < graph[u].size(); j++) {
            uint32_t v = graph[u][j];
            new_graph[i].push_back(offset[v]);
            new_wlist[i].push_back(wlist[u][j]);
        }
    }
    graph = new_graph;
    wlist = new_wlist;

    // reconstruct degree vector
    for (uint32_t i = 0; i < nsize; i++) {
        ndegree[i] = graph[i].size();
    }

    cout << "Finished sorting by degree." << endl;
    // print_graph();
}

void WGraph::tree_decomp() {
    clock_t start_total, end_total;
    start_total = clock();

    // initialize hgraph as graph
    for (auto i = 0; i < nsize; i++) {
        set<uint32_t> neighbors;
        for (auto v : graph[i]) {
            neighbors.insert(v);
        }
        hgraph.push_back(neighbors);
    }

    uint32_t count = 0;
    vector<uint32_t> old_deg(nsize), new_deg(nsize);

    auto comp = [&old_deg](const uint32_t &a, const uint32_t &b) -> bool {
        return old_deg[a] == old_deg[b] ? a < b : old_deg[a] < old_deg[b];
    };

    set<uint32_t, decltype(comp)> ndsQ(comp); // set of vID, ordered by comp (degree then ID)

    for (uint32_t s = 0; s < nsize; s++) {
	    new_deg[s] = old_deg[s] = graph[s].size();
	    ndsQ.insert(s);
	}

    vector<int32_t> vertexOrder(nsize, -1);
    uint32_t vorder = 0;

    cout << "\nStart tree decomposition...\n" << endl;

    while (ndsQ.size() > 0) {
        clock_t start, end;
        start = clock();

        count++;
        uint32_t vid = *ndsQ.begin();
        ndsQ.erase(*ndsQ.begin());

        // re-order ndsQ if vid's new_deg > old_deg
        while (old_deg[vid] != new_deg[vid]) {
		    old_deg[vid] = new_deg[vid];
		    ndsQ.insert(vid);

		    vid = *ndsQ.begin();
		    ndsQ.erase(ndsQ.begin());
        }

        vertexOrder[vid] = vorder++;
        auto &v_adj = hgraph[vid];

        vector<uint32_t> valid_neighbor_index;
        for (auto u : v_adj) {
            if (vertexOrder[u] == -1) {
                valid_neighbor_index.push_back(u);
                hgraph[u].erase(vid);
            }
        }
        vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), 0);

        // iterate through neighbor-pairs of vid
        for (auto i = 0; i < valid_neighbor_index.size(); i++) {
            uint32_t u = valid_neighbor_index[i];
            neighbor_degree_increase_cnt[i]--;

            for (auto j = i+1; j < valid_neighbor_index.size(); j++) {
                uint32_t w = valid_neighbor_index[j];

                if (hgraph[u].find(w) == hgraph[u].end()) {
                    hgraph[u].insert(w);
                    hgraph[w].insert(u);
                    neighbor_degree_increase_cnt[i]++;
                    neighbor_degree_increase_cnt[j]++;
                }
            }

            new_deg[u] += neighbor_degree_increase_cnt[i];

            if (neighbor_degree_increase_cnt[i] < 0) {
		        ndsQ.erase(u);
		        old_deg[u] = new_deg[u];
		        ndsQ.insert(u);
		    }
        }

        end = clock();
        if (count % 5000 == 0) cout << "decomposing vertex " << count << "time: " << 5000*(float)(end - start) / CLOCKS_PER_SEC << " s" << endl;
    }

    // construct new graph
    cout << "Constructing new graph..." << endl;

    AdjList new_graph(nsize);
    AdjList new_wlist(nsize);
    for (uint32_t i = 0; i < nsize; i++) {
        uint32_t u = nsize - vertexOrder[i] - 1;
        for (uint32_t j = 0; j < graph[i].size(); j++) {
            uint32_t v = nsize - vertexOrder[graph[i][j]] - 1;
            new_graph[u].push_back(v);
            new_wlist[u].push_back(wlist[i][j]);
        }
    }
    graph = new_graph;
    wlist = new_wlist;

    // reconstruct degree vector
    for (uint32_t i = 0; i < nsize; i++) {
        ndegree[i] = graph[i].size();
    }

    end_total = clock();
    cout << "Finished tree decomposition, total time = " << (float)(end_total - start_total) / CLOCKS_PER_SEC << " s" << endl;
}

void WGraph::tree_decomp_partial(float threshold) {
    cout << "Perform tree decomposition for the first " << (int)(threshold*100) << "\% of verties." << endl;
    clock_t start_total, end_total;
    start_total = clock();

    // initialize hgraph as graph
    for (auto i = 0; i < nsize; i++) {
        set<uint32_t> neighbors;
        for (auto v : graph[i]) {
            neighbors.insert(v);
        }
        hgraph.push_back(neighbors);
    }

    uint32_t count = 0;
    vector<uint32_t> old_deg(nsize), new_deg(nsize);

    auto comp = [&old_deg](const uint32_t &a, const uint32_t &b) -> bool {
        return old_deg[a] == old_deg[b] ? a < b : old_deg[a] < old_deg[b];
    };

    set<uint32_t, decltype(comp)> ndsQ(comp); // set of vID, ordered by comp (degree then ID)

    for (uint32_t s = 0; s < nsize; s++) {
	    new_deg[s] = old_deg[s] = graph[s].size();
	    ndsQ.insert(s);
	}

    vector<int32_t> vertexOrder(nsize, -1);
    uint32_t vorder = 0;

    cout << "\nStart tree decomposition...\n" << endl;

    while (ndsQ.size() > 0) {
        clock_t start, end;
        start = clock();

        count++;
        uint32_t vid = *ndsQ.begin();
        ndsQ.erase(*ndsQ.begin());

        // re-order ndsQ if vid's new_deg > old_deg
        while (old_deg[vid] != new_deg[vid]) {
		    old_deg[vid] = new_deg[vid];
		    ndsQ.insert(vid);

		    vid = *ndsQ.begin();
		    ndsQ.erase(ndsQ.begin());
	    }

        vertexOrder[vid] = vorder++;
        auto &v_adj = hgraph[vid];

        vector<uint32_t> valid_neighbor_index;
        for (auto u : v_adj) {
            if (vertexOrder[u] == -1) {
                valid_neighbor_index.push_back(u);
                hgraph[u].erase(vid);
            }
        }
        vector<int> neighbor_degree_increase_cnt(valid_neighbor_index.size(), 0);

        // iterate through neighbor-pairs of vid
        for (auto i = 0; i < valid_neighbor_index.size(); i++) {
            uint32_t u = valid_neighbor_index[i];
            neighbor_degree_increase_cnt[i]--;

            for (auto j = i+1; j < valid_neighbor_index.size(); j++) {
                uint32_t w = valid_neighbor_index[j];

                if (hgraph[u].find(w) == hgraph[u].end()) {
                    hgraph[u].insert(w);
                    hgraph[w].insert(u);
                    neighbor_degree_increase_cnt[i]++;
                    neighbor_degree_increase_cnt[j]++;
                }
            }

            new_deg[u] += neighbor_degree_increase_cnt[i];

            if (neighbor_degree_increase_cnt[i] < 0) {
		        ndsQ.erase(u);
		        old_deg[u] = new_deg[u];
		        ndsQ.insert(u);
		    }
        }
        end = clock();
        if (count % 5000 == 0) cout << "decomposing vertex " << count << ", time: " << 5000*(float)(end - start) / CLOCKS_PER_SEC << " s" << endl;

        if (count > threshold * nsize) {
            cout << "Finishing tree decomposition at count = " << count << endl;
            break;
        }
    }

    // find remaining vertices not sorted
    vector<uint32_t> v_remain;
    for (uint32_t i = 0; i < nsize; i++) {
        if (vertexOrder[i] == -1) {
            v_remain.push_back(i);
        }
    }

    // sort all remaining vertices by degree
    vector<pair<uint32_t, uint32_t>> v_degree_ordered(v_remain.size());
    for (uint32_t i = 0; i < v_remain.size(); i++) {
        // v_degree_ordered[i] = make_pair(v_remain[i], ndegree[v_remain[i]]);
        v_degree_ordered[i] = make_pair(v_remain[i], hgraph[v_remain[i]].size());
    }
    std::sort(v_degree_ordered.begin(), v_degree_ordered.end(), [](auto &left, auto &right) {return left.second > right.second;});
    std::reverse(v_degree_ordered.begin(),v_degree_ordered.end());

    // store order in vertexOrder
    for (auto v : v_degree_ordered) {
        vertexOrder[v.first] = count;
        count++;
    }

    // construct new graph
    cout << "Constructing new graph..." << endl;

    AdjList new_graph(nsize);
    AdjList new_wlist(nsize);
    for (uint32_t i = 0; i < nsize; i++) {
        uint32_t u = nsize - vertexOrder[i] - 1;
        for (uint32_t j = 0; j < graph[i].size(); j++) {
            uint32_t v = nsize - vertexOrder[graph[i][j]] - 1;
            new_graph[u].push_back(v);
            new_wlist[u].push_back(wlist[i][j]);
        }
    }
    graph = new_graph;
    wlist = new_wlist;

    // reconstruct degree vector
    for (uint32_t i = 0; i < nsize; i++) {
        ndegree[i] = graph[i].size();
    }

    end_total = clock();
    cout << "Finished tree decomposition, total time = " << (float)(end_total - start_total) / CLOCKS_PER_SEC << " s" << endl << endl;
}