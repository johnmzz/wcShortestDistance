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