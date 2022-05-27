#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>

#include "WGraph.h"

using namespace std;

int main(int argc, char* argv[]) {
    string data_path = "/import/vldb01/2/scratch/mazhuo/wcShortestDistance/data/";
    string input_graph = argv[1];
    string type = argv[2];
    int threshold = stoi(argv[3]);

    bool isBin = true;

    WGraph g(data_path, input_graph, isBin);
    //g.print_graph();
    //g.save_graph_to_bin(data_path, input_graph);

    //g.sort_by_degree();
    g.tree_decomp();
    //g.tree_decomp_partial(0.7);
    //g.print_graph();
    g.get_graph_statistics();
    //g.check_degree_order();

    //if (!isBin) g.save_graph_to_bin(data_path, input_graph);

    g.build_index(type, threshold);

    // if (input_graph == "mz_test") g.print_index();

    g.get_index_size();

    cout << "Complete!" << endl;

    /*
    // Test correctness
    random_device rd;
    std::mt19937 gen(rd());
    uniform_int_distribution<uint32_t> gen_num(0, g.max_w);

    double q1 = 0, q2 = 0, q3 =0, q4 = 0;

    uint32_t cnt = 0;
    for (uint32_t s = 0; s < g.nsize - 1; s++) {
        if (s % 200 == 0) cout << "testing s = " << s << endl;
        for (uint32_t t = s + 1; t < g.nsize; t++) {
            cnt++;
            uint32_t r = gen_num(gen);

            uint32_t d1 = g.constrained_shortest_distance_naive(s, t, r, q1);
            uint32_t d2 = g.constrained_shortest_distance_dijkstra(s, t, r, q2);
            uint32_t d3 = g.constrained_shortest_distance_plus(s, t, r, q3);
            uint32_t d4 = g.query(s, t, r, q4);
            if (d1 != d2 || d1 != d3 || d2 != d3 || d1 != d4 || d2 != d4 || d3 != d4) {
                cout << "ERROR: s, t, r: " << s << " " << t << " " << r << "; dist: " << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
            }

        }
    }

   /*
    uniform_int_distribution<uint32_t> gen_s(0, g.nsize-1);
    uniform_int_distribution<uint32_t> gen_d(0, g.nsize-1);
    int cnt = 100;
    for (int i = 0; i < cnt; i++) {
        uint32_t r = gen_num(gen);
        uint32_t s = gen_s(gen);
        uint32_t t = gen_d(gen);

        uint32_t d1 = g.constrained_shortest_distance_naive(s, t, r, q1);
        uint32_t d2 = g.constrained_shortest_distance_dijkstra(s, t, r, q2);
        uint32_t d3 = g.constrained_shortest_distance_plus(s, t, r, q3);
        uint32_t d4 = g.query(s, t, r, q4);
        if (d1 != d2 || d1 != d3 || d2 != d3 || d1 != d4 || d2 != d4 || d3 != d4) {
            cout << "ERROR: s, t, r: " << s << " " << t << " " << r << "; dist: " << d1 << " " << d2 << " " << d3 << " " << d4 << endl;
        }
    }
    */
    /*
    cout << "total quries: " << cnt << endl;
    cout << "total query time (in ms): " << q1 << ", " << q2 << ", " << q3 << ", " << q4 << endl;
    cout << "average query time (in ms): " << q1/cnt << ", " << q2/cnt << ", " << q3/cnt << ", " << q4/cnt << endl;
    //
    g.check_minimality();
    */
    return 0;
}
