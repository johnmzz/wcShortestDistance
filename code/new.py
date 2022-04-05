import random

def parse_dblp(path, name, delim, num_weights):
    f = open(path + name + "/graph.txt", 'r')
    fgraph = open(path + name + "/new_graph.txt", 'w')

    vertex_map = {}
    vertex_degree = {}
    adj_list = {}
    edges = []

    ncnt = 0
    mcnt = 0

    line = f.readline()
    for line in f.readlines():
        lst = line.split(delim)

        u = lst[0]
        v = lst[1]

        print(f"{u} {v}")

        if u not in vertex_degree.keys():
            vertex_degree[u] = 1
        else:
            vertex_degree[u] += 1

        if v not in vertex_degree.keys():
            vertex_degree[v] = 1
        else:
            vertex_degree[v] += 1

        if u == v:
            continue

        if u not in adj_list:
            adj_list[u] = [v]
        else:
            adj_list[u].append(v)

        if v not in adj_list:
            adj_list[v] = [u]
        else:
            adj_list[v].append(u)


    sorted_vertex_degree = dict(sorted(vertex_degree.items(), key=lambda item: item[1], reverse=True))
    for k,v in sorted_vertex_degree.items():
        vertex_map[k] = ncnt
        ncnt += 1

    for u in sorted_vertex_degree.keys():
        for v in adj_list[u]:
            from_id = vertex_map[u]
            to_id = vertex_map[v]

            e = (min(from_id,to_id),max(from_id,to_id))

            if e not in edges:
                edges.append(e)

    print("#nodes: " + str(len(vertex_map)))
    print("#edges: " + str(len(edges)))

    fgraph.write(str(len(vertex_map)) + " " + str(len(edges)) + " " + str(num_weights) + "\n")

    for e in edges:
        fgraph.write(str(e[0]) + " " + str(e[1]) + " " + str(random.randint(1,3)) + "\n")

path = "./data/"
graph = "DBLP"

parse_dblp(path,graph,",",3)