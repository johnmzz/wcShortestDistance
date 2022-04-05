import random

def parse_dblp(path, name, delim, num_weights):
    f = open(path + name + "/old_graph.txt", 'r')
    fgraph = open(path + name + "/graph.txt", 'w')

    vertex_map = {}
    vertex_degree = {}
    adj_list = {}
    edges = {}

    ncnt = 0
    mcnt = 0

    line = f.readline()
    for line in f.readlines():
        lst = line.split(delim)

        u = lst[0]
        v = lst[1]

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
        print(f"{k} mapped to {vertex_map[k]}")

    for u in sorted_vertex_degree.keys():
        for v in adj_list[u]:
            from_id = vertex_map[u]
            to_id = vertex_map[v]

            e = (min(from_id,to_id),max(from_id,to_id))

            if e not in edges.keys():
                edges[e] = 1
            else:
                edges[e] += 1

    #num_weights = max(edges.values())

    print("#nodes: " + str(len(vertex_map)))
    print("#edges: " + str(len(edges)))
    print(f"max weight: {num_weights}")

    fgraph.write(str(len(vertex_map)) + " " + str(len(edges)) + " " + str(num_weights) + "\n")

    for e in edges.keys():
        #fgraph.write(str(e[0]) + " " + str(e[1]) + " " + str(edges[e]//2) + "\n")
        fgraph.write(str(e[0]) + " " + str(e[1]) + " " + str(random.randint(1,num_weights)) + "\n")

    print("Finished!")

path = "./data/"
graph = "DBLP"

parse_dblp(path,graph,",",5)