import os

graphs=["E","W"]
types=["V4"]
thresholds = ["10"]
result_path = "./result_road_degree/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 100h ./wcsd road_networks/"+ graph + "_degree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_degree_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)