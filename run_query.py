import os

graphs=["NY","COL","CAL"]
types=["V8"]
thresholds = ["10"]
result_path = "./result_road_query/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 50h ./wcsd_66 road_networks/"+ graph + "_tree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_tree_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)