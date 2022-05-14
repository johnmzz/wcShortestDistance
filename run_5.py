import os

graphs=["eu-2005", "movielens_10m"]
types=["V4"]
thresholds = ["10"]
result_path = "./result_social_tree/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 100h ./wcsd_51 "+ graph  + "_degree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_tree_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)