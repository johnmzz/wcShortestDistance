import os

graphs=["frwiki"]
types=["V4","V5"]
thresholds = ["10"]
result_path = "./results/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 30h ./wcsd_3_24 "+ graph + "_degree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_" + type_ + "_" + threshold + " 2>&1"
            os.system(command)