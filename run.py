import os

graphs=["eu-2005","uk-2007"]
types=["V4"]
thresholds = ["10"]
result_path = "./query_time/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 30h ./wcsd "+ graph + "_degree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_" + type_ + "_" + threshold + " 2>&1"
            os.system(command)