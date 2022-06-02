import os

graphs=["movielens_10m", "eu-2005","eswiki-2013","uk-2007","movielens_25m","frwiki"]
types=["V9"]
thresholds = ["10"]
result_path = "./result_social_tree/"

for graph in graphs:
    for type_ in types:
        for threshold in thresholds:
            command = "/usr/bin/time -v timeout 50h ./wcsd_tree_decomp_partial "+ graph + "_degree " + type_ + " " + threshold + " 1 > " + result_path + graph + "_tree_" + type_ + "_" + threshold + " 2>&1"
            print(command)
            os.system(command)