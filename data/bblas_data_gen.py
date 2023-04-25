import sys
import csv
import os.path
import shutil

file_name = sys.argv[1]
func_name = sys.argv[2]

with open(file_name) as f:
    data = f.read()

data = data.replace("_batchf,", "_batchf_" + func_name + ",")
data = data.replace("_batch,", "_batch_" + func_name + ",")

with open(file_name.replace(".csv","") + "_" + func_name + ".csv", mode="w") as f:
    f.write(data)

