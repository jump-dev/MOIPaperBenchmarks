import gzip
import json
import os
import time

def read_mof(filename):
    start = time.time()
    with gzip.GzipFile(filename, 'r') as io:
        data = json.loads(io.read().decode('utf-8'))
    return time.time() - start

data = {}
for filename in os.listdir("miplib_mof"):
    print(filename)
    try:
        data[filename] = read_mof("miplib_mof/" + filename)
    except:
        pass
    with open('timing_python.json', 'w') as io:
        json.dump(data, io)
