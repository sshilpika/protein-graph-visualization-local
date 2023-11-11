import time

all_dot = "./src/visg/static/data/interactions_full_run_static.dot"
master_file = "./src/visg/static/data/interactions_full_run.dot"

with open(master_file, 'w') as f:
    f.write("digraph G {\n")
    f.write("}")

time.sleep(5)

linesall = []
with open(all_dot, 'r') as f:
    linesall = f.readlines()[1:-1]
print(len(linesall))
step = 500
for i, j in zip(range(0,len(linesall), step), range(step,len(linesall), step)):
    time.sleep(5)

    with open(master_file, 'r') as f:
        smlines = f.readlines()

    writeL = smlines[:-1] + linesall[i:j] + ["}"]
    with open(master_file, 'w') as f:
        f.write("".join(writeL))
