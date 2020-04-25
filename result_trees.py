import json
import sys

with open(sys.argv[2],"r") as f:
    evar_info = json.load(f)

if len(sys.argv)>4 and sys.argv[4]=="-i":
    evar = evar_info["0"]
else:
    evar = evar_info["1"]

Result = []
cnt = []
with open(sys.argv[1],"r") as fin:
    for line in fin:
        result = line.split()
        if (len(result)<2): continue
        cnt.append(int(result[0]))
        Result.append(result[1:])

with open(sys.argv[3],"w") as fout:
    fout.write("#%d trees sampled\n"%(sum(cnt)))
    n_unique = len(Result)
    for i,result in enumerate(Result):
        edges = [int(e_var) for e_var in result ]
        sel_edges = [edge for edge in edges if edge>0]
        fout.write("#%d edges, tree %d, %d samples\n"%(len(sel_edges), i, cnt[i]) )
        for ed in sel_edges:
            if str(ed) not in evar.keys(): continue
            fout.write("%s %s\n"%(evar[str(ed)][0],evar[str(ed)][1]) )
        
    