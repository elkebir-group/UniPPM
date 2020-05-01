import sys
import os

# ugen_path = sys.argv[1]

file = sys.argv[1]

out_dir = sys.argv[2]

out_file = sys.argv[3]

n_samples = int(sys.argv[4])

lamb = None
if len(sys.argv) >= 6:
    lamb = float(sys.argv[5])


command = "python PPM2SAT.py --file=%s --output=%s"%(file,out_dir+out_file+".cnf")
if lamb!=None:
    command+= ("--lambda=%.2f"%lamb)
print(command)
os.system(command)

command = "python UniGen2.py -runIndex=0 -samples=%d %s %s"%(n_samples, out_dir+out_file+".cnf", out_dir)
print(command)
os.system(command)

command = "python result_trees.py %s_0.txt %s.evar.json %s"%(out_dir+out_file,out_dir+out_file,out_dir+out_file+".trees")
print(command)
os.system(command)

command = "rm %s.cnf %s_0.txt %s_0.count %s.evar.json"%(out_dir+out_file,out_dir+out_file,out_dir+out_file,out_dir+out_file)
os.system(command)