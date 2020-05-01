import pandas as pd
import CNF
import sys
import getopt
import math
import scipy.stats
import json
import statsmodels.stats.proportion 


try:
    opts, args = getopt.getopt(sys.argv[1:], "h", ["file=", "output=", "lambda=","scs=", "N="])# , "upper=", "lower="])
except Exception as Err:
    print(Err)
    sys.exit(2)

filename = None
scs_file = None
N=10

# range_upper=range_lower=0
#print(args)
lambda_=0.8
output=None

for opt, arg in opts:
    if opt=="--file":
        filename = arg
    if opt=="--lambda":
        lambda_ = float(arg)
    if opt=="--output":
        output=arg
    if opt=="--scs":
        scs_file = arg
    if opt=="--N":
        N=int(arg)
#     if opt=="--upper":
#         range_upper=int(arg)
#     if opt=="--lower":
#         range_lower=int(arg)

df=pd.read_csv(filename,sep="\t")
df=df.astype({"sample_index":int, "sample_label":"object", "anatomical_site_index":int, "anatomical_site_label":"object",  "character_index":int, "character_label":"object"})
#print(df.dtypes)

if (output==None): sys.exit(2)
evar_json=output.split(".")[0]+".evar.json"


def calculate_f_range(x,alpha):
    # Identify the binomal confindence interval
    # Jeffreys interval, percentile at alpha/2 and (1-alpha/2)
    dis=scipy.stats.beta(x["var"]+.5,x["ref"]+.5)
    if x["var"] == 0:
        l = 0
    else :
        l = dis.ppf((1-alpha)/2)
    
    if x["ref"] == 0:
        r = 1
    else :
        r = dis.ppf((1+alpha)/2)
    return l , r


def calculate_f_prob(df):
    # Calculate the worst possible log likelihood given the ranges for all the entries
    probs1=df.apply( lambda x: scipy.stats.binom.pmf(x["var"],x["ref"]+x["var"],x["f-"]) , axis=1)
    probs2=df.apply( lambda x: scipy.stats.binom.pmf(x["var"],x["ref"]+x["var"],x["f+"]) , axis=1)
    re=0.0
    for p1,p2 in zip(list(probs1),list(probs2)):
        re+=math.log(min(p1,p2))
    return re

def find_alpha(df,lambda_):
    # Find the alpha that would satisfy the approximation factor lambda 
    # binary search
    mid=0
    df["f-"]=df.apply( (lambda x: 1.0*x["var"]/(x["var"]+x["ref"])), axis = 1)
    df["f+"]=df["f-"]
    #df["f-"],df["f+"]=tmp.apply(lambda x:x[0]),tmp.apply(lambda x:x[1])
    upper=calculate_f_prob(df)
    target=upper/(lambda_)
    
    l,r=0,1
    ans=-1e300
    while (abs (ans-target) > 0.000001 and abs(l-r) > 0.00000001):
        mid=(l+r)/2
        tmp=df.apply(lambda x: calculate_f_range(x,mid),axis=1)
        #tmp = df.apply(lambda x: statsmodels.stats.proportion.proportion_confint(x["var"],x["var"]+x["ref"],
        #                                                                             mid,method="beta"),axis=1)
        df["f-"],df["f+"]=tmp.apply(lambda x:x[0]),tmp.apply(lambda x:x[1])
        ans=calculate_f_prob(df)
        #print(mid,ans)
        if (ans>target):
            l=mid
        else :
            r=mid
    print(mid)
    print(target)
    print(ans)
    
    if(mid < .0000000075) : 
        print("This lamda is not supportted, please use a smaller lambda"); sys.exit(-1);
    
    tmp=df.apply(lambda x: calculate_f_range(x,mid**(1/len(df))),axis=1)
    df["f-"],df["f+"]=tmp.apply(lambda x:x[0]),tmp.apply(lambda x:x[1])
    return df


# If the input contains only reads data, require parameter lambda to proceed 
if ("f-" not in df.columns) or ("f+" not in df.columns):
    if lambda_==None: sys.exit(2)
    find_alpha(df,lambda_)

#print(df)


# m,n
m=len(set(df["sample_index"])) #samples
n=len(set(df["character_index"])) #mutations

# check if it is the normal PPM problem
df["single"] = df.apply( lambda x: (abs(x["f+"]-x["f-"])<=1e-6), axis=1)

if ( True in set(df["single"]) ) and (False not in set(df["single"]) ):
    single=True
else:
    single=False


# discretize the frequencies
df["f+int"]=df["f+"].apply( lambda x: int(math.floor(x*(2**N-1))) )
#df["f+bin"]=df["f+int"].apply( lambda x: bin(x)[2:] )

df["f-int"]=df["f-"].apply( lambda x: int(math.floor(x*(2**N-1))) )
#df["f-bin"]=df["f-int"].apply( lambda x: bin(x)[2:] )

# Initialize the CNF class
F=CNF.CNF()


def to_bin_list(val,n_bits=N):
    # Convert a integer to a vector of binary values
    lb=bin(val)[2:][::-1]
    ans=[]
    for i in range(N):
        if i<len(lb):
            ans.append(1 if lb[i]=="1" else -1)
        else: ans.append(-1)
    return ans

r_var=[F.new_var() for p in range(n)] # variables indicating root node
F.only_one_in_all(r_var) # only one node can be the root

edge_set=[] # all possible edges
e_var={} # map the edges to variables

#identify all possible edges in ancestry graph
for p in range(n):
    for q in range(n):
        if(p==q): continue
        pq_flag=1
        for i in range(m):
            if (df[ (df["sample_index"]==i) & (df["character_index"]==q) ]["f-int"].iloc[0] > 
               df[ (df["sample_index"]==i) & (df["character_index"]==p) ]["f+int"].iloc[0]):
                pq_flag=0
                break
        if(pq_flag):
            edge_set.append( (p,q) )
            e_var[(p,q)]=F.new_var(ind=True)

def locate_label(index):
    label=df[df["character_index"]==index]["character_label"].iloc[0]
    return label

with open(evar_json,"w") as fjson:
    # save the mapping between the edges and the variables for further reference
    var2label={}
    var2index={}
    for edge in e_var.keys():
        character_labels = (locate_label(edge[0]),locate_label(edge[1]))
        var2index[ e_var[edge] ] = edge
        var2label[ e_var[edge] ] = character_labels
    json.dump({0:var2index,1:var2label},fjson,sort_keys=True,indent=4)

# For every node other than root, only one parent
for p in range(n):
    pp_var = [r_var[p]]
    for (pp,qq) in edge_set:
        if qq==p: pp_var.append(e_var[(pp,qq)])
    F.only_one_in_all(pp_var)

if (single): 
    # For normal PPM, we have frequencies decided
    f_var = [ [ to_bin_list(df[ (df["sample_index"]==i) & (df["character_index"]==p)]["f-int"].iloc[0]) for p in range(n)] for i in range(m)]
else :
    # For I-PPM, frequencies can vary
    f_var = [ [ [ F.new_var() for k in range(N)]  for p in range(n)] for i in range(m)]


# Sum condition
for p in range(n):
    children = [qq for (pp,qq) in edge_set if pp==p]
    for i in range(m):
        Sum = [F.false() for k in range(N) ]
        for q in children:
            tmp = [F.AND( e_var[(p,q)], f_var[i][q][k])  for k in range(N)  ]
            Sum = F.add(Sum,tmp)
        if (single):
            F.leq(Sum,f_var[i][p])
        else:
            f_upper_var=to_bin_list(df[ (df["sample_index"]==i) & (df["character_index"]==p)]["f+int"].iloc[0])
            F.leq(Sum,f_upper_var)
            f_lower_var=to_bin_list(df[ (df["sample_index"]==i) & (df["character_index"]==p)]["f-int"].iloc[0])
            F.eq(f_var[i][p],F.max_(Sum,f_lower_var) )


# cycle prevention
relation = [ [F.new_var() for _ in range(n)] for _ in range(n) ] # indicator variable: true if i is ancestor of j
for i in range(n):
    F.set_true(relation[i][i])

for i in range(n):
    for j in range(n):
        if (i==j): continue
        inj = []
        for ij in range(n):
            if (  ((ij,j) not in edge_set) or 
              ij == j
           ): continue
            #if (ij == j): continue
            inj.append(F.AND(relation[i][ij],e_var[(ij,j)]))
        F.ORList(inj,relation[i][j])

for i in range(n):
    for j in range(i):
        F.add_clause([-relation[i][j], -relation[j][i]])


if scs_file != None:
    label2index = {}
    def find_lab2idx(x):
        label2index[x["character_label"]] = x["character_index"]
    df.apply(find_lab2idx,axis = 1)
    scs_df = pd.read_csv(scs_file,sep="\t")
    scs = SCS.SCS(scs_df,F,label2index,relation)
#     if(range_upper==0):
#     else:
#         scs = SCS.SCS(scs_df,F,label2index,relation,"type1")
#         scs.FN_range(range_lower,range_upper)
    scs.get_clauses()

F.to_cnf_file(output)
