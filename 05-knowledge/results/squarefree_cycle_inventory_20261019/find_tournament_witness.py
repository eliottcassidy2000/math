import itertools,json
import numpy as np
pairs=list(itertools.combinations(range(6),2))
def matrix(mask):
 a=np.zeros((6,6),dtype=np.int64)
 for j,(u,v) in enumerate(pairs):
  if mask>>j&1:a[u,v]=1
  else:a[v,u]=1
 return a
def cyc(a,k):
 ans=0
 for vs in itertools.combinations(range(len(a)),k):
  for tail in itertools.permutations(vs[1:]):
   w=(vs[0],)+tail
   if all(a[w[i],w[(i+1)%k]] for i in range(k)):ans+=1
 return ans
want=(0,0,12,12,10,48)
found={}
for mask in range(1<<15):
 a=matrix(mask); p=np.eye(6,dtype=np.int64); traces=[]
 for k in range(1,7):
  p=p@a;traces.append(int(np.trace(p)))
 if tuple(traces)!=want:continue
 c6=cyc(a,6)
 if c6 not in found:found[c6]={'mask':mask,'adjacency':a.tolist(),'traces':traces,'c6':c6,'c3':cyc(a,3),'c4':cyc(a,4),'c5':cyc(a,5)}
 if len(found)>1:break
out={'encoding':'lexicographic pairs (u,v), u<v; bit1 means u->v, bit0 v->u','witnesses':list(found.values())}
for w in out['witnesses']:
 a=np.array(w['adjacency']);b=np.zeros((7,7),dtype=np.int64);b[:6,:6]=a;b[6,:6]=1;p=np.eye(7,dtype=np.int64);ts=[]
 for k in range(1,8):p=p@b;ts.append(int(np.trace(p)))
 w['source_extension']={'traces':ts,'c6':cyc(b,6),'c7':cyc(b,7)}
print(json.dumps(out,indent=2))
