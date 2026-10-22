"""Exhaustive almostregular members of fixed transitive switching orbit."""
import numpy as np,itertools,json,pathlib,subprocess,collections
from experiment import bits,pred,brute
ROOT=pathlib.Path(__file__).resolve().parent
p=subprocess.Popen([str(ROOT/'count')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
out=[]
for n in [2,4,6,8,10,12,14]:
 m=n//2;S=np.triu(-np.ones((n,n),dtype=int),1)+np.tril(np.ones((n,n),dtype=int),-1);hist=collections.Counter();rows=[]
 for tail in itertools.product([-1,1],repeat=m-1):
  for mode in [0,2]:
   d=[1,-1] if mode==0 else [1]
   for e in tail:d.extend([e,-e])
   if mode==2:d.append(1)
   d=np.array(d);A=d[:,None]*S*d[None,:];assert np.all(abs(A.sum(1))==1)
   p.stdin.write(f'{n} {bits(A)}\n');p.stdin.flush();H=int(p.stdout.readline());hist[H]+=1
   if n<=6:assert H==brute(A)
   rows.append(dict(n=n,signs=d.tolist(),bits=bits(A),H=str(H),**pred(A,H)))
 assert len(rows)==2**m;out.extend(rows)
 print(n,dict(sorted(hist.items())),flush=True)
with open(ROOT/'switching_family.jsonl','w') as f:
 for r in out:f.write(json.dumps(r)+'\n')
p.stdin.close();p.wait()
