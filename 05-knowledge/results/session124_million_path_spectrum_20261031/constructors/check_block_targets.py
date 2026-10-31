from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
import json
from itertools import combinations
from array import array
for line in open(DATA_DIR / 'block_target_witnesses.jsonl'):
 r=json.loads(line);n,code,cut=r['seed'];out=[0]*n;bits=[(cut>>i)&1 for i in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j)if code>>k&1 else(j,i);out[a]|=1<<b
 for block in r['blocks']:
  old=len(out);size=3 if block.startswith('C') else 1
  out.extend([0]*size)
  for i in range(old):
   for j in range(old,old+size):out[i]|=1<<j
  if size==3:
   out[old]|=1<<(old+1);out[old+1]|=1<<(old+2);out[old+2]|=1<<old
   bits.extend([int(i<int(block[1]))for i in range(3)])
  else:bits.append(int(block))
 N=len(out)+1;out.append(0)
 for i,b in enumerate(bits):
  a,c=(N-1,i)if b else(i,N-1);out[a]|=1<<c
 d=array('Q',[0])*(N*(1<<N))
 for i in range(N):d[N*(1<<i)+i]=1
 for S in range(1,1<<N):
  base=N*S;vs=S
  while vs:
   vb=vs&-vs;vs-=vb;v=vb.bit_length()-1;c=d[base+v]
   if c:
    rem=out[v]&~S
    while rem:
     b=rem&-rem;rem-=b;u=b.bit_length()-1;d[N*(S|b)+u]+=c
 actual=sum(d[-N:])
 if actual!=r['target']:raise RuntimeError((r,actual))
 print('independent_full_adjacency_DP',N,actual,'PASS',flush=True)
