import json
from itertools import combinations
from pathlib import Path

def block(mask,n=5):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(mask>>k)&1;a[j][i]=1-a[i][j]
 return a

def H(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    if t:d[s][v]=sum(d[t][u] for u in range(n) if t>>u&1 and a[u][v])
 return sum(d[-1])

def join(bits,seams):
 n=5*len(bits);a=[[int(i<j) for j in range(n)] for i in range(n)]
 for k,b in enumerate(bits):
  z=block(8+2*b)
  for i in range(5):
   for j in range(5):a[5*k+i][5*k+j]=z[i][j]
 for k,p,q in seams:a[5*k+p][5*(k+1)+q]=0;a[5*(k+1)+q][5*k+p]=1
 return a

if __name__=='__main__':
 rows=[]
 for p in range(5):
  for q in range(5):
   hs=[H(join(bits,[(0,p,q)])) for bits in [(0,0),(0,1),(1,0),(1,1)]]
   rows.append(dict(ports=[p,q],H=hs,mixed=hs[3]-hs[2]-hs[1]+hs[0]))
 out=dict(seed_H=[H(block(8)),H(block(10))],seams=rows)
 Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
