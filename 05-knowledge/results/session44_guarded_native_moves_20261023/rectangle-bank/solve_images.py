#!/usr/bin/env python3
"""Exact dyadic row elimination with integer bit planes, three RHS at once."""
import argparse,csv,gzip,json
from pathlib import Path
from hashlib import sha256
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--out',type=Path,required=True);a=P.parse_args()
N=3528;n=N-1;charges=[(1238,2246),(881,2762),(2650,2843)];path=a.repo/'05-knowledge/results/metagraph_order8_20261023/merged_edges.tsv.gz'
raw=path.read_bytes();edges=[tuple(int(row[k]) for k in ['source','target','blue','black']) for row in csv.DictReader(gzip.decompress(raw).decode().splitlines(),delimiter='\t')]
rows=[[0,0,0] for _ in range(n)];degree=[0]*N
for i,j,b,k in edges:
 if i==j:continue
 w=b+k;degree[i]+=w;degree[j]+=w;v=(-w)%8
 for c in range(3):
  if v>>c&1:
   if i<n and j<n:rows[i][c]|=1<<j;rows[j][c]|=1<<i
for i in range(n):
 for c in range(3):
  if degree[i]>>c&1:rows[i][c]|=1<<i
for rhs,(i,j) in enumerate(charges):rows[i][2]|=1<<(n+rhs);rows[j][2]|=1<<(n+rhs)
coeffmask=(1<<n)-1;fullmask=(1<<(n+len(charges)))-1

def plus(x,y):
 carry=0;out=[]
 for l,r in zip(x,y):
  z=l^r;out.append(z^carry);carry=(l&r)|(z&carry)
 return tuple(out)
def times(x,c):
 k=len(x);out=(0,)*k
 for j in range(k):
  if c>>j&1:out=plus(out,(0,)*j+tuple(x[:k-j]))
 return out
pivots=[];ranks=[];ops=0
for bits in [3,2,1]:
 modulus=1<<bits;remaining=[];start_count=len(rows);i=0
 while i<len(rows):
  p=next((j for j in range(i,len(rows)) if rows[j][0]&coeffmask),None)
  if p is None:break
  rows[i],rows[p]=rows[p],rows[i]
  col=((rows[i][0]&coeffmask)&-(rows[i][0]&coeffmask)).bit_length()-1
  val=sum(((row>>col)&1)<<bit for bit,row in enumerate(rows[i]))
  pivot=times(rows[i],pow(val,-1,modulus));rows[i]=pivot
  scaled=[times(pivot,c) for c in range(modulus)]
  for j in range(i+1,len(rows)):
   val=sum(((row>>col)&1)<<bit for bit,row in enumerate(rows[j]))
   if val:rows[j]=plus(rows[j],scaled[-val%modulus]);ops+=1
  pivots.append((bits,col,pivot));i+=1
 ranks.append(i)
 print('layer',bits,'unit pivots',i,'remaining',len(rows)-i,'rowops',ops,flush=True)
 for row in rows[i:]:
  if row[0]:raise RuntimeError('odd residual right hand side: inconsistent system')
 rows=[tuple(row[1:]) for row in rows[i:]]
solutions=[]
for rhs in range(len(charges)):
 x=[0]*n
 for bits,col,row in reversed(pivots):
  modulus=1<<bits
  b=sum(((plane>>(n+rhs))&1)<<bit for bit,plane in enumerate(row))
  total=0
  for bit,plane in enumerate(row):
   scan=plane&coeffmask&~(1<<col)
   while scan:
    low=scan&-scan;scan-=low;j=low.bit_length()-1;total+=(1<<bit)*x[j]
  x[col]=(b-total)%modulus
 x.append(0);image=[0]*N
 for i,j,b,k in edges:
  d=(b+k)*(x[i]-x[j]);image[i]+=d;image[j]-=d
 expected=[0]*N;u,v=charges[rhs];expected[u]=4;expected[v]=-4
 if any((v-b)%8 for v,b in zip(image,expected)):raise RuntimeError('independent full Laplacian residue check failed')
 solutions.append({'charge':charges[rhs],'coefficient':4,'potential_mod8':x,'image_mod8_nonzero':[[i,y%8] for i,y in enumerate(image) if y%8],'image_equals_rhs_mod8':True,'integer_correction_div8_nonzero':[[i,(y-expected[i])//8] for i,y in enumerate(image) if y!=expected[i]]})
out={'status':'EXACT IMAGE WITNESSES','input_sha256':sha256(raw).hexdigest(),'matrix_order':n,'unit_pivot_counts_by_modulus':dict(zip([8,4,2],ranks)),'row_operations':ops,'wholly_unpivoted_columns':n-len(pivots),'solutions':solutions}
a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('DONE',flush=True)
