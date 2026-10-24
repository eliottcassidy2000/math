#!/usr/bin/env python3
"""Regenerate exact image witnesses for seven sparse native charges (adapted from Session45 multi-source solver) by bit-plane elimination."""
import argparse,csv,gzip,json
from pathlib import Path
from hashlib import sha256
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--out',type=Path,required=True);P.add_argument('--modulus',type=int,default=256);a=P.parse_args()
M=a.modulus
if M<2 or M&(M-1):raise RuntimeError('power of two modulus required')
K=M.bit_length()-1;N=3528;n=N-1
charge_records=json.loads(Path(__file__).with_name('charges.json').read_text());charges=[{int(k):v for k,v in row['charge'].items()} for row in charge_records]
path=a.repo/'05-knowledge/results/metagraph_order8_20261023/merged_edges.tsv.gz';raw=path.read_bytes();edges=[tuple(int(row[k]) for k in ['source','target','blue','black']) for row in csv.DictReader(gzip.decompress(raw).decode().splitlines(),delimiter='\t')]
rows=[[0]*K for _ in range(n)];degree=[0]*N
for i,j,b,k in edges:
 if i==j:continue
 w=b+k;degree[i]+=w;degree[j]+=w;v=(-w)%M
 for c in range(K):
  if v>>c&1:
   if i<n and j<n:rows[i][c]|=1<<j;rows[j][c]|=1<<i
for i in range(n):
 for c in range(K):
  if degree[i]>>c&1:rows[i][c]|=1<<i
for rhs,charge in enumerate(charges):
 for i,v in charge.items():
  for c in range(K):
   if (v%M)>>c&1:rows[i][c]|=1<<(n+rhs)
coeffmask=(1<<n)-1

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
for bits in range(K,0,-1):
 modulus=1<<bits;i=0
 while i<len(rows):
  p=next((j for j in range(i,len(rows)) if rows[j][0]&coeffmask),None)
  if p is None:break
  rows[i],rows[p]=rows[p],rows[i]
  col=((rows[i][0]&coeffmask)&-(rows[i][0]&coeffmask)).bit_length()-1
  val=sum(((row>>col)&1)<<bit for bit,row in enumerate(rows[i]))
  pivot=times(rows[i],pow(val,-1,modulus));rows[i]=pivot
  # Only generate multipliers actually requested in this column.
  scaled={0:(0,)*bits}
  for j in range(i+1,len(rows)):
   val=sum(((row>>col)&1)<<bit for bit,row in enumerate(rows[j]))
   if val:
    coeff=-val%modulus
    if coeff not in scaled:scaled[coeff]=times(pivot,coeff)
    rows[j]=plus(rows[j],scaled[coeff]);ops+=1
  pivots.append((bits,col,pivot));i+=1
 ranks.append(i);print('layer',bits,'unit pivots',i,'remaining',len(rows)-i,'rowops',ops,flush=True)
 for row in rows[i:]:
  if row[0]:raise RuntimeError('inconsistent odd right hand side')
 rows=[tuple(row[1:]) for row in rows[i:]]
solutions=[]
for rhs,charge in enumerate(charges):
 x=[0]*n
 for bits,col,row in reversed(pivots):
  modulus=1<<bits;b=sum(((plane>>(n+rhs))&1)<<bit for bit,plane in enumerate(row));total=0
  for bit,plane in enumerate(row):
   scan=plane&coeffmask&~(1<<col)
   while scan:
    low=scan&-scan;scan-=low;j=low.bit_length()-1;total+=(1<<bit)*x[j]
  x[col]=(b-total)%modulus
 x.append(0);image=[0]*N
 for i,j,b,k in edges:
  d=(b+k)*(x[i]-x[j]);image[i]+=d;image[j]-=d
 expected=[charge.get(i,0) for i in range(N)]
 if any((v-b)%M for v,b in zip(image,expected)):raise RuntimeError('independent full Laplacian residue check failed')
 solutions.append({'charge':sorted(charge.items()),'potential':x,'integer_correction_nonzero':[[i,(y-expected[i])//M] for i,y in enumerate(image) if y!=expected[i]]})
out={'status':'EXACT IMAGE WITNESSES','input_sha256':sha256(raw).hexdigest(),'modulus':M,'matrix_order':n,'unit_pivot_counts':dict(zip([1<<b for b in range(K,0,-1)],ranks)),'row_operations':ops,'wholly_unpivoted_columns':n-len(pivots),'solutions':solutions}
a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('DONE',flush=True)
