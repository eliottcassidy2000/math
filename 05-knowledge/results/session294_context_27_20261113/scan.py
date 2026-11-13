from pathlib import Path
from math import factorial,gcd
from collections import Counter
import json,csv
import sys
BASE=Path(sys.argv[1] if len(sys.argv)>1 else '.').resolve()/'05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities'
def profile(rows):
 n=len(rows); M=1<<n
 D=[[0]*n for _ in range(M)];H=[0]*M;H[0]=1
 for i in range(n):D[1<<i][i]=1
 for S in range(1,M):
  H[S]=sum(D[S])
  for i,q in enumerate(D[S]):
   if q:
    R=rows[i]&~S
    while R:
     b=R&-R;R^=b;D[S|b][b.bit_length()-1]+=q
 if H[-1]%27:return H[-1],None
 F=[1]+[0]*(M-1);out=[]
 for k in range(1,n+1):
  G=[0]*M
  for S in range(1,M):
   T=S
   while T:
    G[S]+=H[T]*F[S^T];T=(T-1)&S
  out.append(G[-1]);F=G
 return H[-1],out
def v(x):
 a=0
 while x%3==0:x//=3;a+=1
 return a
def scan(items):
 cnt=Counter(); found=[]
 for rows,known,aut in items:
  if known is not None and known%27: continue
  h,F=profile(rows)
  if known is not None and known!=h:raise ValueError('H mismatch')
  if F is not None:
   raw=gcd(*F); ctx=gcd(*(f*factorial(i) for i,f in enumerate(F)))
   cnt[v(raw),v(ctx)]+=1
   if v(ctx)>=3:found.append(dict(rows=rows,H=h,raw=raw,context=ctx,F=F,aut=aut))
 return dict(eligible={str(k):z for k,z in sorted(cnt.items())},found=found)
if __name__=='__main__':
 a=((list(map(int,L.split())),None,None) for L in (BASE/'n7_representatives.txt').read_text().splitlines())
 b=(([(int(r['code'])>>(8*i))&255 for i in range(8)],int(r['H']),int(r['aut'])) for r in csv.DictReader((BASE/'n8_classes.csv').open()))
 print(json.dumps({'order7':scan(a),'order8':scan(b)},indent=2))
