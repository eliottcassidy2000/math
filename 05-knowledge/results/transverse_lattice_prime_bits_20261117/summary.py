"""Count distinct realizable jets; use the archived complete class catalogues."""
import argparse,csv,gzip,io,subprocess
from collections import Counter,defaultdict
from itertools import permutations,product
from math import comb,factorial
from pathlib import Path
N8='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities'
N9='05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz'

def need(x):
    if not x:raise ArithmeticError('census control failed')

def canonical(a):
 n=len(a)
 cells=[tuple(i for i,r in enumerate(a) if r.bit_count()==k) for k in range(n)]
 return min(tuple(sum((a[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
       for z in product(*(permutations(c) for c in cells if c))
       for p in [sum(z,())])

def small():
 a={(0,)}
 for n in range(1,7):
  yield a
  if n<6:
   a={canonical(tuple(row+((1-(v>>i&1))<<n) for i,row in enumerate(T))+(v,))
      for T in a for v in range(1<<n)}

def rows9(x):
 out=[0]*9;k=35
 for i in range(9):
  for j in range(i+1,9):
   if x>>k&1:out[i]|=1<<j
   else:out[j]|=1<<i
   k-=1
 return out

def jets(n,G):
 p=0; cols=[[] for _ in range(n//3)] # h=0,1,..., maximal nontrivial exponent
 for j in range(1,(n-1)//2+1):
  a=min(j,n-2*j); B=G[p:p+a];p+=a
  for k in range(len(cols)):
   cols[k].append(sum(comb(c,k)*2**(c-k)*B[c-1] for c in range(max(1,k),a+1)))
 return tuple(tuple(v) for v in cols)

def report(n,seen,mode):
 full={G:i for i,G in seen}; T=list(full)
 J=[jets(n,G) for G in T]
 maps=[defaultdict(list) for _ in range(max(1,n//3))]
 for i,v in enumerate(J):
  for d,B in enumerate(maps,1):B[v[:d]].append(i)
 counts=tuple(len(B) for B in maps)
 largest=max(map(len,maps[0].values()))
 print(n,mode,len(seen),counts,largest,(largest-1).bit_length(),flush=True)
 if n==9 and mode=='all':
  for d,B in enumerate(maps,1):
   if d>=2 and max(map(len,B.values()))>1:
    V=next(v for v in B.values() if len(v)>1)
    print('deeper pair n9',d,[(full[T[i]],rows9(full[T[i]]),T[i]) for i in V[:2]],flush=True)
  if largest>1:
   V=max(maps[0].values(),key=len)
   print('largest n9 diagonal fibre',[(full[T[i]],T[i]) for i in V],flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('repo');p.add_argument('binary');p.add_argument('--ref',default='HEAD');args=p.parse_args()
 def rd(f):return subprocess.check_output(['git','-C',args.repo,'show',args.ref+':'+f])
 smalls=list(small())
 for n in range(1,10):
  H9={};H8={}
  if n<=6:
   inp=''.join(' '.join(map(str,T))+'\n' for T in smalls[n-1]).encode()
  elif n==7:inp=rd(N8+'/n7_representatives.txt')
  elif n==8:
   tab=list(csv.DictReader(io.StringIO(rd(N8+'/n8_classes.csv').decode())))
   inp=''.join(' '.join(str(int(r['code'])>>(8*i)&255) for i in range(8))+'\n' for r in tab).encode()
   H8={i+1:int(r['H']) for i,r in enumerate(tab)}
  else:
   inp=gzip.decompress(rd(N9))
   H9={int(row['code']):int(row['Aut']) for row in csv.DictReader(io.StringIO(inp.decode()),delimiter='\t')}
  out=subprocess.check_output([args.binary,str(n),'code' if n==9 else 'rows'],input=inp).decode()
  all_= [];strong=[];mx=0
  for l in out.splitlines():
   id,b,*a=map(int,l.split());G=tuple(a);all_.append((id,G))
   if b:strong.append((id,G))
   # Compute H from all entries, check archive or global labelled mass.
   if n>=8:
    h=1+sum(jets(n,G)[0])
    if n==8:need(H8[id]==h)
    else:mx+=factorial(n)//H9[id]*h
  need(len(all_)==(1,1,2,4,12,56,456,6880,191536)[n-1])
  if n==9:need(mx==factorial(n)*(1<<comb(n-1,2)))
  report(n,all_,'all')
  if strong:report(n,strong,'strong')

if __name__=='__main__':main()
