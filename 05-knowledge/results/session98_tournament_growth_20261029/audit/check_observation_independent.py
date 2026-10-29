#!/usr/bin/env python3
"""Independent evaluation-matrix lower bounds for claimed rational observation ranks."""
import json,random
P=1000000007

def need(b,s):
 if not b:raise RuntimeError(s)
def oddparts(n,lo=1):
 if n==0:yield ()
 for j in range(lo,n+1,2):
  for q in oddparts(n-j,j):yield (j,)+q
def mul(a,b,h):return [sum(a[j]*b[i-j] for j in range(i+1))%P for i in range(h+1)]
def col(part,b,h):
 a=[1]+[0]*h
 for d in part:
  v=[0]*(h+1)
  for i in range(1,h//d+1):v[d*i]=b[i]
  a=mul(a,v,h)
 return a
def rank(rows):
 basis={}
 for row in rows:
  row=row[:]
  for j in range(len(row)):
   if not row[j]:continue
   if j not in basis:
    inv=pow(row[j],P-2,P);basis[j]=[v*inv%P for v in row];break
   f=row[j];row=[(x-f*y)%P for x,y in zip(row,basis[j])]
 return len(basis)
def main():
 rng=random.Random(9801);out=[]
 # Claimed ranks are explicit regression targets, not imported algorithms.
 expected=[(1,1),(2,2),(3,3),(5,4),(7,6),(10,8),(13,10),(18,13),(23,16),(30,19),(37,23),(47,28),(57,33)]
 caps=[0,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168,48542114686912]
 for h,(formal_rank,ordinary_rank) in enumerate(expected):
  parts=[p for n in range(h+1) for p in oddparts(n)];actual=[]
  for bounded in (False,True):
   rows=[]
   for _ in range(40):
    b=[0]+[rng.randrange(caps[i]+1) if bounded else rng.randrange(P) for i in range(1,h+1)]
    columns=[col(p,b,h) for p in parts];rows.extend([list(r) for r in zip(*columns)])
   actual.append(rank(rows))
  need(actual==[formal_rank,ordinary_rank],f'rank lower bound mismatch h={h}')
  out.append({'h':h,'evaluation_ranks_mod_prime':actual})
 print(json.dumps({'status':'PASS','prime':P,'interpretation':'Independent lower bounds; exact symbolic kernels supply upper bounds','ranks':out},indent=2))
if __name__=='__main__':main()
