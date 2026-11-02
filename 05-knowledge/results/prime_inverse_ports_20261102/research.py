#!/usr/bin/env python3
"""Exact marked-module response, cyclic orbit controls and prime filters."""
import argparse,itertools,json,subprocess,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
GEN='/opt/homebrew/Cellar/nauty/2.9.3/bin/gentourng'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def decode(code,n):
 rows=[0]*n
 for bit,(a,b) in zip(code,itertools.combinations(range(n),2)):
  if bit=='1':rows[a]|=1<<b
  else:rows[b]|=1<<a
 return rows
def all_h(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)];hs=[1]+[0]*((1<<n)-1)
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1:
    R=S^(1<<v)
    dp[S][v]=1 if not R else sum(dp[R][u] for u in range(n) if rows[u]>>v&1)
  hs[S]=sum(dp[S])
 return hs
def strong(rows):
 n=len(rows)
 for r in (rows,[sum((1<<j) for j in range(n) if rows[j]>>i&1) for i in range(n)]):
  S=1
  while True:
   new=S
   for u in range(n):
    if S>>u&1:new|=r[u]
   if new==S:break
   S=new
  if S!=(1<<n)-1:return False
 return True
def ordered_covers(rows):
 hs=all_h(rows);n=len(rows);full=(1<<n)-1;prev=[1]+[0]*full;cs=[0]
 for k in range(1,n+1):
  now=[0]*(full+1)
  for S in range(1,full+1):
   U=S
   while U:
    now[S]+=hs[U]*prev[S^U];U=(U-1)&S
  cs.append(now[full]);prev=now
 return cs

def skeleton(rows,v):
 """b[k]: quotient path words with k nonadjacent copies of marked v."""
 ext=[u for u in range(len(rows)) if u!=v];m=len(ext);n=m+1
 dp=[[[0]*(n+1) for _ in range(n)] for _ in range(1<<m)]
 dp[0][m][1]=1
 for S in range(1<<m):
  if S:
   for j,u in enumerate(ext):
    if not(S>>j&1):continue
    R=S^(1<<j)
    if R==0:dp[S][j][0]=1
    for z,w in enumerate(ext+[v]):
     if rows[w]>>u&1:
      for k,a in enumerate(dp[R][z]):dp[S][j][k]+=a
   for j,u in enumerate(ext):
    if S>>j&1 and rows[u]>>v&1:
     for k in range(n):dp[S][m][k+1]+=dp[S][j][k]
 return [sum(x[k] for x in dp[-1]) for k in range(n+1)]

def substitute(q,v,B):
 ext=[u for u in range(len(q)) if u!=v];rows=list(B);b=len(B);rows += [0]*len(ext)
 for i,u in enumerate(ext):
  for j,w in enumerate(ext):
   if q[u]>>w&1:rows[b+i]|=1<<(b+j)
  if q[v]>>u&1:
   for t in range(b):rows[t]|=1<<(b+i)
  else:rows[b+i]|=(1<<b)-1
 return rows

def circ(p,bits):
 D=[i if not t else p-i for i,t in enumerate(bits,1)]
 return [sum(1<<((u+d)%p) for d in D) for u in range(p)]
def prime(n):return n>1 and all(n%d for d in range(2,math.isqrt(n)+1))
def census(maxn):
 blocks=[]
 for p in [3,5,7]:
  seen=set()
  for bits in itertools.product([0,1],repeat=(p-1)//2):
   rows=circ(p,bits);cs=ordered_covers(rows)
   if tuple(cs) in seen:continue
   seen.add(tuple(cs));blocks.append({'p':p,'rows':rows,'covers':cs,'H':cs[1],'E':2*cs[1]+cs[2]})
 out={'blocks':blocks,'by_order':{},'threshold_30625_symmetric_7p':{},'checks':0}
 for n in range(3,maxn+1):
  source=HERE/'inputs'/f'quotients_n{n}.txt'
  lines=source.read_text().split() if source.exists() else subprocess.check_output([GEN,'-q',str(n)],text=True).split()
  totals=0;front=[{} for _ in blocks];polys=set()
  for code in lines:
   q=decode(code,n)
   if not strong(q):continue
   totals+=1
   for v in range(n):
    b=skeleton(q,v);polys.add(tuple(b))
    for i,B in enumerate(blocks):
     value=sum(x*y for x,y in zip(b,B['covers']))
     f=front[i]
     if 'min' not in f or value<f['min']:f.update(min=value,quotient_rows=q,port=v,polynomial=b)
     if B['p']==7 and value<30625 and prime(value//7):
      need(value%7==0,'7 divides module path count')
      key=str(value);out['threshold_30625_symmetric_7p'].setdefault(key,dict(p=value//7,quotient_rows=q,port=v,block_rows=B['rows'],polynomial=b))
     # Separate full-adjacency consumer, on every marked quotient through n5 for C3,
     # and every marked quotient through n4 for both seven-cycles.
     if (n<=5 and B['p']==3) or (n<=4 and B['p']==7):
      full=substitute(q,v,B['rows']);need(all_h(full)[-1]==value,'module formula');out['checks']+=1
  out['by_order'][str(n)]={'strong_classes':totals,'marked_classes':totals*n,'distinct_polynomials':len(polys),'minima':front}
  print('order',n,'strong',totals,'minima',[x['min'] for x in front],flush=True)
 (HERE/'module_census.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print('seven-prime targets',sorted(map(int,out['threshold_30625_symmetric_7p'])))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--max-order',type=int,default=7);a=ap.parse_args();census(a.max_order)
