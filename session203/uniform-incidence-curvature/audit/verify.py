#!/usr/bin/env python3
"""Independent starting DP plus literal DFS; independent mandatory-gap census."""
from pathlib import Path
from itertools import permutations
from math import comb,factorial
import hashlib,json,random,subprocess,sys
BASE=Path(__file__).resolve().parent

def need(b,s):
 if not b:raise RuntimeError(s)

def core(n,rng):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if rng.randrange(2):a[i]|=1<<j
   else:a[j]|=1<<i
 return a

def graph(a,k,m,p,u,word):
 c=len(a);b=a[:]+[0]*m
 for j in range(m):
  for i in range(c):
   outgoing=bool((word>>j&1)^(p>>i&1)) if i<k else bool(u>>(i-k)&1)
   if outgoing:b[i]|=1<<(c+j)
   else:b[c+j]|=1<<i
  for i in range(j):b[c+i]|=1<<(c+j)
 return b

def gaps(rows,mask):
 n=len(rows);ans=[0]*(n+2)
 for p in permutations(range(n)):
  allowed=mandatory=0;valid=True
  for j in range(n+1):
   left=p[j-1] if j else None;right=p[j] if j<n else None
   e=left is None or right is None or bool(rows[left]>>right&1)
   f=(left is None or bool(mask>>left&1)) and (right is None or not(mask>>right&1))
   if not e and not f:valid=False;break
   allowed+=bool(f);mandatory+=not e
  if valid:
   for r in range(mandatory,allowed+1):ans[r]+=comb(allowed-mandatory,r-mandatory)
 return ans

def v2(n):
 n=abs(n);need(n>0,'zero valuation');return (n&-n).bit_length()-1

def main():
 binary=sys.argv[1];rng=random.Random(203110608);cases=[]
 for k in range(5):
  for extra in range(4 if k<3 else 2):
   for outside in range(3 if k<3 else 2):
    reps=2 if k<3 else 1
    for _ in range(reps):
     m=2*k+outside if k else 3; a=core(k+extra,rng)
     cases.append({'k':k,'m':m,'core':a,'polarity':rng.randrange(1<<k),'uniform':rng.randrange(1<<extra)})
 # Explicit all polarities on one three-sensitive context, plus arbitrary exterior bits.
 for p in range(8):cases.append({'k':3,'m':7,'core':[2,4,1,7],'polarity':p,'uniform':p%2})
 inp=str(len(cases))+'\n'+''.join(f"{len(c['core'])} {c['k']} {c['m']} {c['polarity']} {c['uniform']} "+' '.join(map(str,c['core']))+'\n' for c in cases)
 raw=subprocess.check_output([binary],input=inp.encode()).decode();lines=raw.splitlines();need(len(lines)==len(cases),'output length')
 report=[];graphs=direct=blocks=0
 for ix,(c,line) in enumerate(zip(cases,lines)):
  vals=list(map(int,line.split()));need(vals[0]==ix,'wrong case');vals=vals[1:];k,m=c['k'],c['m'];need(len(vals)==1<<m,'cube length');graphs+=len(vals)
  if len(c['core'])+m<=8:direct+=len(vals)
  pol=vals[:]
  for j in range(m):
   for b in range(1<<m):
    if b>>j&1:pol[b]-=pol[b^(1<<j)]
  degree=max(b.bit_count() for b,z in enumerate(pol) if z)
  need(degree==2*k,'degree mismatch')
  results=[]
  if k:
   for start in range(m-2*k+1):
    free=list(range(start,start+2*k));ext=[j for j in range(m) if j not in free]
    G=graph(c['core'],k,m,c['polarity'],c['uniform'],0);cc=len(c['core']);W=list(range(k,cc))+[cc+j for j in ext];rows=[sum(1<<j for j,w in enumerate(W) if G[v]>>w&1) for v in W]
    representative=cc+free[0];mask=sum(1<<j for j,v in enumerate(W) if G[v]>>representative&1)
    # Every W vertex is genuinely uniform toward the selected consecutive block.
    need(all(all(bool(G[v]>>(cc+j)&1)==bool(mask>>i&1) for j in free) for i,v in enumerate(W)),'nonuniform W')
    ar=gaps(rows,mask);Q=sum(ar[r]*comb(k-1,r-1)*2**(r-1) for r in range(1,min(k,len(ar)-1)+1));C=factorial(k)*factorial(2*k)//2**(k-1);pred=(-1)**k*C*Q
    with_x=[a|((1<<len(rows)) if mask>>i&1 else 0) for i,a in enumerate(rows)]+[((1<<len(rows))-1)^mask]
    literal_a1=sum(all(with_x[p[i]]>>p[i+1]&1 for i in range(len(rows))) for p in permutations(range(len(rows)+1)))
    need(ar[1]==literal_a1 and literal_a1%2==1,'a1 singleton identity')
    need(Q>0 and Q%2==1,'nonpositive/even Q')
    for eb in range(1<<len(ext)):
     base=sum((eb>>i&1)<<j for i,j in enumerate(ext));mixed=0
     for b in range(1<<(2*k)):
      word=base|sum((b>>i&1)<<j for i,j in enumerate(free));mixed+=(-1)**b.bit_count()*vals[word]
     need(mixed==pred,'mixed difference mismatch');need(v2(mixed)==1+2*v2(factorial(k)),'valuation mismatch');blocks+=1
    results.append({'block_start':start,'gap_coefficients':ar,'Q':Q,'coefficient':pred,'valuation':v2(pred),'exterior_assignments':1<<len(ext)})
  report.append({**c,'degree':degree,'graphs':len(vals),'blocks':results})
 # Load-bearing consecutivity hostile: one sensitive vertex, four spine vertices.
 target=next((c for c in report if c['k']==1 and c['m']==4 and len(c['core'])==1),None);need(target is not None,'missing hostile context')
 idx=next(i for i,c in enumerate(cases) if c['k']==1 and c['m']==4 and len(c['core'])==1)
 values=list(map(int,lines[idx].split()))[1:];nonconsecutive=values[0]-values[1]-values[8]+values[9]
 need(nonconsecutive==-4,'unexpected nonconsecutive control')
 
 # Position-dependent polarity hostile, outside the theorem's context class.
 nonuniform=[]
 for word in range(16):
  a=graph([2,0],2,4,0,0,word);a[1]^=1<<2;a[2]^=1<<1
  nonuniform.append(sum(all(a[p[i]]>>p[i+1]&1 for i in range(5)) for p in permutations(range(6))))
 npoly=nonuniform[:]
 for j in range(4):
  for b in range(16):
   if b>>j&1:npoly[b]-=npoly[b^(1<<j)]
 need(npoly[-1]==0,'nonuniform top coefficient expected zero')
 odd_control=[]
 for word in range(4):
  p=graph([0],1,2,0,0,word);n=len(p);c=[a<<3 for a in p]
  triangle=[2,4,1];full=triangle+c
  for j in range(3):full[j]|=((1<<n)-1)<<3
  hp=sum(all(p[a[i]]>>a[i+1]&1 for i in range(n-1)) for a in permutations(range(n)))
  hc=sum(all(full[a[i]]>>a[i+1]&1 for i in range(n+2)) for a in permutations(range(n+3)))
  need(hc==3*hp,'odd multiplier control')
  odd_control.append({'bits':word,'parent_H':hp,'child_H':hc})
 ndegree=max(b.bit_count() for b,z in enumerate(npoly) if z)
 print(json.dumps({'status':'PASS','odd_multiplier_positive_control':odd_control,'nonuniform_hostile':{'core':[2,0],'k':2,'m':4,'changed_incidence':'sensitive vertex1 complements only spine position0','H_values':nonuniform,'degree':ndegree,'top_coefficient':npoly[-1]},'contexts':len(cases),'literal_graphs':graphs,'DFS_checked_graphs':direct,'consecutive_block_checks':blocks,'nonconsecutive_hostile':{'k':1,'free_positions':[0,3],'middle_bits':[0,0],'coefficient':nonconsecutive,'valuation':v2(nonconsecutive)},'contexts_detail':report},indent=2,sort_keys=True))
if __name__=='__main__':main()
