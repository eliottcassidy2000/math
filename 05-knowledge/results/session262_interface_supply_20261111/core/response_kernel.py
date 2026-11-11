#!/usr/bin/env python3
"""Exact compressed-run coefficients for a literal two-module response."""
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import json,itertools
Q=(198,196,136,211,167,207,20,64)
PORTS=(5,7)

def skeleton(rows,end=None):
    n=len(rows); ext=[v for v in range(n) if v not in PORTS]; index={v:i for i,v in enumerate(ext)}
    # A repeated module cannot follow itself; each of its runs needs a separating
    # exterior vertex, except at most one separator in the other module.
    # Thus r,s <= number of exterior vertices+2 is a safe finite bound.
    bound=len(ext)+2; full=(1<<len(ext))-1
    states={(0,0,0,-1):1}; out=defaultdict(int)
    for length in range(len(ext)+2*bound+1):
      nxt=defaultdict(int)
      for (mask,a,b,last),num in states.items():
        if mask==full and a and b and (end is None or last==end):out[a,b]+=num
        for v in range(n):
          if last>=0 and (v==last or not rows[last]>>v&1):continue
          if v==5:
            if a<bound:nxt[mask,a+1,b,v]+=num
          elif v==7:
            if b<bound:nxt[mask,a,b+1,v]+=num
          else:
            bit=1<<index[v]
            if not mask&bit:nxt[mask|bit,a,b,v]+=num
      states=nxt
    return dict(out)

def born(rows,root=0):
    n=len(rows)
    return tuple(row|(0 if i==root else 1<<n) for i,row in enumerate(rows))+(1<<root,)

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def substitute(A,B):
    modules=[A if i==5 else B if i==7 else (0,) for i in range(8)]
    starts=[];pos=0
    for M in modules:starts.append(pos);pos+=len(M)
    out=[]
    for i,M in enumerate(modules):
      outside=sum(((1<<len(modules[j]))-1)<<starts[j] for j in range(8) if Q[i]>>j&1)
      out.extend(outside|(x<<starts[i]) for x in M)
    return tuple(out)

def ham(rows):
    n=len(rows);F=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
      for v in range(n):
        num=dp[mask][v]
        if not num:continue
        remain=rows[v]&~mask&F
        while remain:
          bit=remain&-remain;remain-=bit;w=bit.bit_length()-1
          dp[mask|bit][w]+=num
    return sum(dp[-1])

def covers(rows):
    n=len(rows);H={0:0}
    for m in range(1,1<<n):
      verts=[i for i in range(n) if m>>i&1]
      R=tuple(sum(1<<j for j,w in enumerate(verts) if rows[v]>>w&1) for v in verts)
      H[m]=ham(R)
    @lru_cache(None)
    def c(mask,r):
      if r==0:return int(mask==0)
      if mask.bit_count()<r:return 0
      total=0;s=mask
      while s:
        total+=H[s]*c(mask^s,r-1);s=(s-1)&mask
      return total
    return {r:c((1<<n)-1,r) for r in range(1,n+1)}

def contract(coeff,A,B):
    ca,cb=covers(A),covers(B)
    return sum(v*ca.get(r,0)*cb.get(s,0) for (r,s),v in coeff.items())

def main():
    h=skeleton(Q);f=skeleton(born(Q));end=skeleton(Q,end=0);u={key:f.get(key,0)-h.get(key,0)+end.get(key,0) for key in f.keys()|h.keys()|end.keys()};u={k:v for k,v in u.items() if v}
    if any(v<0 for v in u.values()):raise RuntimeError('negative U coefficient')
    modules=[(0,),tt(2),tt(3),(2,4,1)]
    controls=[]
    for A,B in itertools.product(modules,repeat=2):
      T=substitute(A,B);F=born(T)
      hh,ff=ham(T),ham(F)
      if hh!=contract(h,A,B) or ff!=contract(f,A,B):raise RuntimeError('literal coefficient disagreement')
      controls.append({'A':A,'B':B,'H':hh,'U0':contract(u,A,B),'E0':contract(end,A,B),'H_after_root_birth':ff})
    receipt={'status':'PASS','scope':'Exact finite skeleton coefficients give all nonempty-module formulas; 16 literal controls only','Q':Q,'ports':PORTS,'root':0,'H_coefficients':[[*k,v] for k,v in sorted(h.items())],'U_coefficients':[[*k,v] for k,v in sorted(u.items())],'E_coefficients':[[*k,v] for k,v in sorted(end.items())],'F_coefficients':[[*k,v] for k,v in sorted(f.items())],'literal_controls':controls}
    Path(__file__).with_name('response_kernel_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
