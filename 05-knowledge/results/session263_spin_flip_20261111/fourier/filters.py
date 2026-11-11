#!/usr/bin/env python3
"""Exact new-operation controls. No external dependencies or old-result audit."""
import itertools, json
from pathlib import Path

C3=(2,4,1)
F8=(234,236,161,228,111,192,132,16)
PARTIAL8=(42,120,99,68,13,216,17,95)
def valid(a):
 n=len(a)
 return all(not(a[i]>>i&1) and a[i]>>n==0 for i in range(n)) and all(((a[i]>>j)^(a[j]>>i))&1 for i in range(n) for j in range(i))
def join(a,b):
 n=len(a);right=((1<<len(b))-1)<<n
 return tuple(x|right for x in a)+tuple(x<<n for x in b)
def conv(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def hp_all(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*((1<<n)-1)
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not s>>v&1:continue
   p=s^(1<<v)
   if p:d[s][v]=sum(d[p][u] for u in range(n) if p>>u&1 and a[u]>>v&1)
  h[s]=sum(d[s])
 return h

def hp_all_permutations(a):
 # Independent engine: every ordered distinct-vertex tuple, not a DP recurrence.
 n=len(a);h=[0]*(1<<n);h[0]=1
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):h[sum(1<<v for v in p)]+=1
 return h

def split(a,h=None):
 h=hp_all(a) if h is None else h;n=len(a);full=(1<<n)-1;c=[0]*(n+1)
 for s in range(full+1):c[s.bit_count()]+=h[s]*h[full^s]
 return c

def rooted(a,y):
 n=len(a);h=hp_all(a);full=(1<<n)-1;d=[[0]*n for _ in range(1<<n)];d[1<<y][y]=1;p=[0]*(n+1)
 for s in range(1,1<<n):
  if not s>>y&1:continue
  for v in range(n):
   if v==y or not s>>v&1:continue
   prev=s^(1<<v)
   d[s][v]=sum(d[prev][u] for u in range(n) if prev>>u&1 and a[u]>>v&1)
  p[s.bit_count()]+=sum(d[s])*h[full^s]
 return p

def residues(p,m):return [sum(p[r::m]) for r in range(m)]
def cycles(a):
 n=len(a);out=[]
 for k in range(3,n+1,2):
  for sub in itertools.combinations(range(n),k):
   for tail in itertools.permutations(sub[1:]):
    p=(sub[0],)+tail
    if all(a[p[i]]>>p[(i+1)%k]&1 for i in range(k)):out.append((sum(1<<v for v in p),k))
 return out

def jet_packing(a):
 cs=cycles(a);n=len(a);best=n;lead=1
 def rec(i,used,order,weight):
  nonlocal best,lead
  for j in range(i,len(cs)):
   s,k=cs[j]
   if used&s:continue
   o=order-k+1;w=weight*2*k
   if o<best:best,lead=o,w
   elif o==best:lead+=w
   rec(j+1,used|s,o,w)
 rec(0,0,n,1)
 return best,lead,len(cs)

def divide_xplus1(p):
 q=[0]*(len(p)-1);q[0]=p[0]
 for i in range(1,len(q)):q[i]=p[i]-q[i-1]
 assert p[-1]==q[-1]
 return q

def jet(p):
 t=0
 while sum(x*((-1)**i) for i,x in enumerate(p))==0:
  p=divide_xplus1(p);t+=1
 return t,sum(x*((-1)**i) for i,x in enumerate(p))

def close(a,y,t):
 # B => TT_t => z; z beats y and loses to all other old vertices.
 n=len(a);b=join(a,tuple(sum(1<<j for j in range(i+1,t)) for i in range(t)));z=n+t
 b=list(x|(1<<z) for x in b);b[y]^=1<<z;b.append(1<<y)
 return tuple(b)

def main():
 report={'status':'PASS','blocks':{},'rooted_controls':[]}
 for name,a in [('triangle',C3),('eight_filter',F8),('highest_channel_only',PARTIAL8)]:
  assert valid(a)
  h=hp_all(a);ind=hp_all_permutations(a);assert h==ind
  p=split(a,h);j=jet(p);jp=jet_packing(a);assert j==jp[:2]
  report['blocks'][name]={'rows':a,'all_induced_subsets':len(h),'independent_counts_equal':True,'H':h[-1],'K':p,'K_at_1':sum(p),'residues_mod8':residues(p,8),'parity_jet':j,'odd_cycles':jp[2]}
 assert split(F8)==conv([231,231],[1]*8)
 assert split(C3)==[3]*4
 assert len(set(residues(split(PARTIAL8),8)))>1
 assert all(sum(split(PARTIAL8)[r::8])==sum(split(PARTIAL8)[r+4::8]) for r in range(4))
 bases=[((0,),0),(C3,0),((6,4,0),1),((2,12,9,1),2)]
 for a,y in bases:
  for label,c,m in [('C3',C3,4),('F8',F8,8)]:
   expected=conv(rooted(a,y),split(c));actual=rooted(join(a,c),y);assert expected==actual
   rr=residues(actual,m);assert len(set(rr))==1
   report['rooted_controls'].append({'base':a,'root':y,'block':label,'channels':m,'each_count':rr[0]})
 a=join((0,),F8)
 heads=[hp_all(close(a,0,t))[-1] for t in (0,1)]
 assert heads==[3927,7623]
 report['singleton_F8_heads']=heads
 report['scope']='Explicit new filters and new composition/jet controls; no inherited census audit or asymptotic extrapolation.'
 path=Path(__file__).with_name('receipt.json');path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
