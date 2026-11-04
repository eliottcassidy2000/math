#!/usr/bin/env python3
"""Independent finite controls of the all-order one-bit profile reduction."""
from pathlib import Path
from itertools import combinations,product
import hashlib,json
ROOT=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,c):
 T=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j) if c>>k&1 else(j,i);T[a]|=1<<b
 return T
def kernel(Q,a,c,caps=None):
 m=len(Q);caps=[5]*m if caps is None else caps;R=list(Q)+[1<<(m+1),0]
 for i in range(m):
  for mask,p in((a,m),(c,m+1)):
   if mask>>i&1:R[i]|=1<<p
   else:R[p]|=1<<i
 cap=sum(caps)+2;levels=[{} for _ in range(cap+1)];out={}
 for v in range(m+2):
  e=[0]*m;ports=0
  if v<m:e[v]=1
  else:ports=1<<(v-m)
  levels[1][tuple(e)+(ports,v)]=1
 for length in range(1,cap+1):
  for state,z in levels[length].items():
   e=state[:m];ports,last=state[m:]
   if ports==3:out[e]=out.get(e,0)+z
   if length==cap:continue
   for v in range(m+2):
    if not R[last]>>v&1:continue
    f=list(e);pp=ports
    if v<m:
     if f[v]==caps[v]:continue
     f[v]+=1
    else:
     bit=1<<(v-m)
     if pp&bit:continue
     pp|=bit
    key=tuple(f)+(pp,v);levels[length+1][key]=levels[length+1].get(key,0)+z
 return out

def literal(Q,sizes,a,c):
 off=[sum(sizes[:i]) for i in range(len(sizes))];n=sum(sizes);R=[0]*(n+2);R[n]=1<<(n+1)
 for i,s in enumerate(sizes):
  for u in range(s):
   v=off[i]+u;R[v]=sum(1<<(off[i]+w) for w in range(u+1,s))
   for j,t in enumerate(sizes):
    if Q[i]>>j&1:R[v]|=((1<<t)-1)<<off[j]
   for mask,p in((a,n),(c,n+1)):
    if mask>>i&1:R[v]|=1<<p
    else:R[p]|=1<<v
 N=1<<(n+2);d=[[0]*(n+2) for _ in range(N)]
 for i in range(n+2):d[1<<i][i]=1
 for M in range(1,N):
  for v in range(n+2):
   if not M>>v&1 or M==(1<<v):continue
   d[M][v]=sum(d[M^(1<<v)][u] for u in range(n+2) if R[u]>>v&1)
 return sum(d[-1])

def main():
 types=[(0,1,b,c,8,d) for b,c,d in product((6,14),(4,12),(0,8))];base=(0,1,14,4,8,0);records=[];formal=physical=0
 for m in range(1,4):
  for code in range(1<<(m*(m-1)//2)):
   Q=rows(m,code);qs=[]
   for a,c in product(range(1<<m),repeat=2):
    K=kernel(Q,a,c);coeffs=[(e,v) for e,v in K.items() if min(e)>0];baseline=sum(v*__import__('math').prod(base[k] for k in e) for e,v in coeffs);D=[];higher=[];single=literal(Q,[1]*m,a,c);physical+=1;pair_values=[]
    for i in range(m):
     e=[1]*m;e[i]=2;D.append(K.get(tuple(e),0));sizes=[1]*m;sizes[i]=2;pair=literal(Q,sizes,a,c);physical+=1;pair_values.append(pair);need(pair-single==2*D[-1],'exact small physical derivative')
     hi=[]
     for k in(3,5):e[i]=k;value=K.get(tuple(e),0);need(value%2==0,'first-order higher-run coefficient parity');hi.append(value)
     higher.append(hi)
    if m<=2:need(literal(Q,[4]*m,a,c)%16==baseline%16,'actual TT4 baseline');physical+=int(m<=2)
    for ids in product(range(8),repeat=m):
     P=[types[i] for i in ids];actual=sum(v*__import__('math').prod(P[i][e[i]] for i in range(m)) for e,v in coeffs)%16;pred=(baseline+8*sum(D[i] for i in range(m) if P[i][2]==6))%16;need(actual==pred,'every complete formal response matches one-bit formula');formal+=1
    qs.append(dict(a=a,c=c,complete_integer_run_coefficients=[[*e,v] for e,v in sorted(K.items())],TT4_baseline_functional=baseline,derivative_coefficients=D,higher_run_coefficients=higher,actual_singleton_response=single,actual_pair_responses=pair_values))
   records.append(dict(Q=Q,full_block_query_controls=qs))
 # Nonzero five-run coefficient, beyond the small-quotient pigeonhole boundary.
 Q7=[14,124,120,112,97,65,1];higher7=[]
 for k in(3,5):
  caps=[k]+[1]*6;K=kernel(Q7,1,126,caps);v=K.get(tuple(caps),0);need(v>0 and v%2==0,'nonzero higher-order return-swap coefficient');higher7.append([k,v])
 witness=[0,1,4,0,2,5,0,3,6,0,7,8,0];swapped=[0,2,5,0,1,4,0,3,6,0,7,8,0];R=Q7+[256,0]
 for i in range(7):
  for mask,p in((1,7),(126,8)):
   if mask>>i&1:R[i]|=1<<p
   else:R[p]|=1<<i
 need(witness!=swapped and all(R[a]>>b&1 for w in(witness,swapped) for a,b in zip(w,w[1:])),'explicit five-run free swap')
 paths=[ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'even_cyclic_mod16'/'PROOF.md',ROOT.parent/'parity_profile_transfer'/'PROOF.md']
 return dict(status='PASS',higher_order_return_control=dict(Q=Q7,query=[1,126],coefficients=higher7,witness=witness,swapped=swapped),quotients=len(records),formal_response_checks=formal,literal_physical_counts=physical,records=records,formal_profile_domain=types,source_sha256={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},scope='All-order one-bit normalized coarse-response identity with complete controls on every labelled quotient through3. No universal noncollision claim or ordered-cover analogue.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
