#!/usr/bin/env python3
"""Independent literal-path controls for Session286; no producer imports."""
import itertools,json
from collections import defaultdict,deque
from pathlib import Path

def must(x,m):
 if not x: raise RuntimeError(m)

def attach(D,c):
 return (c<<1,)+tuple((r<<1)|int(not(c>>i&1))for i,r in enumerate(D))

def paths(T):
 n=len(T); H=[0]*(1<<n); s=[[0]*n for _ in H];t=[[0]*n for _ in H]; H[0]=1
 for length in range(1,n+1):
  for p in itertools.permutations(range(n),length):
   if all(T[a]>>b&1 for a,b in zip(p,p[1:])):
    S=sum(1<<v for v in p);H[S]+=1;s[S][p[0]]+=1;t[S][p[-1]]+=1
 return H,s,t

def packet(T):
 H,s,t=paths(T);F=(1<<len(T))-1;D=F^1
 sums=[]
 for f in(s,t):
  sums.extend([sum(((-1)**(S.bit_count()-1)if signed else 1)*f[S][0]*H[F^S] for S in range(1,F+1,2))for signed in(0,1)])
 U,d,V,e=sums
 return (H[F],H[D],sum(H[S]*H[D^S]for S in range(F+1)if not S&1),U,V,d,e)

def graph_check(baseline=9,P=None):
 D=(18,8,3,5,14);P=P or [packet(attach(D,c))for c in range(32)];A=P[baseline]
 h,Q=A[1:3];must((h,Q)==(13,104),'core')
 a0=A[3]-h;c0=A[4]-h
 edges=[]
 for x in range(-12,13):
  for y in range(-12,13):
   for c,p in enumerate(P):
    dc=p[4]-h-c0;da=p[3]-h-a0
    if (x+dc)%8:continue
    xx=(x+dc)//8;yy=8*y-da
    if -12<=xx<=12 and -12<=yy<=12:edges.append(((x,y),(xx,yy),c))
 # First/last endpoint rules represented by sentinel start/end states.
 g=defaultdict(list);r=defaultdict(list)
 for u,v,c in edges:
  g[u].append((v,c));r[v].append((u,c))
 def closure(seeds,G):
  seen=set(seeds);q=deque(seeds)
  while q:
   u=q.popleft()
   for v,_ in G[u]:
    if v not in seen:seen.add(v);q.append(v)
  return seen
 first=[v for u,v,c in edges if u==(0,0) and P[c][6]==A[6]]
 last=[u for u,v,c in edges if v==(0,0) and P[c][5]==A[5]]
 core=closure(first,g)&closure(last,r)
 live=[(u,v,c)for u,v,c in edges if u in core and v in core]
 return dict(packets=P,states=sorted(core),edges=live,first=[(v,c)for u,v,c in edges if u==(0,0)and v in core and P[c][6]==A[6]],last=[(u,c)for u,v,c in edges if v==(0,0)and u in core and P[c][5]==A[5]])


def spin_check():
 from fractions import Fraction
 checks=[]
 for k in range(2,17):
  vals=[]
  for s in itertools.product((0,1),repeat=k-1):vals.append(sum(s)-9*sum(a*b for a,b in zip(s,s[1:])))
  js={0}|{m-8*r for m in range(1,k//2+1) for r in range(k-2*m+1)}
  must(set(vals)==js,'exact spin support')
  if k>=3:
   mean=Fraction(sum(vals),len(vals));var=sum((Fraction(x)-mean)**2 for x in vals)/len(vals)
   must(mean==4-Fraction(7*k,4),'mean');must(var==Fraction(337*k-832,16),'variance')
  checks.append([k,len(vals),len(js),min(js),max(js)])
 holes=[t for t in range(112) if t<17*(t%8)]
 must(len(holes)==56 and max(holes)==111,'Apery holes')
 for k in range(17,151):
  js={0}|{m-8*r for m in range(1,k//2+1) for r in range(k-2*m+1)}
  lo=17-8*k;hi=k//2
  must(set(range(lo,hi+1))-js=={lo+t for t in holes},'persistent gaps')
  must(len(js)==8*k+k//2-72,'size')
 D=(18,8,3,5,14)
 aut=[p for p in itertools.permutations(range(5))if all(((D[i]>>j)&1)==((D[p[i]]>>p[j])&1)for i in range(5)for j in range(5))]
 return dict(spin_checks=checks,semigroup_orders=[17,150],holes=holes,core_automorphisms=aut)

if __name__=='__main__':
 base=graph_check();packets=base['packets'];allgraphs={c:graph_check(c,packets)for c in range(32)}
 exceptional=[c for c,g in allgraphs.items()if len(g['edges'])>1]
 must(exceptional==[9,13,18,22],'all32 exceptional baselines')
 expected={9:[3,9,18,24],13:[7,13,22,28],18:[3,9,18,24],22:[7,13,22,28]}
 for c,g in allgraphs.items():
  must(sorted(e[2]for e in g['edges'])==expected.get(c,[c]),'trimmed alphabet')
 out={'literal_packet_and_carry':base,'all_constant_baselines':{c:{k:v for k,v in g.items()if k!='packets'}for c,g in allgraphs.items()},'spin_semigroup':spin_check()}
 print(json.dumps(out,indent=2))
