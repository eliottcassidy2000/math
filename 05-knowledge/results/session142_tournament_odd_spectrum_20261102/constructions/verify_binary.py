#!/usr/bin/env python3
"""Portable exact controls for the specified S0/R1 binary language."""
import json
from itertools import product

def need(ok,msg):
 if not ok: raise RuntimeError(msg)
def S(v):
 C,L,g,t=v;return (2*C+L+g,L,0,g+t)
def R(v):
 C,L,g,t=v;return (C+L+g,0,t,2*C+g)
def N(v):
 C,L,g,t=v;return (C,C+2*L+g,g+t,0)
def Z(v):return (0,0,0,F(v))
def F(v):C,L,g,t=v;return 2*C+g+t
def dead(v):return I(v)>=F(v) or v[:3]==(0,0,0)
def I(v):C,L,g,t=v;return 2*L+g-t
def eligible(v):return I(v) in (-1,1)
def paths(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    z=s^(1<<v)
    dp[s][v]=1 if not z else sum(dp[z][u] for u in range(n) if z>>u&1 and a[u]>>v&1)
  H[s]=sum(dp[s])
 return H,dp

def actual(a,cut):
 n=len(a);H,end=paths(a)
 rev=tuple(sum(1<<j for j in range(n) if j!=i and not a[i]>>j&1) for i in range(n))
 _,start=paths(rev);full=(1<<n)-1
 q=sum(H[full^s]*sum(start[s][v] for v in range(n) if cut>>v&1) for s in range(1,1<<n))
 g=sum(end[full][v] for v in range(n) if cut>>v&1)
 b=list(a)+[cut]
 for v in range(n):
  if not cut>>v&1:b[v]|=1<<n
 f=paths(b)[0][-1];h=H[-1]
 need((f-h)%2==0 and (q-f+h-g)%2==0,'nonintegral actual source')
 return ((f-h)//2,(q-f+h-g)//2,g,h-g),(h,f,q,g),tuple(b)

def extend(a,cut,op):
 # Full extension has old parent A, old root x, fresh z, A->z->x.
 n=len(a);b=list(a)+[cut,1<<n]
 for v in range(n):
  b[v]|=1<<(n+1)
  if not cut>>v&1:b[v]|=1<<n
 root=n if op=='S' else n+1
 keep=[v for v in range(n+2) if v!=root]
 out=tuple(sum(1<<j for j,w in enumerate(keep) if b[v]>>w&1) for v in keep)
 newcut=sum(1<<j for j,w in enumerate(keep) if b[root]>>w&1)
 return out,newcut

def main():
 # Solve all sign equations, not an unbounded-coordinate sample.
 rs=set()
 for d,e,f in product((-1,1),repeat=3):
  t=(e-f)//2
  if t<0:continue
  rhs=t-e
  if rhs<0:continue
  for C in range(rhs//2+1):
   g=rhs-2*C;num=d-g+t
   if num<0 or num%2:continue
   v=(C,num//2,g,t)
   if g+t>=1:
    need((I(v),I(R(v)),I(S(R(v))))==(d,e,f),'RS equation failure');rs.add(v)
 rr=set()
 for g,t in [(0,1),(1,0)]:
  for e in (-1,1):
   num=t-g-e
   if num<0 or num%2:continue
   v=(num//2,0,g,t)
   need(eligible(v) and I(R(v))==e and I(R(R(v)))==-1,'RR equation failure');rr.add(v)
 need(rs=={(0,0,0,1),(0,0,1,0),(0,1,0,1)},'RS list')
 need(rr=={(0,0,0,1),(0,0,1,0),(1,0,0,1)},'RR list')
 # Finite formal-state control, independently track histories and operation words.
 seeds=histories=0
 for v in product(range(7),repeat=4):
  if v[2]+v[3]<1:continue
  need(F(S(v))==F(R(v))==2*F(v)+I(v),'count identity')
  need(I(S(v))==I(v)-2*v[2] and I(R(v))==2*v[1]-2*v[0]-I(v),'I identities')
  need(F(N(v))==F(Z(v))==F(v),'neutral output identity')
  need(I(N(v))==F(v)+4*v[1]+2*v[2],'neutral I identity')
  if F(v)>=5:
   need(dead(N(v)) and dead(Z(v)),'neutral absorption')
   if dead(v):
    need(not eligible(v),'dead source eligible')
    for op in (S,R):
     if F(op(v))==F(v):need(dead(op(v)),'neutral deadset escape')
  if F(v)<5 or not eligible(v):continue
  seeds+=1;level=[(v,'')];pure=v
  for depth in range(1,9):
   nxt=[]
   for u,w in level:
    if eligible(u):
     for name,op in [('S',S),('R',R)]:
      z=op(u);need(abs(F(z)-2*F(u))==1,'binary mismatch');nxt.append((z,w+name));histories+=1
   pure=S(pure)
   if depth>=4:
    need(all(w[:-1]=='S'*(depth-1) for u,w in nxt),'word obstruction failed')
    need(all(F(u)==F(pure) for u,w in nxt),'scalar thinning failed')
   level=nxt
 # Literal graph checks include outputs of forbidden third steps, not only accepted histories.
 witnesses=[];A=(6,4,0)
 for cut,expected in [(1,(2,0,0,1)),(3,(2,1,0,1))]:
  v,native,full=actual(A,cut);need(v==expected,'actual seed mismatch')
  stages=[]
  for word in ('S','R','SS','SR','RS','RR','RSS','RRS','SSSS'):
   a=A;s=cut;w=v
   for letter in word:
    a,s=extend(a,s,letter);w=(S if letter=='S' else R)(w)
    got,raw,_=actual(a,s);need(got==w,'actual transfer mismatch')
   stages.append({'word':word,'positive_state':w,'H':F(w),'I':I(w),'eligible_for_next':eligible(w)})
  nparent=tuple(x|(1<<len(A)) for x in A)+(0,)
  ncut=cut|(1<<len(A));gotn,rawn,_=actual(nparent,ncut)
  gotz,rawz,_=actual(full,0)
  need(gotn==N(v) and gotz==Z(v),'actual neutral gadget mismatch')
  witnesses.append({'parent_masks':A,'cut':cut,'positive_state':v,'native_state':native,'full_masks':full,'continuations':stages,
   'neutral_S1':{'parent_masks':nparent,'cut':ncut,'positive_state':gotn,'native_state':rawn},
   'neutral_reset':{'parent_masks':full,'cut':0,'positive_state':gotz,'native_state':rawz}})
 # Hostiles to scalar branching and illegal long continuations.
 need(F(S(expected))==F(R(expected)),'mark choice wrongly treated as a digit')
 v=(2,0,0,1)
 need(F(R(v))==9 and not eligible(R(v)),'R hostile changed')
 # Finite seed bank is essential: arbitrary fresh seeds can trivially contain every target.
 receipt={'status':'PASS','scope':'Two-operation S0/R1 binary-growth language only; unrestricted spectrum OPEN',
  'symbolic_RS_cases':[list(v) for v in sorted(rs)],'symbolic_RR_cases':[list(v) for v in sorted(rr)],
  'formal_coordinate_box':'0..6 each, h>=1','formal_eligible_seeds_F_ge5':seeds,
  'formal_checked_history_edges_through_depth8':histories,'actual_witnesses':witnesses,
  'claims':['At F>=5 no R can start three binary steps','At depth>=4 first k-1 operations all S','Finite seed banks yield O(log X) outputs','Neutral S1/reset interleavings add no outputs']}
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
