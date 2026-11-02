#!/usr/bin/env python3
"""Independent algebraic controls for S0/R1 binary-history obstruction."""
from itertools import product, permutations
from functools import lru_cache
import json

def need(c,m):
 if not c:raise RuntimeError(m)
def S(v):
 C,L,g,t=v
 return (2*C+L+g,L,0,g+t)
def R(v):
 C,L,g,t=v
 return (C+L+g,0,t,2*C+g)
def N(v):
 C,L,g,t=v
 return (C,C+2*L+g,g+t,0)
def Z(v):return (0,0,0,F(v))
def dead(v):return I(v)>=F(v) or v==(0,0,0,F(v))
def F(v):return 2*v[0]+v[2]+v[3]
def I(v):return 2*v[1]+v[2]-v[3]
def binary(v,ops):
 for op in ops:
  if I(v) not in (-1,1):return False
  v=op(v)
 return True

def path_words(a,vertices):
 return [p for p in permutations(vertices) if all(a[u][v] for u,v in zip(p,p[1:]))]
def native(a,cut):
 n=len(a);vertices=set(range(n));words=path_words(a,range(n));h=len(words)
 g=sum(p[-1] in cut for p in words)
 q=0
 for mask in range(1,1<<n):
  subset={v for v in range(n) if mask>>v&1}
  q+=len(path_words(a,vertices-subset))*sum(p[0] in cut for p in path_words(a,subset))
 full=[row+[int(i not in cut)] for i,row in enumerate(a)]+[[int(i in cut) for i in range(n)]+[0]]
 f=len(path_words(full,range(n+1)))
 need((f-h)%2==0 and (q-f+h-g)%2==0,'actual source parity')
 v=((f-h)//2,(q-f+h-g)//2,g,h-g)
 need(min(v)>=0,'actual source positivity')
 return v

def actual_checks():
 checks=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   a=[[0]*n for _ in range(n)];k=0
   for i in range(n):
    for j in range(i+1,n):a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j];k+=1
   for cutcode in range(1<<n):
    cut={i for i in range(n) if cutcode>>i&1};v=native(a,cut)
    b=[row+[int(i not in cut),1] for i,row in enumerate(a)]
    b += [[int(i in cut) for i in range(n)]+[0,0],[0]*n+[1,0]]
    for root,op in ((n,S),(n+1,R)):
     keep=[i for i in range(n+2) if i!=root]
     parent=[[b[i][j] for j in keep] for i in keep]
     newcut={j for j,i in enumerate(keep) if b[root][i]}
     need(native(parent,newcut)==op(v),'literal actual graph transfer')
     checks+=1
 return {'parent_orders':'1..4','S_R_graph_transfers':checks,'counter':'literal permutations and subset covers'}

def main():
 # Linear identities hold on the basis, therefore identically on all states.
 for v in ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)):
  need(F(S(v))==F(R(v))==2*F(v)+I(v),'response identity')
  need(I(S(R(v)))==-F(v),'RS identity')
  need(I(R(R(v)))==-(2*v[1]+v[2]+v[3]),'RR identity')
  need(F(N(v))==F(v) and I(N(v))==F(v)+4*v[1]+2*v[2],'neutral N identity')
  need(F(Z(v))==F(v) and I(Z(v))==-F(v),'reset identity')
 checked=0;rs=set();rr=set()
 for v in product(range(11),repeat=4):
  if v[2]+v[3]<1:continue
  checked+=1
  if F(v)>=5:
   need(dead(N(v)) and dead(Z(v)),'neutral entry to dead set')
   if dead(v):
    need(I(v) not in (-1,1),'dead set admits binary growth')
    for op in (S,R,N,Z):
     w=op(v)
     if F(w)==F(v):need(dead(w),'dead set escaped by neutral move')
  for final in (S,R):
   if binary(v,(R,S,final)):
    rs.add(F(v));need(F(v)==1,'RS triple bound')
   if binary(v,(R,R,final)):
    rr.add(F(v));need(F(v) in (1,3),'RR triple bound')
 need(rs=={1} and rr=={1,3},'sharp exception coverage')
 need(binary((0,0,1,0),(R,S,S)),'F1 exception positive control')
 need(binary((1,0,0,1),(R,R,S)),'F3 exception positive control')
 # Every binary history from F>=5 has all R moves in its last two places.
 seeds=[v for v in product(range(4),repeat=4) if v[2]+v[3]>=1 and F(v)>=5]
 valid_histories=0
 for seed in seeds:
  frontier=[(seed,'')]
  for length in range(1,9):
   nxt=[]
   for v,word in frontier:
    if I(v) in (-1,1):
     for op,name in ((S,'S'),(R,'R')):
      w=op(v);text=word+name
      need('R' not in text[:-2],'early reset escaped')
      if length>=4:need(text[:-1]=='S'*(length-1),'penultimate reset escaped')
      need(F(w)>=2**length*(F(seed)-1)+1,'binary growth floor')
      nxt.append((w,text));valid_histories+=1
   need(len(nxt)<=4,'more than four histories at one length')
   if length>=4:need(len({F(w) for w,_ in nxt})<=1,'more than one deep scalar output')
   frontier=nxt
 print(json.dumps({'status':'PASS','linear_identities_checked_on_basis':4,'bounded_nonnegative_states':checked,
  'RS_allowed_initial_F':sorted(rs),'RR_allowed_initial_F':sorted(rr),'history_seeds':len(seeds),
  'dead_set_neutral_extension':True,'actual_source_checks':actual_checks(), 'histories_checked':valid_histories,'max_depth':8,'scope':'Identities plus elementary proof give all-state no-go; census is only a control. O(log X) image size does not imply failure to meet every residue modulo a fixed modulus.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
