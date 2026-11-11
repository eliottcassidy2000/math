#!/usr/bin/env python3
"""Exact controls of NEW operators and endpoint consumers, not an inherited audit."""
from collections import Counter,defaultdict
from itertools import combinations,permutations,product
from math import comb
from random import Random
import json
from explore import cycles,wrap as affine_wrap
from interface import wrapper_core,skeletons,compile as compile_interface
from coalescence import iterate,two_colour

def check(ok,label):
 if not ok:raise RuntimeError(label)

def apex(w):
 n=len(w);a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
 for i,b in enumerate(w):
  if b=='1':a[n]|=1<<i
  else:a[i]|=1<<n
 return tuple(a)
def join(a,b):
 n=len(a);right=((1<<len(b))-1)<<n
 return tuple(x|right for x in a)+tuple(x<<n for x in b)
def W(a):
 n=len(a);return tuple(x|(1<<(n+1)) for x in a)+((1<<n)-1,1<<n)
def converse(a):
 n=len(a);return tuple(((1<<n)-1)^x^(1<<i) for i,x in enumerate(a))
def hp_all(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*((1<<n)-1)
 incoming=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)]
 for s in range(1,1<<n):
  bits=s
  while bits:
   bit=bits&-bits;v=bit.bit_length()-1;bits-=bit;prev=s^bit
   if not prev:d[s][v]=1
   else:
    cand=prev&incoming[v];total=0
    while cand:
     b=cand&-cand;cand-=b;total+=d[prev][b.bit_length()-1]
    d[s][v]=total
  h[s]=sum(d[s])
 return h

def K(a):
 h=hp_all(a);n=len(a);full=(1<<n)-1;p=[0]*(n+1)
 for s in range(full+1):p[s.bit_count()]+=h[s]*h[full^s]
 return p

def signature(a):
 n=len(a);cs=[]
 for k in range(3,n+1,2):
  for ss in combinations(range(n),k):
   count=0
   for tail in permutations(ss[1:]):
    p=(ss[0],)+tail
    count+=all(a[p[i]]>>p[(i+1)%k]&1 for i in range(k))
   if count:cs.append((sum(1<<v for v in ss),k,count))
 out=Counter()
 def rec(index,used,parts,coef):
  out[tuple(sorted(parts+(1,)*(n-used.bit_count())))]+=coef
  for j in range(index,len(cs)):
   s,k,c=cs[j]
   if not s&used:rec(j+1,used|s,parts+(k,),coef*2*c)
 rec(0,0,(),1);return out

def endpoint_gate(a,y,gate):
 """Direct full-path DP records the gate position, without a wrap decomposition."""
 n=len(a);full=(1<<n)-1;dp=[defaultdict(int) for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v,1 if v==gate else 0]=1
 for s in range(1,full+1):
  for (last,pos),c in list(dp[s].items()):
   opts=a[last]&~s&full
   while opts:
    b=opts&-opts;opts-=b;v=b.bit_length()-1
    dp[s|b][v,s.bit_count()+1 if v==gate else pos]+=c
 p=[0]*(n+1)
 for (last,pos),c in dp[full].items():
  if last==y:p[n-pos]+=c
 return p,sum(dp[full].values())

def main():
 gates=0
 for n in range(5):
  for bits in product('01',repeat=n):
   w=''.join(bits);a=apex(w);cs=cycles(w)
   F=Counter({(1,)*(n+1):1})
   for j,c in cs.items():
    if c:F[tuple(sorted((j+1,)+(1,)*(n-j)))]+=2*c
   for r in range(3):
    actual=K(a);op=two_colour(iterate(F,r));ic=compile_interface(n,cs,skeletons(*wrapper_core(r)),2*r)
    check(actual==op==ic,('word',w,'wraps',r));gates+=1
    a=W(a)
   check(affine_wrap(n,cs)==two_colour(iterate(F,1)),('affine',w));gates+=1
 rng=Random(264092)
 for n in range(3,8):
  for sample in range(4):
   a=[0]*n
   for i in range(n):
    for j in range(i):
     if rng.randrange(2):a[i]|=1<<j
     else:a[j]|=1<<i
   F=signature(a)
   for r in range(3):
    check(K(a)==two_colour(iterate(F,r)),('random',n,sample,r));gates+=1;a=W(a)
 F8=(234,236,161,228,111,192,132,16)
 source=join(F8,(0,));y=8;wrapped=W(source);p,H=endpoint_gate(wrapped,y,9)
 residue=[sum(p[r::8]) for r in range(8)]
 check(H==7623 and sum(p)==3696 and residue==[462]*8,'11-vertex consumer');gates+=1
 # Starting endpoint and the other gate are the converse geometry.
 source2=join((0,),F8);wrapped2=W(source2)
 p2,H2=endpoint_gate(converse(wrapped2),0,10)
 check(H2==7623 and sum(p2)==3696 and [sum(p2[r::8]) for r in range(8)]==[462]*8,'start consumer');gates+=1
 # A wrap is not an ordinal multiplier: it can destroy filter zeros.
 KF=K(W(F8));check(KF==[3927,7548,7545,7752,7584,7512,7584,7752,7545,7548,3927],'wrapped F8');gates+=1
 check(len(set(sum(KF[r::8]) for r in range(8)))>1,'nonpreservation hostile');gates+=1
 dipoles=[]
 for a,expected in [((50,44,57,17,34,8),(52,41)),((42,52,57,18,33,8),(50,43))]:
  e,HH=endpoint_gate(W(a),0,6);st,HH2=endpoint_gate(converse(W(a)),0,7)
  check(HH==HH2==349 and (sum(st),sum(e))==expected,'cross-lane endpoint dipole');gates+=1
  dipoles.append({'module_rows':a,'wrapped_H':HH,'start':sum(st),'end':sum(e),'wrapped_K':K(W(a))})
 print(json.dumps({'status':'PASS','new_operator_gates':gates,'word_sources':31,'arbitrary_sources':20,'consumer_rows':wrapped,'consumer_H':H,'old_y_endpoint_count':sum(p),'suffix_length_counts':p,'suffix_mod8':residue,'wrapped_F8_K':KF,'cross_lane_dipoles':dipoles},indent=2))
if __name__=='__main__':main()
