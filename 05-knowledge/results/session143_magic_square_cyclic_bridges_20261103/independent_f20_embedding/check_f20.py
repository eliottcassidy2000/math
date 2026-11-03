#!/usr/bin/env python3
"""Independent F20/class-selection audit. B facts remain cited inputs.
Enumerates five-point permutations and exact Gaussian-integer characters.
"""
from collections import deque
from itertools import product
import argparse,json
p=argparse.ArgumentParser();p.add_argument('--mutate',choices=['centralizer','character','permutation-copies']);args=p.parse_args()
checks=0
def need(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)
I=tuple(range(5));g=tuple((i+1)%5 for i in range(5));h=tuple(2*i%5 for i in range(5))
def mul(s,t):return tuple(s[t[i]] for i in range(5))
def inv(s):return tuple(s.index(i) for i in range(5))
def power(s,k):
 v=I
 for _ in range(k):v=mul(s,v)
 return v
G={I};queue=deque([I])
while queue:
 s=queue.popleft()
 for t in (g,h,inv(g),inv(h)):
  u=mul(t,s)
  if u not in G:G.add(u);queue.append(u)
need(len(G)==20,'faithful affine group order twenty')
need(power(g,5)==I and power(h,4)==I and power(h,2)!=I,'generator orders')
need(mul(mul(h,g),inv(h))==power(g,2),'semidirect exponent two')
C={power(g,i) for i in range(5)}
need(all(mul(mul(s,c),inv(s)) in C for s,c in product(G,C)),'normal C5')
need(sum(mul(s,g)==mul(g,s) for s in G)==5,'centralizer in F20')
need(len({mul(mul(s,g),inv(s)) for s in G})==4,'all four nonidentity translations conjugate')
classes=[];unseen=set(G)
while unseen:
 s=min(unseen);cl={mul(mul(t,s),inv(t)) for t in G};classes.append(cl);unseen-=cl
need(sorted(map(len,classes))==[1,4,5,5,5],'five conjugacy classes')
# The exponent b of the multiplier 2^b gives the quotient C4.
bexp={s:next(b for b in range(4) if (s[1]-s[0])%5==pow(2,b,5)) for s in G}
for s,t in product(G,repeat=2):need(bexp[mul(s,t)]==(bexp[s]+bexp[t])%4,'quotient exponent homomorphism')
i_powers=[(1,0),(0,1),(-1,0),(0,-1)]
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def times(x,y):return (x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def conjugate(x):return (x[0],-x[1])
chars=[{s:i_powers[(j*bexp[s])%4] for s in G} for j in range(4)]
W={s:(sum(s[i]==i for i in range(5))-1,0) for s in G}
if args.mutate=='character':W[g]=(0,0)
chars.append(W)
gram=[]
for a,b in product(range(5),repeat=2):
 total=(0,0)
 for s in G:total=add(total,times(chars[a][s],conjugate(chars[b][s])))
 need(total==(20*int(a==b),0),'exact irreducible-character orthogonality')
 gram.append([a,b,total[0]//20])
degrees=[v[I][0] for v in chars]
need(degrees==[1,1,1,1,4] and sum(d*d for d in degrees)==20,'complete irreducible degree list')
need(all(W[c]==((-1,0) if c!=I else (4,0)) for c in C),'W restricts to four nontrivial C5 characters')
need([chars[j][h] for j in range(4)]==[(1,0),(0,1),(-1,0),(0,-1)],'rationality condition a1 equals a3')
# Cited B/HS inputs; only arithmetic consequences are verified here.
HS=44352000;M_order=1774080000;centralA=443520000;centralB=6000000
if args.mutate=='centralizer':centralB=5*HS+1
lower=5*HS
need(lower>centralB,'C5 times HS excludes class5B')
need(M_order==40*HS,'listed subgroup order and index-two base')
need(4*centralA==M_order,'normalizer bound attained by cited normalizing subgroup')
d=4371;c=21
need((d-c)%5==0 and (d+4*c)%5==0,'integral C5 multiplicities')
nW=(d-c)//5;linear_total=(d+4*c)//5
need((nW,linear_total)==(870,891),'restriction isotypic multiplicities')
need(4*nW+linear_total==d and linear_total-nW==c,'degree and translation character')
# Two rational F20 modules with identical known C5 data but different trivial multiplicity.
formal=[]
for multiplicities in [(891,0,0,0),(0,0,891,0)]:
 chi={s:tuple(nW*x for x in W[s]) for s in G}
 for j,n in enumerate(multiplicities):
  for s in G:chi[s]=add(chi[s],tuple(n*x for x in chars[j][s]))
 need(chi[I]==(4371,0) and all(chi[x]==(21,0) for x in C if x!=I),'formal compatible C5 data')
 need(all(v[1]==0 for v in chi.values()),'formal modules rational-valued')
 formal.append({'linear_multiplicities':multiplicities,'trivial_F20_multiplicity':multiplicities[0],'h_trace':chi[h][0]})
need(formal[0]['h_trace']!=formal[1]['h_trace'],'order-five data do not determine linear split')
if args.mutate=='permutation-copies':need(formal[1]['trivial_F20_multiplicity']>=870,'whole five-dimensional permutation copies require trivial constituents')
need(formal[1]['trivial_F20_multiplicity']==0,'no inference of870 whole permutation copies')
print(json.dumps({'status':'PASS','checks':checks,'group_order':20,'conjugacy_class_sizes':sorted(map(len,classes)),
 'irreducible_degrees':degrees,'character_gram_entries':gram,'centralizer_lower_bound':lower,'selected_B_class':'5A',
 'normalizer_order':M_order,'W_multiplicity':nW,'total_linear_multiplicity':linear_total,
 'rationality_relation':'a1=a3','formal_noninference_witnesses':formal,
 'scope':'Abstract embedding and restriction conditional on cited subgroup/class/character inputs. No Baby Monster generators, table reconstruction, intertwiner or actual individual linear multiplicities computed.'},indent=2,sort_keys=True))
