"""Independent literal Hamilton-word census, canonicalized by directed cycles."""
from itertools import permutations,combinations
from pathlib import Path
import json
P=Path(__file__).parent

def need(ok,s):
 if not ok:raise RuntimeError(s)
def good(T,p):return all(T[a]>>b&1 for a,b in zip(p,p[1:]))
def count(T,cap=None):
 h=0;cycles=[]
 for p in permutations(range(len(T))):
  if not good(T,p):continue
  h+=1
  if cap is not None and h>cap:return h,[]
  if p[0]==0 and T[p[-1]]&1:cycles.append(p)
 return h,cycles
def reorder(T,p):return tuple(sum(1<<j for j,b in enumerate(p)if T[a]>>b&1)for a in p)
def canonical(T,cycles):
 need(cycles,'strong graph needs a cycle');forms=[]
 for p in cycles:
  for k in range(len(T)):forms.append(reorder(T,p[k:]+p[:k]))
 return min(forms)
def add(T,mask):
 n=len(T);return tuple(x|((not(mask>>i&1))<<n)for i,x in enumerate(T))+(mask,)
def sub(T,S):return reorder(T,S)
bank={3:{canonical((2,4,1),count((2,4,1))[1]):3}};levels=[];tested=0
for n in range(4,9):
 B={};cand=accepted=0
 for T in bank[n-1]:
  for mask in range(1,(1<<(n-1))-1):
   U=add(T,mask);h,cy=count(U,27);cand+=1;tested+=1
   if h>27:continue
   key=canonical(U,cy);B[key]=h;accepted+=1
 bank[n]=B;levels.append(dict(order=n,candidates=cand,accepted=accepted,classes=len(B)))
need(not bank[8],'cap27 frontier did not close at order8')
source=json.loads((P/'final-source/VALIDATION.json').read_text())
for n in range(3,9):
 expected={}
 for z in source['bank']:
  if z['order']==n:
   T=tuple(z['rows']);h,cy=count(T);need(h==z['H'],'literal source count');expected[canonical(T,cy)]=h
 need(expected==bank[n],'complete cycle-canonical bank differs')
low=[];high=[];subsets=0
for n,B in bank.items():
 for T,h in B.items():
  if h==13:low.append(dict(n=n,rows=T))
  if h==27:
   vals=[]
   for S in combinations(range(n),5):vals.append(count(sub(T,S))[0]);subsets+=1
   need(13 not in vals,'induced H13 found');high.append(dict(n=n,rows=T,five_subset_H=sorted(vals)))
need(len(low)==1 and low[0]['n']==5 and len(high)==5,'exact strong core totals')
old=(2,12,9,17,7);new=(34,12,41,49,39,2)
need(count(old)[0]==13 and count(new)[0]==25 and sub(new,range(5))==old,'positive opposite branch')
H15=(2,4,41,51,39,3);H29=(2,4,73,83,103,79,3);need(count(H15)[0]==15 and count(H29)[0]==29 and sub(H29,(0,1,2,3,4,6))==H15,'literal H15 to H29 hostile')
out=dict(H15_H29_hostile=dict(source=H15,target=H29,subset=[0,1,2,3,4,6]),status='PASS',method='Literal full Hamilton words and Hamilton-cycle canonicalization, without subset DP or degree canonicalization',extension_candidates=tested,levels=levels,unique_H13=low,H27_cores=sorted(high,key=lambda x:(x['n'],x['rows'])),literal_five_subset_counts=subsets,positive_induced_branch=[13,25],closed_empty_frontier_order=8)
print(json.dumps(out,sort_keys=True,indent=2))
