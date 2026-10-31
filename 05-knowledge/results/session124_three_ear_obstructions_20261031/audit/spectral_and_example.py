#!/usr/bin/env python3
"""Independent common local flag replay and finite warning-example check."""
import itertools,json
import audit as a
P=[(0,1,-1),(1,-1,0),(0,0,1)]
Q=[(0,1,-1),(0,1,0),(1,0,0)]
S=[(0,1),(1,0)]
def coordinates(kind,v):
 if kind==0:return [v[0]+v[1],v[0],sum(v)]
 if kind in [1,2]:return [-v[2],v[1]+v[2],v[0]]
 return [v[1],v[0]]
diagonal=[[0,0,1],[0,1,0],[0,1,0],[1,0]]
for bits in itertools.product([0,1],repeat=3):
 for kind,basis in enumerate([P,Q,Q,S]):
  target=1 if kind==3 or not bits[kind] else 2
  for j,vec in enumerate(basis):
   weight=sum(vec) if kind==0 else sum(vec[1:])+bits[kind-1]*vec[0]
   out=[0]*len(vec);out[target]=weight;transformed=coordinates(kind,out)
   a.need(transformed[j]==diagonal[kind][j],'constant flag diagonal')
   for i in range(j+1,len(vec)):a.need(transformed[i]==0,'common triangular local flag')
multiplicities=[0]*5
for terms in itertools.product(*diagonal):multiplicities[sum(terms)]+=1
a.need(multiplicities==[8,20,18,7,1],'per-order spectral multiplicities')
counts=[]
for k in range(4):
 word=[3,5,7,5,1]+[7,0]*k;H=a.brute_h(1,word);v=a.INITIAL
 for letter in word:v=a.actvec(v,letter)
 a.need(H==sum(x*y for x,y in zip(a.output(1),v)),'example subsetDP versus transfer');counts.append(H)
a.need(counts==[215,3443,55091,878159],'published warning example values')
a.need(16*counts[0]+3==counts[1] and 16*counts[1]+3==counts[2] and 16*counts[2]+3-counts[3]==3300,'third-iterate deficit')
print(json.dumps({'status':'passed','checks':a.checks,'local_flag_configurations':32,'per_order_multiplicities':multiplicities,'full_multiplicities':[6*x for x in multiplicities],'example_counts':counts,'third_iterate_deficit':3300},indent=2))
