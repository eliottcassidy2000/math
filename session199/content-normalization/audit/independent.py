from pathlib import Path
from math import gcd
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def need(c,m):
 if not c:raise RuntimeError(m)
def step(s,l):
 D,E,U,d=s
 return (D,0,2*U,d) if l=='O' else ((D-E+U,d,2*U-d,d) if l=='Z' else (D-E+U,D-E,D+U,D))
seed=(1,0,2,1);bank={seed:''};queue=[seed]
for s in queue:
 for l in 'OZJ':
  v=step(s,l)
  need(v[2]>s[2],'U strictly grows')
  if v[2]<=39 and v not in bank:bank[v]=bank[s]+l;queue.append(v)
bad=[];nonprimitive=0
for s in bank:
 c=gcd(*s)
 if c>1:
  nonprimitive+=1;p=tuple(x//c for x in s)
  if p not in bank:bad.append((s,c,p))
need(len(bank)==97 and nonprimitive==17 and bad==[((39,9,39,9),3,(13,3,13,3))],'exact minimal normalization failure')
s=seed;rows=[2,0];root=0
for l in 'OOZOJ':
 s=step(s,l);n=len(rows);rows.append(0)
 for i in range(n):
  if l!='O' and i==root:rows[n]|=1<<i
  else:rows[i]|=1<<n
 if l=='J':root=n
need(s==(39,9,39,9),'physical word')
# Literal permutation count, independent of subset DP and transfer formulas.
from itertools import permutations
n=len(rows);allpaths={0:[()]}
for mask in range(1,1<<n):
 V=[i for i in range(n) if mask>>i&1]
 allpaths[mask]=[w for w in permutations(V) if all(rows[a]>>b&1 for a,b in zip(w,w[1:]))]
full=(1<<n)-1
H=len(allpaths[full]);E=sum(w[-1]==root for w in allpaths[full]);d=len(allpaths[full^(1<<root)])
U=sum(len(allpaths[mask])*sum(w and w[0]==root for w in allpaths[full^mask]) for mask in range(full+1) if not(mask>>root&1))
need((H,E,U,d)==s and sum(w[0]==root for w in allpaths[full])==1,'literal full profile')
report={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),'complete_bank_U_cap':39,'source_profiles':len(bank),'nonprimitive_profiles':nonprimitive,'only_failure':bad,'independent_literal_permutation_graph_order':n,'manual_all_order_audit':'Balanced inverse forces the absent D3,U5 source. Dyadic weak-tail uniqueness propagates nonphysicality to every padded quotient. Complete U-height minimality is valid. Quotients cannot be the omitted singleton because D-E>=d is homogeneous.','scope':'No unrestricted scalar gap or graph-order minimality claimed.'}
(Path(__file__).resolve().parent/'INDEPENDENT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
