"""All-word bounded-state certificate for later return targets in S156 bank65."""
from pathlib import Path
from collections import deque
import hashlib,json

def need(x,m):
 if not x:raise AssertionError(m)
p=Path(__file__).with_name('bank65.json');raw=p.read_bytes()
need(hashlib.sha256(raw).hexdigest()=='26fad58d82f2428be3bba1e7b72e93e6494983eb8e33d1a3dba1f10c2c62ca20','pinned bank checksum')
D=json.loads(raw);A=[v['matrix'] for v in D.values()];need(len(A)==65,'bank size')
need(all(len(a)==3 and all(len(r)==3 and all(isinstance(z,int) and z>=0 for z in r) for r in a) and a[0][0]>=1 for a in A),'nonnegative h-monotone matrices')
def mul(v,a,mod=None):
 w=tuple(sum(v[i]*a[i][j] for i in range(3)) for j in range(3))
 return w if mod is None else tuple(z%mod for z in w)
T=((855,369,675),(857,368,675),(859,367,675),(873,360,675),(875,359,675),(877,358,675))
cap=877;S={(1,0,0)};q=deque(S)
while q:
 v=q.popleft()
 for a in A:
  w=mul(v,a)
  if w[0]<=cap and w not in S:S.add(w);q.append(w)
need(len(S)==7372,'exact reachable-state count')
need(all(t not in S for t in T),'six targets excluded')
need(all(mul(v,a)[0]>cap or mul(v,a) in S for v in S for a in A),'closed at height 877')
print('exact h<=877 reachable states:',len(S),'six targets absent')
# Show finite modular reachability rather than misreading the exact gap as a
# small-modulus obstruction. This does not prove absence of any modulus gate.
for mod in range(2,14):
 S={(1,0,0)};q=deque(S)
 while q:
  v=q.popleft()
  for a in A:
   w=mul(v,a,mod)
   if w not in S:S.add(w);q.append(w)
 need(all(tuple(x%mod for x in t) in S for t in T),'target residues mod '+str(mod))
 print('modulus',mod,'orbit',len(S),'all six target residues reached')
print('PASS')
