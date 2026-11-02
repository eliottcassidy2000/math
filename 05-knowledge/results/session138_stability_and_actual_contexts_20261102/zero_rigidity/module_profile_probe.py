import sys,itertools,math,json
from reference_paths import tournament,c3_blowup,extend,H

def induce(a,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def profile(a):
 n=len(a);F=[0]*(n+1)
 for p in itertools.permutations(range(n)):
  d=sum(not(a[p[i]]>>p[i+1]&1) for i in range(n-1))
  for r in range(d+1,n+1):F[r]+=math.comb(n-1-d,r-1-d)
 return F[1:]
T=c3_blowup(tournament(2,1));B=induce(T,list(range(1,6)));f=profile(T)
for S in [14,22,26,29]:
 U=extend(B,S);g=profile(U)
 print(json.dumps({'cut':S,'profile':g,'difference':[b-a for a,b in zip(f,g)],'masks':U}))
print('original_profile',f)
