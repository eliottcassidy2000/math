from bridge_blocks import *
def packing(out,x):
 n=len(out);cyc={}
 for k in range(3,n+1,2):
  for S in itertools.combinations(range(n),k):
   mask=sum(1<<i for i in S);c=0
   for t in itertools.permutations(S[1:]):
    p=(S[0],)+t
    c+=all(out[a]>>b&1 for a,b in zip(p,p[1:]+p[:1]))
   if c:cyc[mask]=c
 @lru_cache(None)
 def dp(S):
  if not S:return 1
  bit=S&-S
  return dp(S^bit)+sum(x*c*dp(S^M) for M,c in cyc.items() if M&bit and M&S==M)
 return dp((1<<n)-1)
if "--fail" in sys.argv: raise RuntimeError("deliberate sentinel")
checked=0
for n in range(1,8):
 for c in compositions(n):
  for x in (0,1,2,3,5):
   if count(c,x)!=packing(tournament(c),x): raise RuntimeError(("packing mismatch",c,x))
   checked+=1
print(json.dumps({'independent_odd_cycle_packing_checks':checked,'all_compositions_max_n':7,'bases':[0,1,2,3,5]}))
