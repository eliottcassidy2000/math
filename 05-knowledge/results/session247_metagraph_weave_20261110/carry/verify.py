from functools import lru_cache
import json
@lru_cache(None)
def g(v):
    if v==0:return 1
    w=v//2
    return g(w) if v%2 else g(w)+w+w%2

def phi(a,b,c,v):return a*g(v)-b*(v%2==0)+c*(v+1)
def solve(a,b,c,N):
    def rec(b,c,T):
        if T==0:return 0
        if T<0 or c+b>T:return None
        w=rec(0,2*c,T-c-b)
        if w is not None:return 2*w+1
        w=rec(a,2*c+a,T)
        return None if w is None else 2*w
    return rec(b,c,N+b-a-c)
checks=0
for a in range(1,8):
 for b in range(a+1):
  for c in range(a,a+5):
   for N in range(1,180):
    v=solve(a,b,c,N)
    # Phi(v)>=a-b+c(v+1) gives a finite independent brute-force bound.
    vmax=(N+b-a)//c-1
    good=any(phi(a,b,c,w)==N for w in range(max(0,vmax+1)))
    assert (v is not None)==good,(a,b,c,N,v)
    if v is not None:assert phi(a,b,c,v)==N
    checks+=1
print(json.dumps({'status':'PASS','target_cases':checks,'scope':'new arithmetic decoder only'}))
