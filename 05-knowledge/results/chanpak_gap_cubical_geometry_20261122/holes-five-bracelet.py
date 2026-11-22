from itertools import product,permutations,combinations
def comp(T):
 N=len(T)
 return tuple(x^((1<<N)-1)^(1<<i) for i,x in enumerate(T))
def rows(bits):
 T=[0]*5
 for i in range(5):T[i]|=1<<((i+1)%5)
 for i,b in enumerate(bits):
  j=(i+2)%5
  if b=="1":T[i]|=1<<j
  else:T[j]|=1<<i
 return tuple(T)
def right(T,i,j):return T[i]>>j&1
def cs(T,k):
 ans=0
 for W in combinations(range(len(T)),k):
  for V in permutations(W[1:]):
   if all(right(T,i,j) for i,j in zip(W[:1]+V,V+W[:1])):ans+=1
 return ans
def w(T,i,j):
 S=set(range(5))-{i,j}
 a,b,c=tuple(S)
 H=3 if ((right(T,a,b) and right(T,b,c) and right(T,c,a)) or (right(T,b,a) and right(T,c,b) and right(T,a,c))) else 1
 return H+sum(right(T,a,u)*right(T,u,v)*right(T,v,b)
              for a,b in ((i,j),(j,i)) for u,v in permutations(S,2))
B=[]
for z in map("".join,product("01",repeat=5)):
 if z!=min(s[i:]+s[:i] for s in (z,z[::-1]) for i in range(5)):continue
 T=rows(z);c3,c5=cs(T,3),cs(T,5);H=1+2*(c3+c5)
 pairs=list(combinations(range(5),2)); W={ (i,j):w(T,i,j) for i,j in pairs }
 sets=[]
 for R in (T,comp(T)):
  a=[sum(W[min(i,j),max(i,j)] for j in range(5) if right(R,i,j)) for i in range(5)]
  b=[a[i]+a[j]-W[i,j] for i,j in pairs]
  sets.extend((a,b))
 goal=(21-H)//2
 print(z,c3,c5,H,goal,'w=',list(W.values()),'low=',[sorted(set(x for x in t if x<=goal+1)) for t in sets])
 assert all(goal not in t for t in sets)

from functools import lru_cache

def subsets(T):
 n=len(T)
 @lru_cache(None)
 def dp(A,v):
  S=A^(1<<v)
  if S==0:return 1
  return sum(dp(S,u) for u in range(n) if S>>u&1 and right(T,u,v))
 return [1]+[sum(dp(A,v) for v in range(n) if A>>v&1) for A in range(1,1<<n)]
falsifiers=[(35,(24,117,121,114,96,65,1)),(77,(52,93,120,113,96,66,1)),(91,(116,237,248,241,226,192,128,1))]
for p,T in falsifiers:
 hs=subsets(T);n=len(T);M=(1<<n)-1
 assert hs[-1]==p
 for u in range(n):
  x=1<<u
  while True:
   xx=x
   for v in range(n):
    if x>>v&1:xx|=T[v]
   if xx==x:break
   x=xx
  assert x==M
 assert all(tuple(z.bit_count() for z in T).count(d)<min(q for q in range(2,p+1) if p%q==0)
            for d in range(n))
print("Eight bracelet weight tables and three rigid sevenfold witnesses pass")
