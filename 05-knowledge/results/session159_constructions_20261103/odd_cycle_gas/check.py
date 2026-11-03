from itertools import combinations, permutations
from functools import lru_cache
from math import comb, factorial

@lru_cache(None)
def stirling(n,k):
    if n==k==0:return 1
    if not n or not k:return 0
    return stirling(n-1,k-1)+k*stirling(n-1,k)

def predicted_cycles(a,b,c):
    if (a+b+c)%2==0:return 0
    return sum(factorial(k)**3*stirling(a,k)*stirling(b,k)*stirling(c,k)//k for k in range(1,min(a,b,c)+1))

@lru_cache(None)
def predicted_poly(a,b,c):
    if a==b==c==0:return (1,)
    # Recurrence is symmetric in block sizes; permute so a>0.
    if a==0:return predicted_poly(b,c,a)
    old=predicted_poly(a-1,b,c)
    result=list(old)+[0]*(min(a,b,c)+1-len(old))
    for i in range(1,a+1):
      for j in range(1,b+1):
       for k in range(1,c+1):
        number=comb(a-1,i-1)*comb(b,j)*comb(c,k)*predicted_cycles(i,j,k)
        if not number:continue
        for q,value in enumerate(predicted_poly(a-i,b-j,c-k)):
            result[q+1]+=number*value
    while len(result)>1 and not result[-1]:result.pop()
    return tuple(result)

def direct(a,b,c):
    n=a+b+c
    if n==0:return (1,),1
    blocks=[0]*a+[1]*b+[2]*c
    def edge(i,j):return (blocks[i]==blocks[j] and i<j) or ((blocks[j]-blocks[i])%3==1)
    cycles=[]
    for size in range(3,n+1,2):
      for S in combinations(range(n),size):
       for p in permutations(S[1:]):
        seq=(S[0],)+p
        if all(edge(seq[i],seq[(i+1)%size]) for i in range(size)):
          cycles.append(sum(1<<x for x in S))
    @lru_cache(None)
    def pack(mask):
      if not mask:return (1,)
      first=mask & -mask
      r=list(pack(mask^first))+[0]
      for cycle in cycles:
       if cycle&first and cycle&mask==cycle:
        for i,v in enumerate(pack(mask^cycle)):r[i+1]+=v
      while len(r)>1 and r[-1]==0:r.pop()
      return tuple(r)
    dp={(1<<i,i):1 for i in range(n)}
    for mask in range(1,1<<n):
      for i in range(n):
       v=dp.get((mask,i),0)
       if not v:continue
       for j in range(n):
        if not (mask>>j)&1 and edge(i,j):
         key=(mask|1<<j,j);dp[key]=dp.get(key,0)+v
    return pack((1<<n)-1),sum(dp.get(((1<<n)-1,i),0) for i in range(n))

def path_formula(a,b,c):
  return sum(factorial(k)**3*stirling(a+1,k+1)*stirling(b+1,k+1)*stirling(c+1,k+1) for k in range(min(a,b,c)+1))+sum(factorial(k)**3*stirling(a,k)*stirling(b,k)*stirling(c,k) for k in range(1,min(a,b,c)+1))

def p2(a,b):
    A=lambda n:3**n-2**(n+1)+1
    return (A(a)*A(b)+2*(2-2**a)*(2-2**b)+(1+(-1)**a)*(1+(-1)**b))//4

count=0
for n in range(0,9):
 for a in range(0,n+1):
  for b in range(0,n-a+1):
   c=n-a-b
   actual,h=direct(a,b,c); predicted=predicted_poly(a,b,c)
   if actual!=predicted or h!=sum(x*2**i for i,x in enumerate(actual)) or h!=path_formula(a,b,c):
      raise ValueError((a,b,c,actual,predicted,h,path_formula(a,b,c)))
   if a>=1 and b>=1 and c==2 and (actual[2] if len(actual)>2 else 0)!=p2(a,b):raise ValueError((a,b,c,actual))
   count+=1
print('Verified odd-cycle packing polynomial and path count for',count,'nonnegative triples of total order 0..8.')
print('J(2,3,2)=',predicted_poly(2,3,2),'H=',path_formula(2,3,2))
