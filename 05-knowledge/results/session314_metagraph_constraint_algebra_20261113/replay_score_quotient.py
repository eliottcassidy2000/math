from itertools import combinations
def ctl(n):
 edges=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
 free=[(i,j) for i,j in edges if (i,j)<=(n-1-j,n-1-i)]
 def sc(S):return tuple(sorted(abs(2*x-n+1) for x in S))[n%2::2]
 def strong(S):return all(sum(S[:i])>i*(i-1)//2 for i in range(1,n))
 V={};A=set()
 for w in range(1<<len(free)):
  S=list(range(n-1,-1,-1))
  for k,(i,j) in enumerate(free):
   if w>>k&1:
    S[i]-=1;S[j]+=1
    if i+j!=n-1:S[n-1-j]-=1;S[n-1-i]+=1
  z=sc(S);V[z]=strong(sorted(S))
  S[0]-=1;S[-1]+=1
  a=sc(S);A.add(tuple(sorted((a,z))))
 P={tuple(sorted((u,v))) for u,v in combinations(V,2) if sum(abs(a-b) for a,b in zip(u,v))==2 and (V[u] or V[v])}
 print(n,len(V),len(A),len(P),sorted(P-A)[:10],flush=True)
for n in (3,5,7,9):ctl(n)
