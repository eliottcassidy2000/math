import itertools,math
n=3
A=[[0]*n for _ in range(n)]
for i in range(n):A[i][(i+1)%n]=1
# s[i]=1 if i->u
s=[1,0,0]
P=[0]*(n+2)
def mul(p,q):
 r=[0]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):r[i+j]+=a*b
 return r
for pi in itertools.permutations(range(n)):
 gaps=[(1,1-s[pi[0]])]+[(A[pi[i-1]][pi[i]],s[pi[i-1]]*(1-s[pi[i]])) for i in range(1,n)]+[(1,s[pi[-1]])]
 p=[1]
 for a,b in gaps:p=mul(p,[a-b,b])
 for j,c in enumerate(p):P[j]+=c
print('P=',P)
for k in range(6):
 N=n+k; count=0
 def edge(a,b):
  if a<n and b<n:return A[a][b]
  if a<n:return s[a]
  if b<n:return 1-s[b]
  return a<b
 for p in itertools.permutations(range(N)):
  if all(edge(a,b) for a,b in zip(p,p[1:])):count+=1
 f=sum(c*j**k for j,c in enumerate(P))
 print(k,count,f)
 if count!=f: raise AssertionError((k,count,f))
