from itertools import permutations,combinations
from collections import Counter

def canon(T):
 n=len(T)
 return min(sum(T[p[i]][p[j]]<<k for k,(i,j) in enumerate(combinations(range(n),2))) for p in permutations(range(n)))

def con(T):return tuple(tuple(T[j][i] for j in range(len(T))) for i in range(len(T)))

def rev(T,P):
 U=list(map(list,T))
 for a,b in zip(P,P[1:]): U[a][b],U[b][a]=0,1
 return tuple(map(tuple,U))

def blue(T,P):
 R=P[::-1]
 return all(T[P[i]][P[j]]==T[R[j]][R[i]] for i in range(len(P)) for j in range(i))

def stable(T):
 n=len(T)
 return sum(all(T[i][j]==T[p[i]][p[j]] for i,j in combinations(range(n),2)) for p in permutations(range(n)))

def parent(T,U):
 n=len(T);d=list(map(sum,T));e=list(map(sum,U));ans=Counter()
 for a,b in permutations(range(n),2):
  score=[s-(i==a)+(i==b) for i,s in enumerate(d)]
  if sorted(score)!=sorted(e):continue
  inp=[1-(i==a) for i in range(n)];out=[1-(i==b) for i in range(n)]
  ix=sorted(range(n),key=lambda i:e.count(score[i]))
  def loop(j, f, used, ind, outd, succ):
   if j==n:
    if ind!=inp or outd!=out:return
    P=[a]
    while P[-1] in succ:P.append(succ[P[-1]])
    assert len(P)==n and P[-1]==b
    ans[blue(T,P)]+=1
    return
   u=ix[j]
   for v in range(n):
    if v in used or e[v]!=score[u]:continue
    I,O,S=ind[:],outd[:],succ.copy()
    for c,w in f.items():
     if T[u][c]==U[v][w]:continue
     x,y=(u,c) if T[u][c] else (c,u)
     I[y]+=1;O[x]+=1
     if I[y]>inp[y] or O[x]>out[x]:break
     S[x]=y
     z=y
     while z in S and z!=x:z=S[z]
     if z==x:break
    else:loop(j+1,f|{u:v},used|{v},I,O,S)
  loop(0,{},set(),[0]*n,[0]*n,{})
 return ans

def make(n,w):
 T=[[0]*n for _ in range(n)]
 F=[(i,j) for i,j in combinations(range(n),2) if j>i+1]
 for i in range(n-1):T[i][i+1]=1
 for k,(i,j) in enumerate(F):
  T[i][j]=w>>k&1;T[j][i]=1-T[i][j]
 return tuple(map(tuple,T))

def GF(n):
 # Bivariate score words, path reversal convention, all and reflection fixed
 E=[(i,j) for i,j in combinations(range(n),2) if j>i+1];r=lambda i:n-1-i
 s=tuple(int(i<n-1) for i in range(n))
 def product(factors):
  D=Counter({s:1})
  for W in factors:
   Z=Counter()
   for a,m in D.items():
    for w in W:
     b=list(a)
     for i in w:b[i]+=1
     Z[tuple(b)]+=m
   D=Z
  return D
 allD=product([[(i,),(j,)] for i,j in E])
 bf=[]
 for i,j in E:
  f=(r(j),r(i))
  if f<(i,j):continue
  if f==(i,j):bf.append([(i,),(j,)])
  else:bf.append([(i,r(j)),(j,r(i))])
 return allD,product(bf)

for n in range(3,6):
 D={};canmemo={}
 def cc(T):
  if T not in canmemo:canmemo[T]=canon(T)
  return canmemo[T]
 def mk(T):return min(cc(T),cc(con(T)))
 expected=Counter(); scores=Counter();bscores=Counter()
 for w in range(1<<((n-1)*(n-2)//2)):
  T=make(n,w);u=mk(T);P=list(range(n));v=mk(rev(T,P));B=blue(T,P)
  D[u]=T;expected[u,v,B]+=1;s=tuple(map(sum,T));scores[s]+=1
  if B:bscores[s]+=1
 assert GF(n)==(scores,bscores)
 found=Counter()
 for u,T in D.items():
  for v,U in D.items():
   allu=stable(T);allv=stable(U); assert allv==stable(con(U))
   c=1+(cc(T)!=cc(con(T)));variants=[U] if cc(U)==cc(con(U)) else [U,con(U)]
   A=Counter()
   for V in variants:A+=parent(T,V)
   for b,k in A.items():
    assert (c*k)%(allu*allv)==0
    found[u,v,b]=c*k//(allu*allv)
 assert found==expected
 print(n,len(D),sum(expected.values()),sum(k for (u,v,b),k in expected.items() if b),len({tuple(sorted((u,v))) for u,v,b in expected}))
