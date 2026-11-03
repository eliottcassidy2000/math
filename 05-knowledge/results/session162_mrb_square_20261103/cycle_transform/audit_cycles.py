from itertools import combinations,permutations
rows=['0111100','0011010','0001110','0000011','0101011','1000000','1110010']
# matrix reconstructed from NOTE's listed arcs, with 6 ->5 and 2->3
n=7
arcs=[(0,1),(0,2),(0,3),(0,4),(5,0),(6,0),(1,2),(1,3),(4,1),(1,5),(6,1),(2,3),(2,4),(2,5),(6,2),(4,3),(3,5),(3,6),(4,5),(4,6),(6,5)]
A=[[0]*n for _ in range(n)]
for i,j in arcs:A[i][j]=1
B=[r[:] for r in A];B[0][3],B[3][0]=0,1

def paths(M):
 N=len(M);D=[[0]*N for _ in range(1<<N)]
 for i in range(N):D[1<<i][i]=1
 for s in range(1,1<<N):
  for i in range(N):
   if D[s][i]:
    for j in range(N):
     if not s>>j&1 and M[i][j]:D[s|1<<j][j]+=D[s][i]
 h=[1]+[sum(D[s]) for s in range(1,1<<N)]
 return h

def covers(M,r):
 N=len(M);h=paths(M);v=[0]*(1<<N);v[0]=1
 for _ in range(r):
  w=[0]*(1<<N)
  for s in range(1<<N):
   t=s
   while True:
    w[s]+=v[s^t]*h[t]
    if t==0:break
    t=(t-1)&s
  v=w
 return v[-1]

def cy(M):
 out=[]
 for k in (3,5,7):
  for s in combinations(range(7),k):
   for q in permutations(s[1:]):
    p=(s[0],)+q
    if all(M[x][y] for x,y in zip(p,p[1:]+p[:1])):out.append(p)
 return set(out)
old,new=cy(A),cy(B)
if (paths(A)[-1],paths(B)[-1],len(new-old),len(old-new)) != (115,111,8,6): raise RuntimeError('unexpected counts')
print('full path counts DP:',paths(A)[-1],paths(B)[-1])
print('full path counts permutations:',sum(all(A[x][y] for x,y in zip(p,p[1:])) for p in permutations(range(7))),sum(all(B[x][y] for x,y in zip(p,p[1:])) for p in permutations(range(7))))
print('created/destroyed cycles:',len(new-old),len(old-new))
for r in (0,1,2,3):
 weighted=0
 for sign,cs in ((1,new-old),(-1,old-new)):
  for C in cs:
   remain=[x for x in range(7) if x not in C]
   sub=[[A[x][y] for y in remain] for x in remain]
   val=1 if not remain else covers(sub,r)
   weighted+=sign*val
 a,b=covers(A,r),covers(B,r)
 if b-a!=2*r*weighted:raise RuntimeError((r,a,b,weighted))
 print('r, F_before, F_after, signed residual sum, RHS:',r,a,b,weighted,2*r*weighted)
