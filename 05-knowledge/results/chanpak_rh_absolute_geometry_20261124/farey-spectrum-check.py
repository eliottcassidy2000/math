from collections import Counter
def certificate(p,t_num,t_den):
 def S(x):return p if x==0 else 0 if x==p else -pow(x,-1,p)%p
 def T(x):return p if x==p-1 else 0 if x==p else -pow(x+1,-1,p)%p
 group=[]; label={}
 for x in range(p+1):
  if x in label:continue
  orb=sorted({x,T(x),T(T(x))})
  assert all(y not in label for y in orb)
  for y in orb:label[y]=len(group)
  group.append(orb)
 n=len(group);A=[[0]*n for _ in range(n)]
 for i,orb in enumerate(group):
  for x in orb:A[i][label[S(x)]]+=t_den
  A[i][i]-=t_num*len(orb)
 assert all(A[i][j]==A[j][i] for i in range(n) for j in range(n))
 prev=1; pos=0
 for k in range(n):
  d=A[k][k];assert d
  pos+= d*prev>0
  for i in range(k+1,n):
   for j in range(i,n):
    v,r=divmod(A[i][j]*d-A[i][k]*A[k][j],prev)
    assert r==0
    A[i][j]=A[j][i]=v
  prev=d
 assert 9*t_num**2>8*t_den**2
 print(p,n,pos)
 return pos
assert certificate(433,18857,20000)==2
p=7
L=[(x+1)%p for x in range(p)]+[p]
R=[x*pow(x+1,-1,p)%p for x in range(p-1)]+[p,1]
A=[[int(L[j]==i)+int(R[j]==i) for j in range(8)] for i in range(8)]
I=[[int(i==j) for j in range(8)] for i in range(8)]
def multiply(U,V):
 return [[sum(a*b for a,b in zip(row,col)) for col in zip(*V)] for row in U]
def polymul(u,v):
 return [sum(u[j]*v[k-j] for j in range(len(u)) if 0<=k-j<len(v))
         for k in range(len(u)+len(v)-1)]
N=[0];tr=[8];B=I
for n in range(1,17):
 B=multiply(A,B);tr.append(sum(B[i][i] for i in range(8)));N.append(B[6][1])
D=[1]
for n in range(1,9):
 v,r=divmod(-sum(tr[j]*D[n-j] for j in range(1,n+1)),n)
 assert r==0
 D.append(v)
assert D==polymul(polymul(polymul([1,-2],[1,-1]),[1,0,2]),[1,1,0,2,4])
P=[0,0,0,4,-4];Q=polymul([1,-2],[1,1,0,2,4])
assert polymul(N,Q)[:len(N)]==P+[0]*(len(N)-len(P))
print('p=7 determinant and entrance series passed')
