from itertools import combinations
from math import comb
C=(24,9,3,4,14)
X=(388,449,202,195,335,351,129,48,204)
Y=(284,33,506,418,138,17,59,355,114)
P=(44,33,2,6,15,28)
U=(766,248,1018,832,264,280,944,312,3,434)
V=(32,909,1,821,295,454,415,541,133,373)
def inflate(R,U):
 m=len(U);q=(0,)*m+tuple(range(1,len(R)))
 return tuple(sum(1<<j for j,b in enumerate(q) if i!=j and ((U[i]>>j&1) if a==b else (R[a]>>b&1))) for i,a in enumerate(q))
def add(A,B):
 n=len(A)
 return tuple(x+((not B>>j&1)<<n)for j,x in enumerate(A))+(B,)
def coeff(A,B,whole=False):
 n=len(A);N=(1<<n)-1;v=min(B.bit_count(),(N^B).bit_count())+1
 D=[[0]*(n+1) for _ in range(N+1)]
 for j in range(n):D[1<<j][j]=1
 ans=[]
 for r in range(v+1):
  E=[[0]*(n+1)for _ in range(N+1)]
  if r==0:E[0][n]=1
  for S,w in enumerate(D):
   for j,d in enumerate(w):
    if not d:continue
    R=(B if j==n else A[j])&~S
    while R:
     b=R&-R;R^=b
     D[S|b][b.bit_length()-1]+=d
    if j!=n and not B>>j&1:E[S][n]+=d
  ans.append([sum(d) for d in D] if whole else sum(D[N]))
  D=E
 return tuple(ans)
def filled(c,m):
 return sum(c[j]*sum((-1)**(j-t)*comb(j,t)*t**m for t in range(j+1))for j in range(1,len(c)))
def two(A,B):
 L=coeff(A,B,True);N=(1<<len(A))-1
 L[0][0]=1 # empty slot
 out=[0]*(2*len(L)-1)
 for a,u in enumerate(L):
  for b,v in enumerate(L):
   out[a+b]+=sum(x*v[N^s] for s,x in enumerate(u) if x)
 return out
def run():
 A=inflate(C,X);B=inflate(C,Y)
 D=inflate(P,U);E=inflate(P,V)
 for a,b,k in ((A,B,7167),(D,E,2048),(D,E,5120),(D,E,6144),(add(D,2048),add(E,2048),38912)):
  x=coeff(a,k);y=coeff(b,k)
  want={
   7167:((532623,1068387,534717),)*2,
   2048:((4591771,7655505,3005811),(4591771,7655505,3006531)),
   5120:((4591771,13775313,13435435,2947888),(4591771,13775313,13439755,2949328)),
   6144:((4591771,9189287,4595601,0),)*2,
   38912:((7655505,15316755,7659335,0,0),)*2,
  }[k]
  if (x,y)!=want:raise ValueError((k,x,y))
  print(len(a),k,x,y,[(filled(x,j),filled(y,j))for j in range(1,5)],flush=True)
  if k in (7167,6144,38912):
   p,q=two(a,k),two(b,k)
   z=tuple(v-w for v,w in zip(p,q))
   want={7167:(-134400,-397920,-394080,-130560),
         6144:(-2053440,-6143040,-6134400,-2044800),
         38912:(-5330880,-15963840,-15949440,-5316480)}[k]
   if z[:4]!=want or any(z[4:]):raise ValueError((k,z))
   print('Q2 layout difference',z,flush=True)
if __name__=='__main__':run()
