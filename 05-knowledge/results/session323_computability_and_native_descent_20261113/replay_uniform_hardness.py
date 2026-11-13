from math import gcd
def need(ok):
 if not ok: raise ValueError('control failed')
def count(i):
 n=i+2
 A=[sum(1<<w for w in range(v+1,n)) for v in range(n)]
 A[0]^=1<<(n-1); A[-1]^=1
 D=[[0]*n for _ in range(1<<n)]
 for j in range(n):D[1<<j][j]=1
 for s,row in enumerate(D):
  for j,x in enumerate(row):
   if not x:continue
   w=A[j]&~s
   while w:
    q=w&-w;w-=q
    D[s|q][q.bit_length()-1]+=x
 return sum(D[-1])
for i in range(1,13):
 need(count(i)==(1<<i)+1)
for k in range(8):
 f=(1<<(1<<k))+1
 for i in range(1,1<<k):need(gcd((1<<i)+1,f)==1)
print('R_i path count and Fermat separation controls passed')
