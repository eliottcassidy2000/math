from itertools import product
for k in range(0,8):
 m=k+2; b=[2**max(0,m-2-j) for j in range(m)]
 S=[[0]*m for _ in range(m)]
 for i in range(m):
  for j in range(i+2,m):S[i][j]=2**(j-i-2);S[j][i]=-S[i][j]
 words=[(1,)+z+(0,) for z in product(range(2),repeat=k)];total=0
 for w in words:
  for u in words:
   J=sum(b[i]*(u[i]-w[i]) for i in range(m))+sum(w[i]*S[i][j]*u[j] for i in range(m) for j in range(m))
   total += J%2
 centers=0 if k%2 else 1 if k%4==0 else 4
 expected=0 if k==0 else (2**k-centers)*2**(k-1)
 if total != expected: raise AssertionError((k,total,expected))
 print(k,total,'of',2**(2*k),'strong presentations')
