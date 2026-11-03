from itertools import product
for m in range(1,9):
 b=[2**max(0,m-2-j) for j in range(m)]
 S=[[0]*m for _ in range(m)]
 for i in range(m):
  for j in range(i+2,m):S[i][j]=2**(j-i-2);S[j][i]=-S[i][j]
 words=list(product(range(2),repeat=m)); total=0
 for w in words:
  for u in words:
   J=sum(b[i]*(u[i]-w[i]) for i in range(m))+sum(w[i]*S[i][j]*u[j] for i in range(m) for j in range(m))
   total += J%2
 expected=2**(2*m-1)-(2**(m-1) if m%4==0 else 0)
 if total != expected: raise AssertionError((m,total,expected))
 print(m,total,'of',2**(2*m),'Boolean pairs')
