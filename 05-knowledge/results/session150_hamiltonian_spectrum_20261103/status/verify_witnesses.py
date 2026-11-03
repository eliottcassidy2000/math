from itertools import permutations
for goal,code in [(49,201564224),(63,104288510),(343,136003753)]:
 n=8; A=[[False]*n for _ in range(n)]; k=0
 for i in range(n):
  for j in range(i+1,n):
   A[i][j]=bool(code & (1<<k)); A[j][i]=not A[i][j]; k+=1
 count=sum(all(A[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))
 R=[row[:] for row in A]
 for k in range(n):
  for i in range(n):
   for j in range(n): R[i][j]|=R[i][k] and R[k][j]
 strong=all(R[i][j] for i in range(n) for j in range(n))
 print(goal,code,count,strong)
 if not (count==goal and strong): raise AssertionError('check failed')
