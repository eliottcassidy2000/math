from functools import lru_cache
rows=[(80265,2181427721976895),(80387,31503056742464543),(80405,22513037423950911),(80875,33617973800943135),(84259,36024248539118623)]
for h,code in rows:
 n=11;A=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   A[i][j]=bool(code>>k&1);A[j][i]=not A[i][j];k+=1
 @lru_cache(None)
 def F(v,s):
  return 1 if not s else sum(F(w,s^(1<<w)) for w in range(n) if s>>w&1 and A[v][w])
 z=sum(F(v,((1<<n)-1)^(1<<v)) for v in range(n))
 for v in range(n):
  seen={v}; before=0
  while len(seen)!=before:
   before=len(seen);seen|={j for i in seen for j in range(n) if A[i][j]}
  if len(seen)!=n:raise ValueError((h,v,'not strong'))
 if z!=h:raise ValueError((h,z))
 print(h,code,z,'strong')
