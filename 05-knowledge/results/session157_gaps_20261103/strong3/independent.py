def require(x):
 if not x: raise AssertionError("strong-family check failed")
from itertools import combinations
for n in range(1,6):
 N=2*n+4; word=[0,1]*n+[2,1]; out=[0]*N; a,b=N-2,N-1
 for i in range(len(word)):
  for j in range(i+1,len(word)):out[i]|=1<<j
  for z,v in enumerate((a,b)):
   if word[i] & (1<<z):out[i]|=1<<v
   else:out[v]|=1<<i
 out[b]|=1<<a
 D=[[0]*N for _ in range(1<<N)]
 for v in range(N):D[1<<v][v]=1
 for S in range(1,1<<N):
  for last in range(N):
   c=D[S][last]
   if not c:continue
   x=out[last]&~S
   while x:
    t=x&-x;x-=t;D[S|t][t.bit_length()-1]+=c
 h=sum(D[-1]); expression=(25+150*n+(668+60*n)*4**n+432*9**n)//225
 require(h==expression)
 print(n,N,h,expression)
