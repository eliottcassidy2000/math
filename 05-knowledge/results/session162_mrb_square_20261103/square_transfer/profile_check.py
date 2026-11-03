def profile(n,back,r):
 out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if (j,i) in back: out[j]|=1<<i
   else: out[i]|=1<<j
 N=1<<n; dp=[[0]*n for _ in range(N)]; path=[0]*N;path[0]=1
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,N):
  path[m]=sum(dp[m]);
  for j in range(n):
   if dp[m][j]:
    x=out[j]&~m
    while x:
     bit=x&-x; i=bit.bit_length()-1; dp[m|bit][i]+=dp[m][j];x-=bit
 f=[0]*N;f[0]=1
 for k in range(r):
  z=[0]*N
  for m in range(N):
   sub=m
   while True:
    z[m]+=f[m^sub]*path[sub]
    if not sub:break
    sub=(sub-1)&m
  f=z
 return f[-1]
for d in range(1,8):
 for r in range(1,5):
  a=profile(d+1,{(d,0)},r)
  b=r**(d+1)+r*((r+1)**(d-1)-(r-1)**(d-1))
  if a != b: raise AssertionError((d,r,a,b))
for j in range(3):
 n=9
 # Vertex labels 1–9 map to indices 0–8.
 edges={(5,i) for i in range(5)}|{(6,5),(8,6)}
 if j>=1:edges.add((2,0))
 if j>=2:edges.add((3,0))
 got=profile(n,edges,1)
 if got != [99,153,207][j]: raise AssertionError((j,got))
print('PASS d=1..7 and r=1..4; direct 9-vertex matched chain 99,153,207')
