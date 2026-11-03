for n in range(5):
 word=[1]+[0,0,1]*n+[0,0];N=len(word)+1;r=N-1;O=[0]*N
 for i,x in enumerate(word):
  for j in range(i+1,len(word)):O[i]|=1<<j
  if x:O[r]|=1<<i
  else:O[i]|=1<<r
 D=[[0]*N for _ in range(1<<N)]
 for i in range(N):D[1<<i][i]=1
 for s in range(1,1<<N):
  for i in range(N):
   c=D[s][i]; t=O[i]&~s
   while t:
    b=t&-t;t-=b;D[s|b][b.bit_length()-1]+=c
 H=sum(D[-1]); F=(192*8**n+28*n+53)//49
 print(n,N,H,F,H==F)
