from pathlib import Path
# independent bitmask subset Hamilton DP
def profile(a):
 n=len(a);N=1<<n;D=[[0]*n for _ in range(N)]
 for i in range(n):D[1<<i][i]=1
 for mask,row in enumerate(D):
  for i,w in enumerate(row):
   if not w:continue
   z=a[i]&~mask
   while z:
    b=z&-z;z-=b;D[mask|b][b.bit_length()-1]+=w
 h=[1]+[sum(row) for row in D[1:]]
 return h[-1],sum(h[z]*h[N-1-z] for z in range(N))
def blow(a,j,t):
 groups=[];c=0
 for i in range(len(a)):
  k=t if i==j else 1;groups.append(list(range(c,c+k)));c+=k
 o=[0]*c
 for x,g in enumerate(groups):
  for j,i in enumerate(g):
   for h in g[j+1:]:o[i]|=1<<h
   for y,hh in enumerate(groups):
    if a[x]>>y&1:
     for h in hh:o[i]|=1<<h
 return tuple(o)
a=(174,124,120,240,65,80,1,118)
b=(174,124,120,240,65,16,33,118)
if __name__=='__main__':
 for j in (0,1,2,3,4,7):
  print(j,[(x:=profile(blow(a,j,t)),tuple(v-u for u,v in zip(x,profile(blow(b,j,t))))) for t in range(1,5)])
