from itertools import permutations
from collections import Counter

def vals(n,m):
 e=list(__import__('itertools').combinations(range(n),2)); a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(e):
  if m>>k&1:a[i][j]=True
  else:a[j][i]=True
 def enc(p):
  z=0
  for k,(i,j) in enumerate(e):z|=a[p[i]][p[j]]<<k
  return z
 can=min(map(enc,permutations(range(n))))
 p=[[0]*n for _ in range(1<<n)]
 for i in range(n):p[1<<i][i]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    prev=mask^(1<<j)
    if prev:p[mask][j]=sum(p[prev][i] for i in range(n) if prev>>i&1 and a[i][j])
 h=sum(p[-1]); reach=[row[:] for row in a]
 for i in range(n):reach[i][i]=True
 for k in range(n):
  for i in range(n):
   for j in range(n):reach[i][j]|=reach[i][k] and reach[k][j]
 strong=all(map(all,reach))
 return can,h,strong

for n in range(1,6):
 seen={}
 for mask in range(1<<(n*(n-1)//2)):
  can,h,strong=vals(n,mask)
  seen[can]=(h,strong)
 allc=Counter(h for h,s in seen.values()); sc=Counter(h for h,s in seen.values() if s and n>1)
 print('n=',n,'classes=',len(seen),'all H=3,5,9:',[allc[x] for x in (3,5,9)],'strong:',[sc[x] for x in (3,5,9)])
 if allc[3] != max(0,n-2):raise RuntimeError('H=3 formula failed')
