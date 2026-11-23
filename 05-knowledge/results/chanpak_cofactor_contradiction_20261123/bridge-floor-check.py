from math import gcd,factorial
def heights(rows):
 n=len(rows); table=[[0]*n for _ in range(1<<n)]
 for i in range(n): table[1<<i][i]=1
 for S,row in enumerate(table):
  for i,v in enumerate(row):
   if not v:continue
   mask=rows[i]&~S
   while mask:
    b=mask&-mask;mask-=b
    table[S|b][b.bit_length()-1]+=v
 return [1]+[sum(r) for r in table[1:]]
def slot(h):
 N=len(h); old=[1]+[0]*(N-1); vals=[]
 for k in range(1,N.bit_length()):
  new=[0]*N
  for S in range(1,N):
   A=S
   while A:
    new[S]+=h[A]*old[S^A]; A=(A-1)&S
  vals.append(new[-1]);old=new
 return vals
def ext(M,r):
 n=len(M); b=sum(1<<(n+2*j) for j in range(r)); rows=[x|b for x in M]
 for j in range(2*r):
  rows.append(sum(1<<(n+i) for i in range(j+1,2*r))|((1<<n)-1 if j%2 else 0))
 return rows
def count(n,p,cut):
 fixed=[1<<(i+1) for i in range(n-1)]+[0]
 free=[(i,j) for i in range(n) for j in range(i+2,n)]
 hist={}
 for t in range(1<<len(free)):
  rows=fixed.copy()
  for k,(i,j) in enumerate(free):
   if (t>>k)&1: rows[i]|=1<<j
   else: rows[j]|=1<<i
  hs=heights(rows); H=hs[-1]
  if H%p or H>cut:continue
  P=slot(hs)
  pp=( (2*P[0]+P[1])%p==0, (3*P[0]+3*P[1]+P[2])%p==0, all(v%p==0 for v in P[:p]))
  hist[H,pp]=hist.get((H,pp),0)+1
 print(n,p,sorted(hist.items()))
 return hist
if __name__=="__main__":
 for M in [(106,52,121,114,97,64,2),(10,116,121,114,97,65,1)]:
  assert all(((M[i]>>j)&1)+((M[j]>>i)&1)==1 for i in range(len(M)) for j in range(i))
  for start in range(len(M)):
   reach=1<<start
   for _ in M:
    for i,out in enumerate(M):
     if reach>>i&1:reach|=out
   assert reach==(1<<len(M))-1
  P=slot(heights(M))
  H=[heights(ext(M,r))[-1] for r in range(3)]
  print(P,gcd(*(factorial(i)*v for i,v in enumerate(P))),H)
 for n in range(3,7):assert not count(n,7,49)
 assert count(7,7,49)=={(35,(False,True,False)):70,(49,(False,False,False)):245,(49,(True,True,True)):98}
 assert count(5,5,25)=={(5,(False,False,False)):10,(15,(False,False,False)):5,(15,(True,True,True)):3}
