from itertools import combinations,permutations
def rows(order):
 r=[0]*10
 for i,j in combinations(range(5),2):r[i]|=1<<j;r[5+j]|=1<<(5+i)
 for i in range(5):
  for j in range(5):
   if (order[i]-order[j])%5 in (0,1,4):r[i]|=1<<(5+j)
   else:r[5+j]|=1<<i
 return tuple(r)
def heights(r):
 n=len(r);f=[[0]*n for _ in range(1<<n)]
 for i in range(n):f[1<<i][i]=1
 for s in range(1,1<<n):
  for i,p in enumerate(f[s]):
   if not p:continue
   x=r[i]&~s
   while x:
    b=x&-x;x^=b;f[s|b][b.bit_length()-1]+=p
 return sum(f[-1]),f
def flip(r,i,j):
 t=list(r);t[i]^=1<<j;t[j]^=1<<i;return tuple(t)
if __name__=='__main__':
 hist={}; notches=[]
 for p in permutations(range(1,5)):
  if p[0]>p[-1]:continue
  p=(0,*p);r=rows(p);h,_=heights(r);hist[h]=hist.get(h,0)+1
  hits=[]
  for i,j in combinations(range(10),2):
   a,_=heights(flip(r,i,j))
   if abs(a-h)==2:hits.append((i,j,a, j==i+5))
  if hits:notches.append((p,h,hits))
 print(hist,notches)

def coefficients(rows,x):
 n=len(rows);W=((1<<n)-1)^(1<<x);from functools import cache
 @cache
 def run(s,i):
  if not s:
   if i!=x and (i<0 or (rows[i]>>x)&1):return (int(i>=0),1)
   return (1,)
  out=[0]*(n+1)
  y=s if i<0 else s&rows[i]
  while y:
   b=y&-y;y^=b
   for k,v in enumerate(run(s^b,b.bit_length()-1)):out[k]+=v
  if i!=x and (i<0 or (rows[i]>>x)&1):
   for k,v in enumerate(run(s,x)):out[k+1]+=v
  while len(out)>1 and not out[-1]:out.pop()
  return tuple(out)
 return run(W,-1)

if __name__=='__main__':
 r=rows((0,1,4,3,2));s=flip(r,3,4)
 for x in range(10):
  a=coefficients(r,x);b=coefficients(s,x)
  print(x,a,b)

def tail_data(T):
 n=len(T);W=(1<<n)-1
 h,end=heights(T)
 opp=tuple(W^(1<<i)^row for i,row in enumerate(T))
 _,start=heights(opp)
 u=[0]*n;v=[[0]*n for _ in range(n)]
 zero=start[W][:]
 for A in range(1,W):
  B=A^W; e=end[A]; s=start[B];hA=sum(e)
  for j,x in enumerate(s):
   if x:
    u[j]+=x*hA
    for i,y in enumerate(e):
     if y:v[i][j]+=x*y
 for j in range(n):u[j]+=zero[j]
 return h,u,v,zero
def tail_params(X,w):
 h,u,v,z=X;n=len(u);J=[j for j in range(n) if w>>j&1];I=list(set(range(n))-set(J))
 U=sum(u[j] for j in J); V=sum(z[j]+sum(v[i][j] for i in I) for j in J)
 return U,h+V-U
if __name__=='__main__':
 r=rows((0,1,4,3,2));s=flip(r,3,4)
 X=tail_data(r);Y=tail_data(s)
 print('affine matches',[(w,tail_params(X,w),tail_params(Y,w)) for w in range(1,1024) if tail_params(X,w)[0]==tail_params(Y,w)[0] and abs(tail_params(X,w)[1]-tail_params(Y,w)[1])==2])
