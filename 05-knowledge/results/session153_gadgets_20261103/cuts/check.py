from itertools import permutations,combinations

def tours(n):
 for code in range(1<<(n*(n-1)//2)):
  a=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(combinations(range(n),2)):
   a[i][j]=code>>k&1;a[j][i]=1-a[i][j]
  yield a

def path(a,p): return all(a[u][v] for u,v in zip(p,p[1:]))
def calc(a,m):
 n=len(a);ps=list(permutations(range(n)));good=[p for p in ps if path(a,p)];h=len(good)
 Eplus=sum(m>>p[-1]&1 for p in good); Pplus=sum(m>>p[0]&1 for p in good)
 X=sum((m>>p[k]&1) for p in ps for k in range(n) if path(a,p[:k]) and path(a,p[k:]))
 Y=sum((not(m>>p[k-1]&1)) for p in ps for k in range(1,n+1) if path(a,p[:k]) and path(a,p[k:]))
 B=sum(sum(not a[u][v] for u,v in zip(p,p[1:]))==1 for p in ps)
 child=[row+[int(not(m>>i&1))] for i,row in enumerate(a)]+[[m>>i&1 for i in range(n)]+[0]]
 d=[[0]*(n+1) for _ in range(1<<(n+1))]
 for i in range(n+1): d[1<<i][i]=1
 for z in range(1,1<<(n+1)):
  for i in range(n+1):
   if z>>i&1:
    d[z][i]+=sum(d[z^(1<<i)][j]*child[j][i] for j in range(n+1) if j!=i and z>>j&1)
 F=sum(d[-1])
 return h,Eplus,Pplus,X,Y,B,F

def child_hp(a,m):
 n=len(a)
 child=[row+[int(not(m>>i&1))] for i,row in enumerate(a)]+[[m>>i&1 for i in range(n)]+[0]]
 d=[[0]*(n+1) for _ in range(1<<(n+1))]
 for i in range(n+1): d[1<<i][i]=1
 for z in range(1,1<<(n+1)):
  for i in range(n+1):
   if z>>i&1:
    d[z][i]+=sum(d[z^(1<<i)][j]*child[j][i] for j in range(n+1) if j!=i and z>>j&1)
 return sum(d[-1])

def join(a,b):
 x=len(a);y=len(b)
 return [a[i]+[1]*y for i in range(x)]+[[0]*x+b[j] for j in range(y)]
count=0
for x,y in [(1,1),(1,2),(2,2),(2,3),(3,3),(1,4)]:
 for a in tours(x):
  for b in tours(y):
   t=join(a,b)
   for ma in range(1<<x):
    ha,ea,pa,xa,ya,ba,fa=calc(a,ma)
    for mb in range(1<<y):
     hb,eb,pb,xb,yb,bb,fb=calc(b,mb)
     f=child_hp(t,ma|(mb<<x))
     rhs=fa*hb+ha*fb-(ha-ea)*(hb-pb)-ea*pb-(ha-ea)*pb+xa*yb
     if rhs!=f: raise AssertionError((x,y,ma,mb,rhs,f))
     if ma==(1<<x)-1 and mb==0:
      if xa!=x*ha+ba or yb!=y*hb+bb: raise AssertionError("factor")
      if f!=2*ha*hb+(x*ha+ba)*(y*hb+bb): raise AssertionError("connector")
     count+=1
 print('sizes',x,y,'checked')
print('PASS',count,'ear cases')
