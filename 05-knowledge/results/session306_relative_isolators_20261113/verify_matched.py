from probe import profile
def req(x):
 if not x:raise ValueError('check failed')
def make(a,b,c):
 p=[0,1+a,2+a+b,3+a+b+c]; n=p[-1]+1
 old={tuple(x) for x in [(p[0],p[2]),(p[1],p[3])]}
 A=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n): A[i][j]=int((i,j) not in old);A[j][i]=1-A[i][j]
 return A,p
def h(a,b,c):
 return 1+2**(a+c+1)*3**b+2*(2**a+2**c)*2**b
for u in range(1,12):
 for b in range(1,12):req(3*h(u+1,b,u+1)-1==2*h(u,b+1,u+1))
for u,b in [(1,1),(1,2),(2,1),(2,2),(3,1)]:
 a,p=make(u+1,b,u+1)
 z,q=make(u,b+1,u+1)
 perm=list(range(len(a)));i=p[1];perm[i-1],perm[i]=perm[i],perm[i-1]
 zz=[[z[i][j] for j in perm] for i in perm]
 req(sum(a[i][j]!=zz[i][j] for i in range(len(a)) for j in range(i+1,len(a)))==1)
 req(profile(a)[1]==h(u+1,b,u+1) and profile(z)[1]==h(u,b+1,u+1))
for a in range(1,11):
 for c in range(a,11):
  for v in (1,2):
   for d in range(1,11):
    e=a+c-v-d
    if d<=e<=10:
     equal=3*(2**a+2**c)==2**(v+1)*(2**d+2**e)
     req(equal==(a==c and a>=2 and v==1 and (d,e)==(a-1,a)))
print('PASS single-flip DP and coefficient-rigidity controls')
for ell in range(2,6):
 for v in range(ell,ell+5):
  for a in range(1,13):
   for c in range(a,13):
    for d in range(1,13):
     e=a+c-v-d
     if d<=e<=12:
      req(3**ell*(2**a+2**c)!=2**(v+ell)*(2**d+2**e))
print('PASS long-macro coefficient obstructions')
