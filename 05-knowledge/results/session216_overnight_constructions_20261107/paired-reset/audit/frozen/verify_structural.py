from inherited import profile,need
import itertools,json

def seven(g):return tuple(profile(g)[i] for i in (0,1,2,3,4,5,7))
def star(p,q):
 h,s,e,L,R,d,b=p;H,S,E,LL,RR,D,B=q
 return (d*H+(R-d)*LL+(e-d)*S+(h-e-R+d)*D,s*D,d*E,b*LL+D*(L-b),B*R+d*(RR-B),d*D,b*B)
def glue(g,f):
 n=len(g);m=len(f);out=list(g)+[0]*(m-1)
 for i in range(1,n):out[i]|=((1<<(m-1))-1)<<n
 for i in range(m):
  ii=0 if i==0 else n+i-1
  for j in range(m):
   if f[i]>>j&1:out[ii]|=1<<(0 if j==0 else n+j-1)
 return out

def graphs(n):
 for bits in range(1<<(n*(n-1)//2)):
  g=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if bits>>k&1:g[i]|=1<<j
    else:g[j]|=1<<i
    k+=1
  yield g
bank=[g for n in range(2,4) for g in graphs(n)];controls=0
for g in bank:
 for f in bank:
  need(seven(glue(g,f))==star(seven(g),seven(f)),'ordered root glue');controls+=1
core=[4,1,10,3];roots=[]
for mask in range(16):roots.append([mask<<1]+[(core[j]<<1)|(not(mask>>j&1)) for j in range(4)])
for g in roots:
 for f in roots:
  need(seven(glue(g,f))==star(seven(g),seven(f)),'double H5 root glue');controls+=1
stable=0
for n in range(2,6):
 for g in graphs(n):
  h,s,e,L,R,d,b=seven(g)
  if s==e==d:
   K=h-L-R+b;need(K>=d and (K-d)%2==0,'defect nonnegative integral');stable+=1
print(json.dumps({'status':'PASS','literal_retained_root_products':controls,'stable_defect_controls':stable},indent=2))
