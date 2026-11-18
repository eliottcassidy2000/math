from itertools import product

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def prime(p):return p>=2 and all(p%d for d in range(2,int(p**.5)+1))
def leg(a,p):
 z=pow(a%p,(p-1)//2,p);return -1 if z==p-1 else z
def branch(x,k,b,p):return ((3*x+b)*pow(pow(2,k,p),-1,p))%p
def word(n,ks,b):
 for k in ks:
  z=3*n+b
  if z<=0 or z%(1<<k) or (z//(1<<k))%2==0:return None
  n=z//(1<<k)
 return n
def inverse_cylinder(ks,b):
 r=1;s=1
 for k in reversed(ks):
  s+=k;m=1<<s;r=((1<<k)*r-b)*pow(3,-1,m)%m
 return r,1<<s
checks=0
for p in range(7,200):
 if not prime(p) or p%4!=3:continue
 table={7:(1,-1),11:(-1,1),19:(-1,-1),23:(1,1)}
 need((leg(2,p),leg(3,p))==table[p%24],('sign table',p))
 for b in [-1,1]:
  for k in range(1,7):
   sg=leg(3,p)*leg(2,p)**k
   for x in range(p):
    for y in range(p):
     if x==y:continue
     need(leg(branch(y,k,b,p)-branch(x,k,b,p),p)==sg*leg(y-x,p),('arc action',p,b,k,x,y))
   checks+=1
for b in [-1,1]:
 for L in range(1,5):
  for ks in product(range(1,4),repeat=L):
   r,m=inverse_cylinder(ks,b)
   need(r%2==1, 'odd cylinder')
   # Shift above zero to avoid low boundary for minus branch.
   for j in [1,2,9]:need(word(r+j*m,ks,b)!=None,('cylinder',b,ks))
   # All other residue classes are rejected, up to a full period at modest size.
   if m<=1024:
    found=[n%m for n in range(m+1,2*m,2) if word(n,ks,b)!=None]
    need(found==[r],('unique cylinder',b,ks,found,r))
   for p in [7,11,19,23]:
    for a in range(p):
     n=r+m*((a-r)*pow(m,-1,p)%p)+m*p
     need(n%p==a and word(n,ks,b)!=None,('CRT lift',p,a))
# Literal Hamilton path counts for two manageable Paley tournaments.
def hp(p):
 r=[sum(1<<j for j in range(p) if j!=i and leg(j-i,p)==1) for i in range(p)]
 d=[[0]*p for _ in range(1<<p)]
 for i in range(p):d[1<<i][i]=1
 for m in range(1,1<<p):
  for i,v in enumerate(d[m]):
   if not v:continue
   a=r[i]&~m
   while a:
    q=a&-a;a-=q;j=q.bit_length()-1;d[m|q][j]+=v
 return sum(d[-1])
for p in [7,11]:
 h=hp(p);need(h%(p*(p-1)//2)==0,('affine divisibility',p,h));print('Paley',p,'H',h,'known affine divisor',p*(p-1)//2)
# Clock subgroup may be smaller than full square subgroup: track it, do not assume equality.
for p in [7,11,19,23]:
 K={1};old=None
 while K!=old:
  old=set(K);K|={2*x%p for x in K}|{3*x%p for x in K}
 print('p',p,'multiplier_order',len(K),'branch affine order',p*len(K),'signs',(leg(2,p),leg(3,p)))
print('PASS',checks,'branch action cases; all finite cylinders and CRT controls')
