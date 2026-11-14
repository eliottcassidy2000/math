import random
from functools import cache
from itertools import combinations

def rows_of(n,code):
 rows=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if code>>k&1:rows[i]|=1<<j
  else:rows[j]|=1<<i
 return rows

def starts_table(R):
 @cache
 def starts(m,i):
  t=m^(1<<i)
  if not t:return 1
  return sum(starts(t,j) for j in range(len(R)) if t&R[i]&(1<<j))
 return starts

def table(R):
 st=starts_table(R)
 return [1]+[sum(st(m,i) for i in range(len(R)) if m>>i&1)
             for m in range(1,1<<len(R))]

def profile(R,v):
 n=len(R)
 full=(1<<n)-1
 H=table(R)
 other=[sum(1<<j for j,x in enumerate(R) if x>>i&1) for i in range(n)]
 P,Q=starts_table(R),starts_table(other)
 a=b=c=d=0
 for m in range(1,full+1):
  if not m>>v&1:continue
  u,w=H[full^m]*P(m,v),H[full^m]*Q(m,v)
  sign=(-1)**m.bit_count()
  a+=u;b+=w;c+=sign*u;d+=sign*w
 rem=full^(1<<v)
 f=sum(H[m]*H[rem^m] for m in range(full+1) if not m>>v&1)
 return H[full],a,b,c,d,f

def mkchain(blocks, roots):
 rows=[]
 offsets=[]
 for b in blocks:
  offsets.append(len(rows))
  rows.extend(x<<offsets[-1] for x in b)
 for i,b in enumerate(blocks):
  root=offsets[i]+roots[i]
  all_later=((1<<len(rows))-1)^((1<<(offsets[i]+len(b)))-1)
  for j in range(offsets[i],offsets[i]+len(b)): rows[j]|=all_later
  if i+1<len(blocks):
   rr=offsets[i+1]+roots[i+1]
   rows[root]^=1<<rr
   rows[rr]|=1<<root
 return rows

def pred(blocks, roots):
 F,A,B=1,0,0
 for R,r in zip(blocks,roots):
  h,p,q,s,t,f=profile(R,r)
  D=(h*F+q*A-t*B,p*F+f*A,s*F-(len(R)==1)*B)
  z,x,y=F+B,A-F,-B
  e,d=-s,-t
  E=(z+x+y,y,z) if len(R)==1 else (
   (h-e+q)*z+q*x+(h-e+q-d)*y,
   (p-h+f-q)*(z+y)+(f-q)*x+d*y,
   e*(z+y))
  assert E==(D[0]+D[2],D[1]-D[0],-D[2])
  F,A,B=D
 return F

def run():
 rng=random.Random(132613)
 for trial in range(1200):
  sizes=[]
  total=0
  for _ in range(rng.randint(1,6)):
   n=rng.randint(1,4)
   if total+n>11:break
   total+=n
   sizes.append(n)
  blocks=[rows_of(n,rng.randrange(1<<(n*(n-1)//2))) for n in sizes]
  roots=[rng.randrange(n) for n in sizes]
  R=mkchain(blocks,roots)
  assert pred(blocks,roots)==table(R)[-1],(R,blocks,roots)
 # A fixed terminal phase repair for several dyadic widths.
 for j in (1,2,3,4):
  f=1<<j
  TT=lambda n:[sum((1<<k) for k in range(i+1,n)) for i in range(n)]
  for lastmark,intercept in [(0,15),(1,17)]:
   blocks=[TT(1),TT(j+1),TT(2),TT(2),TT(1)]
   roots=[0,j,0,lastmark,0]
   assert pred(blocks,roots)==table(mkchain(blocks,roots))[-1]==18*f+intercept
 print('1,200 arbitrary-root multi-block chains checked by independent path subset DP')
if __name__=='__main__':run()
