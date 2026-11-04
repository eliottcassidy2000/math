"""Exact construction checks; standard library, no external data."""
import itertools as it,json,random
from explore import hp,profile,rows,strong,inflate

def need(c,msg):
 if not c:raise RuntimeError(msg)

def streaming(a,b):
 """Vertices 0..b-1 are a transitive backbone; the rest are feedback."""
 n=len(a);total=0
 for p in it.permutations(range(b,n)):
  k=len(p);d={(0,)*(k+1):1}
  for x in range(b):
   nd={}
   for state,count in d.items():
    for j in range(k+1):
     if state[j]==0 and j and not(a[p[j-1]]>>x&1):continue
     t=list(state);t[j]=1 if j==k or (a[x]>>p[j]&1) else 2;t=tuple(t)
     nd[t]=nd.get(t,0)+count
   d=nd
  for s,c in d.items():
   if 2 in s:continue
   if any(s[j]==0 and not(a[p[j-1]]>>p[j]&1) for j in range(1,k)):continue
   total+=c
 return total

def backbone(n,k,mask):
 a=[0]*(n+k)
 for i in range(n):
  for j in range(i+1,n):a[i]|=1<<j
 edges=[(i,j) for i in range(n+k) for j in range(i+1,n+k) if j>=n]
 for b,(i,j) in enumerate(edges):
  if mask>>b&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a

def substitute(seed,a):
 n=len(a);b=list(a)+[0]*4
 for i in range(n):
  for j in range(1,5):
   if seed[0]>>j&1:b[i]|=1<<(n+j-1)
   else:b[n+j-1]|=1<<i
 for j in range(1,5):
  for k in range(1,5):
   if seed[j]>>k&1:b[n+j-1]|=1<<(n+k-1)
 return b

def induced(a,subset):
 ids=[i for i in range(len(a)) if subset>>i&1]
 return [sum(1<<j for j,y in enumerate(ids) if a[x]>>y&1) for x in ids]

def literal(a):
 return sum(all(a[x]>>y&1 for x,y in zip(p,p[1:])) for p in it.permutations(range(len(a))))

def factorial_stirling(r,k):
 return sum((-1)**(k-j)*__import__('math').comb(k,j)*j**r for j in range(k+1))

def run():
 streaming_count=0;literal_count=0
 for n in range(6):
  k=2
  for mask in range(1<<(n*k+1)):
   a=backbone(n,k,mask);h=hp(a)
   need(streaming(a,n)==h,('stream',n,mask));streaming_count+=1
   if n<=4:need(literal(a)==h,('literal',n,mask));literal_count+=1
 rng=random.Random(176)
 for _ in range(64):
  a=backbone(5,3,rng.randrange(1<<18));need(streaming(a,5)==hp(a),'three feedback');streaming_count+=1
 lo=[16,17,3,7,12];hi=[16,25,3,5,12]
 need(profile(lo,0)==(11,6,0,0,0),'low profile')
 need(profile(hi,0)==(13,6,0,0,0),'high profile')
 gate=[]
 for r in range(1,10):
  x=inflate(lo,0,r);y=inflate(hi,0,r);h=hp(x);g=hp(y)
  need(h==6*2**r-1 and g==h+2,'gate count')
  need(strong(x) and strong(y),'gate strong')
  need(sum((u^v).bit_count() for u,v in zip(x,y))==2,'single edge')
  gate.append({'r':r,'low':h,'high':g,'strong':True})
 # General insertion profiles for every labelled seed of order four, every port.
 profile_checks=0
 for mask in range(64):
  a=rows(4,mask)
  for v in range(4):
   p=profile(a,v)
   for r in range(1,5):
    calc=sum(c*factorial_stirling(r,k) for k,c in enumerate(p,1))
    need(hp(inflate(a,v,r))==calc,'profile identity');profile_checks+=1
 module_checks=0
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   a=rows(n,mask);h=hp(a);full=(1<<n)-1
   q=sum(hp(induced(a,s))*hp(induced(a,full^s)) for s in range(1,full))
   x=substitute(lo,a);y=substitute(hi,a)
   need(hp(x)==11*h+6*q and hp(y)==13*h+6*q,'arbitrary module')
   need(strong(x) and strong(y),'arbitrary module strong')
   module_checks+=1
 return {'module_checks':module_checks,'status':'proved identities and finite checks; completeness open','streaming_checks':streaming_count,'literal_checks':literal_count,'profile_checks':profile_checks,'gate_seeds':[lo,hi],'gate_profiles':[profile(lo,0),profile(hi,0)],'gates':gate}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
