BANK={'d': [((186, 92, 73, 80, 36, 78, 145, 62), (2, 5), (401, 6080, 395, 6080)), ((170, 220, 57, 48, 161, 194, 157, 12), (0, 5), (409, 6184, 415, 6184)), ((194, 44, 121, 209, 227, 73, 130, 38), (2, 3), (441, 6576, 435, 6576)), ((74, 4, 233, 146, 103, 203, 138, 19), (2, 4), (369, 5424, 393, 5424)), ((198, 204, 136, 209, 167, 79, 148, 32), (5, 6), (383, 5528, 359, 5528))], 'q': [((182, 156, 88, 17, 32, 206, 155, 28), (0, 6), (357, 5136, 357, 4944)), ((238, 188, 184, 240, 225, 192, 134, 0), (3, 6), (81, 2256, 81, 2232)), ((46, 100, 8, 50, 39, 196, 157, 31), (4, 6), (335, 5408, 335, 5384)), ((234, 244, 105, 82, 37, 200, 144, 28), (5, 6), (461, 6728, 461, 6752)), ((186, 188, 185, 16, 96, 72, 143, 56), (5, 7), (147, 3120, 147, 3096)), ((238, 148, 24, 242, 33, 70, 150, 52), (3, 4), (393, 5592, 393, 5976))]}

from math import isqrt

DYADIC = {
 (19,8):(38,4,8,19,39,14),
 (217,12):(70,60,24,17,33,77,158,63),
 (223,12):(230,28,8,17,37,78,158,62),
 (599,13):(178,100,137,19,38,76,157,58),
}

def profile(A):
 n=len(A)
 assert all((A[i]>>j&1)+(A[j]>>i&1)==1 for i in range(n) for j in range(i+1,n))
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,len(d)):
  for i,c in enumerate(d[S]):
   if c:
    w=A[i]&~S
    while w:
     b=w&-w;w-=b;d[S|b][b.bit_length()-1]+=c
 H=[1]+list(map(sum,d[1:])); F=len(H)-1
 return H[-1],sum(h*H[F^s] for s,h in enumerate(H))

def flip(A,e):
 u,v=e;T=list(A);T[u]^=1<<v;T[v]^=1<<u;return tuple(T)

def prime(p):
 return p>1 and all(p%d for d in range(2,isqrt(p)+1))

def generator(v,p):
 h,q,H,Q=v
 k,a,b=(q,h,H) if q==Q else (h,q,Q)
 if any(x%p==0 for x in (k,a,b)):return False
 z=b*pow(a,-1,p*p)%(p*p)
 if p==3:return z%3==1 and z!=1
 return all(pow(z,(p-1)//d,p)!=1 for d in range(2,p) if (p-1)%d==0 and prime(d)) and pow(z,p-1,p*p)!=1

def inverse(target, even_order=False):
 # Results specify (D,Q) and exponents of (19,217,223,599,singleton).
 if target<=41:return
 for k in range(1,(target//40).bit_length()):
  D=target-40*(1<<k)
  if D<=0:continue
  a=[]; R=D
  for d,w in DYADIC:
   j=0
   while R%d==0:R//=d;j+=1
   a.append(j)
  e=k-sum(a[i]*w for i,(d,w) in enumerate(DYADIC))
  if R==1 and e>=0 and (not even_order or e%2==0):
   yield (D,1<<k,(*a,e))

if __name__=='__main__':
 for (d,w),A in DYADIC.items():assert profile(A)==(d,1<<w)
 for mode,B in BANK.items():
  seen=set()
  for A,e,expected in B:
   assert (*profile(A),*profile(flip(A,e)))==expected
   ps=[p for p in range(3,100,2) if prime(p) and generator(expected,p)]
   seen.update(ps)
   print(mode, A, e, expected, ps)
  assert seen=={p for p in range(3,100,2) if prime(p)}
 F=((198,76,184,81,35,75,148,58),(0,1),(653,8816,659,8816))
 G=((42,12,9,16,39,14),(0,3),(37,352,37,376))
 for A,e,data in (F,G):
  assert (*profile(A),*profile(flip(A,e)))==data and generator(data,5)
 for p in (164057,655583,5243479):
  assert prime(p) and list(inverse(p,True))
  print('exact prime',p,list(inverse(p,True)))
 print('fiber profile, all odd prime clocks under 100, and exact inverse PASS')


from math import gcd
CRT_BANK={'d': [((166, 92, 56, 81, 225, 202, 133, 14), (1, 3), (651, 8808, 657, 8808)), ((90, 20, 249, 18, 160, 75, 154, 43), (1, 4), (369, 5424, 393, 5424)), ((146, 148, 25, 83, 224, 207, 135, 12), (0, 3), (387, 6024, 393, 6024)), ((42, 76, 57, 240, 227, 66, 133, 39), (0, 7), (497, 7232, 503, 7232)), ((106, 92, 105, 112, 229, 66, 128, 47), (3, 7), (383, 5528, 359, 5528)), ((294, 212, 296, 179, 37, 322, 157, 309, 90), (7, 8), (2543, 40432, 2531, 40432))], 'q': [((254, 12, 232, 80, 38, 202, 146, 26), (1, 4), (137, 3272, 137, 3296)), ((70, 44, 136, 17, 39, 77, 158, 59), (2, 4), (469, 6616, 469, 6832)), ((14, 140, 104, 48, 103, 195, 139, 29), (2, 7), (547, 7624, 547, 7648)), ((226, 148, 73, 51, 37, 70, 154, 60), (0, 1), (657, 8808, 657, 8784)), ((2, 76, 41, 17, 231, 75, 141, 47), (4, 6), (281, 4424, 281, 4592))]}
PRIMES=(5,11,13,29,37)
L=[p*(p-1) for p in PRIMES]
Z=[p for p in range(2,max(L)+1) if prime(p) and any(n%p==0 for n in L)]
def full_log(p):
 m=p*p;s=p*(p-1)
 for g in range(2,m):
  if gcd(g,m)==1 and all(pow(g,s//d,m)!=1 for d in Z if s%d==0):break
 else:raise AssertionError
 A={pow(g,e,m):e for e in range(s)}
 assert len(A)==s
 return A
LOG=[full_log(p) for p in PRIMES]
def log_ratio(a,b):
 return [D[b*pow(a,-1,p*p)%(p*p)] for p,D in zip(PRIMES,LOG)]
def ranks(v):
 out=[]
 for q in Z:
  idx=[i for i,m in enumerate(L) if m%q==0]
  r=len(idx); piv={}
  for row in v:
   w=[row[i]%q for i in idx]
   for i in range(r):
    if not w[i]:continue
    if i in piv:
     w=[(a-w[i]*b)%q for a,b in zip(w,piv[i])]
    else:
     v0=pow(w[i],-1,q)
     piv[i]=[(a*v0)%q for a in w];break
  out.append(len(piv))
 return out
GOAL=[sum(m%q==0 for m in L) for q in Z]

if __name__=='__main__':
 for a,Ws in CRT_BANK.items():
  v=[]
  for A,e,d in Ws:
   assert (*profile(A),*profile(flip(A,e)))==d
   assert all(gcd(x,m)==1 for x in d for m in PRIMES)
   d,q,D,Q=d
   assert (q==Q and d!=D) if a=='d' else (d==D and q!=Q)
   v.append(log_ratio(d,D) if a=='d' else log_ratio(q,Q))
  z=ranks(v)
  assert z==GOAL
  print('simultaneous CRT mode',a,'primes',PRIMES,'ranks',dict(zip(Z,z)))
