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
