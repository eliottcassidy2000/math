from verify import dp,prof,append,reorder
from functools import lru_cache

@lru_cache(None)
def f(z):
 if not z:return 1
 u,p=divmod(z,2)
 return f(u)+(0 if p else u+(u&1))
@lru_cache(None)
def inverse(N,a,b,c):
 if N<b:return ()
 hits=[0] if N==b else []
 if N>=2*a+b and (N-b)%2==0:
  hits += [2*u for u in inverse(N,2*a+b,b,b) if u]
 if N>=a+b+c and (N-a-c-b)%2==0:
  hits += [2*u+1 for u in inverse(N-a-c,2*a,b,0)]
 return tuple(sorted(hits))
def deck(H):
 out=[]
 for k in range(H.bit_length()):
  d=1<<k
  for x in inverse(H-6*d-2,8*d+1,6*d-1,2*d-1):out.append((k,x))
 return out
def terminal(E,y):
 N=(1<<len(E))-1
 rev=tuple(N^r^(1<<i) for i,r in enumerate(E))
 h=dp(E);B=[0,0]
 for u in range(len(E)):
  if not(E[u]>>y&1):continue
  C=dp(rev,u)
  for s in range(1,N+1):
   if not(s>>y&1):B[s.bit_count()&1]+=C[s]*h[(N^(1<<y))^s]
 return h[N],h[N^(1<<y)],*B
def graft(E,y,X):
 w=format(X,'b') if X else ''
 n=len(w);m=len(E)
 return tuple((((1<<n)-1)^((1<<(i+1))-1)) |
    ((((1<<m)-1)<<n) & ~(int(w[i])<<(n+y))) for i in range(n))+tuple(
    (v<<n)|(sum((b=='1')<<i for i,b in enumerate(w)) if j==y else 0)
    for j,v in enumerate(E))
def main():
 E=(6,20,8,3,13);x,y=3,4
 assert terminal(reorder(E,(0,1,2,4)),3)==(3,1,1,1)
 T=E
 for target in ((13,5,4,5),(25,11,7,10),(49,23,13,20)):
  assert terminal(T,y)==target
  T=append(T,x)
 for X in range(41):
  G=graft(E,y,X);p=X%2;fx=f(X)
  assert prof(G,len(G)-2)==(5*fx+9*X+8+p,fx+2*X+2,3*fx+4*X+3+p,4*fx+6*X+5+p)
 for H,(k,X) in ((13,(0,0)),(41,(0,2)),(125,(0,8)),(10009,(2,214)),
                (30011,(0,2673)),(30089,(0,2625)),(30091,(0,2332))):
  assert (k,X) in deck(H)
  G=graft(E,y,X);r=len(G)-2
  for _ in range(k):G=append(G,r)
  assert dp(G)[-1]==H
  print(H,deck(H),'checked at order',len(G))
 for H in (10007,30013,30021,30027):assert not deck(H)
if __name__=='__main__':main()
