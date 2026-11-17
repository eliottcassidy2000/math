from collections import Counter,defaultdict
from functools import cache
def G(rows):
 n=len(rows);f=(1<<n)-1;B=defaultdict(list)
 for u in range(n):
  F=f^((1<<(u+1))-1); b=1<<u
  @cache
  def d(s,v):
   z=s^(1<<v)
   if not z:return rows[u]>>v&1
   return sum(d(z,w) for w in range(u+1,n) if z>>w&1 and rows[w]>>v&1)
  k=F
  while k:
   if k.bit_count()>=2 and k.bit_count()%2==0:
    c=sum(d(k,w) for w in range(u+1,n) if k>>w&1 and rows[w]>>u&1)
    if c:B[b].append((k|b,c))
   k=(k-1)&F
 @cache
 def fcn(S):
  if not S:return {(0,0):1}
  b=S&-S;P=Counter({(c,s+1):v for (c,s),v in fcn(S^b).items()})
  for u,k in B[b]:
   if S&u==u:
    for (c,s),v in fcn(S^u).items():P[c+1,s]+=k*v
  return P
 return fcn(f)
from itertools import product
from math import gcd,isqrt,comb,prod
import ast
from pathlib import Path

def need(x):
    if not x:raise ArithmeticError('exact replay failed')

def ev(g,x,y):return sum(a*x**i*y**j for (i,j),a in g.items())

def conv_rows(T,k):
    """Independent of cycle packing. Count path covers by subset convolution."""
    n=len(T);f=(1<<n)-1
    need(all(0<=r<=f and not r>>i&1 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(n) for j in range(i)))
    d=[[0]*n for _ in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for S in range(1,f+1):
        for j,v in enumerate(d[S]):
            if not v:continue
            s=T[j]&~S
            while s:
                b=s&-s;s^=b;d[S|b][b.bit_length()-1]+=v
    h=[1]+[sum(x) for x in d[1:]]
    @cache
    def P(S,q):
        if q==1:return h[S]
        b=S;out=0
        while True:
            out+=h[b]*P(S^b,q-1)
            if not b:return out
            b=(b-1)&S
    return tuple(P(f,q) for q in range(1,k+1))

def convert(g):
    D=Counter()
    for (a,b),v in g.items():
        for k in range(a+1):D[k,a+b-k]+=v*comb(a,k)*2**(a-k)
    return D

def mul(A,B):
    D=Counter()
    for (a,b),x in A.items():
        for (c,d),y in B.items():
            if x and y:D[a+c,b+d]+=x*y
    return D

def shift(A,t):
    D=Counter()
    for (a,b),x in A.items():
        if x:
            for j in range(b+1):D[a,j]+=comb(b,j)*x*t**(b-j)
    return D

def sub(A,B):
    D=Counter(A)
    D.subtract(B)
    return Counter({m:v for m,v in D.items() if v})

def L(A,B):
    D=mul(sub(shift(A,1),A),sub(shift(B,1),B))
    D.update(mul(sub(shift(A,-1),A),sub(shift(B,-1),B)))
    need(all(v%2==0 for v in D.values()))
    out=mul(A,B)
    out.update(D)
    ans=Counter({(a,b+1):v for (a,b),v in out.items() if v})
    ans.update({(a+1,b):v//2 for (a,b),v in D.items() if v})
    return ans

def lift(A,B,C,D):
    return mul(L(A,C),L(B,D)),mul(L(A,D),L(B,C))

def erase(T,S):
    I=[i for i in range(len(T)) if i not in S]
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(I))for i in I)

def prime(p):
    return p==2 or p>2 and p%2==1 and all(p%d for d in range(3,isqrt(p)+1,2))

def gap(p):
    N=7*p; s=N.bit_length()-5
    return (18<<s)<=N<=(28<<s)

def factors(v):
    d=2
    while d*d<=v:
        c=0
        while v%d==0:c+=1;v//=d
        if c:yield d,c
        d=3 if d==2 else d+2
    if v>1:yield v,1

def primitive():
    S=Path(__file__).with_name('primitive.txt').read_text().splitlines()
    V=[[Counter() for _ in range(3)] for _ in (-1,1)]
    Q=[1,1];orders=[0,0];mod=1000000007
    for l in S:
        d,w=l.split('[',1)
        a,*r=map(int,d.split());T=ast.literal_eval('['+w);b=int(a>0);e=abs(a)
        P=conv_rows(T,4)
        need(P[:3]==tuple(r))
        Q[b]=Q[b]*pow(P[3],e,mod)%mod
        orders[b]+=e*len(T)
        for j,p in enumerate(r):
            for k,x in factors(p):V[b][j][k]+=e*x
    need(orders==[270948,271000])
    need(all(v1==v0+Counter({7:1}) for v0,v1 in zip(*V)))
    need(not (V[0][0].keys() & V[0][1].keys()))
    need((Q[1]-7*Q[0])%mod!=0)
    u,R,S=(prod(pow(p,e) for p,e in v.items()) for v in V[0])
    # Common singleton cyclic wrap has H=u+R and P2=4R+2(S-u).
    need(gcd(u+R,S-3*u)==1)
    print('primitive 3-colour scale and primitive first wrap:',orders,'verified')

def main():
    A=(38,8,2,5,15,30);B=(58,4,1,6,14,30)
    C=(260,9,50,53,195,83,15,111,254);D=(36,65,266,131,15,94,157,55,251)
    a,b,c,d=map(lambda t:convert(G(t)),(A,B,C,D))
    need(sub(a,b)=={(1,1):2,(2,0):1})
    need(sub(c,d)=={(2,1):2,(3,0):1})
    for L1,L2,r,N in ((a,b,1,6),(c,d,2,9)):
        need(max(i+j for i,j in L2)==N and min(i for i,j in sub(L1,L2))==r)
    u,v=lift(c,d,c,d)
    need(max(i+j for i,j in v)==38 and min(k for k,d in sub(u,v))==4)
    alpha,beta=3390151449,-3390038347
    for X,Y,remove,profile,answer,r in (
      (A,B,{1},((5,64),(9,84)),(416989193737,539035098309),1),
      (C,D,{1,2},((117,1332),(29,640)),(8634728973537,4241082742657),2)):
        uv=tuple(conv_rows(erase(T,remove),2) for T in (X,Y))
        need(uv==profile)
        need(tuple(alpha*(R*2)+beta*u for u,R in uv)==answer)
        need(tuple((convert(G(T))[r,1]>>1)&1 for T in (X,Y))==(1,0))
        ps=tuple(prime(p) for p in answer)
        need(sum(ps)==1 and all(p%3==0 for p,b in zip(answer,ps) if not b))
        need(gap(answer[ps.index(True)]))
        print('depth',r,'prime quotient',answer[ps.index(True)],'composite mate',answer[ps.index(False)])
    v=alpha*224+beta*17
    need(conv_rows(A,2)==conv_rows(B,2)==(17,224) and v==701763272677 and prime(v))
    primitive()
    print('exact depth, single-bit controls, prime and primitive receipts: PASS')

if __name__=='__main__':main()
