from functools import lru_cache
from itertools import product

def dp(r,root=None):
 n=len(r);N=1<<n
 d=[[0]*n for _ in range(N)]
 for j in range(n) if root is None else (root,):d[1<<j][j]=1
 for i,row in enumerate(d):
  for j,v in enumerate(row):
   if not v:continue
   z=r[j]&~i
   while z:
    b=z&-z;z^=b
    d[i|b][b.bit_length()-1]+=v
 return [1]+[sum(t) for t in d[1:]]

def prof(t,x):
 d=dp(t);z=dp(t,x);N=len(d)-1
 ev,od=(sum(z[s]*d[N^s] for s in range(1,N+1) if s.bit_count()%2==k) for k in (0,1))
 return d[-1],d[N^(1<<x)],ev,od

def tail(t,x,w):
 t=list(t)
 for b in map(int,w):
  n=len(t);m=0 if b else 1<<x
  t=[a|((not(m>>j&1))<<n) for j,a in enumerate(t)]+[m]
 return tuple(t)

def rc(w):return w[::-1].translate(str.maketrans('01','10'))

@lru_cache(None)
def f(z):
 if not z:return 1
 u,p=divmod(z,2)
 return f(u)+(0 if p else u+(u&1))

def response(profile,w):
 assert w and w[-1]=='0'
 h,g,r,s=profile
 z=int(rc(w),2)
 x=int(w,2)
 b=s-g
 hr=h+g*(f(z)-1)+(r+b)*z+(r-b)*(w[0]=='0')
 rr=((1<<len(w))*(r+b)+g*x)//2
 return hr,g,rr,g+rr



D=(
    ((0,1,27,51,99,135,47,95),164,(723,37,378,415),0),
    ((2,72,131,53,39,71,157,59),22,(1395,287,799,1000),-1),
    ((0,33,11,19,39,77,159,63),9,(111,49,119,119),2),
)
def init(rows,q):
    n=len(rows)
    return tuple(r|(((q>>i&1)^1)<<n) for i,r in enumerate(rows))+(q,)

if __name__=='__main__':
    seeds=[]
    for rows,m,expected,t in D:
        T=init(rows,m)
        assert prof(T,8)==expected
        assert prof(T,8)[0]==dp(T)[-1]
        seeds.append((T,expected,t))
    for s in range(2,5):
        E=378*(1<<s)
        vals=[]
        for i,(T,expected,t) in enumerate(seeds):
            A=tail(T,8,'0'*(s+t))
            h,g,e,o=prof(A,8)
            assert (h,g,e,o)==(756*(1<<s)-33+2*i,expected[1],E,E+expected[1])
            for j in range(len(A)):
                R=1<<j
                while True:
                    S=R
                    for k in range(len(A)):
                        if R>>k&1:S|=A[k]
                    if R==S: break
                    R=S
                assert R==(1<<len(A))-1
            for w in ('0','01'+'0','10','100'):
                assert prof(tail(A,8,w),8)==response((h,g,e,o),w)
            vals.append(h)
        print(s,vals,E)

    from math import isqrt
    def block(T,root,B):
        n=len(T);m=len(B)
        new=((1<<m)-1)<<n
        return tuple(t|(new if j!=root else 0) for j,t in enumerate(T))+tuple((b<<n)|(1<<root) for b in B)
    for B,g,v,s,result in [
        ((2,12,9,49,7,23),25,304,6,14707961),
        ((2,12,17,37,11,87,31),49,824,5,19932689),
    ]:
        R=dp(B);F=len(R)-1
        assert (R[F],sum(R[j]*R[F^j] for j in range(F+1)))==(g,v)
        d=seeds[1]
        first=tail(d[0],8,'0') # s=2 middle branch
        assert prof(block(first,8,B),8)==(756*4*v-31*g,287*g,378*4*v,378*4*v+287*g)
        assert result==756*(1<<s)*v -31*g
        prime=result if result%7 else result//7
        assert all(prime%p for p in range(2,isqrt(prime)+1))
        print('composite block',g,'target',result,'prime part',prime)
