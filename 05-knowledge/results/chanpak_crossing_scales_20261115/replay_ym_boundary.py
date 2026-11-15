from fractions import Fraction as F
from itertools import product
from random import Random

def gap(w,t):
    a,c=0,0
    for x in w:
        if x=='B':c+=a
        a=t*a+(x=='A')
    return c

def profile(w):
    n=len(w)
    return 1+gap(w,2)+gap(w,0),int(2**(n+1)+2**(n-1)*(gap(w,F(3,2))+gap(w,F(1,2))))

def feature(w,t):
    a=0
    for x in w:a=t*a+(x=='A')
    return a

def star(w):return w[::-1].translate(str.maketrans('AB','BA'))

def check(ok,value):
    if not ok:raise AssertionError(value)

def literal_dp(w):
    n=len(w)
    outgoing=[((1<<n)-(1<<(i+1)))+((1<<n) if b=='B' else 0) for i,b in enumerate(w)]
    outgoing.append(sum(1<<i for i,b in enumerate(w) if b=='A'))
    f=(1<<(n+1))-1
    dp=[[0]*(n+1) for s in range(f+1)]
    h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n+1):
            if not s>>i&1:continue
            t=s^(1<<i)
            dp[s][i]=1 if not t else sum(dp[t][j] for j in range(n+1) if t&outgoing[i]&(1<<j))
        h[s]=sum(dp[s])
    return h[f],sum(h[s]*h[f^s] for s in range(f+1))

def candidate(h,q,n):
    x=h-1-2**(n-2)
    y=q-2**(n+1)-2*3**(n-2)
    return any(2**(n-2-s)<=x<=6*2**(n-2-s)+n-1 and
               2**(s+1)*3**(n-2-s)<=y<=12*2**(s+1)*3**(n-2-s)+n*2**n
               for s in range(1,n//2+1))

def main():
    controls=0
    for n in range(1,8):
        for bits in product('AB',repeat=n):
            w=''.join(bits)
            h,q=profile(w)
            check((h,q)==literal_dp(w),w)
            if n>=3 and w.startswith('A') and w.endswith('B'):
                check(candidate(h,q,n),w)
            controls+=1
    print('literal paths and covers',controls)
    rng=Random(26)
    for _ in range(280):
        u='A'+''.join(rng.choices('AB',k=rng.randrange(1,20)))
        v='A'+''.join(rng.choices('AB',k=rng.randrange(1,20)))
        def h(w):return profile(w)[0]
        def q(w):return gap(w,F(3,2))+gap(w,F(1,2))
        for fn,shift,p,t in ((h,1,2,0),(q,0,F(3,2),F(1,2))):
            def K(a,b):return fn(a+star(b))-fn(a)-fn(b)+shift
            check(K(u,v)==feature(u,p)*feature(v,p)+feature(u,t)*feature(v,t),(u,v))
            check(K(u,u)*K(v,v)-K(u,v)**2==(feature(u,p)*feature(v,t)-feature(v,p)*feature(u,t))**2,(u,v))
        w=u+star(v);n=len(w);r=n//2
        for t in (F(2),F(3,2)):
            L=sum(t**-i for i in range(r) if w[i]=='A')
            R=sum(t**-i for i in range(r) if w[n-1-i]=='B')
            err=gap(w,t)/(t**(n-2))-L*R
            check(0<=err<=2*t**-r/(1-1/t)**2,w)
        check(candidate(*profile(w),n),w)
    print('exact-rational seam and boundary instances',280)
    for r in range(14,30):
        w='A'+'B'*r+'A'*r+'B';n=len(w)
        h,trueq=profile(w)
        q=6*3**(n-2)
        q+=next(d for d in range(24) if (q+d)%8==0 and (q+d-(-1)**(n+1)*h)%3==0)
        check(q>=(n+2)*h and q>=2**(n+1)+2*3**(n-2)+2,(r,q))
        check(candidate(h,trueq,n) and not candidate(h,q,n),(r,h,q,trueq))
    print('false cross-scale pairs eliminated',16)

if __name__=='__main__':main()
