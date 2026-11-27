from collections import Counter
from itertools import product
from math import comb, gcd, isqrt, lcm


def F(x):
    return (x**3-x)//6


def H(x):
    h=1
    while x:
        x,b=divmod(x,2)
        if not b:
            h+=x+(x%2)
    return h


def lower(u,w):
    return 2*u*w+u+(u%2)+H(u)+H(w)-1


def C(n,k):
    return comb(n,k) if 0<=k<=n else 0


def profile(u,w):
    return tuple(sum(((w>>i&1)-(u>>i&1))*C(i,j)
                 for i in range(max(u,w).bit_length())) for j in range(0,24,2))


def fast_profile(u,w):
    mask=u^w
    ans=[0]*12
    while mask:
        lo=mask & -mask
        i=lo.bit_length()-1
        sign=1 if w&lo else -1
        for j in range(12):
            ans[j]+=sign*C(i,2*j)
        mask^=lo
    return tuple(ans)


def expected(r):
    return ((1-r,)+tuple(C(r+1,j)-C(r+1,j+1) for j in range(2,24,2))
            if r else (1,)+(0,)*11)


def setup(r,Q):
    a=(Q&-Q).bit_length()-1
    m=Q>>a
    ell=max(a,1)
    while pow(2,ell,m)!=1%m:
        ell+=1
    K=r+5
    choices=(4,2) if r==0 else [(2*z+1)<<r for z in range(4)]
    def flank(p,s0):
        s=(p<<K)+s0
        return 2*s-(r!=0),2*s+1
    p=(1<<ell)-1
    bases={lower(*flank(p,s0))%4:(s0,lower(*flank(p,s0))%Q) for s0 in choices}
    assert set(bases)=={1,3}
    D=((1<<ell)-1)*p+1
    step=2*D%Q
    assert gcd(step,Q)==4
    def address(phase,extra=0,short=False):
        s0,h=bases[phase%4]
        b=Q//4
        j=((phase-h)//4)*pow(step//4,-1,b)%b if b>1 else 0
        if short:
            return (h+step*j)%Q
        p=(1<<(ell+extra))-1
        # Extra prefix width must be a multiple of the period.
        assert pow(2,extra,m)==1%m
        for _ in range(j):
            p=(p<<(2*ell))+((1<<ell)-1)
        rest=b-1-j
        p=(p<<(2*ell*rest))+((1<<(2*ell*rest))-1)
        return flank(p,s0)
    return address


def pollock_check():
    primes=(2,3,5,7,11,13,17,19,23,31,37)
    runs=0
    for p in primes:
        for e in (1,2,3):
            q=p**e
            if p<5:
                assert len({F(6*t)%q for t in range(q)})==q
                runs+=q
                continue
            A={F(t)%p:t for t in range(p) if (3*t*t-1)%p}
            chi=1 if pow(-3,(p-1)//2,p)==1 else -1
            assert len(A)==(2*p+chi)//3
            for n in range(q):
                x,y=next((x,A[(n-a)%p]) for a,x in A.items() if (n-a)%p in A)
                for i in range(1,e):
                    step=p**i
                    x=next(x+j*step for j in range(p) if (F(x+j*step)+F(y)-n)%(step*p)==0)
                assert (F(x)+F(y)-n)%q==0
                runs+=1
    print('Pollock prime-power targets:',runs)


def scan(K):
    # All k in [K,2K) and the entire canonical shell remainder range.
    # For each qualifying pair use the smallest possible common 3-term sum, zero.
    lo=K*(K+1)//2
    hi=(2*K-1)*(2*K)//2
    cap=(2*K+1)**2
    fs=[]
    while F(len(fs)+1)<=cap:
        fs.append(F(len(fs)+1))
    found=Counter()
    for i,b in enumerate(fs):
        for a in fs[:i]:
            d=b-a
            if d>hi:
                continue
            if d<lo:
                break
            z=isqrt(8*d+1)
            if z*z==8*d+1:
                k=(z-1)//2
                if a<(k+1)*(k+2)//2:
                    assert b < (k+1)**2
                    found[k]+=1
    return len(found),sum(found.values())


def inverse(k):
    A=3*k*(k+1)
    ds=set()
    for e in range(2,isqrt(A)+1):
        if A%e==0:
            ds.update((e,A//e))
    ds.add(A)
    hits=[]
    for d in sorted(ds):
        q,rem=divmod(4*(A//d)-d*d+4,3)
        if rem or q<0:
            continue
        m=isqrt(q)
        if m*m!=q or m<=d or (m-d)%2:
            continue
        u,v=(m-d)//2,(m+d)//2
        if F(u)<(k+1)*(k+2)//2:
            assert F(v)-F(u)==k*(k+1)//2
            hits.append((u,v))
    return hits


def exact_small():
    def dp(word):
        n=len(word)+1
        out=[0]*n
        for i in range(n-1):
            for j in range(i+1,n-1):
                out[i]|=1<<j
            if word[i]=='1':
                out[-1]|=1<<i
            else:
                out[i]|=1<<(n-1)
        h=[[0]*n for _ in range(1<<n)]
        for i in range(n):
            h[1<<i][i]=1
        for mask,row in enumerate(h):
            for i,num in enumerate(row):
                if not num:
                    continue
                for j in range(n):
                    if out[i] & (1<<j) and not mask & (1<<j):
                        h[mask|(1<<j)][j]+=num
        return sum(h[-1])
    def cycles(word):
        result=[0]*12
        for mask in range(1,1<<len(word)):
            size=mask.bit_count()
            if size>=2 and size%2==0:
                first=(mask&-mask).bit_length()-1
                last=mask.bit_length()-1
                if word[first]=='1' and word[last]=='0':
                    result[size//2-1]+=1
        return result
    cases=0
    for s in range(1,17):
        for i,j in product((0,1),repeat=2):
            u,w=2*s-i,2*s+2-j
            right=''.join('1' if c=='0' else '0' for c in bin(w)[2:][::-1])
            left=bin(u)[2:]
            before,after=left+'0'+right,left+'1'+right
            assert dp(before)==lower(u,w)
            assert dp(after)==lower(u,w)+2
            assert tuple(b-a for a,b in zip(cycles(before),cycles(after)))==profile(u,w)
            cases+=1
    print('Independent Hamilton-path DP and odd-cycle enumerations:',cases)

if __name__=='__main__':
    pollock_check()
    for K in (64,256,1024,4096,16384):
        count,pairs=scan(K)
        if K<=1024:
            chk=[inverse(k) for k in range(K,2*K)]
            assert (count,pairs)==(sum(bool(p) for p in chk),sum(map(len,chk)))
        print(f'Pollock shell [{K},{2*K}): eligible shells={count}; pairs={pairs}')
    total=0
    for s in range(1,4097):
        for i,j in product((0,1),repeat=2):
            u,w=2*s-i,2*s+2-j
            got=profile(u,w)
            assert got[0]<=1
            assert got==expected(1-got[0])
            if got[0]==0:
                assert got[1:]==(1,)+(0,)*10
            total+=1
    print('Exact +2 odd-cycle profiles:',total)
    phases=0
    for M in (*range(1,65),77,121,180,315):
        Q=lcm(4,M)
        for r in range(0,9):
            addresses=setup(r,Q)
            for phase in range(1,Q,2):
                assert addresses(phase,short=True)==phase
                phases+=1
            for phase in {1,3%Q,Q-1}:
                u,w=addresses(phase)
                assert lower(u,w)%Q==phase and fast_profile(u,w)==expected(r)
    print('Congruence and cycle-profile combinations:',phases)
    exact_small()
