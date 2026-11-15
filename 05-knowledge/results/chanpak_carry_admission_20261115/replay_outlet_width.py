from itertools import permutations
from math import factorial, comb
from random import Random

def pick_rows(n,rng):
    r=[0]*n
    for u in range(n):
        for v in range(u+1,n):
            if rng.randrange(2):r[u]|=1<<v
            else:r[v]|=1<<u
    return tuple(r)

def tab(r,kmax):
    n=len(r);f=(1<<n)-1
    dp=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if not (s>>i&1):continue
            t=s^(1<<i)
            dp[s][i]=sum(dp[t][j] for j in range(n) if (t&r[i])>>j&1) if t else 1
        h[s]=sum(dp[s])
    q=[[1]+[0]*f,h]
    for _ in range(2,kmax+1):
        a=[]
        p=q[-1]
        for m in range(f+1):
            v=0;s=m
            while True:
                v+=h[s]*p[m^s]
                if not s:break
                s=(s-1)&m
            a.append(v)
        q.append(a)
    return q[:kmax+1]

def wrap(r,u):
    # Mark vertex zero of r as the one to be replaced, and append U.
    n=len(r)-1
    m=len(u)
    outs=r[0]>>1
    return tuple((v>>1)|((((1<<m)-1)<<n) if v&1 else 0) for v in r[1:])+tuple((v<<n)|outs for v in u)

def gate(u):
    return wrap((2,4,1),u)

def candidates(r,j):
    # Distinct numbered placeholders are kept sorted through path slots.
    n=len(r)-1
    d=r[0].bit_count()
    ans=[0]*(d+j+1)
    ext=tuple(range(1,n+1))
    for k in range(1,d+j+1):
        # List k runs in their mandatory order and then independently
        # insert j-1 indistinguishable separators; this is small-host only.
        from itertools import combinations
        ts=n+k+j-1
        for z in combinations(range(ts),j-1):
            zp=set(z)
            rem=[i for i in range(ts) if i not in zp]
            for w in combinations(rem,k):
                wp=set(w)
                for ex in permutations(ext):
                    e=iter(ex)
                    chars=[-1 if i in zp else 0 if i in wp else next(e) for i in range(ts)]
                    ok=True
                    for a,b in zip(chars,chars[1:]):
                        if a!=-1 and b!=-1 and not ((r[a]>>b)&1):
                            ok=False
                            break
                    if ok:ans[k]+=1
    return ans

def locals_for_outlet(r,j):
    n=len(r)
    assert r[0].bit_count()==1
    q=tab(tuple(v>>1 for v in r[1:]),j)[j]
    A,E,O=q[-1],0,0
    for k in range(2,n):
        for v in permutations(range(1,n),k):
            seq=(0,)+v+(0,)
            if not all((r[a]>>b)&1 for a,b in zip(seq,seq[1:])):continue
            c=((1<<(n-1))-1)^sum(1<<(a-1) for a in v)
            if k%2:E+=q[c]
            else:O+=q[c]
    return A,E,O

def main():
    u=(92,141,32,228,206,19,38,101)
    v=(52,229,208,231,10,148,49,81)
    for w,Q3,G2 in [(u,45225,118305),(v,44505,116865)]:
        qq=[arr[-1] for arr in tab(w,3)][1:]
        hd=tab(gate(w),1)[1][-1]
        hdd=tab(gate(gate(w)),1)[1][-1]
        if qq!=[345,5640,Q3] or hd!=5985 or hdd!=G2:
            raise AssertionError((qq,hd,hdd))
        print('core and nested gate readouts:',qq,hd,hdd)
    def apex(w,mask):
        return tuple(x|((1<<len(w)) if not mask>>i&1 else 0) for i,x in enumerate(w))+(mask,)
    z=tab(apex(u,49),1)[1][-1],tab(apex(v,152),1)[1][-1]
    if z!=(1939,1941):raise AssertionError(z)
    print('selective width-one readouts:',z)
    rng=Random(4829)
    for r in range(160):
        host=pick_rows(rng.randrange(3,6),rng)
        # Force exactly one outlet at random.
        qn=rng.randrange(1,len(host))
        host=(1<<qn,)+tuple((x&~1)|int(i!=qn) for i,x in enumerate(host[1:],1))
        u=pick_rows(rng.randrange(1,5),rng)
        j=rng.randrange(1,4)
        A,E,O=locals_for_outlet(host,j)
        qi=[a[-1] for a in tab(u,j+1)]
        theory=(A-2*j*E)*qi[j]+j*(E+O)*qi[j+1]+j*(E-O)*qi[j-1]
        actual=tab(wrap(host,u),j)[j][-1]
        if theory!=actual:raise AssertionError((host,u,j,theory,actual))
    print('160 literal all-slot one-outlet checks with independently enumerated exterior cycles')
    for _ in range(72):
        r=pick_rows(rng.randrange(1,5),rng)
        u=pick_rows(rng.randrange(1,5),rng)
        j=rng.randrange(1,3)
        a=candidates(r,j)
        q=[x[-1] for x in tab(u,len(a)-1)]
        L=[sum((-1)**(k-i)*comb(k,i)*q[i] for i in range(k+1)) for k in range(len(a))]
        got=sum(x*y for x,y in zip(a,L))
        expected=tab(wrap(r,u),j)[j][-1]
        if got!=expected:raise AssertionError((r,u,j,a,L,got,expected))
        d=r[0].bit_count()
        if j==1 and a[-1]%factorial(d):raise AssertionError((r,a))
    print('72 literal skeleton-transfer and top-coefficient controls')
if __name__=='__main__':main()
