"""Exact keep-root word inverse. Import canonical literal constructions for controls."""
from itertools import product

def step(S,b):
    h,e,t,d=S
    return (h+t-e,d,2*t-d,d) if b else (h,0,2*t,d)

def run(S,bs):
    for b in bs: S=step(S,b)
    return S

def inverse(S,d):
    """Only after the first bit: old e is d times old t parity."""
    H,E,T,D=S
    assert D==d and d%2==1
    b=T%2
    assert E==d*b
    t=(T+b*d)//2
    e=d*(t%2)
    return (H-b*(t-e),e,t,d),b

def bits_lsb(w,q):
    return tuple(w>>i&1 for i in range(q))

def fixed_depth(S,N,k,trace=None):
    p=(k+1)//2
    q=k-p
    assert S[2]>=2*S[3]
    for prefix in product((0,1),repeat=p):
        h,e,t,d=R=run(S,prefix)
        # Here all candidate suffixes obey h+tW-t < N <= h+tW,
        # so the division cannot drop a second W on the same prefix.
        w=-((h-N)//t)
        if not 0<=w<1<<q: continue
        suf=bits_lsb(w,q)
        if not (prefix+suf)[-1]: continue
        actual=run(R,suf)[0]
        if trace is not None: trace.append((prefix+suf,actual-N))
        if actual==N: yield prefix+suf

def words_for(S,N):
    if N==S[0]: yield ()
    if N<S[0]: return
    # Nonempty minimal word ends in 1. Then t_k<=2N-d.
    k=1
    while (S[2]-S[3])*(1<<k)+S[3]<=2*N-S[3]:
        yield from fixed_depth(S,N,k)
        k+=1

def main():
    S=(29,11,31,13)
    collisions={}
    for k in range(1,13):
        bank={}
        for w in product((0,1),repeat=k):
            R=run(S,w)
            for i in range(k,1,-1):
                R,b=inverse(R,S[3])
                assert b==w[i-1] and R==run(S,w[:i-1])
            if not w[-1]: continue
            bank.setdefault(run(S,w)[0],[]).append(w)
        collision=next(((h,v) for h,v in bank.items() if len(v)>1),None)
        if collision: collisions[k]=collision
        # For selected depths, check all positive witnesses and selected nonmembers.
        if k<=8 or k==12:
            tests=list(bank) if k<=8 else sorted(bank)[::51]
            tests+=list(range(31,450,14))
            for n in tests:
                assert set(fixed_depth(S,n,k))==set(bank.get(n,()))
        print('depth',k,'ending-one words',1<<(k-1),'distinct',len(bank),
              'first collision',collision)
    for n in (30011,30013,33749,33751,10**9+7,10**9+9):
        for A in (S,run(S,(0,0,0,1))):
            W=list(words_for(A,n))
            print('seed',A,'target',n,'words',W)
            assert all(run(A,w)[0]==n for w in W)

if __name__=='__main__': main()
