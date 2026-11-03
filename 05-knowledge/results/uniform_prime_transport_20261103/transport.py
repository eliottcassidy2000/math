"""Actual homogeneous modules, native C-R-C rewrites, and cyclic constructors."""
from itertools import combinations, permutations
from math import comb, factorial, gcd
from fractions import Fraction


def require(ok, text):
    if not ok:
        raise RuntimeError(text)


def subset_counts(rows):
    n=len(rows); size=1<<n
    dp=[[0]*n for _ in range(size)]; hs=[0]*size; hs[0]=1
    for v in range(n):
        dp[1<<v][v]=1
    for s in range(1,size):
        hs[s]=sum(dp[s])
        for v,ways in enumerate(dp[s]):
            rem=rows[v]&~s
            while ways and rem:
                bit=rem&-rem;rem-=bit
                dp[s|bit][bit.bit_length()-1]+=ways
    return hs


def H(rows):
    return subset_counts(rows)[-1]


def covers(rows):
    """Unordered path-cover DP anchored at the least vertex, then order runs."""
    n=len(rows);hs=subset_counts(rows);pc=[[0]*(n+1) for _ in hs];pc[0][0]=1
    for s in range(1,len(hs)):
        pivot=s&-s; part=s
        while part:
            if part&pivot:
                rest=s^part;h=hs[part]
                for k in range(rest.bit_count()+1):
                    pc[s][k+1]+=h*pc[rest][k]
            part=(part-1)&s
    return [0]+[factorial(k)*pc[-1][k] for k in range(1,n+1)]


def permutation_covers(rows):
    n=len(rows);ans=[0]*(n+1)
    for p in permutations(range(n)):
        bad=sum(not(rows[a]>>b&1) for a,b in zip(p,p[1:]))
        for k in range(bad+1,n+1):
            ans[k]+=comb(n-1-bad,k-1-bad)
    return ans


def path_ok(rows,p):
    return sorted(p)==list(range(len(rows))) and all(rows[a]>>b&1 for a,b in zip(p,p[1:]))


def tau(rows,p):
    require(path_ok(rows,p),'illegal selected complement path')
    keep={tuple(sorted(e)) for e in zip(p,p[1:])};out=list(rows)
    for a,b in combinations(range(len(rows)),2):
        if (a,b) not in keep:
            out[a]^=1<<b;out[b]^=1<<a
    return tuple(out)


def rewrite(rows,p,i,j):
    """Physical vertex names stay fixed; only the selected path is re-marked."""
    n=len(rows)
    require(0<i<n-1 and 0<j<n-1 and abs(i-j)>=2,'interior relocation required')
    mid=tau(rows,p);q=list(p);q.insert(j,q.pop(i))
    if not path_ok(mid,q):
        return None
    out=tau(mid,q)
    require([x.bit_count() for x in rows]==[x.bit_count() for x in out],'score invariant')
    old={tuple(sorted(e)) for e in zip(p,p[1:])};new={tuple(sorted(e)) for e in zip(q,q[1:])}
    flipped=[(a,b) for a,b in combinations(range(n),2) if (rows[a]^out[a])>>b&1]
    require(set(flipped)==old^new and len(flipped)==6,'six-edge support')
    return {'rows':out,'path_before':tuple(p),'relocation':(i,j),'path_after':tuple(q),
            'middle':mid,'flipped':flipped}


def native_atlas(rows):
    n=len(rows);hist={};examples={}
    for p in permutations(range(n)):
        if not path_ok(rows,p):
            continue
        for i in range(1,n-1):
            for j in range(1,n-1):
                if abs(i-j)<2:
                    continue
                step=rewrite(rows,p,i,j)
                if step is None:
                    continue
                target=H(step['rows']);hist[target]=hist.get(target,0)+1
                examples.setdefault(target,step)
    return {'source':rows,'H':H(rows),'frame_counts':hist,'examples':examples}


def transitive(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))


def ordinal(*blocks):
    total=sum(map(len,blocks));out=[];off=0
    for b in blocks:
        later=((1<<total)-1)^((1<<(off+len(b)))-1)
        out.extend((r<<off)|later for r in b);off+=len(b)
    return tuple(out)


def cyclic(*blocks):
    require(len(blocks)==3 and all(blocks),'three nonempty cyclic blocks')
    sizes=list(map(len,blocks));offset=[0,sizes[0],sizes[0]+sizes[1]];out=[]
    for i,b in enumerate(blocks):
        j=(i+1)%3;nextmask=((1<<sizes[j])-1)<<offset[j]
        out.extend((r<<offset[i])|nextmask for r in b)
    return tuple(out)


def circulant(n,steps):
    return tuple(sum(1<<((i+s)%n) for s in steps) for i in range(n))


def ear(rows,cut):
    n=len(rows)
    return tuple(r|((not(cut>>i&1))<<n) for i,r in enumerate(rows))+(cut,)


def prime(n):
    if n<2:return False
    if n%2==0:return n==2
    d=3
    while d*d<=n:
        if n%d==0:return False
        d+=2
    return True


def profiles():
    p=circulant(7,[1,2,4]);u=list(p)
    for a,b in [(0,1),(1,3),(0,3)]:u[a]^=1<<b;u[b]^=1<<a
    return {'B11':(18,4,9,19,6),'B13':(10,4,9,18,7),
            'B15':(6,8,18,21,3),'R5':(6,12,24,17,3),
            'U7':tuple(u),'R7':circulant(7,[1,2,3]),'P7':p}


def family_value(c,m,b=2):
    require(m>=1 and b in (1,2),'family parameters')
    if b==1:return (c[1]+c[2])*2**m+c[1]-c[2]
    return 2*(c[2]+c[3])*3**m+(3*c[1]+c[2]-4*c[3])*2**m-c[1]-5*c[2]+2*c[3]


def family_rewrite(source,step,m,b=2):
    n=len(source);rows=cyclic(source,transitive(m),transitive(b))
    p=step['path_before']+tuple(range(n,len(rows)))
    answer=rewrite(rows,p,*step['relocation'])
    require(answer is not None,'contextual rewrite guard')
    expected=cyclic(step['rows'],transitive(m),transitive(b))
    require(answer['rows']==expected,'full exterior relation preserved')
    return answer


def slot_polynomial(c):
    """Integer monomial coefficients of sum_j C_j binom(X,j)."""
    n=len(c)-1;answer=[Fraction(0)]*(n+1);basis=[Fraction(1)]
    for j in range(1,n+1):
        new=[Fraction(0)]*(j+1)
        for i,a in enumerate(basis):
            new[i]-=a*(j-1)/j;new[i+1]+=a/j
        basis=new
        for i,a in enumerate(basis):answer[i]+=a*c[j]
    require(all(a.denominator==1 for a in answer),'OCF integral monomial coefficients')
    return [int(a) for a in answer]


def canonical(rows):
    n=len(rows)
    return min(tuple(sum(((rows[p[i]]>>p[j])&1)<<j for j in range(n))
                     for i in range(n)) for p in permutations(range(n)))


def substitute_root(q,block):
    """Replace quotient vertex0 by the supplied homogeneous block."""
    n=len(block);out=list(block)+[0]*(len(q)-1)
    for j in range(1,len(q)):
        v=n+j-1
        if q[0]>>j&1:
            for i in range(n):out[i]|=1<<v
        else:out[v]|=(1<<n)-1
    for i,j in combinations(range(1,len(q)),2):
        a,b=n+i-1,n+j-1
        if q[i]>>j&1:out[a]|=1<<b
        else:out[b]|=1<<a
    return tuple(out)


def strong(rows):
    n=len(rows)
    for v in range(n):
        seen=1<<v
        while True:
            new=seen
            for i in range(n):
                if seen>>i&1:new|=rows[i]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True
