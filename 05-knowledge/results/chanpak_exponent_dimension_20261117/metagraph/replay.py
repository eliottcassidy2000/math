from itertools import combinations,product
from collections import Counter
from math import gcd,factorial

def check(b, info):
    if not b: raise ValueError(info)

def reflections(v):
    if len(v)%2:
        for x in v:
            for d in reflections([y for y in v if y!=x]):
                yield {**d,x:x}
    elif v:
        x,*tail=v
        for y in tail:
            for d in reflections([z for z in tail if z!=y]):
                yield {**d,x:y,y:x}
    else: yield {}

def cycles(h):
    rem=set(h)
    while rem:
        c=[]; x=min(rem)
        while x in rem:
            c.append(x); rem.remove(x); x=h[x]
        yield c

for n in range(2,8):
    v=list(range(n)); k=n//2; r={x:n-1-x for x in v}; edges=list(combinations(v,2))
    rep=lambda a,b:min(tuple(sorted((a,b))),tuple(sorted((r[a],r[b]))))
    reps=sorted(set(rep(*e) for e in edges)); idx={e:i for i,e in enumerate(reps)}
    # Under standard reversal, bit=1 means increasing, and both members of
    # each r-orbit have the same bit.
    def bit(word,a,b):
        return word>>idx[rep(a,b)]&1
    frames=[w for w in range(1<<len(reps)) if all(bit(w,x,x+1) for x in v[:-1])]
    layers=Counter(); minima=[]
    for s in reflections(v):
        h={x:r[s[x]] for x in v}; cs=list(cycles(h)); ms=list(map(len,cs))
        tests=[(idx[rep(a,b)],idx[rep(s[a],s[b])],int(s[a]<s[b])) for a,b in edges]
        def accepts(w):
            return all(((w>>i)^(w>>j))&1==c for i,j,c in tests)
        plain=sum(accepts(w) for w in range(1<<len(reps)))
        framed=sum(accepts(w) for w in frames)
        if any(z%2==0 for z in ms):
            check(not plain and not framed, (n,'even obstruction')); continue
        e=sum((m-1)//2 for m in ms)+sum(gcd(a,b) for a,b in combinations(ms,2))
        exponent=(e+k)//2
        check(plain==1<<exponent, (n,'unframed',s))
        # Frame transport: collect directed constraints on distinct signed orbits.
        def arc_orbit(a,b):
            out=set()
            for aa,bb in ((a,b),(r[b],r[a])):
                x,y=aa,bb
                while (x,y) not in out:
                    out.add((x,y)); x,y=h[x],h[y]
            return out
        occupied=set(); conflict=False
        for a in v[:-1]:
            orb=arc_orbit(a,a+1)
            if (a+1,a) in occupied: conflict=True
            occupied|=orb
        d=len({min(arc_orbit(a,a+1)) for a in v[:-1]})
        check(framed==(0 if conflict else 1<<(exponent-d)), (n,'frame',s))
        fixed_circle=next((c for c in cs if r[c[0]] in c),None)
        j=0 if fixed_circle is None else (len(fixed_circle)-1)//2
        a=Counter(len(c) for c in cs if c is not fixed_circle)
        check(all(x%2==0 for x in a.values()), (n,'cycle partners'))
        t=tuple(sorted((l,z//2) for l,z in a.items()))
        layers[j,t]+=1
        if h!=dict(zip(v,v)): minima.append(len(reps)-exponent)
    for (j,types),count in layers.items():
        mult=factorial(k)*2**j; den=1
        for l,a in types:
            mult*=2**((l-1)*a); den*=l**a*factorial(a)
        check(count*den==mult, (n,'layers',j,types))
    if minima:
        check(min(minima)==(n-2 if n%2 else 2*n-7), (n,'sharp loss'))
    print(n, len(reps), len(frames), len(layers), min(minima, default=None))
