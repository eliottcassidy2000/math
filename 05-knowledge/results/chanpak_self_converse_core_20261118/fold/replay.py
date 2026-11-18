from itertools import combinations, permutations
from functools import cache

def check(b,info):
    if not b: raise ValueError(info)
for n in range(2,7):
    V=tuple(range(n)); m=n//2; c=m; odd=n%2
    r=lambda i:n-1-i
    pairs=tuple((i,r(i)) for i in range(m))
    pair=lambda v:min(v,r(v))
    rep=lambda a,b:min(tuple(sorted((a,b))),tuple(sorted((r(a),r(b)))))
    reps={e:i for i,e in enumerate(sorted({rep(*e) for e in combinations(V,2)}))}
    sets=[tuple(i for i in range(m) if mask>>i&1) for mask in range(1<<m)]
    for w in range(1<<len(reps)):
        def arc(a,b):
            return (a<b)==bool(w>>reps[rep(a,b)]&1)
        @cache
        def signed(S,v,start=-1):
            if len(S)==1: return int(start==-1 or arc(start,v))
            B=tuple(i for i in S if i!=pair(v))
            return sum(signed(B,u,start) for i in B for u in pairs[i] if arc(u,v))
        fixed=sum(signed(sets[-1],v) for v in V if (not odd or v!=c)
                  and arc(v,c if odd else r(v)))
        h=b=0
        for p in permutations(V):
            if all(arc(u,v) for u,v in zip(p,p[1:])):
                h+=1; b+=p==tuple(map(r,reversed(p)))
        check(fixed==b,(n,w,'fixed path'))
        @cache
        def E(S,a):
            if len(S)==1:return 1
            B=tuple(v for v in S if v!=a)
            return sum(E(B,b) for b in B if arc(b,a))
        middle=0
        for x in V if odd else (-1,):
            U=tuple(v for v in V if v!=x)
            for S in combinations(U,m):
                C=tuple(v for v in U if v not in S); R=tuple(sorted(map(r,C)))
                middle+=sum(E(S,a)*E(R,r(b)) for a in S for b in C
                            if (arc(a,x) and arc(x,b) if odd else arc(a,b)))
        check(h==middle,(n,w,'full midpoint'))
        cycles=[]
        for l in range(3,n+1,2):
            for C in combinations(V,l):
                a,*tail=C
                cycles.extend((a,)+p for p in permutations(tail)
                              if all(arc(u,v) for u,v in zip((a,)+p,p+(a,))))
        def mirror(C):
            S=tuple(map(r,reversed(C))); at=S.index(min(S)); return S[at:]+S[:at]
        q={S:sum(len(C)==len(S) and tuple(sorted(map(pair,C)))==S and min(S) in C
                 for C in cycles if not odd or c not in C) for S in sets if S}
        @cache
        def L(S):
            if not S:return (1,0,0)
            out=list(L(S[1:]))
            for A in sets:
                if A and q[A] and A[0]==S[0] and set(A)<=set(S):
                    out[2]+=q[A]*L(tuple(i for i in S if i not in A))[0]
            return tuple(out)
        A=sum(signed(S,v,c) for S in sets if S for i in S for v in pairs[i] if arc(v,r(v))) if odd else 0
        F=(1,A,L(sets[-1])[2])
        direct=(1,sum(C==mirror(C) for C in cycles),
                   sum(not set(C)&set(D) and {C,D}=={mirror(C),mirror(D)}
                       for C,D in combinations(cycles,2)))
        check(F==direct,(n,w,'fixed packing'))
        check((h-1-2*A-4*((len(cycles)-A)//2+F[2]))%8==0,(n,w,'mod 8'))
        if not odd:
            overlapping=sum(bool(set(C)&set(mirror(C))) for C in cycles)//2
            check((h-1-4*overlapping)%8==0,(n,w,'mod 8 overlap'))
    print(n,1<<len(reps))

A,B=(18,4,9,19,6),(2,12,9,17,7)
for rows,expected in ((A,11),(B,13)):
    v=tuple(range(5));mirror=lambda i:4-i
    arc2=lambda i,j:rows[i]>>j&1
    words=list(permutations(v))
    check(all(arc2(i,j)==arc2(mirror(j),mirror(i)) for i in v for j in v),'anti')
    H=[p for p in words if all(arc2(i,j) for i,j in zip(p,p[1:]))]
    check(len(H)==expected and sum(p==tuple(map(mirror,reversed(p))) for p in H)==3,'local hostile paths')
    check(sum(all(arc2(i,j)==arc2(p[i],p[j]) for i in v for j in v) for p in words)==1,'local hostile rigid')
    check(sorted(x.bit_count() for x in rows)==[1,2,2,2,3],'local hostile score')
    check(all(any(arc2(i,j) for i in v[k:] for j in v[:k]) for k in range(1,5)),'local hostile strong')
print('same-order rigid blue-degree and score collision: 11 versus 13')
