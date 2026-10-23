#!/usr/bin/env python3
# Load definitions without executing verifier's main loop.
exec(open(__file__.replace('spectrum.py','verify.py')).read().split('rows=[]')[0])
def graded(n,p):
    remaining=set(combinations(range(n),2));poly=[1]
    while remaining:
        start=next(iter(remaining));cur=start;flip=0;length=0
        while True:
            remaining.remove(cur);a,b=map(p.__getitem__,cur)
            flip^=(a>b);cur=tuple(sorted((a,b)));length+=1
            if cur==start:break
        out=poly+[0]*length
        for k,a in enumerate(poly):out[k+length]+=(-1)**flip*a
        poly=out
    return poly
out=[]
for n in range(1,9):
    pp=parts(n);N=n*(n-1)//2;P={mu:graded(n,perm(mu)) for mu in pp}
    for mu in pp:
        assert sum(P[mu])==fixed(n,perm(mu),False)
        assert sum((-1)**k*a for k,a in enumerate(P[mu]))==fixed(n,perm(mu),True)
    mult={la:[sum(Fraction(chi(la,mu)*P[mu][k],z(mu)) for mu in pp) for k in range(N+1)] for la in pp}
    for la in pp:
        assert all(a.denominator==1 and a>=0 for a in mult[la])
        assert mult[la]==list(reversed(mult[transpose(la)]))
    out.append({'n':n,'N':N,'walsh_degree_multiplicity':{str(la):list(map(int,v)) for la,v in mult.items()},'quotient_adjacency_spectrum':{str(N-2*k):int(a) for k,a in enumerate(mult[(n,)]) if a}})
print(json.dumps(out,indent=2))
