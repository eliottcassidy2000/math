#!/usr/bin/env python3
"""Exact difference carriers and the golden 19-adic phase tower.

Universe and proof boundaries are in the companion research note.
All checks remain active under python -O; no floating arithmetic establishes
any lattice, period, guard, or integer-cycle claim.
"""
from collections import Counter
from fractions import Fraction
from math import gcd
from decimal import Decimal, localcontext


def check(condition, message="check failed"):
    if not condition:
        raise AssertionError(message)


def sign(a, b):
    """Exact sign of a+b*phi, for integer a,b."""
    u, v = 2*a+b, b
    if v == 0:
        return (u > 0)-(u < 0)
    if u == 0 or u*v >= 0:
        return 1 if v > 0 else -1
    return (1 if u > 0 else -1)*(1 if u*u > 5*v*v else -1)


def domain(a, b, q):
    """Open L-shaped fundamental domain; exact q>1 has no boundary points."""
    return (sign(a, b) > 0 and sign(a-q, b) < 0
            and sign(a+b-q, -b) < 0
            and (sign(a+b-q, q-b) > 0
                 or (sign(a+q, b-q) < 0 and sign(a+b, q-b) > 0)))


def representative(a, b, q):
    # x in (0,1), x' in (-phi,1) bounds b/q in (-1,2), a/q in (-3,3).
    candidates = [(a+i*q, b+j*q) for i in range(-3, 3) for j in range(-1, 2)
                  if domain(a+i*q, b+j*q, q)]
    check(len(candidates) == 1, (a, b, q, candidates))
    return candidates[0]


def beta(a, b, q):
    digit = int(sign(b-q, a+b) > 0)
    return (b-digit*q, a+b), digit


def mmul(M, N, modulus=None):
    a,b,c,d=M; e,f,g,h=N
    result=(a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)
    return tuple(x % modulus for x in result) if modulus else result


def mpow(M, exponent, modulus=None):
    result=(1,0,0,1)
    while exponent:
        if exponent & 1:
            result=mmul(result,M,modulus)
        M=mmul(M,M,modulus); exponent//=2
    return result


G=(0,1,1,1)


def j2(q):
    return sum(gcd(gcd(a,b),q)==1 for a in range(q) for b in range(q))


def phase_cycles(q, decode=False, profile=False):
    visited=bytearray(q*q); lengths=Counter(); roots=[]; total=0
    odd_counts=Counter()
    for a in range(q):
        for b in range(q):
            start=a*q+b
            if visited[start] or gcd(gcd(a,b),q)!=1:
                continue
            u,v=a,b; length=0
            while not visited[u*q+v]:
                visited[u*q+v]=1; length+=1
                u,v=v,(u+v)%q
            check((u,v)==(a,b),"module map is not a permutation")
            lengths[length]+=1; total+=length
            if decode:
                x,y=representative(a,b,q); initial=(x,y)
                A=p=B=0
                for _ in range(length):
                    (x,y),digit=beta(x,y,q)
                    if digit:
                        B=3*B+2**A; p+=1
                    else:
                        A+=1
                check((x,y)==initial,"golden lift does not close")
                odd_counts[p]+=1
                denominator=2**A-3**p
                if B%denominator==0:
                    root=B//denominator; current=root; values=[]
                    x,y=initial
                    for _ in range(length):
                        (x,y),digit=beta(x,y,q)
                        check(current%2==digit,"integer realization has wrong parity")
                        values.append(current)
                        current=3*current+1 if digit else current//2
                    check(current==root)
                    roots.append(min(values,key=abs))
    if profile:
        check(len(lengths)==1,'profile sign threshold needs a common period')
        L=next(iter(lengths))
        negative=sum(number for p,number in odd_counts.items() if 3**p>2**(L-p))
        if q==1444:check(min(odd_counts)==83 and max(odd_counts)==106 and negative==0)
        print('q',q,'odd-count profile',dict(sorted(odd_counts.items())),
              'negative rational cycles',negative)
    return total,dict(sorted(lengths.items())),sorted(roots)


def pair_add(P,Q):
    return (P[0]+Q[0],P[1]+Q[1])


def pair_mul(P,Q):
    a,b=P; c,d=Q
    return (a*c+b*d,a*d+b*c)


def pair_collatz(P):
    a,b=P
    return (3*a+1,3*b) if (a-b)%2 else (a//2,b//2)


def value(P):
    return P[0]-P[1]


def owner_identities():
    triangular=lambda n:n*(n+1)//2
    splice=lambda a,b:a+b+1
    cross=lambda a,b:(a+1)*(b+1)
    for a in range(20):
        for b in range(20):
            N=a+b+2
            check(triangular(N-1)==triangular(a)+triangular(b)+cross(a,b))
            for c in range(20):
                check(cross(a,b)+cross(splice(a,b),c)==cross(b,c)+cross(a,splice(b,c)))
    check(triangular(9)+triangular(10)==100)
    check(triangular(19)==2*triangular(9)+100)
    for n in range(-10,11):
        for radix in range(1,31):
            for k in range(100):
                quotient,digit=divmod(k,radix)
                P=(max(n,0)+k,max(-n,0)+k)
                Q=(max(n,0)+quotient,max(-n,0)+quotient)
                check(value(P)==value(Q)==n and radix*quotient+digit==k)
    check(value((11,7))==value((12,8))==4)
    print('Triangular split: 400 pairs; cocycle: 8000 triples; fibre radix bijection: 63000 instances')

    for c in (Fraction(-7,4),Fraction(-29,16)):
        for i in range(-9,10):
            for j in range(-9,10):
                x,y=Fraction(i,4),Fraction(j,4)
                delta,s=x-y,x+y
                fx,fy=x*x+c,y*y+c
                check(fx-fy==s*delta)
                check(fx+fy==(s*s+delta*delta)/2+2*c)
    orbit=[Fraction(-7,4),Fraction(5,4),Fraction(-1,4)]
    check(all(orbit[i]**2-Fraction(29,16)==orbit[(i+1)%3] for i in range(3)))
    check(8*orbit[0]*orbit[1]*orbit[2]==Fraction(35,8))
    print('Quadratic difference and sum: 722 pairs; -29/16 cycle and multiplier 35/8')
    n=-53; path=[]
    while n not in path:
        path.append(n); n=3*n+1
        while n%2==0:n//=2
    check(path==[-53,-79,-59,-11,-1] and n==-1)
    print('Negative-root extrapolation hostile:',path,'then repeats -1')

    # Exact coefficient-pair checks for the sixth powers, separate from
    # decimal displays of the roots and the illustrative tree mass.
    phi_minus4=(Fraction(5),Fraction(-3))
    def mul(X,Y):
        a,b=X;c,d=Y
        return (a*c+b*d,a*d+b*c+b*d)
    for scale in (5,31250):
        t=tuple(v/scale for v in phi_minus4); t2=mul(t,t)
        check((scale**2*t2[0]-7*scale*t[0]+1,
               scale**2*t2[1]-7*scale*t[1])==(0,0))
    with localcontext() as ctx:
        ctx.prec=55
        phi=(1+Decimal(5).sqrt())/2
        literal=(5*phi**4)**(-Decimal(1)/6)
        historical=1/(5*(2*phi**4)**(Decimal(1)/6))
        for label,lam in [('literal',literal),('historical',historical)]:
            mass=Decimal('246.22')*(2*lam).sqrt()
            print(label,'lambda',format(lam,'.18f'),'illustrative mass GeV',format(mass,'.12f'))


def carry_normalization():
    cases=maximum=0
    for q in (2,11,76):
        for a in range(-40,41):
            for b in range(-40,41):
                if gcd(gcd(a,b),q)!=1 or sign(a,b)<=0 or sign(a-q,b)>=0:
                    continue
                P=(a,b); R=representative(a%q,b%q,q); steps=0
                while P!=R:
                    u,v=(P[0]-R[0])//q,(P[1]-R[1])//q
                    NP,d=beta(*P,q);NR,rd=beta(*R,q)
                    check(((NP[0]-NR[0])//q,(NP[1]-NR[1])//q)==(v-(d-rd),u+v))
                    P,R=NP,NR;steps+=1
                    check(steps<1000,'normalization control exceeded its bound')
                cases+=1;maximum=max(maximum,steps)
    print('Golden lattice-carry normalization:',cases,'field inputs; max entry time',maximum)
    for start,target,q,root in [(151,1,2,(0,1)),(-9,-5,11,(-1,7))]:
        n=start;bits=[]
        while n!=target:
            bits.append(n%2);n=3*n+1 if n%2 else n//2
        a,b=root
        for d in bits[::-1]:a,b=b-a-d*q,a+d*q
        P=(a,b);R=representative(a%q,b%q,q)
        offset=((a-R[0])//q,(b-R[1])//q)
        print('Certified field input:',start,'theta numerator',P,'q',q,'lattice offset',offset)
        for d in bits:
            P,digit=beta(*P,q);R,rd=beta(*R,q)
            check(digit==d)
        check(P==R==root)

    # The inherited rational anchor supplies a concrete 11/19/29/76 bridge.
    P=(-4,20);q=29;n=Fraction(-19,11);initial=n;bits=[]
    for _ in range(7):
        digit=n.numerator%2
        NP,d=beta(*P,q);check(d==digit)
        n=3*n+1 if digit else n/2
        P=NP;bits.append(digit)
    check(P==(-4,20) and n==initial and bits==[1,0,1,0,1,0,0])
    check(3*29-11==76==4*19)
    print('Rational -19/11 anchor: raw word 1010100; theta=(-4+20*phi)/29; 3*29-11=76=4*19')


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--large',action='store_true',help='decode all 4560 cycles at exact denominator 1444')
    args=parser.parse_args()

    # Direct comparison with the earlier *trap* enumeration, a different
    # domain and algorithm, on a complete stated finite universe.
    from collatz_golden_carriers_20261004 import trap,cycles_in
    checks=0
    for q in list(range(2,31))+[38,76,100]:
        trap_points=trap(q); cycles=cycles_in(trap_points,q)
        periodic=set(sum(cycles,[]))
        window={s for s in trap_points if domain(*s,q)}
        primitive={(a,b) for a in range(q) for b in range(q) if gcd(gcd(a,b),q)==1}
        check(periodic==window)
        check({(a%q,b%q) for a,b in window}==primitive)
        check(len(window)==len(primitive))
        for a,b in primitive:
            check(representative(a,b,q) in window)
        checks+=1
    print('L-domain / independent old trap / every primitive phase:',checks,'complete denominators')
    for q in (2,3,4,5,10,11,19,29,38,76,100):
        total,lengths,roots=phase_cycles(q,decode=True,profile=(q==76))
        check(total==j2(q))
        print('q',q,'J2',total,'periods',lengths,'integer roots',roots)

    # The tower proof uses this exact identity and p-adic binomial lifting.
    M=mpow(G,18); K=mpow(G,9)
    check(tuple(M[i]-(1 if i in (0,3) else 0) for i in range(4))==tuple(76*x for x in K))
    for k in range(1,9):
        q=4*19**k; length=18*19**(k-1)
        check(mpow(G,length,q)==(1,0,0,1))
        for divisor in (2,3,19):
            if length%divisor==0:
                check(mpow(G,length//divisor,q)!=(1,0,0,1))
        print('tower k',k,'q',q,'primitive phases',4320*19**(2*k-2),
              'common period',length,'cycles',240*19**(k-1))
    if args.large:
        total,lengths,roots=phase_cycles(1444,decode=True,profile=True)
        check(total==1559520 and lengths=={342:4560} and roots==[])
        print('COMPLETE q=1444 integer filter:',total,'phases;',lengths,'cycles; roots',roots)

    pairs=[(a,b) for a in range(13) for b in range(13)]
    for P in pairs:
        for Q in pairs:
            check(value(pair_add(P,Q))==value(P)+value(Q))
            check(value(pair_mul(P,Q))==value(P)*value(Q))
            check(pair_add(P,Q)==pair_add(Q,P) and pair_mul(P,Q)==pair_mul(Q,P))
        n=value(P); expected=3*n+1 if n%2 else n//2
        check(value(pair_collatz(P))==expected)
    print('Difference-pair arithmetic:',len(pairs)**2,'ordered pairs; Collatz lifts',len(pairs))
    # An extra fibre coordinate is not automatically a descending coordinate.
    P=(1,2); records=[]
    for _ in range(8):
        records.append((value(P),min(P)))
        P=pair_collatz(pair_collatz(P))
    check(all(n==-1 for n,k in records) and all(records[i+1][1]>records[i][1] for i in range(7)))
    print('Hostile over the -1 cycle (value, common part):',records)
    owner_identities()
    carry_normalization()
    print('ALL CHECKS PASSED')


if __name__=='__main__':
    main()
