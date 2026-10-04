"""Exact difference diagonals, primitive lattice points, and golden norms.

Decimal lambda values are illustrations only. All mathematical acceptance
checks use integers or rational pairs; the physical interpretation is cited
and scoped in the companion note.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def load_clock():
    path = Path(__file__).with_name('golden_prime_clocks_20261003.py')
    spec = spec_from_file_location('golden_clock', path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def arithmetic_tables(limit):
    phi = list(range(limit+1))
    mu = [1]*(limit+1)
    mu[0] = 0
    primes = []
    for p in range(2, limit+1):
        if phi[p] == p:
            primes.append(p)
            for k in range(p, limit+1, p):
                phi[k] -= phi[k]//p
                mu[k] = -mu[k]
            for k in range(p*p, limit+1, p*p):
                mu[k] = 0
    return phi, mu, primes


def polynomial_add(a,b):
    out=[Fraction(0)]*max(len(a),len(b))
    for i,x in enumerate(a): out[i]+=x
    for i,x in enumerate(b): out[i]+=x
    while len(out)>1 and out[-1]==0: out.pop()
    return out


def polynomial_mul_mod(a,b,modulus):
    out=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    while len(out)>=len(modulus):
        factor=out[-1]
        offset=len(out)-len(modulus)
        for j,x in enumerate(modulus): out[offset+j]-=factor*x
        out.pop()
    while len(out)>1 and out[-1]==0: out.pop()
    return out


def main():
    ring = load_clock()
    phi, mu, primes = arithmetic_tables(10000)
    squares = []
    for n in (1, 2, 5, 10, 30, 100, 300, 1000, 10000):
        by_totient = 2*sum(phi[1:n+1])-1
        by_mobius = sum(mu[k]*(n//k)**2 for k in range(1, n+1))
        need(by_totient == by_mobius, 'independent primitive-square formulas')
        if n <= 300:
            direct = sum(gcd(a,b)==1 for a in range(1,n+1) for b in range(1,n+1))
            by_gaps = 1+2*sum(sum(gcd(a,d)==1 for a in range(1,n-d+1))
                              for d in range(1,n))
            need(direct == by_gaps == by_totient, 'literal square / gap diagonals')
        squares.append((n, by_totient, str(Fraction(by_totient,n*n))))

    gap_checks = 0
    for d in range(1, 101):
        for blocks in range(1, 9):
            count = sum(gcd(a,a+d)==1 for a in range(1,blocks*d+1))
            need(count == blocks*phi[d], 'fixed difference primitive period')
            gap_checks += 1

    ring_checks = 0
    for a in range(-50,51):
        for d in range(-50,51):
            alpha = a, a+d
            need(ring.pair_norm(alpha) == a*a-a*d-d*d, 'gap norm')
            need(ring.pair_mul(ring.phi_power(2), (a-d,d)) == alpha, 'unit gap chart')
            shifted = ring.pair_mul((0,1), alpha)
            need(shifted == (a+d,2*a+d), 'Fibonacci gap dynamics')
            need(gcd(*alpha) == gcd(*shifted) == gcd(a,d), 'primitive invariant')
            need(ring.pair_norm(shifted) == -ring.pair_norm(alpha), 'norm sign under unit shift')
            ring_checks += 1
    need(ring.pair_mul((-1,2),ring.phi_power(5)) == (7,11), 'sqrt5 phi5 pair')
    need(ring.pair_mul((4,0),ring.phi_power(4)) == (8,12), 'four phi4 pair')
    need(ring.pair_norm((7,11)) == 5 and ring.pair_norm((8,12)) == 16, 'same-gap ideal hostile')

    unit_counts = []
    norm_root_checks = 0
    for p in (q for q in primes if q <= 97):
        primitive = sum(gcd(gcd(a,b),p)==1 for a in range(p) for b in range(p))
        units = sum(ring.is_unit((a,b),p) for a in range(p) for b in range(p))
        split = p != 2 and p != 5 and pow(5,(p-1)//2,p)==1
        expected = p*(p-1) if p==5 else ((p-1)**2 if split else p*p-1)
        need(primitive == p*p-1 and units == expected, 'primitive versus unit density')
        if p in (2,3,5,7,11,19):
            unit_counts.append((p,primitive,units,str(Fraction(units,primitive))))
        for d in range(1,31):
            roots = sum((a*a-a*d-d*d)%p==0 for a in range(p))
            expected_roots = 1 if d%p==0 or p==5 else (2 if split else 0)
            need(roots == expected_roots, 'golden norm roots on a fixed difference')
            norm_root_checks += 1
    need(gcd(7,11)==1 and gcd(8,12)==4, 'same difference does not preserve primitive content')
    need(gcd(1,3)==1 and not ring.is_unit((1,3),5), 'primitive is not modular unit')

    # Retain the sum when transporting differences through a quadratic map.
    quadratic_checks=0
    for a in range(-20,21):
        for b in range(-20,21):
            x,y=Fraction(a,4),Fraction(b,4)
            s,d=x+y,y-x
            for c in (Fraction(-7,4),Fraction(-29,16)):
                fx,fy=x*x+c,y*y+c
                need(fy-fx==s*d, 'secant difference transport')
                need(fx+fy==(s*s+d*d)/2+2*c, 'retained sum transport')
                quadratic_checks+=1
    cycle=[Fraction(-7,4),Fraction(5,4),Fraction(-1,4)]
    secant=derivative=Fraction(1)
    for i,x in enumerate(cycle):
        y=cycle[(i+1)%3]
        need(x*x-Fraction(29,16)==y,'rational quadratic three-cycle')
        secant*=x+y
        derivative*=2*x
    need(secant==1 and derivative==Fraction(35,8),'off-diagonal versus diagonal multiplier')
    P=[Fraction(-1,8),Fraction(-9,4),Fraction(1,2),Fraction(1)]
    x=[Fraction(0),Fraction(1)]
    F=lambda z: polynomial_add(polynomial_mul_mod(z,z,P),[Fraction(-7,4)])
    y,z=F(x),F(F(x))
    need(F(z)==x and y!=x,'parabolic cubic cycle scheme')
    dprod=polynomial_mul_mod([8],polynomial_mul_mod(x,polynomial_mul_mod(y,z,P),P),P)
    sprod=polynomial_mul_mod(polynomial_add(x,y),polynomial_mul_mod(polynomial_add(y,z),polynomial_add(z,x),P),P)
    need(dprod==sprod==[1],'parabolic secant and derivative controls')
    fixed=[Fraction(-7,4),Fraction(-1),Fraction(1)]
    need(polynomial_mul_mod(P,[1],fixed)==[Fraction(5,2),Fraction(1)],'fixed-point gcd remainder')
    need(Fraction(-5,2)**2+Fraction(5,2)-Fraction(7,4)==7,'fixed-point gcd exclusion')
    b,c,d=P[2],P[1],P[0]
    need(b*b*c*c-4*c**3-4*b**3*d-27*d*d+18*b*c*d==49,'parabolic cubic discriminant')

    # A literal bridge to the incoming rational-anchor chart, and its quotient loss.
    def word_chart(word):
        P,Q,S=1,1,0
        for valuation in word:
            P,S,Q=3*P,3*S+Q,Q*2**valuation
        u,v=ring.phi_power(len(word)+sum(word))
        parameter=Fraction(ring.pair_norm((u-1,v)),Q)
        return Fraction(-S,P-Q),parameter
    anchor,parameter=word_chart((1,1,2))
    need(anchor==Fraction(-19,11) and parameter==Fraction(-29,16),'rational anchor / golden norm parameter')
    need(ring.pair_mul((Fraction(-4,29),Fraction(20,29)),(7,13))==(8,12),
         'q29 golden parity fraction')
    need((7+13*24)%29==0 and (24*24-24-1)%29==0,'index29 ideal evaluation')
    first,second=word_chart((1,1,1,3)),word_chart((1,1,2,2))
    need(first==(Fraction(-65,17),Fraction(-121,64)) and
         second==(Fraction(-73,17),Fraction(-121,64)),'norm parameter loses ordered carry')
    def odd_cycle(root,word):
        orbit=[]
        for valuation in word:
            need(root.denominator%2==1 and root.numerator%2==1,'exact rational odd valuation')
            orbit.append(root)
            root=(3*root+1)/2**valuation
        need(root==orbit[0],'rational shadow cycle closes')
        return set(orbit)
    need(odd_cycle(first[0],(1,1,1,3)).isdisjoint(odd_cycle(second[0],(1,1,2,2))),
         'same norm parameter can represent different cyclic arithmetic orbits')

    # Three parenthesizations have different exact algebraic parameters.
    # x^6=phi^-4/C implies C^2*x^12-7C*x^6+1=0.
    u = ring.phi_power(-4)
    need(u == (5,-3), 'phi minus4')
    for C in (5, 15625, 31250):
        x6 = tuple(Fraction(v,C) for v in u)
        square = ring.pair_mul(x6,x6)
        polynomial = (C*C*square[0]-7*C*x6[0]+1,
                      C*C*square[1]-7*C*x6[1])
        need(polynomial == (0,0), 'lambda exact annihilating equation')
    with localcontext() as ctx:
        ctx.prec=45
        golden=(Decimal(1)+Decimal(5).sqrt())/2
        sixth=Decimal(1)/6
        lambdas = {
            'literal (5*phi^4)^(-1/6)': (5*golden**4)**(-sixth),
            'previous 1/[5*(2*phi^4)^(1/6)]': 1/(5*(2*golden**4)**sixth),
            'alternative 1/[5*(phi^4)^(1/6)]': 1/(5*(golden**4)**sixth),
        }
    print('FINITE-EXACT primitive square (N,count,exact density):', squares)
    print(f'Fixed-gap periods: {gap_checks} tests, d1..100 and1..8 complete periods')
    print(f'Golden difference chart: {ring_checks} signed pairs, a,d=-50..50')
    print('Same difference4: (7,11)=sqrt5*phi^5, norm5,gcd1; (8,12)=4phi^4,norm16,gcd4')
    print('Prime (p,primitive pairs,ring units,conditional fraction):',unit_counts)
    print(f'Norm-root controls: {norm_root_checks}, primes<=97 and gaps1..30')
    print(f'Quadratic pair transport: {quadratic_checks} exact controls, a,b=-20..20 divided by4, two parameters')
    print('c=-29/16 cycle(-7/4,5/4,-1/4): secant product1, diagonal derivative product35/8')
    print('c=-7/4 cubic P=x^3+x^2/2-9x/4-1/8: F^3=x, both products1; fixed-point gcd1, discriminant49')
    print('Word112: rational anchor -19/11; golden(-4+20phi)/29, denominator ideal(phi^7-1)=(29,phi-24); N(phi^7-1)/16=-29/16. Norm parameter loses carry:1113 and1122 give -121/64 but disjoint cycles -65/17 and-73/17.')
    print('Exact lambda equation: C^2*x^12-7*C*x^6+1=0 for C=5,15625,31250')
    for label, value in lambdas.items():
        print(label,'=',value)
    print('No numerical Higgs identity is asserted; approximations above only distinguish written expressions.')
    print('PASS: every mathematical acceptance check remains active under python -O')


if __name__=='__main__':
    main()
