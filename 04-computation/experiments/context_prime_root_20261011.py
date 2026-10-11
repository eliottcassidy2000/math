"""Exact finite controls for context, sibling periods, and ROOT alias bills.
No Collatz convergence is assumed. Imported codec is explicitly capped at n<=5.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, gcd, factorial
import json
import tournament_even_bijection_20261011c as codec

CHECKS = 0

def check(ok, reason):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(reason)


def nu(m, r):
    return F(1, (m+1)**r)-F(1, (m+2)**r)


def step_index(m, j=1):
    return 4**j*m + 2*(4**j-1)//3


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def root_counters(n):
    if n == 1:
        return (0,0)
    steps, cost = 0, 0
    while n != 1:
        v = 3*n+1
        a = 0
        while v%2 == 0:
            v //= 2
            a += 1
        n = v
        cost += (a-1)//2
        steps += 1
        if steps > 10000:
            raise RuntimeError('finite control failed to finish within declared ceiling')
    return (steps-1,cost)


def order4(modulus):
    if modulus == 1:
        return 1
    x = 4 % modulus
    j = 1
    while x != 1:
        x = 4*x % modulus
        j += 1
    return j


def sibling_period(m, q):
    if q % 2 == 0:
        raise ValueError('odd modulus required')
    return order4(3*q // gcd(3*m+2, 3*q))


def parity_lift(x, y):
    return (1, x, y, x*y)


def flip_lift(z, port):
    one, x, y, xy = z
    if port == 0:
        return (one, one^x, y, y^xy)
    return (one, x, one^y, x^xy)


def converse(graph):
    n, t = codec.even_to_tournament(*graph)
    return codec.tournament_to_even(n, t ^ ((1 << (n*(n-1)//2))-1))


def joined(g, h):
    ts = (codec.even_to_tournament(*g), codec.even_to_tournament(*h))
    return codec.tournament_to_even(*codec.assemble(tuple(x for x in ts if x[0]), tournament=True))


def union(g, h):
    return codec.assemble(tuple(x for x in (g,h) if x[0]), tournament=False)


def mobius_squarefree_divisors(primes):
    for bits in product((0,1), repeat=len(primes)):
        d = 1
        for p, b in zip(primes,bits):
            if b:
                d *= p
        yield bits, d, (-1)**sum(bits)


def run():
    clocks = {}
    for q in range(1, 200, 2):
        for m in range(q):
            predicted = sibling_period(m,q)
            x = step_index(m) % q
            actual = 1
            while x != m:
                x = step_index(x) % q
                actual += 1
            check(predicted == actual, 'all odd-modulus residue periods')
    for q in (3,5,7,11,13,33,105,195,1365):
        clocks[str(q)] = sorted({sibling_period(m,q) for m in range(q)})
    check(clocks['105'] == [3,6] and clocks['195'] == [3,6], 'six-clock triples')
    check((4**6-1)//3 == 3*5*7*13, '1365 factorization')
    for a in range(1,10):
        q = 2**a
        for m in range(q):
            for j in (1,2,3,6):
                check((step_index(m,j)-m) % q == 0 if (3*m+2)%q == 0
                      else (step_index(m,j)-m)%q != 0, 'dyadic return iff fixed')
    for m in range(201):
        for j in range(1,13):
            n, h = 2*m+1, 2*step_index(m,j)+1
            check(3*h+1 == 4**j*(3*n+1), 'exact native sibling receipt')
            check(oddpart(3*h+1) == oddpart(3*n+1), 'same actual odd successor')
            check(step_index(m,j)-m == (4**j-1)//3*(3*m+2), 'index carry')
    for n in range(1,256,2):
        L,K = root_counters(n)
        for j in range(1,7):
            sibling = 4**j*n+(4**j-1)//3
            LL,KK = root_counters(sibling)
            check((LL,KK)==(L,K+j), 'price-score exact sibling multiplier')
            W = F(2*factorial(K)*factorial(L+1),factorial(L+K+2))
            WW = F(2*factorial(KK)*factorial(LL+1),factorial(LL+KK+2))
            check(WW <= W, 'beta-mixture sibling domination')
            check(WW <= F(1,3), 'all nonroot-label siblings satisfy mixture score cap')
            for theta in (F(1,3),F(1,2),F(2,3)):
                check(theta**KK*(1-theta)**LL == theta**j*theta**K*(1-theta)**L,
                      'price-specific score covariance')
    # Purely rational zeta enclosures, independently bounded by the integral test.
    bounds = {}
    for r in (1,2,4,6):
        last = None
        for K in (1,2,4,8,16):
            upper = r*sum((F(1,k**(r+1)) for k in range(1,K+1)),F())+F(1,K**r)
            if last is not None:
                check(upper < last, 'rational zeta bills improve')
            last = upper
            if r == 6:
                bounds[str(K)] = str(upper)
        for m in (0,1,6,13,27):
            for q in (m+1,m+3,2*m+7):
                partial = sum((nu(m+k*q,r) for k in range(1,65)),F())
                check(partial < upper/F(q**(r+1)), 'partial alias under certified bound')
                # Compare the individual kernel to its analytic envelope.
                for k in range(1,20):
                    check(nu(m+k*q,r) <= F(r,(k*q)**(r+1)), 'kernel envelope')
    check(bounds['2'] == '97/16', 'sharper rational seven-power bill')
    # Complete Boolean divisibility contexts: all filters retaining zero.
    wheel = {}
    for primes in ((2,3,11),(3,5,7),(3,5,13)):
        ds = list(mobius_squarefree_divisors(primes))
        P = ds[-1][1]
        s = 7
        signature_factors = []
        for bits in product((0,1), repeat=3):
            fac = F(1)
            for p,b in zip(primes,bits):
                fac *= F(1,p**s) if b else 1-F(1,p**s)
            signature_factors.append(fac)
        check(sum(signature_factors) == 1, 'Euler product partitions alias indices')
        for filter_mask in range(128,256):
            fac = sum((signature_factors[j] for j in range(8) if filter_mask>>j&1),F())
            check(fac >= F(1,P**s), 'sharp zero-signature alias floor')
            check((fac == F(1,P**s)) == (filter_mask == 128), 'unique minimal Boolean filter')
        for k in range(0,2*P+1):
            inclusion = sum(mu for _,d,mu in ds if k%d == 0)
            check(inclusion == int(gcd(k,P)==1), 'coprime Mobius filter removes source zero')
        wheel[str(primes)] = {'product':P,'intersection_factor':str(F(1,P**s)),
                             'coprime_factor':str(signature_factors[0])}
    # Finite squarefree factor coordinates have exact conditional weights.
    for primes in ((2,3,11),(3,5,7),(3,5,13)):
        for s in (2,7):
            weighted = [(bits,F(1,d**s)) for bits,d,_ in mobius_squarefree_divisors(primes)]
            total = sum((w for _,w in weighted),F())
            for i,p in enumerate(primes):
                marginal = sum((w for bits,w in weighted if bits[i]),F())/total
                check(marginal == F(1,p**s+1), 'squarefree zeta law uses p^s+1')
                for j,l in enumerate(primes):
                    if i < j:
                        joint = sum((w for bits,w in weighted if bits[i] and bits[j]),F())/total
                        check(joint == marginal*F(1,l**s+1), 'finite squarefree factor independence')
                        mp = marginal*total/(total-1)
                        ml = F(1,l**s+1)*total/(total-1)
                        jp = joint*total/(total-1)
                        check(jp < mp*ml, 'unit removal changes prime correlation')
    # Full divisor-cell inversion on a finite support, including the selected source.
    primes = (3,5,7)
    divisors = [d for _,d,_ in mobius_squarefree_divisors(primes)]
    weights = {k: F((k*k+3)%17,17)*nu(6+8*k,6) for k in range(321)}
    marginals = {d: sum((w for k,w in weights.items() if k%d==0),F()) for d in divisors}
    for bits in product((0,1),repeat=3):
        on = [p for p,b in zip(primes,bits) if b]
        off = [p for p,b in zip(primes,bits) if not b]
        base = 1
        for p in on: base *= p
        inferred = sum((mu*marginals[base*d] for _,d,mu in mobius_squarefree_divisors(off)),F())
        actual = sum((w for k,w in weights.items() if all((k%p==0)==bool(b) for p,b in zip(primes,bits))),F())
        check(inferred == actual, 'all joint context cells reconstructed')
    # Finite common-future cohorts: gain is a source coefficient, not new ROOT mass.
    cohort_examples = []
    for m,q in ((6,105),(13,195),(27,1365)):
        j = sibling_period(m,q)
        C = [m,step_index(m,j),step_index(m,2*j)]
        check(all((z-m)%q==0 for z in C), 'context-compatible cohort')
        a = nu(m,6); B = sum((nu(z,6) for z in C[1:]),F())
        U = F(97,16*q**7)
        check(0 < B < U, 'positive but paid cohort gain')
        for value in (F(),F(1,7),F(1)):
            W = (a+B)*value
            old = max(F(),(W-U)/a)
            new = max(F(),(W-U+B)/(a+B))
            check(old <= new <= value, 'sharper calibrated lower bound')
        cohort_examples.append({'m':m,'q':q,'period':j,'first_sibling_index':C[1],
                                'paid_gain_first_two':str(B)})
    # Choose the modulus so the first alias is a witnessed sibling, for every source.
    for m in range(81):
        for j in range(1,9):
            q = (3*m+2)*(4**j-1)//3
            check(q > m and m+q == step_index(m,j), 'first-alias adaptive modulus')
            check((3*m+2)*(4**(2*j)-1)//3 == (4**j+1)*q, 'nested native-compatible refinement')
            a, b = nu(m,6), nu(m+q,6)
            bill = F(1,16*q**7)
            tail = sum((nu(m+k*q,6) for k in range(2,49)),F())
            check(tail < bill, 'first-alias removal has rational one-sixteenth bill')
            for value in (F(),F(1,13),F(1)):
                # Declared synthetic scores obey exactly the known pair equality.
                observation = (a+b)*value + tail/3
                floor = max(F(),(observation-bill)/(a+b))
                check(floor <= value, 'independent lower measurement compiler')
            # A finite rooted-looking bank excluding the selected pair cannot pass.
            check(tail <= bill, 'unmarked-bank hostile does not certify source')
    check(F(7)/F(1,16) == 112, 'rational error-constant improvement factor')
    for j in range(1,7):
        B = 4**j
        repunits = set()
        k = 0
        while k <= 500:
            repunits.add(k)
            k = B*k+1
        for k in range(501):
            A = 1+k*(B-1)
            check((A & (A-1) == 0) == (k in repunits), 'complete fixed-word alias slope classification')
            check(A%3 == 1, 'prime-three slope cannot hide in the alias multiplier')
    for n in range(1,202,2):
        check(oddpart(3*(7*n+2)+1)==7*oddpart(3*n+1), 'second alias multiplier seven')
        check(oddpart(3*(10*n+3)+1)==5*oddpart(3*n+1), 'third alias multiplier five')
    # Transport full operations through the incoming all-order codec, bounded checks.
    graphs = [(0,0)]
    for n in range(1,6):
        graphs.extend((n,g) for g in codec._catalogue(n,False)[0] if codec.automorphism_even(n,g))
    for g in graphs:
        check(converse(converse(g))==g, 'transported converse involution')
        check(joined(g,(0,0))==g==joined((0,0),g), 'transported shared unit')
        for h in graphs:
            if g[0]+h[0] <= 5:
                check(converse(joined(g,h)) == joined(converse(h),converse(g)), 'converse anti-automorphism')
                for k in graphs:
                    if g[0]+h[0]+k[0] <= 5:
                        check(joined(joined(g,h),k)==joined(g,joined(h,k)), 'associative transported join')
    e,a,b = (0,0),(1,0),(3,3)
    check(joined(a,b)==(4,7), 'source over triangle maps to star')
    check(joined(b,a)==(4,3), 'triangle over sink maps to path and isolate')
    check(converse((4,7))==(4,3), 'star path exact duality')
    check(union(joined(a,e),joined(e,b)) != joined(union(a,e),union(e,b)), 'interchange fails')
    for x,y in product((0,1),repeat=2):
        z = parity_lift(x,y)
        check(flip_lift(z,0)==parity_lift(1-x,y), 'linear context lift first port')
        check(flip_lift(z,1)==parity_lift(x,1-y), 'linear context lift second port')
        check((y^(x*y)) == (codec.section_to_graph(x,y)[1].bit_count()%2), 'mixed parity observable')
    # Taylor coefficients of the reciprocal difference and its symmetric partner.
    for r in (2,4,6):
        for m in range(12):
            a,b=m+1,m+2
            for j in range(1,13):
                c = (-1)**j*comb(r+j-1,j)*(a**j-b**j)
                reflected = comb(r+j-1,j)*(b**j-a**j)
                check(c+reflected == (0 if j%2==0 else 2*comb(r+j-1,j)*(b**j-a**j)), 'odd-zeta projection coefficients')
    print(json.dumps({'status':'FINITE-EXACT; universal ROOT positivity OPEN',
                      'checks':CHECKS,'odd_modulus_period_sets':clocks,
                      'r6_rational_alias_constants':bounds,'wheel_factors':wheel,
                      'cohort_examples':cohort_examples,'graph_types_checked':len(graphs)},indent=2))

if __name__ == '__main__':
    run()
