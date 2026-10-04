#!/usr/bin/env python3
"""Exact controls for sixth_clock_branches_20261004.md; Python stdlib only.

All checks survive python -O. Integer computations are proofs of the stated
finite controls; all-height conclusions require the arguments in the note.
Optional --pari independently checks the number field's maximal order.
"""
import argparse
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
import shutil
import subprocess


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def valuation(n, p):
    check(n != 0, "valuation of zero")
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def qmul(x, y, modulus=None):
    a, b = x
    c, d = y
    z = (a*c+b*d, a*d+b*c+b*d)
    return tuple(v % modulus for v in z) if modulus else z


def power(x, exponent, mul, identity):
    y = identity
    while exponent:
        if exponent & 1:
            y = mul(y, x)
        x = mul(x, x)
        exponent //= 2
    return y


def qp(k, q=None):
    return power((0, 1), k, lambda a, b: qmul(a, b, q), (1, 0))


def norm(x):
    a, b = x
    return a*a+a*b-b*b


ONE = (1, 0, 0, 0, 0, 0)
T = (0, 1, 0, 0, 0, 0)
BETA = (1, 1, 0, 0, 0, 0)


def smul(a, b, q=None):
    """Q[t], t^6 = 5*t^3+5, in the power basis."""
    c = [0]*11
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    for j in range(10, 5, -1):
        c[j-3] += 5*c[j]
        c[j-6] += 5*c[j]
    return tuple(v % q for v in c[:6]) if q else tuple(c[:6])


def sp(a, k, q=None):
    return power(a, k, lambda x, y: smul(x, y, q), ONE)


def sadd(a, b):
    return tuple(x+y for x, y in zip(a, b))


def det(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        d = a[j][j]
        result *= d
        for i in range(j+1, len(a)):
            c = a[i][j]/d
            a[i] = [x-c*y for x, y in zip(a[i], a[j])]
    check(result.denominator == 1, "noninteger determinant")
    return int(result)


def trace(a):
    basis = [tuple(int(i == j) for i in range(6)) for j in range(6)]
    return sum(smul(a, b)[j] for j, b in enumerate(basis))


def check_order(x, order, prime_factors, mul, identity):
    check(power(x, order, mul, identity) == identity, "not a period")
    for p in prime_factors:
        check(power(x, order//p, mul, identity) != identity, "not exact period")


def mersenne_and_rays():
    past = 1
    exceptions = []
    for n in range(1, 61):
        term = 2**n-1
        remainder = term
        while (d := gcd(remainder, past)) > 1:
            remainder //= d
        if remainder == 1:
            exceptions.append(n)
        past *= term
    check(exceptions == [1, 6], "bounded primitive-prime exceptions")
    print("Mersenne primitive-support exceptions, n=1..60:", exceptions)
    check(63 == 3*19+6 and 63 != 2*19+6, "63 arithmetic")
    print("Corrected arithmetic: 63=3*19+6; 2*19+6=44")
    for residue, h0 in ((3, 5), (5, 3), (1, 7)):
        seeds = []
        for j in range(30):
            h = h0+6*j
            n = (2**(h+1)-1)//3
            check(n % 6 == residue and (3*n+1)//2 == 2**h, "ray")
            if j < 3:
                seeds.append(n)
        print(f"Plus single-halving row {residue}: h={h0}+6j; seeds={seeds}")
    for j in range(100):
        n = (4**(j+1)-1)//3
        check(valuation(n, 3) == valuation(j+1, 3), "root valuation")
        check(64*n+21 == (4**(j+4)-1)//3, "three-ray shift")
    for s in range(1, 8):
        q = 3**s
        orbit = []
        x = 1
        for _ in range(q):
            orbit.append(x)
            x = (4*x+1) % (2*q)
        check(x == 1 and len(set(orbit)) == q, "triadic odometer")
        first_indices = [q if r == 0 else r for r in range(q)]
        check(sum(a != r for r, a in enumerate(first_indices)) == 1,
              "one deleted-root offset")
    print("Root rays: 100 valuation checks; odometer and one offset through 3^7")
    check((3*3+1)//2 == 5, "all-row power-of-two hostile")
    # M_6 has no new prime in the full sequence, but 7 is new among M_2,M_4,M_6.
    check(gcd(7, 3*15) == 1 and (2**3-1) % 7 == 0, "subsequence novelty")
    print("Hostiles: F(3)=5; 7 is new in the even-index subsequence at M_6")


def golden_and_19():
    check(qp(3) == (1, 2) and qp(6) == (5, 8), "six golden identity")
    check((qp(18)[0]-1, qp(18)[1]) == tuple(76*x for x in qp(9)), "76")
    c9, c18 = (7, 10), (5, 6)
    check(norm(c9) == norm(c18) == 19, "split norm19")
    check(qmul(c9, c18) == tuple(19*x for x in qp(6)), "split product")
    check(2**6-2**3+1 == 3*19, "ordinary Phi18")
    check((7+10*5) % 19 == (5+6*15) % 19 == 0, "split roots")
    check(pow(2, 16, 19) == 5 and pow(2, 11, 19) == 15, "split logs")
    print("Phi_9(phi)=(7,10), Phi_18(phi)=(5,6), both norm19")
    print("phi^18-1=76*phi^9; 5=2^16 and 15=2^11 modulo19")
    for k in range(1, 13):
        q = 2**k
        check_order((0, 1), 3*2**(k-1), {3} | ({2} if k > 1 else set()),
                    lambda x, y: qmul(x, y, q), (1, 0))
        r = 3**k
        check_order(2, 2*3**(k-1), {2} | ({3} if k > 1 else set()),
                    lambda x, y: x*y % r, 1)
    print("Dual clocks checked through prime-power depth12")
    check(valuation(2**18-1, 19) == 1, "ordinary lift seed")
    for k in range(1, 9):
        q, order = 19**k, 18*19**(k-1)
        primes = {2, 3} | ({19} if k > 1 else set())
        check_order(2, order, primes, lambda x, y: x*y % q, 1)
        check_order((0, 1), order, primes,
                    lambda x, y: qmul(x, y, 4*q), (1, 0))
    print("Synchronized base2/mod19^k and golden/mod(4*19^k) orders checked k=1..8")
    lengths = []
    unseen = {(a, b) for a in range(76) for b in range(76) if gcd(gcd(a, b), 76) == 1}
    while unseen:
        x = min(unseen)
        start = x
        length = 0
        while x in unseen:
            unseen.remove(x)
            x = (x[1], (x[0]+x[1]) % 76)
            length += 1
        check(x == start, "golden permutation cycle")
        lengths.append(length)
    check(lengths == [18]*240, "complete primitive q76 universe")
    print("Complete primitive q76 control: 4320 phases, 240 cycles, all length18")


def sextic():
    alpha = (1, 3)
    check(qmul(alpha, alpha) == tuple(5*x for x in qp(4)), "lambda cube")
    check(norm(alpha) == -5, "cube obstruction norm")
    phi = (Fraction(-1, 3), 0, 0, Fraction(1, 3), 0, 0)
    check(smul(phi, phi) == sadd(phi, ONE), "quadratic subfield")
    power_basis = [sp(T, j) for j in range(6)]
    integral_basis = [sp(T, j) for j in range(3)] + [smul(phi, sp(T, j)) for j in range(3)]
    d_power = det([[trace(smul(a, b)) for b in power_basis] for a in power_basis])
    d_integral = det([[trace(smul(a, b)) for b in integral_basis] for a in integral_basis])
    check(d_power == 3**12*5**5 and d_integral == 3**6*5**5, "trace discriminants")
    check(d_power == 27**2*d_integral, "order index")
    print(f"Sextic discriminants: power order={d_power}; integral basis={d_integral}; index27")
    # At 2 the index27 is a unit: this power basis gives the full residue field.
    check_order(T, 9, {3}, lambda x, y: smul(x, y, 2), ONE)
    check_order(BETA, 63, {3, 7}, lambda x, y: smul(x, y, 2), ONE)
    orbit = {sp(BETA, j, 2) for j in range(63)}
    check(len(orbit) == 63 and (0,)*6 not in orbit, "F64 all units")
    check(sp(BETA, 21, 2) == (1, 0, 0, 1, 0, 0), "beta21=phi")
    check(sp(BETA, 56, 2) == T, "beta56=t")
    check(sp(BETA, 63, 4) == (1, 0, 0, 0, 2, 0), "binary lift seed1")
    check(sp(BETA, 126, 8) == (1, 0, 4, 0, 4, 4), "binary lift seed2")
    # Independent six-bit recurrence in the beta power basis.
    start = state = (1, 0, 0, 0, 0, 0)
    states = set()
    while state not in states:
        states.add(state)
        state = state[1:]+(state[4] ^ state[3] ^ state[1] ^ state[0],)
    check(state == start and len(states) == 63 and (0,)*6 not in states,
          "six-bit recurrence")
    print("F64: t order9; beta=t+1 order63; beta^21=phi; beta^56=t")
    print("Independent binary recurrence: all63 nonzero six-bit states")
    for k in range(1, 17):
        q, order = 2**k, 63*2**(k-1)
        check_order(BETA, order, {3, 7} | ({2} if k > 1 else set()),
                    lambda x, y: smul(x, y, q), ONE)
    # A complete second-level census is small; higher levels use exact-order controls.
    unseen = {v for v in product(range(4), repeat=6) if any(x % 2 for x in v)}
    lengths = []
    while unseen:
        x = start = min(unseen)
        length = 0
        while x in unseen:
            unseen.remove(x)
            x = smul(BETA, x, 4)
            length += 1
        check(x == start, "sextic permutation cycle")
        lengths.append(length)
    check(lengths == [126]*32, "all sextic q4 primitive phases")
    print("Sextic binary tower: exact orders checked k=1..16")
    print("Complete q4 control: 4032 primitive phases, 32 cycles, all length126")
    # Binary return translation c gives a reproducible five-bit child selector.
    for k in range(1, 6):
        q, order = 2**k, 63*2**(k-1)
        a = sp(BETA, order, 2*q)
        c = tuple(((x-int(i == 0))//q) % 2 for i, x in enumerate(a))
        pivot = next(i for i, x in enumerate(c) if x)
        labels = {}
        for w in product(range(2), repeat=6):
            label = tuple(x ^ (w[pivot]*y) for x, y in zip(w, c))
            labels.setdefault(label, []).append(w)
        check(len(labels) == 32 and all(len(v) == 2 for v in labels.values()), "child selectors")
    print("Five-bit child selectors checked at depths1..5")


def finite_fourier():
    def evaluate(word, root, q):
        answer = (0,)*6
        for a in reversed(word):
            answer = sadd(smul(root, answer, q), (a, 0, 0, 0, 0, 0))
            answer = tuple(x % q for x in answer)
        return answer
    count = 0
    previous = None
    for k in range(1, 9):
        q = 2**k
        e = 2**(k-1)*pow(2**(k-1), -1, 63)
        tau = sp(BETA, e, q)
        if previous is not None:
            check(tuple(x % (q//2) for x in tau) == previous, "compatible odd clock")
        previous = tau
        check_order(tau, 63, {3, 7}, lambda x, y: smul(x, y, q), ONE)
        for length in (3, 7, 9, 21, 63):
            zeta = sp(tau, 63//length, q)
            for p in range(1, length):
                if length % p:
                    continue
                word = tuple(1+(i*i+3*i) % 4 for i in range(p))*(length//p)
                check(evaluate(word, zeta, q) == (0,)*6, "repetition annihilation")
                for delta in (1, 2, 3, q):
                    changed = word[:-1]+(word[-1]+delta,)
                    check((evaluate(changed, zeta, q) == (0,)*6) == (delta % q == 0),
                          "last-entry detection bound")
                    count += 1
        # Inherited rational anchor112 and negative integer anchor1112114.
        root21 = sp(tau, 3, q)
        for word in ((1, 1, 2)*7, (1, 1, 1, 2, 1, 1, 4)*3):
            check(evaluate(word, root21, q) == (0,)*6, "anchor21 Fourier zero")
    print(f"Compatible odd clock and finite Fourier controls: {count} last-entry tests, depths1..8")
    print("Period21 controls: anchor112 repeated7; anchor1112114 repeated3")
    print("Fourier hostile: a last-entry difference2 vanishes modulo2, survives modulo4")


def divisors():
    def direct(n):
        ds = [d for d in range(2, n) if n % d == 0]
        squarefree = lambda d: all(d % (p*p) for p in range(2, isqrt(d)+1))
        prime = lambda d: all(d % p for p in range(2, isqrt(d)+1))
        return len(ds), sum(map(squarefree, ds)), sum(map(prime, ds))
    check(direct(63) == (4, 3, 2), "63 unbalanced")
    check(direct(126) == (10, 7, 3), "126 balanced")
    for k in range(2, 9):
        # Enumerate exponent boxes for 2^k*3*5; S_k means k-free divisors.
        exponents = list(product(range(k+1), range(2), range(2)))
        ds = [v for v in exponents if v != (0, 0, 0) and v != (k, 1, 1)]
        f, s, u = len(ds), sum(all(a < k for a in v) for v in ds), sum(sum(v) == 1 for v in ds)
        check((f, s, u) == (4*k+2, 4*k-1, 3), "k-free inherited extension")
    check(2**(9-1)-9-1 == 246, "nine axes hostile")
    print("Divisor controls: 63 gives(4,3,2); 126 gives(10,7,3); k-free k=2..8")
    print("Nine-axis marked profile defect=246; fixed squarefree balance does not recurse")


def history_and_macros():
    histories = {(0,)*6: ()}
    for _ in range(12):
        new = {}
        for state, word in histories.items():
            for bit in (0, 1):
                z = sadd(smul(BETA, state), (bit, 0, 0, 0, 0, 0))
                check(z not in new, "history collision")
                new[z] = word+(bit,)
        histories = new
    check(len(histories) == 4096, "all binary histories")
    for source in range(-200, 201):
        n, z, word = source, (0,)*6, ()
        for _ in range(12):
            bit = n % 2
            z = sadd(smul(BETA, z), (bit, 0, 0, 0, 0, 0))
            n = (3*n+1)//2 if bit else n//2
            word += (bit,)
        check(histories[z] == word, "signed history decode")
    for u in (1, 5, -1, -5, -17):
        a0 = 1 if u % 3 == 2 else 2
        for j in range(26):
            n, remainder = divmod(2**a0*4**j*u-1, 3)
            check(remainder == 0 and n % 2 == 1, "integer odd inverse")
            a = valuation(3*n+1, 2)
            check(a == a0+2*j and (3*n+1)//2**a == u, "signed macro")
    print("History controls: 4096 words of length12, 401 signed source histories")
    print("Signed inverse macros: 130 exact ray members for targets1,5,-1,-5,-17")


def pari_audit():
    gp = shutil.which("gp") or "/opt/homebrew/bin/gp"
    script = '''f=x^6-5*x^3-5;
if(polisirreducible(f)!=1,error("reducible"));
n=nfinit(f);
if(poldisc(f)!=1660753125,error("polynomial discriminant"));
if(nfdisc(f)!=2278125,error("independent field discriminant"));
if(n.disc!=2278125 || n.index!=27,error("field discriminant or index"));
forprime(p=2,5,d=idealprimedec(n,p);print(p," prime decomposition [e,f]: ",vector(#d,i,[d[i][3],d[i][4]])));
print("PARI independent number-field audit PASS");
'''
    run = subprocess.run([gp, "-q", "-f"], input=script, text=True, capture_output=True, check=True)
    check("PARI independent number-field audit PASS" in run.stdout and not run.stderr.strip(), "PARI audit")
    print(run.stdout.strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pari", action="store_true")
    args = parser.parse_args()
    mersenne_and_rays()
    golden_and_19()
    sextic()
    finite_fourier()
    divisors()
    history_and_macros()
    if args.pari:
        pari_audit()
    print("PASS: exact controls; universal signed Collatz coverage remains OPEN")
