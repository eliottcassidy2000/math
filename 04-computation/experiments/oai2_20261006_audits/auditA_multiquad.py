"""Audit A (independent of s4_multiquad.py): decomposition of primes in multiquadratic L = Q(sqrt a_1..a_k) in C,
relative to L+ = L cap R, via Frobenius/inertia computed from quadratic subfields; c = complex conjugation.
Bound from Lemma R: c in I_p -> 2 (p=2) / 3; c in D_p \ I_p -> kappa(q), q = p^(f/2) (f = |D/I| = 2 here, so q = p)."""
import itertools
from sympy import factorint, primerange

KAPPA = {2: 4, 3: 3, 4: 4, 5: 4, 7: 4, 8: 4, 9: 3, 11: 5, 16: 4}   # audit A SAT values
def sqf(n):
    s = -1 if n < 0 else 1
    r = 1
    for p, e in factorint(abs(n)).items():
        if e % 2: r *= p
    return s * r

def disc(D):
    return D if D % 4 == 1 else 4 * D

def leg(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def behaviour(D, p):   # in Q(sqrt D), D squarefree != 1
    if disc(D) % p == 0: return 'R'
    if p == 2: return 'S' if D % 8 == 1 else 'I'
    return 'S' if leg(D, p) == 1 else 'I'

def analyse(avec, pmax=400):
    k = len(avec)
    T = [t for t in itertools.product((0, 1), repeat=k) if any(t)]
    Dt = {}
    for t in T:
        prod = 1
        for ti, a in zip(t, avec):
            if ti: prod *= a
        Dt[t] = sqf(prod)
    assert all(v != 1 for v in Dt.values()), "dependent generators"
    G = list(itertools.product((0, 1), repeat=k))
    dot = lambda e, t: sum(x * y for x, y in zip(e, t)) % 2
    c = tuple(1 if a < 0 else 0 for a in avec)
    out = []
    for p in primerange(2, pmax + 1):
        unr = [t for t in T if behaviour(Dt[t], p) != 'R']
        spl = [t for t in T if behaviour(Dt[t], p) == 'S']
        I = [e for e in G if all(dot(e, t) == 0 for t in unr)]
        D = [e for e in G if all(dot(e, t) == 0 for t in spl)]
        f = len(D) // len(I)
        if c in I:
            out.append((p, 'ram', 2 if p == 2 else 3, len(I), f))
        elif c in D:
            q = p ** (f // 2)
            out.append((p, 'inert', KAPPA.get(q, f'kappa({q})'), len(I), f))
    return out

fields = {
    'Q(sqrt-1) [Q^2]': [-1], 'Q(sqrt-2)': [-2], 'Q(sqrt-3)': [-3], 'Q(sqrt-5)': [-5], 'Q(sqrt-7)': [-7], 'Q(sqrt-6)': [-6],
    'Q(zeta8)=Q(sqrt2)^2': [-1, 2], 'Q(zeta12)=Q(sqrt3)^2': [-1, 3], 'Q(sqrt7)^2': [-1, 7],
    'Moser M=Q(sqrt-3,sqrt-11)': [-3, -11], 'M(i)=Q(sqrt3,sqrt11)^2': [-1, 3, 11],
    'Q(sqrt2,sqrt3)^2': [-1, 2, 3],
    'Polymath Q(sqrt-3,sqrt-11,sqrt-15)': [-3, -11, -15],
    'EI Q(sqrt-3,sqrt-11,sqrt-247)': [-3, -11, -247],
    'Heegner compositum': [-3, -11, -19, -43, -67, -163],
    'Frac Z[w1,w2,w3]': [-3, -7, -11],
    'Frac Z[w1,w3,w4,w5]': [-3, -11, -15, -19],
    'Frac Z[w1,w2,w3,w4]': [-3, -7, -11, -15],
    'Q(sqrt-3,sqrt-11,sqrt-23) [w6]': [-3, -11, -23],
    'Q(sqrt-3,sqrt-7,sqrt-15)': [-3, -7, -15],
    'Q(sqrt-3,sqrt-15) [w1,w4]': [-3, -15],
}
for N in range(1, 50):
    m = sqf(-(4 * N - 1))
    if m != -3:
        fields[f'L_{N}=Q(sqrt-3,sqrt{m}) (N mod 3={N%3})'] = [-3, m]
for name, av in fields.items():
    res = analyse(av)
    good = [r for r in res if isinstance(r[2], int)]
    best = min(good, key=lambda r: r[2]) if good else None
    first = res[:4]
    unram = all(r[1] != 'ram' for r in res)
    print(f"{name:45s} best: {('chi<=%d via p=%d (%s)' % (best[2], best[0], best[1])) if best else 'none with known kappa'} | first c-stable primes: {[(r[0], r[1], r[2]) for r in first]}")
