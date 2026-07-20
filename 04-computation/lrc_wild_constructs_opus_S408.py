"""
opus-2026-07-19-S408 (HYP-8010): four NEW constructs, defined and probed.

W1  THETA-INTERPOLATED LONELINESS: M_theta(V) := max_t min_i v_i^theta * ||v_i t||.
    theta = 0 is LRC; theta = 1 weights by speed (the Kakeya-tube flavor of kp's
    lattice-tube thread).  Probe: does the AP remain the covering-corpus minimizer
    at theta = 1?  (Exact at theta in {0,1}.)
W2  DISPLACEMENT PARTITION: lambda(V) := sorted multiset {q*||v t*||} at the
    (theta=0) maximizer t* = a/q.  Observation target: AP13 gives the DOUBLED
    STAIRCASE 1,1,2,2,...,6,6,7 -- a Young-diagram bridge to the home mandate's
    staircase delta_{n-2}.
W3  EXCLUSION SIGNATURE: sigma(V) := {(q_j, delta_j) : convergent denominators of
    t* with delta_j < qM} -- provably disjoint from V (they would beat M); the
    S407 ghost loop says extremal ghost families have sigma = {(g, 1)} with the
    duty-payer K*g in V, while D=1 families have sigma = EMPTY.
W4  REGULARITY SIGNAL: c3(t) over rationals; hypothesis "the tournament sees the
    witnesses": witness times of tight families lie in argmax(c3).
"""
from math import gcd, comb
from fractions import Fraction

def scan(V, qmax, theta):
    best, wit = Fraction(0), None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = None
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                val = Fraction(r, q) * (v if theta else 1)
                if m is None or val < m: m = val
            if m > best: best, wit = m, (a, q)
    return best, wit

def convs(a, q):
    num, den, quots = q, a, []
    while den: quots.append(num//den); num, den = den, num % den
    hp, h = 1, quots[0]
    out = [1, h]
    for c in quots[1:]:
        hp, h = h, c*h + hp
        out.append(h)
    return sorted(set(out))

def c3_at(V, a, q):
    pos = [Fraction(0)] + [Fraction((v*a) % q, q) for v in V]
    n = len(pos); out = [0]*n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            d = (pos[j]-pos[i]) % 1
            if (0 < d < Fraction(1,2)) or ((d == 0 or d == Fraction(1,2)) and i < j):
                out[i] += 1
    return comb(n,3) - sum(comb(o,2) for o in out)

FAMS = [("AP13", list(range(1,14))), ("GW", list(range(1,12))+[13,24]),
        ("DW", list(range(1,13))+[182]), ("ladder3", list(range(1,12))+[13,36]),
        ("{1..12,14}", list(range(1,13))+[14]), ("{1..12,15}", list(range(1,13))+[15]),
        ("primes13", [2,3,5,7,11,13,17,19,23,29,31,37,41]),
        ("cluster20", list(range(20,33)))]

print("W1: theta = 0 vs theta = 1 (exact):")
r0, r1 = [], []
for name, V in FAMS:
    M0, w0 = scan(V, 220, 0)
    M1, w1 = scan(V, 220, 1)
    r0.append((M0, name)); r1.append((M1, name))
    print(f"  {name}: M_0 = {M0} at {w0}; M_1 = {M1} = {float(M1):.3f} at {w1}")
print("  theta=0 ranking:", [n for _, n in sorted(r0)])
print("  theta=1 ranking:", [n for _, n in sorted(r1)])

print("\nW2: displacement partitions at the theta=0 maximizer:")
for name, V in FAMS[:4]:
    M0, (a, q) = scan(V, 220, 0)
    lam = sorted(min((v*a) % q, q-(v*a) % q) for v in V)
    print(f"  {name} (t*={a}/{q}): lambda = {lam}")

print("\nW3: exclusion signatures:")
for name, V in FAMS[:6]:
    M0, (a, q) = scan(V, 220, 0)
    mt = int(M0*q) if (M0*q).denominator == 1 else None
    sig = []
    for cd in convs(a, q):
        if cd >= q: continue
        d = min((cd*a) % q, q-(cd*a) % q)
        if mt and d < mt:
            duty = [v for v in V if v % cd == 0]
            sig.append((cd, d, f"in V: {cd in V}", f"duty-payers: {duty}"))
    print(f"  {name}: sigma = {sig if sig else 'EMPTY'}")

print("\nW4: does the tournament see the witnesses?  argmax c3 vs witness set:")
for name, V, qm in [("AP13", list(range(1,14)), 60), ("GW", list(range(1,12))+[13,24], 60),
                    ("DW", list(range(1,13))+[182], 200)]:
    best_c3, argmax, wit_c3 = -1, [], []
    thr = Fraction(1, 14)
    for q in range(2, qm+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            c = c3_at(V, a, q)
            if c > best_c3: best_c3, argmax = c, [(a,q)]
            elif c == best_c3 and len(argmax) < 30: argmax.append((a,q))
            if all(Fraction(min((v*a)%q, q-(v*a)%q), q) >= thr for v in V):
                wit_c3.append(((a,q), c))
    wset = set(t for t,_ in wit_c3)
    inmax = [t for t in wset if t in set(argmax)]
    print(f"  {name}: max c3 = {best_c3}; #argmax(<=30 shown) = {len(argmax)}; "
          f"witnesses found = {len(wit_c3)} with c3 values {sorted(set(c for _,c in wit_c3))}; "
          f"witnesses in argmax: {len(inmax)}/{len(wset)}")
