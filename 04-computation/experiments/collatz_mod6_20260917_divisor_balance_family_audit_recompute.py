#!/usr/bin/env python3
"""Independent recomputation for the divisor_balance_family lane (audit, recompute lens).

Does NOT import the explorer's script.  Different algorithms throughout:
  * primes: bytearray sieve (not numpy boolean slicing on the same array shape);
  * p^2qr numbers: generated DIRECTLY by enumerating triples (p,q,r), not via omega/Omega sieves;
  * F, S, U on 2..200000: computed literally from divisor LISTS built by a multiples sieve,
    with squarefree/prime indicators from independent sieves (no profile formulas);
  * profile window: own enumerator, own definitions of F,S,U,S_k from the profile;
  * every check raises RuntimeError.
RAM well below 1 GB, runtime about 1-2 minutes.
"""
import math
import sys
import time
from itertools import combinations_with_replacement
from math import isqrt, prod

import numpy as np

T0 = time.time()


def req(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def say(*a):
    print(*a)
    sys.stdout.flush()


X = 10_000_000

# ---------------------------------------------------------------- primes (bytearray sieve)
sieve = bytearray([1]) * (X + 1)
sieve[0] = sieve[1] = 0
for i in range(2, isqrt(X) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(range(i * i, X + 1, i)))
is_prime = np.frombuffer(bytes(sieve), dtype=np.uint8).astype(bool)
primes = np.flatnonzero(is_prime)
plist = primes.tolist()
cum_prime = np.cumsum(is_prime, dtype=np.int64)
say(f"pi(10^7) = {len(plist)}   [t={time.time()-T0:.1f}s]")
req(len(plist) == 664579, "pi(10^7)")


def pi(y):
    return int(cum_prime[y]) if y >= 2 else 0


# ---------------------------------------------------------------- p^2 q r by direct triple enumeration
p2qr = np.zeros(X + 1, dtype=bool)
n_marks = 0
for p in plist:
    if p * p * 6 > X:
        break
    lim_q = isqrt(X // (p * p))          # need q < r, so q^2 < X/p^2
    for q in plist:
        if q > lim_q:
            break
        if q == p:
            continue
        rmax = X // (p * p * q)
        if rmax <= q:
            continue
        lo = np.searchsorted(primes, q, side="right")
        hi = np.searchsorted(primes, rmax, side="right")
        rs = primes[lo:hi]
        if q < p <= rmax:
            rs = rs[rs != p]
        if len(rs) == 0:
            continue
        vals = (p * p * q) * rs
        req(vals.max() <= X, "p2qr overflow")
        req(not p2qr[vals].any(), "duplicate p^2qr representation (should be unique)")
        p2qr[vals] = True
        n_marks += len(rs)
cum_p2qr = np.cumsum(p2qr, dtype=np.int64)
say(f"p^2qr numbers <= 10^7 generated directly: {n_marks} marks, {int(cum_p2qr[X])} distinct   [t={time.time()-T0:.1f}s]")
req(n_marks == int(cum_p2qr[X]), "each p^2qr has a unique representation")

cube = np.zeros(X + 1, dtype=bool)
for p in plist:
    if p ** 3 > X:
        break
    cube[p ** 3] = True
cum_cube = np.cumsum(cube, dtype=np.int64)

say()
say("x        primes   p^3   p^2qr   ratio    asym(P2)  asym/exact")
P2 = 0.45224742004106549850  # literature prime zeta P(2)
expected = {10**5: (9592, 14, 9346), 10**6: (78498, 25, 87338), 10**7: (664579, 47, 804249)}
for x in (10**5, 10**6, 10**7):
    a, b, c = int(cum_prime[x]), int(cum_cube[x]), int(cum_p2qr[x])
    asym = P2 * x * math.log(math.log(x)) / math.log(x)
    say(f"{x:<8} {a:>7} {b:>5} {c:>7}  {c/a:.4f}  {asym:>9.0f}  {asym/c:.4f}")
    req((a, b, c) == expected[x], f"density counts at {x}: {(a,b,c)}")
    cr = int(round(x ** (1 / 3)))
    while cr ** 3 > x:
        cr -= 1
    while (cr + 1) ** 3 <= x:
        cr += 1
    req(b == pi(cr), f"cube count = pi(floor(x^(1/3))) at {x}")
say("density table CONFIRMED against explorer (independent generation)")

# partial prime zeta
P2_partial = float(np.sum(1.0 / primes.astype(np.float64) ** 2))
say(f"P(2) partial over p<=10^7: {P2_partial:.10f}; literature {P2:.10f}; diff {P2-P2_partial:.2e} (tail bound 1/X = 1e-7)")
req(abs(P2 - P2_partial) < 1e-7, "P(2) partial sum")
naive = math.exp(math.exp(1 / P2))
say(f"naive crossover exp(exp(1/P2)) = {naive:.1f}")

# ---------------------------------------------------------------- crossover
diff = cum_p2qr - cum_prime
lead = diff > 0
first = int(np.argmax(lead))
req(lead[first], "no crossover")
say(f"first n with #p2qr > #primes: {first}, counts primes={int(cum_prime[first])} p2qr={int(cum_p2qr[first])}")
req(first == 145119 and int(cum_prime[first]) == 13432 and int(cum_p2qr[first]) == 13433, "crossover")
req(145119 == 3 * 13 * 61 ** 2 and p2qr[145119], "145119 factorization / shape")
flips = np.flatnonzero(lead[1:] != lead[:-1]) + 1
say(f"flip points of the predicate (#p2qr > #primes): {flips.tolist()} ({len(flips)} flips)")
req(flips.tolist() == [145119, 145121, 145132, 145133, 145138, 145139, 145148, 147709, 147725, 147727, 147908, 147919, 147925], "flip list")
last_nonlead = int(np.flatnonzero(~lead)[-1])
say(f"last n with #primes >= #p2qr: {last_nonlead}; lead at 10^7 = {int(diff[X])}")
req(last_nonlead == 147924 and int(diff[X]) == 139670, "tail lead")
# also: does the lead ever TIE (diff==0) after 147925? sign changes counted on strict lead; report ties/zero crossings too
zeros_after = int(np.count_nonzero(diff[147925:] == 0))
say(f"ties (diff==0) at or after 147925: {zeros_after}; min diff on [147925,10^7] = {int(diff[147925:].min())}")
say(f"factor ratio first/naive = {first/naive:.2f}   [t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- sandwich endpoint census, K = 10^6
K = 1_000_000
ks = np.arange(1, K + 1, dtype=np.int64)
L = 6 * ks - 1
R = 6 * ks + 1
say()
say("endpoint census over 6k, k<=10^6:")
res = {}
for name, arr in (("prime", is_prime), ("p^3", cube), ("p^2qr", p2qr)):
    l, r = int(arr[L].sum()), int(arr[R].sum())
    res[name] = (l, r)
    say(f"  {name:<6} left {l:>7}  right {r:>7}")
req(res == {"prime": (206502, 206345), "p^3": (21, 19), "p^2qr": (35061, 34876)}, f"census {res}")
sol = is_prime | cube | p2qr
both = int((sol[L] & sol[R]).sum())
twin = int((is_prime[L] & is_prime[R]).sum())
say(f"  both endpoints solutions: {both}; twin primes: {twin}; mixed: {both-twin}")
req((both, twin) == (55613, 37915), "both/twin")
# twin primes independent: pairs (p,p+2) with p+2 <= 6K+1, minus (3,5)
tp = int(np.count_nonzero(is_prime[3:6 * K] & is_prime[5:6 * K + 2]))  # p in [3, 6K-1], p+2 <= 6K+1
say(f"  twin prime pairs (p,p+2), 3<=p, p+2<=6*10^6+1: {tp} (includes (3,5), which is not a 6k+-1 pair)")
req(tp - 1 == twin, "twin count vs generic twin-prime count")
# cube side law by direct residues
cl = sorted(int(p) ** 3 for p in plist if p ** 3 <= 6 * K + 1 and (p ** 3) % 6 == 5)
cr_ = sorted(int(p) ** 3 for p in plist if p ** 3 <= 6 * K + 1 and (p ** 3) % 6 == 1)
req(all(round(c ** (1/3)) % 6 == 5 for c in cl) and all(round(c ** (1/3)) % 6 == 1 for c in cr_), "cube side law")
say(f"  cube endpoints left {len(cl)} first {cl[:4]}; right {len(cr_)} first {cr_[:4]}")
req(cl[:4] == [125, 1331, 4913, 12167] and cr_[:4] == [343, 2197, 6859, 29791], "first cubes")
say(f"  [t={time.time()-T0:.1f}s]")

# ---------------------------------------------------------------- literal definitions on 2..M via divisor lists
M = 200_000
divs = [[] for _ in range(M + 1)]
for d in range(2, M + 1):
    for m in range(d, M + 1, d):
        divs[m].append(d)
sqf = bytearray([1]) * (M + 1)
for p in plist:
    if p * p > M:
        break
    sqf[p * p::p * p] = bytearray(len(range(p * p, M + 1, p * p)))
isp_small = [bool(is_prime[n]) for n in range(M + 1)]
# exponent profile by trial division against the prime list (independent of divs)
def profile(n):
    out = []
    m = n
    for p in plist:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            out.append(e)
    if m > 1:
        out.append(1)
    return tuple(sorted(out, reverse=True))


def kfree(d, k):
    # d is k-free iff no p^k | d
    m = d
    for p in plist:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            if e >= k:
                return False
    return True


cells = [(a, b) for a in range(4) for b in range(4)]
cell_counts = {c: 0 for c in cells}
shape_counts = {"p": 0, "p^3": 0, "p^2qr": 0}
sk_counts = {k: 0 for k in range(1, 7)}
sk_profiles = {k: set() for k in range(1, 7)}
cell_profiles = {c: set() for c in cells}
for n in range(2, M + 1):
    D = divs[n]           # all divisors >= 2 including n
    proper = D[:-1]       # 1 < d < n
    F = len(proper)
    S = sum(1 for d in proper if sqf[d])
    U = sum(1 for d in proper if isp_small[d])
    prof = profile(n)
    if F == S + U:
        if isp_small[n]:
            shape_counts["p"] += 1
        elif cube[n]:
            shape_counts["p^3"] += 1
        elif p2qr[n]:
            shape_counts["p^2qr"] += 1
        else:
            raise RuntimeError(f"F=S+U solution of unexpected shape: {n}")
    else:
        req(not (isp_small[n] or cube[n] or p2qr[n]), f"shape but not solution at {n}")
    for (a, b) in cells:
        if F == a * S + b * U:
            cell_counts[(a, b)] += 1
            cell_profiles[(a, b)].add(prof)
    # k-free: only need to recompute S_k for k != 2 when n is not squarefree-ish; do all k for safety but cheaply
    for k in range(1, 7):
        if k == 2:
            Sk = S
        elif k == 1:
            Sk = 0
        else:
            Sk = sum(1 for d in proper if kfree(d, k))
        if F == Sk + U:
            sk_counts[k] += 1
            sk_profiles[k].add(prof)
say()
say(f"literal-definition audit on 2..{M} done   [t={time.time()-T0:.1f}s]")
say(f"F=S+U shape counts: {shape_counts}")
req(shape_counts == {"p": 17984, "p^3": pi(58), "p^2qr": int(cum_p2qr[M])}, "shape counts <= 2e5")
say("cell counts #N<=2e5 (literal F,S,U):")
expected_cells = {(0,0):17984,(0,1):63214,(0,2):68867,(0,3):22140,(1,0):121666,(1,1):36372,(1,2):18122,(1,3):21629,
                  (2,0):32827,(2,1):21840,(2,2):20334,(2,3):18493,(3,0):17992,(3,1):17989,(3,2):18042,(3,3):18565}
for c in cells:
    say(f"  {c}: {cell_counts[c]}  profiles seen: {sorted(cell_profiles[c], key=lambda p:(len(p),p))}")
req(cell_counts == expected_cells, f"cell counts {cell_counts}")
say("k-free counts #N<=2e5 (literal S_k):")
for k in range(1, 7):
    say(f"  k={k}: {sk_counts[k]}  profiles seen: {sorted(sk_profiles[k], key=lambda p:(len(p),p))}")
req([sk_counts[k] for k in range(1, 7)] == [63214, 36372, 25045, 21096, 19442, 18684], "S_k counts")

# ---------------------------------------------------------------- profile window: own enumerator, exps<=15, r<=9
def pF(pr):
    return prod(a + 1 for a in pr) - 2


def pS(pr):
    return 2 ** len(pr) - 1 - (1 if all(a == 1 for a in pr) else 0)


def pU(pr):
    return len(pr) - (1 if pr == (1,) else 0)


def pSk(pr, k):
    return prod(min(a + 1, k) for a in pr) - 1 - (1 if all(a < k for a in pr) else 0)


predicted = {
    (0,0): {(1,)}, (0,1): {(1,),(2,),(1,1)}, (0,2): {(1,),(3,),(2,1),(1,1,1)}, (0,3): {(1,),(4,),(3,1)},
    (1,1): {(1,),(3,),(2,1,1)}, (1,2): {(1,),(4,),(2,2)}, (1,3): {(1,),(5,),(2,2,1),(2,1,1,1,1)},
    (2,1): {(1,),(4,),(4,1),(2,2,1,1)}, (2,2): {(1,),(5,),(3,2),(5,1),(4,1,1,1)}, (2,3): {(1,),(6,),(6,1)},
    (3,0): {(1,),(4,)}, (3,1): {(1,),(5,)}, (3,2): {(1,),(6,),(4,2)}, (3,3): {(1,),(7,),(3,3,1),(7,1,1)},
}
def pred_member(c, pr):
    a, b = c
    if c == (1, 0):
        return all(e == 1 for e in pr) or pr == (2,)
    if c == (2, 0):
        return pr == (1,) or pr == tuple([3] + [1] * (len(pr) - 1))
    return pr in predicted[c]

hits = {c: set() for c in cells}
skhits = {k: set() for k in range(1, 7)}
nprof = 0
for r in range(1, 10):
    for comb in combinations_with_replacement(range(15, 0, -1), r):
        pr = tuple(comb)
        nprof += 1
        F, S, U = pF(pr), pS(pr), pU(pr)
        for c in cells:
            if F == c[0] * S + c[1] * U:
                hits[c].add(pr)
            req((F == c[0] * S + c[1] * U) == pred_member(c, pr), f"window cell {c} at {pr}")
        for k in range(1, 7):
            if F == pSk(pr, k) + U:
                skhits[k].add(pr)
say()
say(f"own window search: {nprof} profiles (r<=9, exponents<=15): all 16 cells match the claimed classification")
for k in range(1, 7):
    exp_k = {(1,), (2,), (1, 1)} if k == 1 else ({(1,), (k + 1,), (k, 1, 1)} | ({(k, 2)} if k >= 3 else set()))
    req(skhits[k] == exp_k, f"S_{k} window: {skhits[k]}")
say("k-free window: k=1..6 match {p,p^2,pq} / {p,p^(k+1),p^k qr} + {p^k q^2 for k>=3}")
# (2,0) via own factorization reasoning: check T(r)=2^(r+1) solutions r<=9 are exactly (3,1^(r-1))
for r in range(1, 10):
    fam = sorted(pr for pr in hits[(2, 0)] if len(pr) == r and pr != (1,))
    req(fam == [tuple([3] + [1] * (r - 1))], f"(2,0) r={r}: {fam}")
# Omega values per cell -- HOSTILE: is (2,0) Omega-injective?
say()
say("Omega values of solution profiles per cell (window):")
for c in cells:
    om = sorted(sum(p) for p in hits[c])
    inj = len(set(om)) == len(om)
    say(f"  {c}: {om} {'injective' if inj else 'COLLIDES'}")
om20 = [sum(p) for p in hits[(2, 0)]]
req(len(set(om20)) == len(om20), "(2,0) Omega values are pairwise distinct (Omega = r+2 determines the profile)")
say("HOSTILE: (2,0) F=2S solutions p^3 q_1..q_(r-1) have Omega = r+2, pairwise DISTINCT; so (1,1) is NOT the unique")
say("         prime-cube cell with pairwise distinct Omega -- it is the unique FINITE one. Explorer's .out labels (2,0) 'collides' only")
say("         because its injectivity flag is defined as finite AND injective.")

# ---------------------------------------------------------------- partitions per Omega
def partitions(n, m=None):
    m = n if m is None else m
    if n == 0:
        yield ()
        return
    for f in range(min(n, m), 0, -1):
        for rest in partitions(n - f, f):
            yield (f,) + rest

say()
for k in range(1, 9):
    parts = list(partitions(k))
    sols = [p for p in parts if pF(p) == pS(p) + pU(p)]
    say(f"  Omega={k}: p(k)={len(parts)}, solutions {sols}")
    req(len(parts) == [1,2,3,5,7,11,15,22][k-1], "partition count")
    req(sols == ({1:[(1,)],3:[(3,)],4:[(2,1,1)]}.get(k, [])), f"solutions at Omega={k}")
d22 = pF((2,2)) - pS((2,2)) - pU((2,2)); d31 = pF((3,1)) - pS((3,1)) - pU((3,1))
req((d22, d31) == (2, 1), "D at p^2q^2 / p^3q")
say(f"  D(p^2q^2)={d22}, D(p^3q)={d31}")

# ---------------------------------------------------------------- lattice identity: |Q| = 2^(r-1) for (2,1^(r-1)); F-S-U = |Q|-(r+1)
for r in range(1, 10):
    pr = tuple([2] + [1] * (r - 1))
    Q = prod(a + 1 for a in pr) - 2 ** r
    req(Q == 2 ** (r - 1), "Q size")
    req(pF(pr) - pS(pr) - pU(pr) == Q - (r + 1), "F-S-U = |Q|-(r+1)")
    req((Q == r + 1) == (r == 3), "iff r=3")
say("lattice identity F-S-U = |Q|-(r+1), |Q|=2^(r-1) at tight profile, equality iff r=3: confirmed r<=9")
# general nonsquarefree profiles: F-S-U = |Q| - (r+1) where |Q| = prod(a_i+1) - 2^r
for r in range(1, 7):
    for comb in combinations_with_replacement(range(6, 0, -1), r):
        pr = tuple(comb)
        if max(pr) >= 2:
            req(pF(pr) - pS(pr) - pU(pr) == (prod(a + 1 for a in pr) - 2 ** r) - (r + 1), f"lattice identity at {pr}")
say("lattice identity holds for every nonsquarefree profile r<=6, exps<=6")

# ---------------------------------------------------------------- support bound sanity for alpha in {2,3} (Lemma A numerics)
def T(a, b, r):
    return 2 + a * (2 ** r - 1) + b * r
for a in (2, 3):
    for b in range(4):
        if (a, b) == (2, 0):
            continue
        for r in range(4, 40):
            c = b * r + 2 - a
            req(c != 0 and abs(c) < 2 ** r, "c_r nonzero and small")
            v2T = (T(a,b,r) & -T(a,b,r)).bit_length() - 1
            v2c = (abs(c) & -abs(c)).bit_length() - 1
            req(v2T == v2c, f"v2(T)=v2(c_r) at {(a,b,r)}")
            req(v2c <= math.log2(3 * r + 1), "v2 bound")
say("2-adic step of the support bound verified numerically for alpha in {2,3}, r=4..39")
# any factorization of T(a,b,r) into r factors >=2 with a factor >=3 for r=8..14?  (independent multiplicative-partition search)
def mult_parts(T, r, lo=2):
    if r == 1:
        return [[T]] if T >= lo else []
    out = []
    f = lo
    while f ** r <= T:
        if T % f == 0:
            for rest in mult_parts(T // f, r - 1, f):
                out.append([f] + rest)
        f += 1
    return out
for a in (2, 3):
    for b in range(4):
        if (a, b) == (2, 0):
            continue
        for r in range(8, 15):
            req(all(max(f) < 3 for f in mult_parts(T(a,b,r), r)), f"unexpected nonsquarefree solution {(a,b,r)}")
say("no nonsquarefree solutions for alpha in {2,3} at r=8..14 (own multiplicative-partition search)")

say()
say(f"ALL AUDIT CHECKS PASSED   [t={time.time()-T0:.1f}s]")
