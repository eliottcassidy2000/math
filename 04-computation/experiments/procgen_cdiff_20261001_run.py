#!/usr/bin/env python3
"""procgen cdiff lane, 2026-10-01: drop multiplicities of the Syracuse map.

Re-verifies every claim of
  05-knowledge/results/procgen_cdiff_20261001_drop_multiplicities.md
Prints to stdout only; ends with ALL CHECKS PASSED.

Sections
  A  independent re-verification of THM-4527 (opus S15, twentieth note, Theorems 1-2)
  B  the divisor problem over 2^v - 3 (direction D74): exact densities, lattice, maximal order
  C  the summatory function of K
  D  k-step drops, consecutive drops, cycles
  E  other q (direction D75): the drop law of x -> oddpart(qx +- 1)
  F  the drop identity mod 3^k

Resources: about 1 minute on one core; peak RSS is printed (about 300 MB: one numpy int8 array of 1e8 + 1 entries).
Independence: the main density computation (section B) factors only the small pairwise gcds of the
moduli by trial division (no sympy); sympy is used only for a second, cross-checking structure
computation and for multiplicative orders.
"""
import sys
import time
from math import gcd, comb, log2
from fractions import Fraction as Fr
from collections import Counter, defaultdict
from itertools import product

import numpy as np
import mpmath as mp

T0 = time.time()
NCHECK = 0


def check(cond, msg):
    global NCHECK
    if not cond:
        print("FAIL:", msg)
        sys.exit(1)
    NCHECK += 1
    print("PASS", msg)


def info(msg):
    print("   ", msg)


def peak_mb():
    import resource
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2 ** 20 if sys.platform == "darwin" else r / 2 ** 10


def v2(x):
    return (x & -x).bit_length() - 1


def syr(A):
    """Syracuse map on an odd integer A of either sign: (oddpart(3A+1), v_2(3A+1))."""
    x = 3 * A + 1
    v = v2(x)
    return x >> v, v


def Tq(q, s, A):
    """oddpart(qA + s) and its valuation."""
    x = q * A + s
    v = v2(x)
    return x >> v, v


def Mv(v):
    return (1 << v) - 3


def polymul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r


def trial_factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def structure(mods, factor_small=trial_factor):
    """Coprime decomposition of a list of positive odd moduli: the primes shared by two moduli
    (found from pairwise gcds only), each modulus' exponents at those primes, and its private part."""
    n = len(mods)
    shared = set()
    for i in range(n):
        for j in range(i + 1, n):
            g = gcd(mods[i], mods[j])
            if g > 1:
                shared |= set(factor_small(g))
    shared = sorted(shared)
    exps, priv = [], []
    for m_ in mods:
        e, r = {}, m_
        for p in shared:
            while r % p == 0:
                e[p] = e.get(p, 0) + 1
                r //= p
        exps.append(e)
        priv.append(r)
    for i in range(n):
        for j in range(i + 1, n):
            assert gcd(priv[i], priv[j]) == 1
    return shared, exps, priv


def exact_law(shared, exps, priv):
    """Exact distribution of #{i : mods[i] | omega} for omega Haar-random in Zhat.
    Enumerates the states of the shared primes; given the state, the private events are independent."""
    Amax = {p: max(e.get(p, 0) for e in exps) for p in shared}
    free = [i for i, e in enumerate(exps) if not e]
    dep = [i for i, e in enumerate(exps) if e]
    base = [Fr(1)]
    for i in free:
        q_ = Fr(1, priv[i])
        base = polymul(base, [1 - q_, q_])
    total = [Fr(0)]
    nstates = 0
    for x in product(*[range(Amax[p] + 1) for p in shared]):
        nstates += 1
        pr = Fr(1)
        for p, xp in zip(shared, x):
            pr *= Fr(1, p ** xp) * ((1 - Fr(1, p)) if xp < Amax[p] else 1)
        st = dict(zip(shared, x))
        poly = [pr]
        for i in dep:
            if all(st[p] >= a for p, a in exps[i].items()):
                q_ = Fr(1, priv[i])
                poly = polymul(poly, [1 - q_, q_])
        if len(poly) > len(total):
            total += [Fr(0)] * (len(poly) - len(total))
        for j, c in enumerate(poly):
            total[j] += c
    dist = polymul(base, total)
    while len(dist) > 1 and dist[-1] == 0:
        dist.pop()
    return dist, nstates


MAXCELLS = 1_200_000  # memory guard for the variable-elimination engine


def polymul_arr(A, B, D):
    """Product of polynomial-valued arrays (last axis = coefficients), truncated at degree D."""
    shp = np.broadcast_shapes(A.shape[:-1], B.shape[:-1]) + (D + 1,)
    n = 1
    for s_ in shp:
        n *= s_
    assert n <= MAXCELLS * (D + 1), ("table too large", shp)
    R = np.zeros(shp)
    for i in range(D + 1):
        a = A[..., i:i + 1]
        if not np.any(a):
            continue
        R[..., i:] += a * B[..., :D + 1 - i]
    return R


def law_events(events, D=12):
    """Law (float64) of #{i : x = r mod M_i for some r in R_i, counted with multiplicity} for Haar-random x.
    events: list of (M_i, {residue: multiplicity}).  Conditions on the l-adic digits of x at the primes l
    shared by two moduli; given them the events are independent; the digits are summed out by variable
    elimination (one digit value at a time)."""
    from sympy import factorint
    facs = [factorint(M) for M, _ in events]
    occ = Counter()
    for f in facs:
        for l in f:
            occ[l] += 1
    shared_ = sorted(l for l, c in occ.items() if c >= 2)
    const = np.zeros(D + 1)
    const[0] = 1.0
    factors = []
    for (M, R), f in zip(events, facs):
        S = [l for l in shared_ if l in f]
        sp = 1
        for l in S:
            sp *= l ** f[l]
        pv = M // sp
        if not S:
            poly = np.zeros(D + 1)
            poly[0] = 1.0
            for r, mlt in R.items():
                poly[0] -= 1.0 / M
                if mlt <= D:
                    poly[mlt] += 1.0 / M
            const = polymul_arr(const, poly, D)
            continue
        scope = tuple((l, i) for l in S for i in range(f[l]))
        tab = np.zeros(tuple(l for (l, i) in scope) + (D + 1,))
        tab[..., 0] = 1.0
        for r, mlt in R.items():
            idx = []
            for l in S:
                y = r % (l ** f[l])
                for i in range(f[l]):
                    idx.append(y % l)
                    y //= l
            idx = tuple(idx)
            tab[idx + (0,)] -= 1.0 / pv
            if mlt <= D:
                tab[idx + (mlt,)] += 1.0 / pv
        factors.append((scope, tab))
    variables = sorted({var for S, _ in factors for var in S})
    while variables:
        best_ = None
        for var in variables:
            sc = set()
            for S, _ in factors:
                if var in S:
                    sc |= set(S)
            size = 1
            for (l, i) in sc:
                size *= l
            if best_ is None or size < best_[0]:
                best_ = (size, var, sorted(sc))
        size, var, sc = best_
        sc = tuple(sc)
        keep, parts = [], []
        for S, A in factors:
            if var in S:
                shape = [w[0] if w in S else 1 for w in sc] + [D + 1]
                perm = [S.index(w) for w in sc if w in S] + [len(S)]
                parts.append(np.transpose(A, perm).reshape(shape))
            else:
                keep.append((S, A))
        ax = sc.index(var)
        newA = None
        for val in range(var[0]):
            prod = None
            for A2 in parts:
                sl = A2[(slice(None),) * ax + ((val if A2.shape[ax] > 1 else 0),)]
                prod = sl if prod is None else polymul_arr(prod, sl, D)
            newA = prod.copy() if newA is None else newA + prod
        newA = newA / var[0]
        newS = tuple(w for w in sc if w != var)
        if newS:
            keep.append((newS, np.broadcast_to(newA, tuple(w[0] for w in newS) + (D + 1,)).copy()))
        else:
            const = polymul_arr(const, newA.reshape(D + 1), D)
        factors = keep
        variables.remove(var)
    return const


def fl(c, n=12):
    return mp.nstr(mp.mpf(c.numerator) / c.denominator, n)


mp.mp.dps = 40

# ---- heavy engine computations first (keeps the peak RSS low); their checks are printed in section D10 ----
def cvec(word, q=3):
    k = len(word)
    c, ps = 0, 0
    for i, v in enumerate(word):
        c += q ** (k - 1 - i) * (1 << ps)
        ps += v
    return c, ps


def words_k(k, p):
    from itertools import combinations
    for cuts in combinations(range(1, p), k - 1):
        parts, prev = [], 0
        for c_ in cuts:
            parts.append(c_ - prev)
            prev = c_
        parts.append(p - prev)
        yield tuple(parts)


def kstep_events(k, P, side=+1):
    """Descending side (side=+1, D >= 0): x = 2*3^k*D, events x = -c(v) mod (2^p - 3^k) for words with
    2^p > 3^k, p <= P.  Ascending side (side=-1, D < 0): the words with 2^p < 3^k, modulus 3^k - 2^p."""
    import math
    plo = int(math.floor(k * math.log2(3)))
    ps = range(plo + 1, P + 1) if side > 0 else range(k, plo + 1)
    ev = []
    for p in ps:
        Mk = abs(2 ** p - 3 ** k)
        R = {}
        for w in words_k(k, p):
            r = (-cvec(w)[0]) % Mk
            R[r] = R.get(r, 0) + 1
        ev.append((Mk, R))
    return ev


def pair_moment(ev):
    """Exact E[mu(mu-1)] = sum over ordered pairs of distinct words of P(both hit)."""
    tot = Fr(0)
    for i, (M1, R1) in enumerate(ev):
        tot += Fr(sum(m * (m - 1) for m in R1.values()), M1)
        for (M2, R2) in ev[i + 1:]:
            g = gcd(M1, M2)
            c1, c2 = Counter(), Counter()
            for r, m in R1.items():
                c1[r % g] += m
            for r, m in R2.items():
                c2[r % g] += m
            tot += Fr(2 * sum(c1[c] * c2[c] for c in c1 if c in c2), M1 // g * M2)
    return tot


def engine_payload():
    t_ = time.time()
    out = {"k1": law_events([((1 << v) - 3, {0: 1}) for v in range(3, 65)]).tolist()}
    for k_ in (2, 3):
        ev_ = kstep_events(k_, 50)
        pm_ = pair_moment(ev_)
        out[str(k_)] = [law_events(ev_, D=12).tolist(), "%d/%d" % (pm_.numerator, pm_.denominator)]
    out["asc2"] = law_events(kstep_events(2, 0, side=-1), D=6).tolist()
    out["asc3"] = law_events(kstep_events(3, 0, side=-1), D=6).tolist()
    out["time"] = time.time() - t_
    out["rss"] = peak_mb()
    return out


if "--engine" in sys.argv:  # child process: heavy elimination-engine computations, isolated for memory
    import json
    print(json.dumps(engine_payload()))
    sys.exit(0)
else:
    import json
    import subprocess
    _res = subprocess.run([sys.executable, __file__, "--engine"], capture_output=True, text=True, check=True)
    _pl = json.loads(_res.stdout)
    ENG = {"k1": np.array(_pl["k1"]), "asc2": np.array(_pl["asc2"]), "asc3": np.array(_pl["asc3"]),
           "time": _pl["time"], "rss": _pl["rss"]}
    for _k in (2, 3):
        ENG[_k] = (np.array(_pl[str(_k)][0]), Fr(_pl[str(_k)][1]))

# ======================================================================================
print("=== A. independent re-verification of THM-4527 (S15 Theorems 1-2) ===")
LIM = 10 ** 6


def F(M_):
    return (syr(2 * M_ - 1)[0] + 1) // 2


def KN(N):
    A = 4 * N - 3
    return (A - syr(A)[0]) // 2


check(all(F(2 * N) == 3 * N for N in range(1, LIM // 2)), "A1 F(2N) = 3N for 1 <= N < 5e5")
check(all(F(4 * j + 1) == 3 * j + 1 for j in range(LIM // 4)), "A1 F(4j+1) = 3j+1 for 0 <= j < 2.5e5")
check(all(F(4 * n - 1) == F(n) for n in range(1, LIM // 4)), "A1 F(4n-1) = F(n) (the microcosm) for 1 <= n < 2.5e5")
check([KN(N) for N in range(1, 12)] == [0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5],
      "A2 owner's K_1..K_11 = 0,2,1,4,2,10,3,9,4,15,5 (K_N = (A - S(A))/2 at A = 4N-3)")
check([syr(a)[0] for a in (1, 5, 9, 13)] == [1, 1, 7, 5], "A2 owner's examples S(1)=1, S(5)=1, S(9)=7, S(13)=5")
check(all(F(2 * N - 1) == 2 * N - 1 - KN(N) for N in range(1, LIM // 2)), "A2 F(2N-1) = 2N-1-K_N for N < 5e5")
check(all(KN(2 * j + 1) == j for j in range(LIM // 4))
      and all(KN(4 * m) == 5 * m - 1 and KN(4 * m - 2) == KN(m) + 6 * m - 4 for m in range(1, LIM // 8)),
      "A2 K_{2j+1} = j, K_{4m} = 5m-1, K_{4m-2} = K_m + 6m - 4 (2-regular recursion)")
ok = True
for A in range(-2 * 10 ** 6 + 1, 2 * 10 ** 6, 2):
    s, v = syr(A)
    if 6 * ((A - s) // 2) + 1 != ((1 << v) - 3) * s:
        ok = False
        break
check(ok, "A3 6K(A) + 1 = (2^v - 3) S(A) for every odd A with |A| < 2e6 (both signs)")

D1 = 10 ** 6
cnt = np.zeros(D1 + 1, dtype=np.int32)
for A in range(1, 8 * D1 + 2, 4):  # A = 1 mod 4 <=> v >= 2 <=> K(A) >= 0 ; K >= (A-1)/8
    s, v = syr(A)
    K = (A - s) // 2
    if K <= D1:
        cnt[K] += 1
mf = np.ones(D1 + 1, dtype=np.int32)  # v = 2 always
v = 3
while Mv(v) <= 6 * D1 + 1:
    m_ = Mv(v)
    mf[(-pow(6, -1, m_)) % m_::m_] += 1
    v += 1
check(np.array_equal(cnt, mf), "A4 multiplicity of m in the owner's K = #{v >= 2 : (2^v-3) | 6m+1}, all 0 <= m <= 1e6 (brute force over A <= 8e6+1)")
h5 = Counter(cnt[:10 ** 5 + 1].tolist())
check(dict(h5) == {1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6}, "A4 S15's histogram on [0,1e5] {1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6} reproduced")
neg = Counter()
for A in range(3, 4 * 10 ** 5, 4):
    s, v = syr(A)
    neg[(A - s) // 2] += 1
check(sorted(neg) == list(range(-10 ** 5, 0)) and set(neg.values()) == {1},
      "A4 every negative m in [-1e5,-1] is a drop exactly once (the ascents A = 4N-1, K = -N)")
# two-sided count over all odd A in Z \ {0}: every d occurs 2 + #{v>=3 : M_v | 6d+1} times
Y = 20000
two = Counter()
for A in range(-8 * Y - 9, 8 * Y + 10, 2):
    s, v = syr(A)
    K = (A - s) // 2
    if -Y <= K <= Y:
        two[K] += 1


def Ncount(n):
    n = abs(n)
    c, v = 0, 3
    while Mv(v) <= n:
        if n % Mv(v) == 0:
            c += 1
        v += 1
    return c


check(all(two[d] == 2 + Ncount(6 * d + 1) for d in range(-Y, Y + 1)),
      "A5 over all odd A in Z, every d in [-2e4,2e4] is a drop exactly 2 + #{v>=3 : (2^v-3) | 6d+1} times (unit branches: A = 8d+1 and A = -4d-1)")
check(all(syr(8 * d + 1) == (6 * d + 1, 2) and syr(-4 * d - 1) == (-6 * d - 1, 1)
          and (8 * d + 1 - (6 * d + 1)) // 2 == d and (-4 * d - 1 - (-6 * d - 1)) // 2 == d for d in range(-10 ** 5, 10 ** 5)),
      "A5 the two unit preimages of every d: S(8d+1) = 6d+1 (v=2), S(-4d-1) = -6d-1 (v=1), both with drop d")
# the same in the owner's labels M = (A+1)/2 over all integers M: differences M - F(M)
lab = {0: Counter(), 1: Counter(), 3: Counter()}
Yl = 5000
for M_ in range(-16 * Yl, 16 * Yl + 1):
    A = 2 * M_ - 1
    dlt = M_ - (syr(A)[0] + 1) // 2
    if -Yl <= dlt <= Yl:
        lab[0 if M_ % 2 == 0 else M_ % 4][dlt] += 1
check(all(lab[0][d] == 1 and lab[1][d] == 1 and lab[3][d] == Ncount(6 * d + 1) for d in range(-Yl, Yl + 1)),
      "A5 in labels M (all integers): every d = M - F(M) occurs exactly once with M even, once with M = 1 mod 4, and #{v>=3 : (2^v-3) | 6d+1} times with M = 3 mod 4 (the microcosm), |d| <= 5000")
rho1 = mp.nsum(lambda v: 1 / (mp.mpf(2) ** v - 3), [2, mp.inf])
info("mean multiplicity sum_{v>=2} 1/(2^v-3) = %s" % mp.nstr(rho1, 25))
check(abs(rho1 - mp.mpf("1.3436734331817690185444830")) < mp.mpf(10) ** -24, "A6 mean multiplicity = 1.34367343318176901854448...")

# ======================================================================================
info("peak RSS so far: %.0f MB, elapsed %.0f s" % (peak_mb(), time.time() - T0))
print("=== B. the divisor problem over 2^v - 3 (D74) ===")
V = 64
vs = list(range(3, V + 1))
mods = [Mv(v) for v in vs]
shared, exps, priv = structure(mods)  # trial division of pairwise gcds only
gvals = sorted({gcd(mods[i], mods[j]) for i in range(len(mods)) for j in range(i + 1, len(mods))} - {1})
info("pairwise gcds > 1 among 2^v-3, 3 <= v <= 64: %s" % gvals)
check(shared == [5, 11, 13, 19, 23, 29, 37, 47, 71, 431],
      "B1 for v <= 64 exactly ten primes divide two of the 2^v - 3: 5, 11, 13, 19, 23, 29, 37, 47, 71, 431 (private parts pairwise coprime)")
Lall = 1
for m_ in mods:
    Lall = Lall * m_ // gcd(Lall, m_)
prod_bits = sum(log2(m_) for m_ in mods)
info("overlap: log2 prod M_v - log2 lcm M_v (3 <= v <= 64) = %.1f of %.1f bits (%.1f%%)" % (prod_bits - log2(Lall), prod_bits, 100 * (prod_bits - log2(Lall)) / prod_bits))
dep_v = [vs[i] for i, e in enumerate(exps) if e]
info("moduli involving shared primes (%d of %d): %s" % (len(dep_v), len(vs), dep_v))
t = time.time()
law64, ns = exact_law(shared, exps, priv)
info("exact law of N_64 computed over %d shared-prime states in %.1fs" % (ns, time.time() - t))
check(sum(law64) == 1 and all(c >= 0 for c in law64), "B2 the law of N_64 is a probability vector (exact rationals)")
EN = sum(j * c for j, c in enumerate(law64))
check(EN == sum(Fr(1, m_) for m_ in mods), "B2 exact first moment: E[N_64] = sum_{v=3}^{64} 1/(2^v-3)")
EN2 = sum(Fr(j * (j - 1), 2) * c for j, c in enumerate(law64))
S2 = sum(Fr(gcd(mods[i], mods[j]), mods[i] * mods[j]) for i in range(len(mods)) for j in range(i + 1, len(mods)))
check(EN2 == S2, "B2 exact second binomial moment: E[C(N_64,2)] = sum_{v<w} 1/lcm(M_v, M_w)")
# cross-check the structure with sympy factorizations (independent route to the shared primes)
try:
    from sympy import factorint, n_order
    occ = defaultdict(int)
    for m_ in mods:
        for p in factorint(m_):
            occ[p] += 1
    check(sorted(p for p, c in occ.items() if c >= 2) == shared, "B2 sympy factorizations give the same ten shared primes")
except ImportError:
    n_order = None
eps64 = Fr(1, 2 ** 63)  # sum_{v>64} 1/(2^v-3) < sum 2^(1-v) = 2^-63
check(sum(Fr(1, Mv(v)) for v in range(65, 400)) < eps64, "B2 truncation error sum_{v>64} 1/(2^v-3) < 2^-63 = 1.08e-19")
delta = [law64[j] for j in range(10)]
for j in range(10):
    info("delta_%d = density{d : m(d) = %d} = %s  (+- 1.1e-19)" % (j + 1, j + 1, fl(delta[j], 22)))
# empirical counts on [0, 1e8]
X = 10 ** 8
m8 = np.ones(X + 1, dtype=np.int8)
v = 3
while Mv(v) <= 6 * X + 1:
    m_ = Mv(v)
    m8[(-pow(6, -1, m_)) % m_::m_] += 1
    v += 1
mmax = int(m8.max())
h8 = [0] * (mmax + 1)
first = {}
CH = 10 ** 7
for c0 in range(0, X + 1, CH):  # chunked, so that no full-size temporary is created
    ch = m8[c0:c0 + CH]
    for j in range(1, mmax + 1):
        h8[j] += int(np.count_nonzero(ch == j))
        if j >= 2 and j not in first and np.any(ch >= j):
            first[j] = c0 + int(np.argmax(ch >= j))
info("counts on [0,1e8]: %s" % {j: int(c) for j, c in enumerate(h8) if c})
check(all(abs(h8[j] / (X + 1) - float(delta[j - 1])) < 3e-7 for j in range(1, 7)),
      "B3 empirical frequencies of m(d) = 1..6 on [0,1e8] agree with the exact densities to 3e-7")
info("first d with m(d) >= j: %s" % first)
del m8
# independent model for comparison
pind = [mp.mpf(1)]
for m_ in mods:
    q_ = mp.mpf(1) / m_
    new = [mp.mpf(0)] * (len(pind) + 1)
    for i, c in enumerate(pind):
        new[i] += c * (1 - q_)
        new[i + 1] += c * q_
    pind = new
info("independent-events model P(m=1..5): %s" % [mp.nstr(c, 8) for c in pind[:5]])
check(float(delta[0]) - float(pind[0]) > 0.005 and float(delta[3]) / float(pind[3]) > 1.5,
      "B3 the gcd structure matters: P(m=1) exceeds the independent model by 0.0059, P(m=4) by a factor 1.59")
EN2full = sum(j * j * c for j, c in enumerate(law64))
info("Var(m) = Var(N_64) = %s (mean %s)" % (fl(EN2full - EN * EN, 15), fl(EN, 15)))
EU = {u: sum(c * u ** j for j, c in enumerate(law64)) for u in (2, 3)}
info("E[2^N_64] = %s, E[3^N_64] = %s (finite for u < 4 in the limit, Lemma 2.2)" % (fl(EU[2], 10), fl(EU[3], 10)))
# divisibility lattice
ok1 = all(gcd(Mv(v), Mv(w)) % 3 != 0 and ((1 << (w - v)) - 1) % gcd(Mv(v), Mv(w)) == 0
          for v in range(3, 201) for w in range(v + 1, 201))
check(ok1, "B4 gcd(2^v-3, 2^w-3) divides 2^(w-v) - 1 and is prime to 6, 3 <= v < w <= 200")
if n_order is not None:
    ords = {v: n_order(2, Mv(v)) for v in range(3, 41)}
    ok2 = all((Mv(w) % Mv(v) == 0) == ((w - v) % ords[v] == 0) for v in range(3, 41) for w in range(v + 1, 401))
    check(ok2, "B4 M_v | M_w <=> ord_{M_v}(2) | w - v, for 3 <= v <= 40 and v < w <= 400 (E_w inside E_v)")
    info("ord_{M_v}(2) for v = 3..8: %s (so M_3 | M_v iff v = 3 mod 4, M_4 iff v = 4 mod 12, ...)" % [ords[v] for v in range(3, 9)])
    check(all((Mv(v) % 25 == 0) == (v % 20 == 7) and (Mv(v) % 125 == 0) == (v % 100 == 7) for v in range(3, 401))
          and Mv(16) == 13 * 71 ** 2,
          "B4 5^2 | M_v iff v = 7 mod 20, 5^3 | M_v iff v = 7 mod 100 (v <= 400); M_16 = 13*71^2")
    prog = {}
    for p in shared:
        o = n_order(2, p)
        hits = [v for v in range(3, 201) if Mv(v) % p == 0]
        prog[p] = (hits[0] % o, o)
        check(hits == [v for v in range(3, 201) if v % o == hits[0] % o], "B4 p = %d divides 2^v - 3 exactly for v = %d mod %d (v <= 200)" % (p, hits[0] % o, o))
# maximal order: records n_k = least n = 1 mod 6 with #{v>=3 : M_v | n} >= k
B = 1 << 100
best = {}
t = time.time()
stack, nodes = [(3, 1, 0)], 0
while stack:  # exhaustive over all T with lcm(M_T) <= B
    start, L, size = stack.pop()
    v = start
    while Mv(v) <= B:
        L2 = L * Mv(v) // gcd(L, Mv(v))
        if L2 <= B:
            nodes += 1
            n = L2 if L2 % 6 == 1 else 5 * L2
            if n <= B and (size + 1 not in best or n < best[size + 1]):
                best[size + 1] = n
            stack.append((v + 1, L2, size + 1))
        v += 1
info("record search: %d sets T with lcm(M_T) <= 2^100" % nodes)
rec = {}
for k in sorted(best):
    rec[k] = min(best[kk] for kk in best if kk >= k)
info("records (k, n_k, d_k = (n_k-1)/6, log2 n_k) [%.1fs]:" % (time.time() - t))
for k in sorted(rec):
    info("  k=%d  n_k=%d  d_k=%d  log2 n_k=%.2f  sqrt(2 log2 n_k)=%.2f" % (k, rec[k], (rec[k] - 1) // 6, log2(rec[k]), (2 * log2(rec[k])) ** 0.5))
check([(rec[k] - 1) // 6 for k in range(1, 6)] == [first[j] for j in range(2, 7)] == [2, 24, 314, 7854, 479104],
      "B5 record d_k (first d with m(d) >= k+1) = 2, 24, 314, 7854, 479104 for k <= 5, matching the 1e8 sieve")
check(rec[9] == 690259514980267375 and rec[10] == 704754964794852989875 and rec[11] == 2884562070905333287558375
      and rec[12] == 14541077399433785102581768375,
      "B5 exhaustive below 2^100: n_9 = 690259514980267375, n_10 = 704754964794852989875, n_11 = 2884562070905333287558375, n_12 = 14541077399433785102581768375 (13 copies first at d = 2423512899905630850430294729)")
check([v for v in range(3, 101) if rec[12] % Mv(v) == 0] == [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 19],
      "B5 n_12 is divisible exactly by M_v for v = 3..12, 16, 19 (M_16 = 13*71^2 and M_19 = 5*23*47*97 reuse shared primes)")


def prodM(k):
    r = 1
    for v in range(3, k + 3):
        r *= Mv(v)
    return r


check(all(k + 1 <= 1 + log2(rec[k]) / 2 and rec[k] <= 5 * prodM(k) for k in rec),
      "B5 records respect the proved bounds: k+1 <= 1 + log2(n_k)/2 and n_k <= 5 prod_{v=3}^{k+2} M_v")
# tail lower bound P(N >= k) >= 1/lcm(M_3..M_{k+2})
ok3 = True
for k in range(1, 9):
    L = 1
    for v in range(3, k + 3):
        L = L * Mv(v) // gcd(L, Mv(v))
    if sum(law64[k:]) < Fr(1, L):
        ok3 = False
check(ok3, "B6 P(N >= k) >= 1/lcm(M_3..M_{k+2}) for k <= 8 (lower tail bound)")
info("log2 P(m = j) for j = 2..9: %s" % [round(float(mp.log(mp.mpf(delta[j].numerator) / delta[j].denominator, 2)), 2) for j in range(1, 9)])

# ======================================================================================
info("peak RSS so far: %.0f MB, elapsed %.0f s" % (peak_mb(), time.time() - T0))
print("=== C. the summatory function of K ===")
X = 1 << 20
Kl = np.zeros(X + 1, dtype=np.int64)
for N in range(1, X + 1):
    Kl[N] = KN(N)
Ssum = np.cumsum(Kl)  # int64 is ample: S(2^20) < 2^40
memo = {}


def Srec(x):
    if x <= 0:
        return 0
    if x in memo:
        return memo[x]
    J, U, W = (x - 1) // 2, x // 4, (x + 2) // 4
    r = J * (J + 1) // 2 + 5 * U * (U + 1) // 2 - U + 3 * W * (W + 1) - 4 * W + Srec(W)
    memo[x] = r
    return r


check(all(Srec(x) == int(Ssum[x]) for x in range(0, 20001)) and all(Srec(x) == int(Ssum[x]) for x in range(X - 2000, X + 1)),
      "C1 exact recursion S(X) = J(J+1)/2 + 5U(U+1)/2 - U + 3W(W+1) - 4W + S(W), J=floor((X-1)/2), U=floor(X/4), W=floor((X+2)/4)")
check(all(6 * int(Ssum[4 ** n]) == 3 * 16 ** n - 4 ** n - 2 for n in range(0, 11)), "C2 S(4^n) = (3*16^n - 4^n - 2)/6 for n <= 10")


def ell(x):
    t, r = divmod(x, 4)
    return [Fr(-t, 2), Fr(-(5 * t + 1), 2), Fr(t + 1, 2), Fr(-(3 * t + 2), 2)][r]


def Sdigit(x):
    tot, y = Fr(x * x, 2), x
    while y > 0:
        tot += ell(y)
        y = (y + 2) // 4
    return tot


check(all(Sdigit(x) == int(Ssum[x]) for x in list(range(0, 5001)) + list(range(X - 500, X + 1))),
      "C3 digit formula S(X) = X^2/2 + sum_i ell(X_i), X_{i+1} = floor((X_i+2)/4), ell(4t+r) = -t/2, -(5t+1)/2, (t+1)/2, -(3t+2)/2")
xs = np.arange(1000, X + 1, dtype=np.float64)
phis = (Ssum[1000:].astype(np.float64) - xs * xs / 2) / xs
del xs
info("E(X)/X = (S(X) - X^2/2)/X on [1e3, 2^20]: min %.6f, max %.6f, mean %.6f" % (float(phis.min()), float(phis.max()), float(phis.mean())))
up = [float((Fr(Srec((4 ** n + 2) // 3)) - Fr(((4 ** n + 2) // 3) ** 2, 2)) / ((4 ** n + 2) // 3)) for n in (10, 20, 40)]
dn = [float((Fr(Srec((4 ** n - 1) // 3)) - Fr(((4 ** n - 1) // 3) ** 2, 2)) / ((4 ** n - 1) // 3)) for n in (10, 20, 40)]
check(abs(up[-1] - 1 / 6) < 1e-9 and abs(dn[-1] + 5 / 6) < 1e-9 and -5 / 6 - 1e-3 < float(phis.min()) and float(phis.max()) < 1 / 6 + 1e-3,
      "C4 E(X)/X lies in [-5/6, 1/6] + o(1); extremes attained along X = (4^n+2)/3 (-> 1/6) and X = (4^n-1)/3 (-> -5/6)")
Yc = 10 ** 6
cs = 0
for v in range(2, 80):
    m_ = Mv(v)
    if m_ > 6 * Yc + 1:
        break
    d0 = (-pow(6, -1, m_)) % m_ if m_ > 1 else 0
    if d0 <= Yc:
        cs += (Yc - d0) // m_ + 1
del Kl, Ssum, phis, memo
check(cs == int(mf.sum()), "C5 counting function: sum_{d<=Y} m(d) = sum_v (floor((Y - d_v)/M_v) + 1) exactly (Y = 1e6)")
info("sum_{d<=1e6} m(d) = %d ~ 1.34367 * 1e6" % cs)
lam_cnt = np.zeros(D1 + 1, dtype=np.int32)
for v in range(2, 80):
    m_ = Mv(v)
    if m_ > 6 * D1 + 1:
        break
    r_ = 1 if v % 2 == 0 else 5  # s = M_v^(-1) = M_v mod 6
    for n in range(m_ * r_, 6 * D1 + 2, 6 * m_):
        lam_cnt[(n - 1) // 6] += 1
check(np.array_equal(lam_cnt, mf), "C6 Lambert series: sum_d m(d) q^(6d+1) = sum_{v even} q^(M_v)/(1-q^(6M_v)) + sum_{v odd>=3} q^(5M_v)/(1-q^(6M_v)) (coefficients d <= 1e6)")
del lam_cnt
# all one-step drops (ascents and descents): the arithmetic mean drop is zero
tot, ok, rmin, rmax = 0, True, 9.0, -9.0
for A in range(1, 1 << 22, 2):
    tot += (A - syr(A)[0]) // 2
    if (A + 1) & A == 0:
        n_ = (A + 1).bit_length() - 1
        if tot != -((2 ** (n_ - 1) - (-1) ** (n_ - 1)) // 3):
            ok = False
    if A >= 10 ** 4:
        rmin, rmax = min(rmin, tot / A), max(rmax, tot / A)
check(ok, "C7 sum_{A odd < 2^n} K(A) = -J_{n-1} (Jacobsthal 0,1,1,3,5,11,21,...) exactly for n <= 22: up and down steps cancel to O(X)")
info("sum_{A odd <= X} K(A) / X ranges over [%.4f, %.4f] for 1e4 <= X < 2^22 (extremes next to the trunk (4^j-1)/3)" % (rmin, rmax))

# ======================================================================================
info("peak RSS so far: %.0f MB, elapsed %.0f s" % (peak_mb(), time.time() - T0))
print("=== D. k-step drops, consecutive drops, cycles ===")


ok = True
for k in range(1, 8):
    for A in list(range(1, 100001, 2)) + list(range(-100001, 0, 2)):
        a, word = A, []
        for _ in range(k):
            a, v = syr(a)
            word.append(v)
        c, p = cvec(word)
        if (2 ** p - 3 ** k) * a != 2 * 3 ** k * ((A - a) // 2) + c:
            ok = False
check(ok, "D1 (2^p - 3^k) Syr^k(A) = 2*3^k*D + c(v) for k <= 7 and all odd |A| <= 1e5 (D = (A - Syr^k A)/2)")


def words(k, pmax):
    def rec(pref, s):
        if len(pref) == k:
            yield tuple(pref)
            return
        rem = k - len(pref) - 1
        for v in range(1, pmax - s - rem + 1):
            yield from rec(pref + [v], s + v)
    yield from rec([], 0)


ok = True
for k in range(1, 5):
    for w in words(k, 11):
        c, p = cvec(w)
        s0 = next(s for s in range(1, 2 * 3 ** k + 2, 2) if (2 ** p * s - c) % 3 ** k == 0)
        for tt in range(4):
            s = s0 + 2 * 3 ** k * tt
            A = (2 ** p * s - c) // 3 ** k
            a, ww = A, []
            for _ in range(k):
                a, v = syr(a)
                ww.append(v)
            if not (A > 0 and A % 2 == 1 and tuple(ww) == w and a == s):
                ok = False
check(ok, "D2 bijection: for every word v (k <= 4, p <= 11) and odd s >= 1 with 2^p s = c(v) mod 3^k, A = (2^p s - c)/3^k has word v and Syr^k(A) = s")
for k in (2, 3):
    Yk, Bk = 3000, 400000
    brute = Counter()
    for A in range(1, Bk, 2):
        a = A
        for _ in range(k):
            a = syr(a)[0]
        D = (A - a) // 2
        if -Yk <= D <= Yk:
            brute[D] += 1
    pmax = int(log2(4 * 3 ** k * (Yk + 1))) + k + 8
    form, maxA = Counter(), 0
    for w in words(k, pmax):
        c, p = cvec(w)
        Mk = 2 ** p - 3 ** k
        inv = pow(2 * 3 ** k, -1, abs(Mk))
        D0 = (-c * inv) % abs(Mk)
        for D in range(D0 - abs(Mk) * ((D0 + Yk) // abs(Mk)), Yk + 1, abs(Mk)):
            if D < -Yk:
                continue
            s = (2 * 3 ** k * D + c) // Mk
            if s > 0:
                form[D] += 1
                maxA = max(maxA, (2 ** p * s - c) // 3 ** k)
    check(maxA < Bk and all(brute[D] == form[D] for D in range(-Yk, Yk + 1)),
          "D3 k=%d: multiplicity of a total drop D = #{v : (2^p-3^k) | 2*3^k*D + c(v), quotient > 0}, all |D| <= 3000" % k)
    info("k=%d histogram D in [0,3000]: %s ; D in [-3000,-1]: %s" % (k, sorted(Counter(form[D] for D in range(0, Yk + 1)).items()), sorted(Counter(form[D] for D in range(-Yk, 0)).items())))
lam = mp.log(3) / mp.log(2)


def rho_k(k, q=3):
    L = mp.log(q) / mp.log(2)
    plo = int(mp.floor(k * L))
    rp = mp.fsum(mp.mpf(comb(p - 1, k - 1)) / (mp.mpf(2) ** p - mp.mpf(q) ** k) for p in range(max(k, plo + 1), 14 * k + 200))
    rm = mp.fsum(mp.mpf(comb(p - 1, k - 1)) / (mp.mpf(q) ** k - mp.mpf(2) ** p) for p in range(k, plo + 1)) if plo >= k else mp.mpf(0)
    return rp, rm


tab = {k: rho_k(k) for k in (1, 2, 3, 4, 5, 7, 12, 17, 29, 41, 53, 100, 250)}
for k, (rp, rm) in tab.items():
    info("k=%3d rho_k^+ = %s  rho_k^- = %s  frac(k log2 3) = %s" % (k, mp.nstr(rp, 10), mp.nstr(rm, 10), mp.nstr(k * lam - mp.floor(k * lam), 4)))
check(abs(tab[1][0] - rho1) < 1e-30 and tab[1][1] == 1, "D4 rho_1^+ = sum_{v>=2} 1/(2^v-3) = 1.3437..., rho_1^- = 1")
check(abs(tab[2][1] - mp.mpf(11) / 5) < 1e-25, "D4 rho_2^- = 1/(9-4) + 2/(9-8) = 11/5 (two unit words (1,2),(2,1): 3^2 - 2^3 = 1)")
check(abs(tab[250][0] - 1) < 1e-4 and tab[250][1] < 1e-6, "D4 rho_k^+ -> 1 and rho_k^- -> 0 (k = 250: %s, %s)" % (mp.nstr(tab[250][0], 8), mp.nstr(tab[250][1], 4)))
Irate = mp.log(3) + (lam - 1) * mp.log(lam - 1) - lam * mp.log(lam)
kk = 4000
pstar = int(mp.floor(kk * lam))
emp = -(mp.log(mp.mpf(comb(pstar - 1, kk - 1))) - pstar * mp.log(2)) / kk
info("large-deviation rate at p = k log2 3: I = log 3 + (l-1)log(l-1) - l log l = %s ; -(1/k) log pi_k(floor(k l)) at k=4000: %s" % (mp.nstr(Irate, 8), mp.nstr(emp, 8)))
check(abs(Irate - mp.mpf("0.0549795")) < 1e-6 and abs(emp - Irate) < 3e-3, "D4 resonance spikes are damped by exp(-I k), I = 0.05498 nats (0.0793 bits) per odd step")
ok = True
for D in range(-100000, 0):
    for A, w in ((-16 * D - 5, (1, 2)), (-16 * D - 7, (2, 1))):
        a1, x1 = syr(A)
        a2, x2 = syr(a1)
        if (x1, x2) != w or a2 != A - 2 * D:
            ok = False
check(ok, "D5 'two copies' at k = 2: every D <= -1 (|D| <= 1e5) is the 2-step drop of A = -16D-5 (word 12) and A = -16D-7 (word 21)")
# consecutive drops determine the point
for sgn in (+1, -1):
    seen = set()
    dup = 0
    for A in range(1, 1000001, 2):
        a1 = Tq(3, sgn, A)[0]
        a2 = Tq(3, sgn, a1)[0]
        key = (((A - a1) // 2) << 40) + (a1 - a2) // 2 + (1 << 39)  # injective encoding of the pair (|drops| < 2^39)
        if key in seen:
            dup += 1
        seen.add(key)
    check(dup == 0, "D6 3x%s1: A -> (K(A), K(T(A))) is injective on odd A < 1e6 (brute force)" % ("+" if sgn > 0 else "-"))
    del seen


ok = True
for q, sgn in ((1, 1), (5, 1), (7, 1), (9, 1), (11, 1), (13, 1), (5, -1)):
    seen = set()
    for A in range(1, 200001, 2):
        a1 = Tq(q, sgn, A)[0]
        a2 = Tq(q, sgn, a1)[0]
        key = (((A - a1) // 2) << 48) + (a1 - a2) // 2 + (1 << 47)
        if key in seen:
            ok = False
        seen.add(key)
    del seen
check(ok, "D6 (EMPIRICAL for other q) no repeated consecutive-drop pair for A < 2e5 for qx+1, q = 1,5,7,9,11,13, and for 5x-1")


def M3(v):
    return (1 << v) - 3


sols, sols_neg = [], []
R = 60
for sgn in (+1, -1):
    for v in range(1, R):
        for vp in range(1, R):
            if v == vp:
                continue
            for u in range(1, R):
                for up in range(u + 1, R):
                    d_ = up - u
                    E = (1 << d_) * M3(u) * M3(vp) - M3(up) * M3(v)
                    if E == 0:
                        continue
                    num = sgn * M3(v) * M3(vp) * ((1 << d_) - 1)
                    if num % E:
                        continue
                    n1 = num // E
                    if n1 <= 0 or n1 % M3(v) or n1 % M3(vp) or n1 % 6 != (1 if sgn > 0 else 5):
                        continue
                    if n1 // M3(v) <= 0 or n1 // M3(vp) <= 0:  # S, S' must be positive (positive sheet)
                        sols_neg.append((sgn, v, vp, u, up, n1))
                        continue
                    sols.append((sgn, v, vp, u, up, n1))
info("quadruple hits with a negative point (not on the positive sheet): %s  [the fixed points 1 and -1 share the drop pair (0,0)]" % sols_neg)
check(sols == [] and sols_neg == [(1, 1, 2, 1, 2, 1)], "D6 no quadruple (v,v',u,u') with exponents < 60 solves n1*E = +-M_v M_v' (2^d - 1) with n1 admissible (both sheets): injectivity holds for all A with such valuations")
# the proof's finite case: only (3,4,2,3) survives the size bounds, with n1 = 65 = 5 mod 6
cases = {+1: [], -1: []}
for v in range(2, 34):
    for vp in range(2, 34):
        if vp == v:
            continue
        g_ = gcd(M3(v), M3(vp))
        for u in range(2, 14):
            for d_ in range(1, 34):
                up = u + d_
                E = (1 << d_) * M3(u) * M3(vp) - M3(up) * M3(v)
                if 0 < E <= g_ * ((1 << d_) - 1):
                    cases[+1].append((v, vp, u, up, E))
                if 0 < -E <= g_ * ((1 << d_) - 1):
                    cases[-1].append((v, vp, u, up, E))
check(cases[+1] == [(3, 4, 2, 3, 1)] and cases[-1] == [],
      "D6 proof check (exponents 2..33): 0 < E <= g(2^d-1) only at (v,v',u,u') = (3,4,2,3), E = 1, where n1 = 65 = 5 mod 6 is inadmissible; 0 < -E <= g(2^d-1) never (3x-1)")


def fiber(k):
    out = []

    def rec(i, ps, c):
        if i == k:
            Mk = 2 ** ps - 3 ** k
            if c % Mk == 0:
                out.append(c // Mk)
            return
        term = 3 ** (k - 1 - i) * 2 ** ps
        for v in range(1, 2 * k - ps - (k - i - 1) + 1):
            rec(i + 1, ps + v, c + term)
    rec(0, 0, 0)
    return out


ok = True
for k in range(1, 11):
    fb = fiber(k)
    expect = [1, -1] + ([-5, -7] if k % 2 == 0 else []) + ([-17, -25, -37, -55, -41, -61, -91] if k % 7 == 0 else [])
    if sorted(fb) != sorted(expect):
        ok = False
    info("k=%d: D = 0 fiber (all words with k <= p <= 2k): %s" % (k, sorted(fb)))
check(ok, "D7 the D = 0 fiber is the set of integer periodic points: {1} on the positive sheet, {-1}, {-5,-7}, {-17,...} on the negative sheet (k <= 10)")
# sheet symmetry at k = 2: 3x+1 and 3x-1 on positive integers have the same drop-multiplicity law
Ys = 10 ** 6
hist = {}
for sgn in (+1, -1):
    mu = np.zeros(Ys + 1, dtype=np.int8)
    for p in range(4, 70):
        Mk = 2 ** p - 9
        if Mk > 18 * Ys + 2 ** p:
            break
        inv = pow(18, -1, Mk)
        for v1 in range(1, p):
            c = 3 + 2 ** v1
            r = (-sgn * c * inv) % Mk  # 3x+1: 18D + c = 0 ; 3x-1: 18D - c = 0 (mod M)
            start = r
            if sgn < 0:
                while 18 * start - c <= 0:  # quotient must be positive
                    start += Mk
            mu[start::Mk] += 1
    hist[sgn] = np.bincount(mu, minlength=8)[:8] / (Ys + 1)
    if sgn < 0:
        bc = Counter()
        for A in range(1, 300000, 2):
            a = Tq(3, -1, Tq(3, -1, A)[0])[0]
            D = (A - a) // 2
            if 0 <= D <= 3000:
                bc[D] += 1
        check(all(bc[D] == mu[D] for D in range(1, 3001)) and bc[0] == 3,
              "D8 3x-1 two-step multiplicity sieve agrees with brute force on [1,3000]; at D = 0 brute force finds the 3x-1 periodic points 1, 5, 7")
info("k=2 multiplicity law on [0,1e6]: 3x+1 %s" % np.round(hist[1], 5).tolist())
info("                                 3x-1 %s" % np.round(hist[-1], 5).tolist())
check(float(np.max(np.abs(hist[1] - hist[-1]))) < 2e-3, "D8 the 3x+1 and 3x-1 sheets have the same 2-step drop-multiplicity law (to 2e-3 on [0,1e6]) although their cycles differ")
cnt2 = Counter()
for A in range(1, 1 << 21, 2):
    a1, x1 = syr(A)
    a2, x2 = syr(a1)
    cnt2[(x1, x2)] += 1
check(all(cnt2[(a, b)] == 2 ** (20 - a - b) for a in range(1, 19) for b in range(1, 20 - a)),
      "D9 valuation pairs are exactly geometric-independent on odd A < 2^21 (#{word (a,b)} = 2^(20-a-b)): normalized drops K(A_i)/A_i -> iid (2^V-3)/2^(V+1)")

# ======================================================================================
# exact k-step multiplicity laws (k = 2, 3), computed at start-up by the variable-elimination engine
info("engine computations (run first, in a child process to keep memory low) took %.1fs, child peak RSS %.0f MB" % (ENG["time"], ENG["rss"]))
check(max(abs(ENG["k1"][j] - float(law64[j])) for j in range(8)) < 1e-13, "D10 the elimination engine reproduces the exact k = 1 law (to 1e-13)")
klaws = {}
for k in (2, 3):
    lw, pm = ENG[k]
    klaws[k] = lw
    rp = rho_k(k)[0]
    tail = sum(comb(p - 1, k - 1) / 2 ** p for p in range(51, 2000))
    info("k=%d exact law of the k-step drop multiplicity on D >= 0 (words p <= 50, truncation error <= 2 P(P_k > 50) = %.1e):" % (k, 2 * tail))
    info("   P(mu_%d = 0..9) = %s" % (k, ["%.12f" % x for x in lw[:10]]))
    check(abs(lw.sum() - 1) < 1e-12 and abs(sum(j * x for j, x in enumerate(lw)) - float(rp)) < 1e-9,
          "D10 k=%d: the law sums to 1 and its mean is rho_%d^+ = %s" % (k, k, mp.nstr(rp, 12)))
    check(abs(sum(j * (j - 1) * x for j, x in enumerate(lw)) - float(pm)) < 1e-11,
          "D10 k=%d: the engine's E[mu(mu-1)] = %.12f equals the exact sum over word pairs (CRT compatibilities), words p <= 50" % (k, float(pm)))
check(float(np.max(np.abs(hist[1][:6] - klaws[2][:6]))) < 5e-5, "D10 k=2 law agrees with the 3x+1 sieve on [0,1e6] (D8) to 5e-5")
Ysv = 10 ** 6
mu3 = np.zeros(Ysv + 1, dtype=np.int16)
for p in range(5, 61):
    Mk = 2 ** p - 27
    inv = pow(54, -1, Mk)
    for w in words_k(3, p):
        D0 = (-cvec(w)[0] * inv) % Mk
        if D0 <= Ysv:
            mu3[D0::Mk] += 1
sv3 = np.array([np.count_nonzero(mu3 == j) for j in range(8)]) / (Ysv + 1)
del mu3
info("k=3 sieve on [0,1e6]: %s" % np.round(sv3, 6).tolist())
check(float(np.max(np.abs(sv3 - klaws[3][:8]))) < 5e-4, "D10 k=3 law agrees with the sieve on [0,1e6] to 5e-4")
la2, la3 = ENG["asc2"], ENG["asc3"]
check(abs(la2[2] - 0.8) < 1e-14 and abs(la2[3] - 0.2) < 1e-14 and abs(la3[0] - 144 / 209) < 1e-14 and abs(la3[1] - 62 / 209) < 1e-14 and abs(la3[2] - 3 / 209) < 1e-14,
      "D10 ascending side (D < 0): k=2 law (0, 0, 4/5, 1/5); k=3 law (144/209, 62/209, 3/209) [mu_3 = [19 | D] + [D mod 11 in {1,7,8}]]")

info("peak RSS so far: %.0f MB, elapsed %.0f s" % (peak_mb(), time.time() - T0))
print("=== E. other q (D75): x -> oddpart(qx + s) ===")
ok = True
for q in (1, 3, 5, 7, 9, 11, 13):
    for sgn in (1, -1):
        for A in range(-20001, 20002, 2):
            if q * A + sgn == 0:
                continue
            T, v = Tq(q, sgn, A)
            if 2 * q * ((A - T) // 2) + sgn != ((1 << v) - q) * T:
                ok = False
check(ok, "E1 2qK + s = (2^v - q) T for T = oddpart(qA + s), s = +-1, q in {1,3,...,13}, odd |A| <= 2e4")


def mult_q(q, sgn, d):
    n = 2 * q * d + sgn
    c, v = 0, 1
    while (1 << v) - q <= abs(n) + q:
        m_ = (1 << v) - q
        if m_ != 0 and n % abs(m_) == 0 and (n // m_ > 0 if m_ > 0 else (-n) // (-m_) > 0):
            c += 1
        v += 1
    return c


ok = True
for q, sgn in ((1, 1), (5, 1), (7, 1), (11, 1), (5, -1), (3, -1), (7, -1)):
    Yq = 1500
    br = Counter()
    for A in range(1, 4 * q * (Yq + 2), 2):
        T, _ = Tq(q, sgn, A)
        if T == 0:
            continue
        K = (A - T) // 2
        if abs(K) <= Yq:
            br[K] += 1
    if any(br[d] != mult_q(q, sgn, d) for d in range(-Yq, Yq + 1)):
        ok = False
check(ok, "E2 multiplicity of drop d for oddpart(qx+s) = #{v : (2^v - q) | 2qd + s, quotient > 0} (q,s) in {(1,+),(5,+),(7,+),(11,+),(5,-),(3,-),(7,-)}, |d| <= 1500")
try:
    from sympy import factorint as _fi

    def fsmall(n):
        return _fi(n)
except ImportError:
    fsmall = trial_factor
laws = {}
for q in (1, 3, 5, 7, 9, 11, 13):
    v0 = 1
    while (1 << v0) <= q:
        v0 += 1
    Vq = 24 if q == 1 else 48
    desc = [(1 << v) - q for v in range(v0, Vq + 1)]
    asc = [q - (1 << v) for v in range(1, v0)]
    dd = exact_law(*structure(desc, fsmall))[0]
    da = exact_law(*structure(asc, fsmall))[0] if asc else [Fr(1)]
    laws[q] = (dd, da)
    info("q=%2d desc units %s asc units %s | desc P(m=0..4) %s | asc P(m=0..3) %s" % (
        q, [m_ for m_ in desc if m_ == 1], [m_ for m_ in asc if m_ == 1],
        [fl(c, 6) for c in (dd + [Fr(0)] * 5)[:5]], [fl(c, 6) for c in (da + [Fr(0)] * 4)[:4]]))
check(laws[3][0][1] == law64[0] or abs(float(laws[3][0][1]) - float(law64[0])) < 1e-13,
      "E3 q = 3 descending law (truncation 48) agrees with section B (shifted by the unit)")
# brute-force hole densities
for q, side, exact in ((5, +1, laws[5][0][0]), (9, +1, laws[9][0][0]), (11, +1, laws[11][0][0]), (13, +1, laws[13][0][0]),
                       (7, -1, laws[7][1][0]), (11, -1, laws[11][1][0]), (13, -1, laws[13][1][0])):
    Yh = 30000
    hit = np.zeros(Yh + 1, dtype=np.int32)
    for A in range(1, 4 * q * (Yh + 2), 2):
        T, _ = Tq(q, 1, A)
        K = (A - T) // 2
        if side > 0 and 0 <= K <= Yh:
            hit[K] += 1
        if side < 0 and -Yh <= K <= -1:
            hit[-K] += 1
    frac = float(np.mean(hit[(0 if side > 0 else 1):] == 0))
    check(abs(frac - float(exact)) < 4e-3, "E3 q=%d %s holes: brute force %.5f vs exact density %s" % (q, "descent" if side > 0 else "ascent", frac, fl(exact, 6)))
ok = True
for q in range(1, 1026, 2):
    v0 = 1
    while (1 << v0) <= q:
        v0 += 1
    desc_unit = (1 << v0) - q == 1
    asc_unit = any(q - (1 << v) == 1 for v in range(1, v0))
    if desc_unit != (((q + 1) & q) == 0):
        ok = False
    if asc_unit != (q > 1 and ((q - 1) & (q - 2)) == 0):
        ok = False
    if desc_unit and asc_unit and q != 3:
        ok = False
    if not desc_unit:
        s = sum(Fr(1, (1 << v) - q) for v in range(v0, v0 + 40)) + Fr(1, 2 ** 38)
        if s >= Fr(62, 100):
            ok = False
    if not asc_unit and v0 > 1:
        s = sum(Fr(1, q - (1 << v)) for v in range(1, v0))
        if s >= Fr(62, 100):
            ok = False
check(ok, "E4 odd q <= 1025: a descending unit iff q = 2^a - 1, an ascending unit iff q = 2^a + 1, both only for q = 3; a side without unit has sum 1/|modulus| < 0.62 (holes of density > 0.38)")
for q in (1, 5, 7):
    tq = {k: rho_k(k, q) for k in (1, 10, 50, 100, 200)}
    info("q=%d rho_k^+ + rho_k^- for k = 1,10,50,100,200: %s" % (q, [mp.nstr(sum(tq[k]), 6) for k in tq]))
    if q == 1:
        check(abs(sum(tq[200]) - 1) < 1e-6, "E5 q=1: rho_k -> 1")
    else:
        check(sum(tq[200]) < 0.01 and sum(tq[200]) < sum(tq[50]), "E5 q=%d: rho_k^+ + rho_k^- -> 0 (expanding drift log(q/4) > 0)" % q)
ok = True
for q, sgn in ((5, 1), (7, 1), (3, -1), (5, -1)):
    for k in range(1, 5):
        for A in range(1, 3001, 2):
            a, word = A, []
            for _ in range(k):
                a, v = Tq(q, sgn, a)
                word.append(v)
            c, p = cvec(word, q)
            if (2 ** p - q ** k) * a != 2 * q ** k * ((A - a) // 2) + sgn * c:
                ok = False
check(ok, "E6 (2^p - q^k) T^k(A) = 2 q^k D + s c_q(v) for (q,s) in {(5,+),(7,+),(3,-),(5,-)}, k <= 4")
fb5 = {}
for k in (1, 2, 3):
    out = []
    for w in words(k, 3 * k):
        c, p = cvec(w, 5)
        Mk = 2 ** p - 5 ** k
        if Mk > 0 and c % Mk == 0:
            out.append(c // Mk)
    fb5[k] = sorted(out)
info("5x+1 positive D=0 fibers k=1,2,3 (5^k < 2^p <= 6^k): %s" % fb5)
check(fb5[2] == [1, 3] and set(fb5[3]) == {13, 17, 27, 33, 43, 83}, "E7 5x+1: the D = 0 fibers are its cycles {1,3}, {13,33,83}, {17,43,27} (drop statistics generic, cycles present)")

# ======================================================================================
info("peak RSS so far: %.0f MB, elapsed %.0f s" % (peak_mb(), time.time() - T0))
print("=== F. the drop identity mod 3^k ===")
check(all(pow(2, 3 ** (k - 1), 3 ** k) != 1 and pow(2, 2 * 3 ** (k - 2), 3 ** k) != 1 if k >= 2 else True for k in range(1, 11)),
      "F1 2 is a primitive root mod 3^k (k <= 10), so 2^v mod 3^k has period 2*3^(k-1)")
tab9 = {}
ok = True
for A in range(1, 200001, 2):
    s, v = syr(A)
    K = (A - s) // 2
    key = (K % 3, v % 6)
    if tab9.setdefault(key, s % 9) != s % 9:
        ok = False
    if (s - (6 * K + 1) * pow((1 << v) - 3, -1, 9)) % 9 != 0:
        ok = False
    if (K - 2 * (A - (-1) ** v)) % 3 != 0:
        ok = False
check(ok and len(tab9) == 18, "F2 S(A) mod 9 = (6K+1)(2^v-3)^(-1) mod 9 is a function of (K mod 3, v mod 6); K = 2(A - (-1)^v) mod 3")
for r in range(3):
    info("K = %d mod 3: S mod 9 for v = 1..6 mod 6: %s" % (r, [tab9[(r, v % 6)] for v in range(1, 7)]))
ok = True
for k in range(1, 7):
    for A in range(1, 30001, 2):
        a, word = A, []
        for _ in range(k):
            a, v = syr(a)
            word.append(v)
        c, p = cvec(word)
        if (a - c * pow(2, -p, 3 ** k)) % 3 ** k:
            ok = False
check(ok, "F3 Syr^k(A) = 2^(-p) c(v) mod 3^k depends only on the valuation word (k <= 6)")
coord = {((-1) ** c_ * pow(4, s_, 9)) % 9: (c_, s_) for c_ in (0, 1) for s_ in range(3)}
ok = True
for A in range(1, 200001, 2):
    if A % 3 == 0:
        continue
    S_, v_ = syr(A)
    c_, s_ = coord[A % 9]
    if coord[S_ % 9] != (v_ % 2, (1 + c_ + v_) % 3):
        ok = False
check(ok, "F5 (cross-check of the tcpc lane's G1) with x = (-1)^c 4^s mod 9, one Syracuse step is (c,s) -> (v mod 2, 1+c+v mod 3), odd A < 2e5 prime to 3")
Xf = 10 ** 7
m7 = np.ones(Xf + 1, dtype=np.int8)
v = 3
while Mv(v) <= 6 * Xf + 1:
    m_ = Mv(v)
    m7[(-pow(6, -1, m_)) % m_::m_] += 1
    v += 1
ok = True
for r in range(9):
    sub = m7[r::9]
    for j in (1, 2, 3):
        if abs(np.mean(sub == j) - float(delta[j - 1])) > 1e-3:
            ok = False
check(ok, "F4 the residue d mod 9 carries no information on m(d): each class has the frequencies delta_1..3 to 1e-3 (exact by CRT)")
del m7

print()
print("%d checks, %.0f s, peak RSS %.0f MB" % (NCHECK, time.time() - T0, peak_mb()))
print("ALL CHECKS PASSED")
