# opus-2026-07-15-S314 -- HYP-6925: (1) THM-864 referee (proved ledger vs exact
# err, must dominate every battery row); (2) SWEEP 1: single-cluster sub-N0
# cores, probe-first exact runner; (3) SWEEP 2: multi-cluster reduction list.
from fractions import Fraction
from math import gcd
import math, itertools, time, sys

DELTA = Fraction(1, 13)

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def mu(ivs): return sum(b - a for (a, b) in ivs)

def comb_teeth_in(x, a, b):
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def pair_measure_on(E, x1, x2):
    tot = Fraction(0)
    for (a, b) in E:
        for (lo, hi) in comb_teeth_in(x1, a, b):
            for (l2, h2) in comb_teeth_in(x2, lo, hi): tot += h2 - l2
    return tot

def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        T += (min(2*a, s - 13*m)) * (1 if m == 0 else 2)
        m += 1
    return T

def rho_pair(A, B):
    g = gcd(A, B)
    a, b = A//g, B//g
    if a > b: a, b = b, a
    return Fraction(T_of(a, b), 13*a*b)

# ---------- (1) THM-864 referee ----------
def proved_bound(A, B, q, p, y, kap):
    rho = rho_pair(A, B)
    M0 = (A + B)//13
    nJ = (2*M0)//y
    if nJ < 1: return None
    J = Fraction(nJ*q, A)
    Enet = Fraction(2*kap, y) * rho / J
    # E2 exact-ish: (2 kap nu + 2y) * lmax ; nu = ceil(J*y)+1 ; lmax = 2/(13B)
    nu = math.ceil(J * y) + 1
    E2 = Fraction(2*(2*kap*nu + 2*y), 13*B)
    E1 = Fraction(8*kap, 13*B)
    E3a = Fraction(4*y, 13*B)
    return Enet + E1 + E2 + E3a

def referee():
    PREFIXES = [[1,2,3,4,5], [2,3,5,6,8], [8,9,10,11,12]]
    total = fails = 0
    worst = 0.0
    for P in PREFIXES:
        E = safe_set(P)
        muE = mu(E); kap = len(E)
        for (q, p) in [(1,1), (1,2), (1,3), (2,3), (1,12), (3,4), (2,5)]:
            for y in (5, 11, 23, 47, 95):
                for A0 in (611, 1200):
                    found = None
                    for A in range(A0, A0 + 6*q + 1):
                        if (p*A + y) % q: continue
                        B = (p*A + y)//q
                        if B <= A or gcd(A, B) != 1: continue
                        # planted relation must be the small-q minimum
                        best = None
                        for qq in range(1, 14):
                            p0 = round(qq*B/A)
                            for pp in (p0-1, p0, p0+1):
                                if pp < 1: continue
                                v = abs(qq*B - pp*A)
                                if best is None or v < best[0]: best = (v, qq, pp)
                        if best[0] != y or best[1] != q: continue
                        found = (A, B); break
                    if not found: continue
                    A, B = found
                    if A < 26*q*y: continue          # standing size hypothesis
                    err = abs(pair_measure_on(E, A, B) - muE*rho_pair(A, B))
                    bnd = proved_bound(A, B, q, p, y, kap)
                    if bnd is None: continue
                    total += 1
                    ratio = float(err/bnd) if bnd > 0 else 999
                    worst = max(worst, ratio)
                    if err > bnd:
                        fails += 1
                        print(f"   DOMINATION FAIL: P={P} (q,p,y)=({q},{p},{y}) "
                              f"A={A} B={B} err={float(err):.6f} bound={float(bnd):.6f}")
    print(f"(1) THM-864 referee: {total} rows, domination failures = {fails}, "
          f"worst err/bound = {worst:.4f}")

# ---------- (2) SWEEP 1: single-cluster sub-N0 cores ----------
def union_measure_at(ds, t):
    W = DELTA
    arcs = []
    for d in ds:
        c = (-d * t) % 1
        lo, hi = c - W, c + W
        if lo < 0: arcs.extend([(lo % 1, Fraction(1)), (Fraction(0), hi)])
        elif hi > 1: arcs.extend([(lo, Fraction(1)), (Fraction(0), hi % 1)])
        else: arcs.append((lo, hi))
    arcs.sort()
    tot, cur = Fraction(0), Fraction(0)
    for (lo, hi) in arcs:
        if hi <= cur: continue
        tot += hi - max(lo, cur)
        cur = hi
    return tot

def pattern_data(E, ds):
    # breakpoints within E + exact integral of a + probe cells (a > 0)
    bps = set()
    for (a, b) in E: bps.update([a, b])
    diffs = {abs(x - yy) for x in ds for yy in ds if x != yy} | set(ds)
    for d in diffs:
        if d == 0: continue
        for k in range(0, d + 1):
            for r in (Fraction(0), Fraction(2,13), Fraction(-2,13),
                      Fraction(1,13), Fraction(-1,13)):
                t = (Fraction(k) + r) / d
                if 0 < t < 1: bps.add(t)
    bps = sorted(bps)
    total = Fraction(0); Bcount = 0; probes = []
    for (lo, hi) in E:
        cells = [lo] + [t for t in bps if lo < t < hi] + [hi]
        Bcount += len(cells) - 2
        for c0, c1 in zip(cells, cells[1:]):
            m0 = c0 + (c1-c0)/1000
            m1 = c1 - (c1-c0)/1000
            a0 = 1 - union_measure_at(ds, m0)
            a1 = 1 - union_measure_at(ds, m1)
            total += (a0 + a1) / 2 * (c1 - c0)
            if a0 > 0 or a1 > 0:
                probes.append((c0 + (c1-c0)/2, max(a0, a1)))
    probes.sort(key=lambda z: -z[1])
    return total, Bcount, [t for t, _ in probes[:40]]

def in_E(E, t):
    for (a, b) in E:
        if a < t < b: return True
    return False

def packet_nontight(E, xs, probes):
    # find exact safe t: ||x t|| > 1/13 for all x, t in E
    for t in probes:
        ok = True
        for x in xs:
            f = (x * t) % 1
            if min(f, 1 - f) <= DELTA: ok = False; break
        if ok: return True
    return False

def sweep_single(P, R, time_budget):
    E = safe_set(P)
    t_start = time.time()
    n_pat = n_pack = n_fallback = n_unresolved = 0
    for base in R:
        cs = sorted(((r - base) % 13) for r in R if r != base)
        choices = [[c + 13*k for k in range(4) if c + 13*k <= 42] for c in cs]
        for combo in itertools.product(*choices):
            ds = sorted(set([0] + list(combo)))
            if len(ds) != 7: continue
            if any(b - a > 7 for a, b in zip(ds, ds[1:])): continue
            if time.time() - t_start > time_budget:
                return n_pat, n_pack, n_fallback, n_unresolved, False
            Ia, Bc, probes = pattern_data(E, ds)
            if Ia <= 0:
                n_unresolved += 1
                continue
            CA = 28*Bc + 4*sum(ds)
            N0 = int(2*CA/Ia) + 1
            n_pat += 1
            # extra probe points: rational grid
            probes = probes + [Fraction(2*i+1, 1024) for i in range(512)
                               if in_E(E, Fraction(2*i+1, 1024))][:60]
            N = base + 13*max(1, (14 - base + 12)//13)   # first proper lift >= 14
            while N <= N0:
                xs = [N + d for d in ds]
                n_pack += 1
                if not packet_nontight(E, xs, probes):
                    n_fallback += 1
                    # exact fallback: full uncovered computation
                    U = list(E)
                    for x in xs:
                        out = []
                        for (a, b) in U:
                            cur = a
                            for (lo, hi) in sorted(comb_teeth_in(x, a, b)):
                                if lo > cur: out.append((cur, min(lo, b)))
                                cur = max(cur, hi)
                                if cur >= b: break
                            if cur < b: out.append((cur, b))
                        U = [iv for iv in out if iv[0] < iv[1]]
                    if not U or mu(U) <= 0:
                        n_unresolved += 1
                        print(f"   !! COVERED PACKET (tight candidate): P={P} ds={ds} N={N}")
                N += 13
    return n_pat, n_pack, n_fallback, n_unresolved, True

def sweeps():
    CASES = [
        ([2, 4, 6, 8, 10], [1, 3, 5, 7, 9, 11, 12]),
        ([8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 6, 7]),
        ([2, 5, 7, 10, 12], [1, 3, 4, 6, 8, 9, 11]),
        ([4, 5, 6, 11, 12], [1, 2, 3, 7, 8, 9, 10]),
        ([1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11, 12]),
    ]
    print("\n(2) SWEEP 1 -- single-cluster sub-N0 cores (probe-first exact):")
    for P, R in CASES:
        np_, npk, nfb, nun, complete = sweep_single(P, R, time_budget=1500)
        print(f"   P={P}: patterns={np_} packets={npk} fallbacks={nfb} "
              f"unresolved={nun} complete={complete}", flush=True)

    # SWEEP 2: the multi-cluster reduction list
    print("\n(3) SWEEP 2 -- multi-cluster reduction: X0 = 188 (singles-excess "
          "154/(169 x) <= muE phi*/2 with muE >= 3/13, phi* = 17/546)")
    n6 = 0
    with open('05-knowledge/results/radius7_multicluster_handoff_opus_S314.txt',
              'w', encoding='utf-8', newline='\n') as f:
        f.write("# radius-7 multi-cluster reduction: 6-speed prefixes P u {x1}\n"
                "# for kind-pasteur's THM-815 radius-6 recursion (m'=6 coercive).\n"
                "# x1 = first lift <= X0 = 188; one line per (P, x1).\n")
        for Pset in itertools.combinations(range(1, 13), 5):
            R7 = [r for r in range(1, 13) if r not in Pset]
            for r in R7:
                x1 = r + 13
                while x1 <= 188:
                    f.write(f"{list(Pset)} + [{x1}]\n")
                    n6 += 1
                    x1 += 13
    print(f"   handoff list written: {n6} six-speed prefixes "
          f"(792 five-prefixes x lifts <= 188) -> radius7_multicluster_handoff_opus_S314.txt")

if __name__ == '__main__':
    referee()
    sweeps()
