#!/usr/bin/env python3
"""Lane G5 (death-star-2026-07-30-coinC2): the construction-gate reading of
certificate (27), taken seriously via ANTICIPATORY DYADIC-CHECKPOINT schemes.

Frame: THM-2966 spine normal form. Depth law d_m = floor(g1*m/g2) + D0.
Cell (m,k) side 0: (z,o) = (m+d_m-k, k+1), cap binom(d_m,k); side 1 mirrored
(z,o) = (k+1, m+d_m-k). Doubled deficits delta: |delta| <= cap,
delta == cap (mod 2). Scheme feasible <=> D_M = (1/2) sum delta p^z q^o -> 0
on (0,1).

CHECKPOINT-CLOSURE PROGRAM (this lane's sufficient reduction, proof in the
findings file):
  If the "base" system (rows 1..2^{r0}-1) and every "epoch" system
  (rows 2^r..2^{r+1}-1, r >= r0) admit doubled deficits with
      sum_cells delta * p^z (1-p)^o == 0   identically in p        (E_r)
  then a fair extractor with deadline T(n) = n+1+d_n <= (1+gamma)n + O(1)
  exists: C <= 1 + gamma.  Parity consistency of (E_r) at dyadic endpoints
  is PROVED (Frobenius; verified here in mode `parity`).

Modes:
  parity    -- verify mod-2 consistency of checkpoint epoch systems (PROVED
               expectation) and inconsistency at non-dyadic cut points.
  dfs       -- exhaustive class-sequential DFS on small base/epoch systems:
               exact FEASIBLE (witness) / INFEASIBLE (exhausted) verdicts.
  solve     -- randomized exact witness search on larger epochs.
  abeldini  -- integer hunts for (896,1285), (2974400,8847357) as
               consecutive partial sums / structured products.
  all       -- quick suite.

All asserted computations exact int/Fraction; floats only for orientation.
"""
import sys, math, random, argparse, itertools
from fractions import Fraction
from math import comb

# ---------------------------------------------------------------- ledger ----

def depth(m, g1, g2, D0):
    return (g1 * m) // g2 + D0

def row_cells(m, g1, g2, D0):
    """Cells of row m, both sides: list of (side, m, k, z, o, cap)."""
    d = depth(m, g1, g2, D0)
    out = []
    for k in range(d + 1):
        cap = comb(d, k)
        out.append((0, m, k, m + d - k, k + 1, cap))
        out.append((1, m, k, k + 1, m + d - k, cap))
    return out

def system_cells(m_lo, m_hi, g1, g2, D0):
    cells = []
    for m in range(m_lo, m_hi + 1):
        cells.extend(row_cells(m, g1, g2, D0))
    return cells

def cell_poly(z, o):
    """Coefficient list of p^z (1-p)^o : coeff[O] for O = z..z+o."""
    return {z + j: ((-1) ** j) * comb(o, j) for j in range(o + 1)}

# ------------------------------------------------------------ mode parity ---

def hom2_boxsum(m_lo, m_hi, g1, g2, D0):
    """Parity vector (mod 2) of sum_cells cap * p^z(1-p)^o, as coeff list.
    Consistency of (E_r): system A x = 0 with x == caps (mod 2) requires
    A*caps == 0 (mod 2), i.e. this vector == 0."""
    acc = {}
    for (_s, _m, _k, z, o, cap) in system_cells(m_lo, m_hi, g1, g2, D0):
        if cap % 2 == 1:
            for O, c in cell_poly(z, o).items():
                acc[O] = (acc.get(O, 0) + c) % 2
    return {O: v for O, v in acc.items() if v}

def mode_parity():
    print("== parity: mod-2 consistency of cut-point systems ==")
    gammas = [(1,1),(1,2),(1,3),(2,5),(3,8),(4,11),(3,5),(3,4),(2457,6592)]
    ok = True
    for (g1, g2) in gammas:
        for D0 in (0, 2):
            # checkpoint epochs: rows [2^r, 2^{r+1}-1]
            for r in range(1, 6):
                bad = hom2_boxsum(2 ** r, 2 ** (r + 1) - 1, g1, g2, D0)
                assert not bad, (g1, g2, D0, r, bad)
            # base systems rows [1, 2^r - 1]
            for r in range(1, 6):
                bad = hom2_boxsum(1, 2 ** r - 1, g1, g2, D0)
                assert not bad, ("base", g1, g2, D0, r, bad)
    print("  PASS: all checkpoint epochs+bases parity-consistent "
          "(9 gammas x 2 D0 x r<6)")
    # hostile control: non-dyadic cut points should generically FAIL
    fails = 0; tries = 0
    for (g1, g2) in [(1,2),(1,3),(3,5)]:
        for (lo, hi) in [(3,5),(4,8),(5,9),(6,11),(9,17),(5,11)]:
            tries += 1
            if hom2_boxsum(lo, hi, g1, g2, 0):
                fails += 1
    print(f"  control: non-dyadic-cut systems parity-INCONSISTENT in "
          f"{fails}/{tries} cases (expected: most)")
    return ok

# ----------------------------------------------- class-sequential machinery -

class ClosureSystem:
    """Exact machinery for  sum_c delta_c * P_c(p) == 0,  |delta_c| <=
    cap_c, delta_c == cap_c (mod 2), with each cell polynomial P_c
    normalized to leading (lowest-degree) coefficient +1. Cells grouped by
    leading valuation; fixing classes in ascending valuation order
    finalizes that coefficient as the class is fixed."""

    def __init__(self, m_lo, m_hi, g1, g2, D0, antisym=False):
        raw = system_cells(m_lo, m_hi, g1, g2, D0)
        self.meta = []      # (side,m,k,z,o)
        self.polys = []     # dict O -> int
        self.caps = []
        self.free_cells = []   # antisym z==o cells (drop out; delta free)
        if not antisym:
            for (s, m, k, z, o, cap) in raw:
                self.meta.append((s, m, k, z, o))
                self.polys.append(cell_poly(z, o))
                self.caps.append(cap)
        else:
            for (s, m, k, z, o, cap) in raw:
                if s == 1:
                    continue
                if z == o:
                    self.free_cells.append((m, k, cap))
                    continue
                P = dict(cell_poly(z, o))
                for O, c in cell_poly(o, z).items():
                    P[O] = P.get(O, 0) - c
                P = {O: c for O, c in P.items() if c}
                lead = min(P)
                if P[lead] < 0:
                    P = {O: -c for O, c in P.items()}
                assert P[min(P)] == 1
                self.meta.append((s, m, k, z, o))
                self.polys.append(P)
                self.caps.append(cap)
        tag = "antisym " if antisym else ""
        self.label = f"{tag}rows[{m_lo}..{m_hi}] g={g1}/{g2} D0={D0}"
        self.antisym = antisym
        self.Amax = max(max(P) for P in self.polys)
        self.zclasses = {}
        for idx, P in enumerate(self.polys):
            self.zclasses.setdefault(min(P), []).append(idx)
        self.zs = sorted(self.zclasses)
        nz = len(self.zs)
        self.B = [dict() for _ in range(nz + 1)]
        for j in range(nz - 1, -1, -1):
            self.B[j] = dict(self.B[j + 1])
            for idx in self.zclasses[self.zs[j]]:
                cap = self.caps[idx]
                for O, c in self.polys[idx].items():
                    self.B[j][O] = self.B[j].get(O, 0) + cap * abs(c)

    @property
    def cells(self):
        # legacy views used by callers: (side,m,k,z,o,cap)
        return [self.meta[i] + (self.caps[i],) for i in range(len(self.caps))]

    def verify(self, delta):
        """Exact verification: returns True iff delta closes the system."""
        assert len(delta) == len(self.caps)
        acc = {}
        for idx, P in enumerate(self.polys):
            d = delta[idx]
            cap = self.caps[idx]
            assert abs(d) <= cap and (d - cap) % 2 == 0, (idx, d, cap)
            for O, c in P.items():
                acc[O] = acc.get(O, 0) + d * c
        return all(v == 0 for v in acc.values())

    # ---- class distribution generation ----
    def count_class_distributions(self, caps, need):
        """DP count of integer points: delta_i in [-cap,cap], parity of cap,
        sum = need."""
        from collections import defaultdict
        cur = {0: 1}
        for cap in caps:
            nxt = defaultdict(int)
            start = -cap
            vals = range(start, cap + 1, 2)
            for s, cnt in cur.items():
                for v in vals:
                    nxt[s + v] += cnt
            cur = nxt
        return cur.get(need, 0)

    def class_distributions(self, zj, need, sample=None, rng=None):
        """Distributions of `need` over the z-class: exhaustive list if
        sample is None (memoized), else up to `sample` random ones."""
        idxs = self.zclasses[zj]
        caps = [self.caps[i] for i in idxs]
        n = len(idxs)
        sufcap = [0] * (n + 1)
        for i in range(n - 1, -1, -1):
            sufcap[i] = sufcap[i + 1] + caps[i]
        if sample is None:
            if not hasattr(self, '_dist_cache'):
                self._dist_cache = {}
            key = (zj, need)
            hit = self._dist_cache.get(key)
            if hit is not None:
                return hit
            out = []
            def rec(i, rem, cur):
                if i == n:
                    if rem == 0:
                        out.append(tuple(cur))
                    return
                if abs(rem) > sufcap[i]:
                    return
                cap = caps[i]
                start = -cap
                for v in range(start, cap + 1, 2):
                    cur.append(v); rec(i + 1, rem - v, cur); cur.pop()
            rec(0, need, [])
            if len(self._dist_cache) < 4000:
                self._dist_cache[key] = out
            return out
        # random sampling biased toward small |delta|
        rng = rng or random
        seen = set()
        out = []
        attempts = 0
        while len(out) < sample and attempts < 30 * sample:
            attempts += 1
            cur = []
            rem = need
            ok = True
            for i in range(n):
                cap = caps[i]
                lo = max(-cap, rem - sufcap[i + 1])
                hi = min(cap, rem + sufcap[i + 1])
                lo += (lo - cap) % 2
                if lo > hi:
                    ok = False
                    break
                nch = (hi - lo) // 2 + 1
                if nch <= 4096:
                    choices = list(range(lo, hi + 1, 2))
                    weights = [1.0 / (1.0 + abs(v)) ** 1.5
                               for v in choices]
                    v = rng.choices(choices, weights=weights)[0]
                else:
                    # heavy-tailed magnitude sample toward small |v|,
                    # clipped into [lo, hi], parity-corrected
                    mag = int((1.0 / max(rng.random(), 1e-12)) ** 1.5) - 1
                    v = mag if rng.random() < 0.5 else -mag
                    v = min(max(v, lo), hi)
                    if (v - lo) % 2:
                        v = v - 1 if v - 1 >= lo else v + 1
                cur.append(v)
                rem -= v
            if ok and rem == 0:
                t = tuple(cur)
                if t not in seen:
                    seen.add(t)
                    out.append(t)
        return out

    def dfs_exhaustive(self, node_cap=20_000_000, shuffle=False,
                       sample=None, rng=None):
        """Search. sample=None: exhaustive (verdict INFEASIBLE is a proof
        if not CAPped). sample=k: randomized incomplete (k distributions
        per class visit) -- only FEASIBLE verdicts meaningful."""
        zs = self.zs; nz = len(zs)
        resid = {}          # running coefficients from fixed classes
        assign = {}         # idx -> delta
        nodes = 0
        capped = False
        def prune(j):
            Bj = self.B[j]
            frontier = zs[j] if j < nz else self.Amax + 1
            for O, v in resid.items():
                if O < frontier:
                    if v != 0:
                        return True   # finalized nonzero coefficient
                else:
                    if abs(v) > Bj.get(O, 0):
                        return True
            return False
        def rec(j):
            nonlocal nodes, capped
            if capped:
                return None
            nodes += 1
            if nodes > node_cap:
                capped = True
                return None
            if prune(j):
                return None
            if j == nz:
                if all(v == 0 for v in resid.values()):
                    return dict(assign)
                return None
            zj = zs[j]
            need = -resid.get(zj, 0)
            dists = self.class_distributions(zj, need, sample=sample,
                                             rng=rng)
            if shuffle and sample is None:
                dists = list(dists)
                random.shuffle(dists)
            idxs = self.zclasses[zj]
            for dist in dists:
                for i, v in zip(idxs, dist):
                    assign[i] = v
                    if v:
                        for O, c in self.polys[i].items():
                            resid[O] = resid.get(O, 0) + v * c
                res = rec(j + 1)
                for i, v in zip(idxs, dist):
                    del assign[i]
                    if v:
                        for O, c in self.polys[i].items():
                            resid[O] = resid.get(O, 0) - v * c
                if res is not None:
                    return res
            return None
        res = rec(0)
        if res is not None:
            delta = [res[i] for i in range(len(self.caps))]
            assert self.verify(delta)
            return ('FEASIBLE', delta, nodes)
        return (('CAP' if capped else 'INFEASIBLE'), None, nodes)

    def max_class_count(self):
        """Max over z-classes of the distribution count for |need| <= 2
        (proxy for exhaustive branching)."""
        worst = 0
        for zj in self.zs:
            caps = [self.caps[i] for i in self.zclasses[zj]]
            worst = max(worst, self.count_class_distributions(caps, 0))
        return worst

def system_size_report(sys_):
    ncells = len(sys_.cells)
    caps = [c[5] for c in sys_.cells]
    nodd = sum(1 for c in caps if c % 2)
    return (f"{sys_.label}: {ncells} cells ({nodd} forced-odd), "
            f"max cap {max(caps)}, levels<= {sys_.Amax}, "
            f"{len(sys_.zs)} z-classes")

# ------------------------------------------------------------- mode dfs -----

def mode_dfs(node_cap=20_000_000):
    print("== dfs: exhaustive verdicts on small closure systems ==")
    gammas = [(1,1),(1,2),(1,3),(2,5),(3,8),(4,11),(3,5),(3,4),(2457,6592)]
    results = {}
    for (g1, g2) in gammas:
        for D0 in (0,):
            for (tag, lo, hi) in [("base r0=1", 1, 1), ("base r0=2", 1, 3),
                                  ("epoch[2,3]", 2, 3),
                                  ("epoch[4,7]", 4, 7),
                                  ("base r0=3", 1, 7)]:
                S = ClosureSystem(lo, hi, g1, g2, D0)
                verdict, delta, nodes = S.dfs_exhaustive(node_cap=node_cap)
                results[(g1, g2, D0, tag)] = (verdict, nodes)
                extra = ""
                if verdict == 'FEASIBLE':
                    nz = sum(1 for d in delta if d)
                    extra = f" (witness support {nz}/{len(delta)})"
                print(f"  g={g1}/{g2} D0={D0} {tag:12s}: {verdict:10s} "
                      f"nodes={nodes}{extra}")
    return results

# ------------------------------------------------------------- mode solve ---

def solve_system(S, tries=40, node_cap=300_000, samples=(4, 12, 40),
                 exhaust_cap=6_000_000, seed=12345, verbose=False):
    """Try exhaustive DFS first if branching is modest; else randomized
    sampled DFS restarts. Returns (verdict, delta_or_None, detail).
    verdict in {FEASIBLE, INFEASIBLE(proof), UNRESOLVED}."""
    rng = random.Random(seed)
    branching = S.max_class_count()
    if branching <= 20000:
        verdict, delta, nodes = S.dfs_exhaustive(node_cap=exhaust_cap)
        if verdict in ('FEASIBLE', 'INFEASIBLE'):
            return verdict, delta, f"exhaustive nodes={nodes}"
        if verbose:
            print(f"    (exhaustive capped at {nodes} nodes)")
    for t in range(tries):
        sample = samples[min(t * len(samples) // max(tries, 1),
                             len(samples) - 1)]
        r2 = random.Random(rng.random())
        verdict, delta, nodes = S.dfs_exhaustive(
            node_cap=node_cap, sample=sample, rng=r2)
        if verdict == 'FEASIBLE':
            return 'FEASIBLE', delta, f"sampled try={t} sample={sample}"
    return 'UNRESOLVED', None, f"{tries} sampled tries"

def mode_solve(spec=None, tries=40, node_cap=300_000):
    print("== solve: exact witness search on larger epochs ==")
    targets = spec or [
        (1, 1, 0, 8, 15), (1, 2, 0, 8, 15), (1, 3, 0, 8, 15),
        (2, 5, 0, 8, 15), (3, 5, 0, 8, 15), (3, 8, 0, 8, 15),
        (4, 11, 0, 8, 15), (2457, 6592, 0, 8, 15),
        (1, 2, 0, 16, 31), (1, 2, 2, 8, 15),
    ]
    out = {}
    for (g1, g2, D0, lo, hi) in targets:
        S = ClosureSystem(lo, hi, g1, g2, D0)
        print(" ", system_size_report(S))
        verdict, delta, detail = solve_system(S, tries=tries,
                                              node_cap=node_cap)
        if verdict == 'FEASIBLE':
            assert S.verify(delta)
            nz = sum(1 for d in delta if d)
            detail += f" support {nz}/{len(delta)}"
        print(f"    -> {verdict} ({detail})")
        out[(g1, g2, D0, lo, hi)] = (verdict, delta)
    return out

# ---------------------------------------------------------- mode depth0 -----

def mode_depth0():
    """Sanity anchor: depth-0 (gamma=0) closure must FAIL at some epoch,
    else C*=1 contra THM-2967. Exhaustive verdicts (caps all 1)."""
    print("== depth0: gamma=0 (all d_m = 0) epoch closures, exhaustive ==")
    for (lo, hi) in [(1, 1), (2, 3), (4, 7), (8, 15), (16, 31)]:
        S = ClosureSystem(lo, hi, 0, 1, 0)
        verdict, delta, nodes = S.dfs_exhaustive(node_cap=50_000_000)
        print(f"  rows[{lo},{hi}] d=0: {verdict} (nodes={nodes})")
        if verdict == 'FEASIBLE':
            print(f"      delta = {delta}")

# --------------------------------------------------------- mode threshold ---

def mode_threshold(r_list=(2, 3, 4), tries=40, node_cap=300_000,
                   exhaust_cap=3_000_000):
    """Measure the construction gate: for epoch [2^r, 2^{r+1}-1], scan
    gamma grid; report FEASIBLE (exact witness) / INFEASIBLE (proof) /
    UNRESOLVED per gamma. Compare threshold with alpha = 2457/6592."""
    print("== threshold: gamma*(r) hunt for epoch closure ==", flush=True)
    grid = [(0,1),(1,20),(1,10),(1,8),(1,6),(1,5),(1,4),(2,7),(3,10),
            (1,3),(4,11),(3,8),(2457,6592),(2,5),(5,12),(3,7),(1,2),
            (3,5),(3,4),(1,1)]
    for r in r_list:
        lo, hi = 2 ** r, 2 ** (r + 1) - 1
        print(f"  epoch [{lo},{hi}]:", flush=True)
        for (g1, g2) in grid:
            S = ClosureSystem(lo, hi, g1, g2, 0)
            verdict, delta, detail = solve_system(
                S, tries=tries, node_cap=node_cap, exhaust_cap=exhaust_cap)
            gap = f"{g1}/{g2}"
            sup = ""
            if verdict == 'FEASIBLE':
                sup = f" support {sum(1 for d in delta if d)}/{len(delta)}"
            print(f"    gamma={gap:10s}: {verdict:10s} ({detail}){sup}",
                  flush=True)

# ---------------------------------------------------------- mode inspect ----

def mode_inspect(specs=None, tries=40, node_cap=300_000):
    """Find witnesses for given systems and print per-row structure +
    telescoping ratios at the certificate biases (exact Fractions)."""
    print("== inspect: witness anatomy + telescoping hunt ==")
    pA = Fraction(1285, 2181); pB = Fraction(8847357, 11821757)
    specs = specs or [(1, 2, 0, 4, 7), (1, 2, 0, 8, 15),
                      (2457, 6592, 0, 4, 7), (1, 3, 0, 8, 15)]
    for (g1, g2, D0, lo, hi) in specs:
        S = ClosureSystem(lo, hi, g1, g2, D0)
        verdict, delta, detail = solve_system(S, tries=tries,
                                              node_cap=node_cap)
        print(f"  {S.label}: {verdict} ({detail})")
        if verdict != 'FEASIBLE':
            continue
        assert S.verify(delta)
        # per-row printout and per-row polynomial value at biases
        rows = {}
        for idx, (s, m, k, z, o) in enumerate(S.meta):
            rows.setdefault(m, []).append((s, k, delta[idx], z, o))
        cumA = Fraction(0); cumB = Fraction(0)
        seqA = []; seqB = []
        for m in sorted(rows):
            ent = sorted(rows[m])
            side0 = [(k, d) for (s, k, d, _z, _o) in ent if s == 0]
            side1 = [(k, d) for (s, k, d, _z, _o) in ent if s == 1]
            anti = all(dict(side1).get(k, None) == -d for k, d in side0)
            vA = sum(d * pA ** z * (1 - pA) ** o
                     for (s, k, d, z, o) in ent)
            vB = sum(d * pB ** z * (1 - pB) ** o
                     for (s, k, d, z, o) in ent)
            cumA += vA; cumB += vB
            seqA.append(cumA); seqB.append(cumB)
            print(f"    m={m:3d} side0={side0} side1={side1} "
                  f"antisym={anti}")
        # telescoping hunt: ratios of consecutive cumulative row-values
        def ratio_report(seq, name):
            for i in range(1, len(seq)):
                if seq[i - 1] != 0 and seq[i] != 0:
                    r = seq[i] / seq[i - 1]
                    if r.denominator < 10 ** 7:
                        print(f"      {name} cum ratio row{i}: {r}")
        ratio_report(seqA, "p_A"); ratio_report(seqB, "p_B")

# ---------------------------------------------------------- mode abeldini ---

TARGETS = {
    "S_A-": 896, "S_A": 1285, "x_A": 389, "b_A": 2181,
    "S_B-": 2974400, "S_B": 8847357, "x_B": 5872957, "b_B": 11821757,
    "alpha_num": 2457, "alpha_den": 6592, "C_num": 9049,
}

def factorize(n):
    try:
        from sympy import factorint
        return dict(factorint(n))
    except Exception:
        f = {}
        d = 2
        while d * d <= n:
            while n % d == 0:
                f[d] = f.get(d, 0) + 1
                n //= d
            d += 1
        if n > 1:
            f[n] = f.get(n, 0) + 1
        return f

def mode_abeldini():
    print("== abeldini: integer hunts for the certificate pairs ==")
    print("  factorizations:")
    for name, n in TARGETS.items():
        print(f"    {name:9s} = {n} = {factorize(n)}")
    tgt_pairs = [(896, 1285), (2974400, 8847357)]
    tgt_set = set(TARGETS.values())

    print("  [A] consecutive checkpoint cumulative budgets "
          "S_r = sum_(m<=2^r-1) 2^(d_m):")
    hitsA = []
    for g2 in range(1, 49):
        for g1 in range(1, g2 + 1):
            if math.gcd(g1, g2) != 1:
                continue
            for D0 in range(0, 7):
                S = 0; m = 0; r = 1
                vals = []
                dead = False
                while True:
                    lim = 2 ** r - 1
                    while m < lim:
                        m += 1
                        d = depth(m, g1, g2, D0)
                        if d > 30 or m > 2 ** 17:
                            dead = True
                            break
                        S += 2 ** d
                        if S > 2 * 10 ** 7:
                            dead = True
                            break
                    vals.append(S)
                    if dead or r > 17:
                        break
                    r += 1
                for a, b in zip(vals, vals[1:]):
                    if (a, b) in tgt_pairs:
                        hitsA.append((g1, g2, D0, a, b))
                for v in vals:
                    if v in tgt_set:
                        hitsA.append((g1, g2, D0, 'single', v))
    print(f"      hits: {hitsA if hitsA else 'NONE'}")

    print("  [B] binomial partial sums sum_(k<=K) binom(n,k):")
    hitsB = []
    for n in range(1, 80):
        S = 0
        for K in range(n + 1):
            S += comb(n, K)
            if S in tgt_set:
                hitsB.append((n, K, S))
    print(f"      hits: {hitsB if hitsB else 'NONE'}")

    print("  [C] structured products a^m (b-a)^j * prod_(i in I)"
          "(a^(2^i)+(b-a)^(2^i)) == target:")
    hitsC = []
    for b in range(2, 25):
        for a in range(1, b):
            if math.gcd(a, b) != 1:
                continue
            c = b - a
            fer = [a ** (2 ** i) + c ** (2 ** i) for i in range(5)]
            for mask in range(32):
                P = 1
                for i in range(5):
                    if mask >> i & 1:
                        P *= fer[i]
                if P > 10 ** 8:
                    continue
                for me in range(0, 31):
                    Pa = P * a ** me
                    if Pa > 10 ** 8:
                        break
                    Q = Pa
                    for j in range(0, 31):
                        if Q in tgt_set:
                            d_bits = [2 ** i for i in range(5)
                                      if mask >> i & 1]
                            hitsC.append((f"p={a}/{b}", f"a^{me}", f"c^{j}",
                                          f"bits{d_bits}", Q))
                        Q *= c
                        if Q > 10 ** 8:
                            break
    # dedupe
    seen = set(); uniq = []
    for h in hitsC:
        if h not in seen:
            seen.add(h); uniq.append(h)
    print(f"      {len(uniq)} hits; targets found: "
          f"{sorted(set(h[-1] for h in uniq))}")
    for h in uniq[:40]:
        print(f"        {h}")
    if len(uniq) > 40:
        print(f"        ... ({len(uniq)-40} more)")

    print("  [D] consecutive cumulative FORCED-ODD cell counts and cell "
          "totals at checkpoints:")
    hitsD = []
    for g2 in range(1, 49):
        for g1 in range(1, g2 + 1):
            if math.gcd(g1, g2) != 1:
                continue
            for D0 in range(0, 5):
                for kind in ('odd', 'cells', 'corners'):
                    S = 0; m = 0; vals = []
                    for r in range(1, 18):
                        lim = 2 ** r - 1
                        dead = False
                        while m < lim:
                            m += 1
                            d = depth(m, g1, g2, D0)
                            if kind == 'odd':
                                S += 2 * 2 ** min(bin(d).count('1'), 30)
                            elif kind == 'cells':
                                S += 2 * (d + 1)
                            else:
                                S += 2
                            if S > 2 * 10 ** 7 or m > 2 ** 17:
                                dead = True
                                break
                        vals.append(S)
                        if dead:
                            break
                    for x, y in zip(vals, vals[1:]):
                        if (x, y) in tgt_pairs:
                            hitsD.append((kind, g1, g2, D0, x, y))
                    for v in vals:
                        if v in tgt_set and v > 2000:
                            hitsD.append((kind, g1, g2, D0, 'single', v))
    print(f"      hits: {hitsD if hitsD else 'NONE'}")
    return dict(A=hitsA, B=hitsB, C=uniq, D=hitsD)

# --------------------------------------------------------------- driver -----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="all")
    ap.add_argument("--tries", type=int, default=60)
    ap.add_argument("--nodecap", type=int, default=20_000_000)
    ap.add_argument("--spec", type=str, default=None,
                    help="g1,g2,D0,lo,hi;g1,g2,D0,lo,hi;...")
    args = ap.parse_args()
    spec = None
    if args.spec:
        spec = [tuple(int(x) for x in s.split(","))
                for s in args.spec.split(";")]
    if args.mode in ("parity", "all"):
        mode_parity()
    if args.mode in ("dfs", "all"):
        mode_dfs(node_cap=args.nodecap)
    if args.mode == "solve":
        mode_solve(spec=spec, tries=args.tries, node_cap=args.nodecap)
    if args.mode == "depth0":
        mode_depth0()
    if args.mode == "inspect":
        mode_inspect(specs=spec, tries=args.tries, node_cap=args.nodecap)
    if args.mode == "threshold":
        rl = tuple(int(x) for x in (args.spec or "2,3,4").split(","))
        mode_threshold(r_list=rl, tries=args.tries, node_cap=args.nodecap)
    if args.mode in ("abeldini", "all"):
        mode_abeldini()
    print("DONE")

if __name__ == "__main__":
    main()
