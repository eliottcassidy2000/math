#!/usr/bin/env python3
"""
THE [43,48] ESCAPE-WINDOW FINITE STATEMENT -- prototype enumerator.
Compute agent, 2026-07-27.  Scratch only; nothing written to canon.

PART A: independent re-verification of klein's corrected 59-row escape atlas
        (S323c, frozen in 05-knowledge/results/lrc14_subgap_census_klein_S319.out):
        filters (covering 2..13, primitivity, spread>12, rung pair-sum) +
        escape at every q in [14,41] + first witness q* in [43,48] == stated q.
PART B: THE WINDOW EXTENSION of kps's HYP-8040 escape atlas: E_q for
        q in [42,48] x 78 two-far complements, with the HYP-8050 three-family
        classification (DIV / IMP / INV / LINE / REST); clean-modulus verdict.
        Cross-check vs kps's frozen [14,41] rows for S = 13\{1,2}.
PART C: the surviving channel, exhibited:
        C1  kps c96 witness {3..13, 1+L41, 2+L41} escapes [14,41], exits at 43 (j=2).
        C2  the window-closing analogue {3..13, 1+L48, 2+L48}, L48 = lcm(14..48):
            escapes ALL q in [14,48]; first exit located in [49, 200].
        C3  EXACT impersonation x-floor for complement {1,2} over [14,48]:
            min x >= 14 with x mod q in {+-1, +-2} for every q in [14,48].
        C4  EXACT divisibility-channel floor: min max(x,y) with every
            q in [14,48] dividing x or y (prime-power linkage argument).
Conventions == kps lrc14_small_witness_law_kps_S128c95.py:
        t_q = ceil(3q/41); witness j in [1,q-1] needs dist_q(v*j) >= t_q for all v.
"""
import sys, time
from math import gcd
from itertools import combinations, product
from functools import reduce

T0 = time.time()

def tq(q):  return -(-3*q // 41)          # ceil(3q/41)
def dist(v, j, q):
    r = (v*j) % q
    return min(r, q-r)

def has_witness(W, q):
    """smallest-j witness at modulus q, or None (j-range [1, q//2] by symmetry)."""
    t = tq(q)
    for j in range(1, q//2 + 1):
        if all(dist(v, j, q) >= t for v in W):
            return j
    return None

def first_exit(W, qlo=14, qhi=200):
    for q in range(qlo, qhi+1):
        j = has_witness(W, q)
        if j is not None:
            return q, j
    return None, None

def lcm(a, b): return a//gcd(a, b)*b
def lcm_range(lo, hi): return reduce(lcm, range(lo, hi+1), 1)

def rung_ok(s):
    """exists integer D with s/14 < D < 3s/41"""
    Dmin = s//14 + 1
    return 41*Dmin < 3*s

# ----------------------------------------------------------------------------
ATLAS = """37 43 1 10 17 19 22 23 24 25 26 27 29 35 37
37 44 2 20 22 23 24 25 26 27 28 31 34 35 37
38 43 1 10 17 22 23 24 25 26 27 29 35 37 38
38 45 1 9 10 13 22 23 24 25 26 34 35 37 38
39 47 3 9 21 23 24 25 26 27 29 30 33 36 39
40 45 3 22 23 24 25 26 27 28 29 31 34 37 40
41 45 1 20 21 22 23 24 25 26 27 28 34 35 41
41 45 3 20 22 23 25 26 29 32 35 36 37 38 41
47 43 1 10 22 23 24 25 26 27 34 35 37 40 47
48 43 1 2 10 22 23 25 26 27 29 35 37 38 48
48 43 1 10 15 22 23 25 26 27 29 35 37 38 48
48 43 1 10 17 22 23 25 26 27 29 35 37 38 48
48 43 1 10 22 23 25 26 27 29 32 35 37 38 48
49 43 1 10 22 23 25 26 27 29 37 38 47 48 49
49 45 1 10 22 23 25 26 27 29 37 38 43 48 49
49 45 1 22 23 24 25 26 27 29 30 31 37 38 49
49 45 1 22 23 25 26 27 29 30 31 37 38 48 49
50 43 1 11 26 27 37 38 41 44 46 47 48 49 50
50 43 1 17 24 26 27 29 31 38 41 42 44 46 50
50 43 1 17 26 27 29 31 38 41 42 44 46 48 50
50 43 1 22 23 24 26 27 29 31 32 37 38 49 50
50 43 1 22 23 24 26 27 31 32 37 38 47 49 50
50 43 1 22 23 26 27 29 31 32 37 38 48 49 50
50 43 2 7 8 22 23 26 29 31 36 37 38 41 50
50 43 2 17 24 26 27 29 31 38 41 42 44 46 50
50 43 2 17 26 27 29 31 38 41 42 44 46 48 50
50 43 3 7 8 22 23 26 29 32 36 37 38 41 50
50 43 3 7 22 23 26 29 32 35 36 37 38 41 50
50 43 3 20 22 23 26 29 32 35 36 37 38 41 50
50 43 3 22 23 26 29 32 35 36 37 38 41 47 50
50 45 3 9 12 20 23 26 29 32 35 38 41 44 50
50 45 3 9 12 22 23 26 27 29 32 35 38 41 50
50 45 3 9 12 23 26 27 29 32 35 38 41 44 50
50 45 3 9 12 23 26 29 32 35 38 41 44 47 50
50 45 3 12 18 20 23 26 29 32 35 38 41 44 50
50 45 3 12 18 22 23 26 27 29 32 35 38 41 50
50 45 3 12 19 22 23 26 27 29 32 35 38 41 50
50 45 3 12 22 23 26 27 29 32 35 37 38 41 50
50 45 3 12 22 23 26 27 29 32 35 38 41 46 50
50 45 3 22 23 24 26 27 29 32 35 37 38 41 50
50 45 3 22 23 26 27 29 32 35 36 37 38 41 50
50 45 3 22 23 26 27 29 32 35 37 38 41 48 50
50 45 4 14 22 23 26 27 29 32 36 37 38 41 50
50 45 4 22 23 26 27 28 29 32 36 37 38 41 50
50 46 3 12 22 23 26 27 28 29 32 35 38 41 50
52 46 1 4 8 12 16 20 28 32 36 40 44 48 52
52 48 4 7 8 12 16 20 24 28 32 36 40 44 52
52 48 4 8 12 16 20 24 27 28 32 36 40 44 52
52 48 4 8 12 16 20 24 28 32 36 40 41 44 52
53 43 1 10 22 23 26 27 37 38 47 48 49 50 53
53 45 1 10 22 23 26 27 37 38 43 48 49 50 53
54 43 1 5 19 20 22 23 24 25 34 35 45 52 54
54 43 1 13 17 25 29 30 31 33 41 42 46 48 54
54 43 1 13 17 29 30 31 33 41 42 46 48 50 54
54 43 1 17 19 26 29 31 33 35 37 46 48 50 54
54 43 1 17 25 26 29 30 31 33 41 42 46 48 54
54 47 2 17 19 29 31 33 42 44 46 48 50 52 54
55 45 3 12 22 23 26 27 29 32 35 38 41 50 55
55 48 4 8 12 16 20 24 28 32 36 40 44 52 55"""

def part_A():
    print("="*78)
    print("PART A -- independent verification of klein's 59-row escape atlas (S323c)")
    print("="*78)
    rows = []
    for line in ATLAS.strip().split("\n"):
        xs = list(map(int, line.split()))
        rows.append((xs[0], xs[1], xs[2:]))
    assert len(rows) == 59, len(rows)
    hist, band, kfar, rungs = {}, {}, {}, {}
    nfail = 0
    for (w, qstated, W) in rows:
        assert len(W) == 13 and len(set(W)) == 13
        ok_w    = (max(W) == w)
        ok_cov  = all(any(v % d == 0 for v in W) for d in range(2, 14))
        ok_prim = reduce(gcd, W) == 1
        ok_spr  = max(W) > 12*min(W)
        sums    = {a+b for a, b in combinations(W, 2)}
        ok_rung = any(rung_ok(s) for s in sums)
        ok_esc  = all(has_witness(W, q) is None for q in range(14, 42))
        qstar, jstar = first_exit(W, 42, 200)
        ok_q    = (qstar == qstated)
        allok = ok_w and ok_cov and ok_prim and ok_spr and ok_rung and ok_esc and ok_q
        if not allok:
            nfail += 1
            print(f"  FAIL w={w} q={qstated} W={W}: w{ok_w} cov{ok_cov} prim{ok_prim} "
                  f"spr{ok_spr} rung{ok_rung} esc{ok_esc} q*={qstar}")
        hist[qstar] = hist.get(qstar, 0) + 1
        band[qstar - w] = band.get(qstar - w, 0) + 1
        k = sum(1 for v in W if v >= 14)
        kfar[k] = kfar.get(k, 0) + 1
        rs = sorted(s for s in sums if rung_ok(s))
        rungs[tuple(rs)] = rungs.get(tuple(rs), 0) + 1
    print(f"  rows verified: {59-nfail}/59 pass ALL of [w=max, covering2..13, primitive,")
    print(f"    spread>12, rung-pair-exists, escape at every q in [14,41], first-witness==stated]")
    print(f"  first-witness histogram (mine): "
          + " ".join(f"{q}:{n}" for q, n in sorted(hist.items())))
    print(f"    klein's claim:               43:27 44:1 45:23 46:2 47:2 48:4")
    print(f"  q - v_max band (mine): "
          + " ".join(f"{d}:{n}" for d, n in sorted(band.items())))
    print(f"    klein's claim:       -11:5 -10:2 -8:1 -7:15 -6:2 -5:18 -4:8 4:2 5:2 6:1 7:2 8:1")
    print(f"  k-far profile of the 59 rows (k = #elements >= 14): "
          + " ".join(f"k={k}:{n}" for k, n in sorted(kfar.items())))
    print(f"  admissible rung sums present per row (top shapes): ")
    for rs, n in sorted(rungs.items(), key=lambda kv: -kv[1])[:5]:
        print(f"    {list(rs)}: {n} rows")
    print(f"  [{time.time()-T0:.1f}s]")
    return rows

# ----------------------------------------------------------------------------
def escape_table(S, q):
    """A_q, E_q for two-far complement family S u {x,y}; returns
    (A, esc_rows) where esc_rows[u] = bitmask of v with (u,v) in E_q;
    A empty -> useless (E_q = everything)."""
    t = tq(q)
    A = [j for j in range(1, q) if all(dist(s, j, q) >= t for s in S)]
    full = (1 << q) - 1
    if not A:
        return A, [full]*q
    Cj = {}
    for j in A:
        m = 0
        for u in range(q):
            if dist(u, j, q) >= t:
                m |= (1 << u)
        Cj[j] = m
    esc_rows = []
    for u in range(q):
        cov = 0
        bu = 1 << u
        for j in A:
            if Cj[j] & bu:
                cov |= Cj[j]
        esc_rows.append(full & ~cov)
    return A, esc_rows

def classify(ab, q, esc_rows, A):
    """count DIV / IMP / INV / LINE / REST over E_q (A nonempty assumed)."""
    a, b = ab
    Ra = {a % q, (q - a) % q}
    Rb = {b % q, (q - b) % q}
    inv_pairs = set()
    for c in (a, b):
        if gcd(c, q) == 1:
            ci = pow(c, -1, q)
            for u in (c % q, (q-c) % q):
                for v in (ci, (q-ci) % q):
                    inv_pairs.add((u, v)); inv_pairs.add((v, u))
    n = {'DIV':0, 'IMP':0, 'INV':0, 'LINE':0, 'REST':0}
    # detect blocking-class lines: rows/cols nearly full
    rowsz = [bin(esc_rows[u]).count("1") for u in range(q)]
    colmask = [0]*q
    for u in range(q):
        m = esc_rows[u]
        while m:
            v = (m & -m).bit_length() - 1
            colmask[v] |= (1 << u)
            m &= m-1
    colsz = [bin(colmask[v]).count("1") for v in range(q)]
    urow_line = {u for u in range(q) if u != 0 and rowsz[u] >= q-6}
    vcol_line = {v for v in range(q) if v != 0 and colsz[v] >= q-6}
    tot = 0
    for u in range(q):
        m = esc_rows[u]
        while m:
            v = (m & -m).bit_length() - 1
            m &= m-1
            tot += 1
            if u == 0 or v == 0:                              n['DIV']  += 1
            elif (u in Ra and v in Rb) or (u in Rb and v in Ra): n['IMP'] += 1
            elif (u, v) in inv_pairs:                          n['INV']  += 1
            elif u in urow_line or v in vcol_line:             n['LINE'] += 1
            else:                                              n['REST'] += 1
    return tot, n

def part_B():
    print()
    print("="*78)
    print("PART B -- THE WINDOW EXTENSION: escape sets E_q, q in [42,48] x 78 complements")
    print("="*78)
    # cross-check my table builder against kps's frozen [14,41] rows for 13\{1,2}
    S12 = [v for v in range(1, 14) if v not in (1, 2)]
    want = {15:52, 16:56, 17:8, 18:64, 19:8, 20:72}
    got = {}
    for q in sorted(want):
        A, er = escape_table(S12, q)
        nondiv = sum(bin(er[u] & ~1).count("1") for u in range(1, q))
        got[q] = nondiv
    A14, _ = escape_table(S12, 14)
    print(f"  cross-check vs kps frozen atlas, S = 13\\{{1,2}}:")
    print(f"    useless at q=14: {'YES' if not A14 else 'NO'} (kps: YES)")
    print(f"    non-div escape counts {got}")
    print(f"    kps frozen           {want}  -> {'MATCH' if got == want else 'MISMATCH'}")

    WINDOW = range(42, 49)
    agg = {q: {'E':0, 'nondiv':0, 'useless':0, 'clean':0,
               'DIV':0,'IMP':0,'INV':0,'LINE':0,'REST':0} for q in WINDOW}
    q43_rows = []
    grand_cells = 0
    for ab in combinations(range(1, 14), 2):
        S = [v for v in range(1, 14) if v not in ab]
        for q in WINDOW:
            grand_cells += q*q
            A, er = escape_table(S, q)
            if not A:
                agg[q]['useless'] += 1
                if q == 43: q43_rows.append((ab, 'USELESS', 0, 0, None))
                continue
            tot, n = classify(ab, q, er, A)
            nondiv = tot - n['DIV']
            agg[q]['E'] += tot; agg[q]['nondiv'] += nondiv
            for k in ('DIV','IMP','INV','LINE','REST'): agg[q][k] += n[k]
            if nondiv == 0: agg[q]['clean'] += 1
            if q == 43:
                q43_rows.append((ab, len(A), tot, nondiv, n))
    print(f"\n  window aggregate over 78 complements (cells scanned: {grand_cells:,}):")
    print(f"  {'q':>3} {'t_q':>3} {'useless':>8} {'clean':>6} {'sum|E_q|':>9} "
          f"{'nondiv':>7} {'DIV':>6} {'IMP':>5} {'INV':>5} {'LINE':>6} {'REST':>6}")
    for q in WINDOW:
        a = agg[q]
        print(f"  {q:>3} {tq(q):>3} {a['useless']:>8} {a['clean']:>6} {a['E']:>9} "
              f"{a['nondiv']:>7} {a['DIV']:>6} {a['IMP']:>5} {a['INV']:>5} "
              f"{a['LINE']:>6} {a['REST']:>6}")
    nclean_win = sum(agg[q]['clean'] for q in WINDOW)
    print(f"\n  CLEAN WINDOW MODULI (E_q inside divisibility): {nclean_win} "
          f"complement-x-modulus pairs out of {78*7}")
    print(f"  q=43 row class, per complement (first 12 + duty-split examples):")
    print(f"    {'removed':>10} {'|A_43|':>7} {'|E_43|':>7} {'nondiv':>7}  classification")
    shown = 0
    for (ab, nA, tot, nondiv, n) in q43_rows:
        if shown < 12 or ab in ((7,11),(12,13),(6,12),(5,10)):
            if nA == 'USELESS':
                print(f"    {str(ab):>10}  USELESS (A_43 empty: every pair escapes)")
            else:
                print(f"    {str(ab):>10} {nA:>7} {tot:>7} {nondiv:>7}  {n}")
            shown += 1
    tot43 = sum(r[2] for r in q43_rows if r[1] != 'USELESS')
    print(f"  q=43 totals: sum|E_43| = {tot43} over 78 complements "
          f"(raw table = 78 x 43^2 = {78*43*43:,} cells)")
    print(f"  [{time.time()-T0:.1f}s]")
    return agg

# ----------------------------------------------------------------------------
def part_C():
    print()
    print("="*78)
    print("PART C -- the surviving channel past the window, exhibited + exact floors")
    print("="*78)
    L41 = lcm_range(14, 41)
    L48 = lcm_range(14, 48)
    print(f"  L41 = lcm(14..41) = {L41}")
    print(f"  L48 = lcm(14..48) = {L48}")
    S = list(range(3, 14))

    # C1: kps c96 witness
    W1 = S + [1 + L41, 2 + L41]
    esc41 = all(has_witness(W1, q) is None for q in range(14, 42))
    q1, j1 = first_exit(W1, 42, 200)
    print(f"\n  C1  kps c96 witness {{3..13, 1+L41, 2+L41}}: escapes [14,41]: {esc41}; "
          f"exit q*={q1}, j={j1}  (kps claim: exits q=43, j=2)")

    # C2: the window-closing analogue
    W2 = S + [1 + L48, 2 + L48]
    esc48 = all(has_witness(W2, q) is None for q in range(14, 49))
    q2, j2 = first_exit(W2, 49, 200)
    print(f"  C2  WINDOW-CLOSER {{3..13, 1+L48, 2+L48}} (fars ~ {float(1+L48):.3e}):")
    print(f"      escapes EVERY q in [14,48]: {esc48}   "
          f"first exit q*={q2}, j={j2}   nextprime(48)=53")
    print(f"      -> the [43,48] window is NOT absolute at all heights; the ladder "
          f"steps to {q2}.")

    # C3: exact impersonation x-floor for {1,2} over [14,48]
    cluster = [32, 27, 25, 7, 11, 13, 17, 19, 23]
    iso     = [29, 31, 37, 41, 43, 47]
    Mc  = reduce(lambda x, y: x*y, cluster)
    Mi  = reduce(lambda x, y: x*y, iso)
    assert Mc * Mi == L48
    comp_qs = [q for q in range(14, 49) if q not in iso]
    # enumerate 4^9 cluster choices by iterative CRT
    combos = [(0, 1)]   # (residue, modulus)
    for m in cluster:
        R = sorted({1 % m, (m-1) % m, 2 % m, (m-2) % m})
        new = []
        for (r, M) in combos:
            inv = pow(M % m, -1, m)
            for c in R:
                k = ((c - r) * inv) % m
                new.append((r + M*k, M*m))
        combos = new
    surv = []
    for (r, M) in combos:
        assert M == Mc
        if all(r % q in (1, 2, q-1, q-2) for q in comp_qs):
            surv.append(r)
    print(f"\n  C3  impersonation x-channel for {{a,b}}={{1,2}}: x mod q in {{+-1,+-2}} "
          f"for all q in [14,48]")
    print(f"      cluster prime powers {cluster} linked by composite q; "
          f"4^9 = {4**9} sign/target patterns -> {len(surv)} coherent residues mod "
          f"Mc = {Mc}")
    ex = sorted(surv)[:8]
    print(f"      surviving cluster residues (sorted, first 8): {ex}")
    # combine with iso primes: minimal x >= 14
    iso_res = [(0, 1)]
    for p in iso:
        R = (1, 2, p-1, p-2)
        new = []
        for (r, M) in iso_res:
            inv = pow(M % p, -1, p)
            for c in R:
                k = ((c - r) * inv) % p
                new.append((r + M*k, M*p))
        iso_res = new
    invMc = pow(Mc % Mi, -1, Mi)
    best = None
    for r in surv:
        for (s, _) in iso_res:
            k = ((s - r) * invMc) % Mi
            x = r + Mc * k
            if x >= 14 and (best is None or x < best):
                best = x
    print(f"      total channel classes mod L48: {len(surv)} x 4^6 = "
          f"{len(surv)*4096:,}")
    print(f"      EXACT x-floor of the channel: min x >= 14 = {best}")
    print(f"        ~ {float(best):.3e}  (vs drift threshold X0 = 339; "
          f"vs L48 = {float(L48):.3e})")
    # sanity: verify best really lands in +-{1,2} everywhere
    assert all(best % q in (1, 2, q-1, q-2) for q in range(14, 49))

    # C4: exact divisibility-channel floor
    ppws = cluster + iso
    # linkage: q in [14,48] composite forces its prime supports to one side
    import collections
    parent = {p: p for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47)}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y): parent[find(x)] = find(y)
    def primes_of(q):
        ps, d, qq = [], 2, q
        while d*d <= qq:
            if qq % d == 0:
                ps.append(d)
                while qq % d == 0: qq //= d
            d += 1
        if qq > 1: ps.append(qq)
        return ps
    for q in range(14, 49):
        ps = primes_of(q)
        for p in ps[1:]:
            union(ps[0], p)
    comps = collections.defaultdict(list)
    ppw_of = {2:32, 3:27, 5:25, 7:7, 11:11, 13:13, 17:17, 19:19, 23:23,
              29:29, 31:31, 37:37, 41:41, 43:43, 47:47}
    for p in parent: comps[find(p)].append(p)
    blocks = [reduce(lambda x, y: x*y, [ppw_of[p] for p in ps]) for ps in comps.values()]
    blocks.sort(reverse=True)
    bestmax = None
    nb = len(blocks)
    for mask in range(1 << (nb-1)):          # fix largest block on side X
        x = blocks[0]; y = 1
        for i in range(1, nb):
            if mask >> (i-1) & 1: x *= blocks[i]
            else:                 y *= blocks[i]
        m = max(x, y)
        if bestmax is None or m < bestmax: bestmax, bx, by = m, x, y
    print(f"\n  C4  divisibility channel (every q in [14,48] divides x or y):")
    print(f"      linkage components (as ppw products): {blocks}")
    print(f"      EXACT floor: min max(x,y) = {bestmax} ~ {float(bestmax):.3e}")
    print(f"        split: x = {bx}, y = {by}")
    print(f"  [{time.time()-T0:.1f}s]")

if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
    print(f"\nTOTAL {time.time()-T0:.1f}s")
