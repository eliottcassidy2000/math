#!/usr/bin/env python3
"""
adversarial_maxH_check_kpswf.py
================================================================
INDEPENDENT adversarial re-check of HYP-2688 H-maximizer claim.

Written FRESH (does NOT import / reuse half_tiling_maxH_kpswf.py or
half_tiling_framework_kps.py).  Implements its own:
  - tile enumeration  (a,b) with a>=b+2
  - tiling -> tournament adjacency
  - grid reflection rho(a,b) = (n+1-b, n+1-a)  [THM-280]
  - Hamiltonian-path counter (independent bitmask DP)
  - converse_relabel for an INDEPENDENT cross-check that rho-fixed
    tilings are exactly the phi-self-converse tournaments.

QUESTIONS:
  Q1. Does some grid-sym tiling attain global max H at n=8 (=661)?
      Exhaust the 2^half(8)=2^12 grid-sym tilings, report max.
  Q2. Re-confirm grid-sym maxima at n=5,6,7 are 15, 45, 189.
  Q3. Full-cube cross-check at n=5,6,7: global max over ALL 2^m tilings,
      and whether grid-sym attains it; count maximizers & how many GS.
  Q4. n=9 grid-sym max (2^16 grid-sym tilings) -> is it 3357?
"""

import itertools
from math import comb, floor
from collections import Counter

# ---------------------------------------------------------------------------
# tile enumeration  (a > b, a - b >= 2).  Order does not matter for correctness;
# we key tilings by the (a,b) pair directly.
# ---------------------------------------------------------------------------
def tiles(n):
    out = []
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            if a - b >= 2:
                out.append((a, b))
    return out

def half_n(n):
    return floor((n - 1) ** 2 / 4)

# ---------------------------------------------------------------------------
# tiling (dict tile->bit) -> directed adjacency matrix (set of arcs u->v)
#   base path P0 = n -> n-1 -> ... -> 1  (consecutive arcs k -> k-1)
#   tile bit 0 => forward a->b ; bit 1 => backward b->a
# ---------------------------------------------------------------------------
def tiling_to_arcs(t, n):
    """Return adj: adj[u] = set of w with arc u->w. Vertices 1..n."""
    adj = {v: set() for v in range(1, n + 1)}
    # base path: k -> k-1
    for k in range(n, 1, -1):
        adj[k].add(k - 1)
    for (a, b), bit in t.items():
        if bit == 0:        # forward a -> b
            adj[a].add(b)
        else:               # backward b -> a
            adj[b].add(a)
    return adj

def count_ham_paths(adj, n):
    """Independent bitmask DP counting directed Hamiltonian paths."""
    # vertices labelled 1..n ; bit (v-1) in mask
    full = (1 << n) - 1
    # adjacency as bitmask of out-neighbours
    out = [0] * (n + 1)
    for u in range(1, n + 1):
        msk = 0
        for w in adj[u]:
            msk |= (1 << (w - 1))
        out[u] = msk
    # dp[mask] : list over last-vertex of count
    dp = [[0] * (n + 1) for _ in range(1 << n)]
    for v in range(1, n + 1):
        dp[1 << (v - 1)][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        # quick skip
        if mask == 0:
            continue
        for last in range(1, n + 1):
            c = row[last]
            if c == 0:
                continue
            if not (mask & (1 << (last - 1))):
                continue
            reachable = out[last] & ~mask
            w = 1
            rm = reachable
            while rm:
                lb = rm & (-rm)
                wv = lb.bit_length()   # vertex label
                dp[mask | lb][wv] += c
                rm ^= lb
    return sum(dp[full])

# ---------------------------------------------------------------------------
# rho reflection (THM-280): (a,b) -> (n+1-b, n+1-a)
# ---------------------------------------------------------------------------
def rho(a, b, n):
    u, v = n + 1 - b, n + 1 - a   # u>v since b<a
    return (u, v)

# ---------------------------------------------------------------------------
# independent converse_relabel for self-converse cross-check
#   reverse all arcs, then relabel i -> n+1-i
# ---------------------------------------------------------------------------
def converse_relabel_arcs(adj, n):
    new = {v: set() for v in range(1, n + 1)}
    for u in range(1, n + 1):
        for w in adj[u]:
            # arc u->w ; reverse -> w->u ; relabel phi(i)=n+1-i
            pu, pw = n + 1 - u, n + 1 - w
            new[pw].add(pu)
    return new

def arcs_equal(a1, a2, n):
    for v in range(1, n + 1):
        if a1[v] != a2[v]:
            return False
    return True

# ---------------------------------------------------------------------------
# grid-sym enumeration: free bit per rho-orbit, propagate to mirror cell.
# ---------------------------------------------------------------------------
def rho_orbits(n):
    tl = tiles(n)
    tset = set(tl)
    seen = set()
    orbits = []   # each orbit is a list of tiles (1 if fixed, 2 if paired)
    for p in tl:
        if p in seen:
            continue
        img = rho(*p, n)
        assert img in tset, f"rho image {img} of {p} not a tile (n={n})"
        if img == p:
            orbits.append([p])
            seen.add(p)
        else:
            orbits.append([p, img])
            seen.add(p)
            seen.add(img)
    return tl, orbits

def grid_sym_max(n):
    """Exhaust grid-sym tilings (one free bit per orbit). Return (hmax, counter, n_gridsym)."""
    tl, orbits = rho_orbits(n)
    half = len(orbits)
    assert half == half_n(n), f"orbit count {half} != half_n {half_n(n)} at n={n}"
    spec = Counter()
    hmax = 0
    for bits in itertools.product([0, 1], repeat=half):
        t = {}
        for i, orbit in enumerate(orbits):
            for cell in orbit:
                t[cell] = bits[i]
        adj = tiling_to_arcs(t, n)
        H = count_ham_paths(adj, n)
        spec[H] += 1
        if H > hmax:
            hmax = H
    return hmax, spec, 2 ** half

def full_cube_analysis(n):
    """Exhaust ALL 2^m tilings. Return (global_max, n_max, n_max_gridsym, n_max_nongs)."""
    tl = tiles(n)
    m = len(tl)
    assert m == comb(n - 1, 2), f"tile count {m} != C(n-1,2)={comb(n-1,2)}"
    rho_map = {p: rho(*p, n) for p in tl}
    gmax = 0
    Hvals = []
    gridsym_flags = []
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        adj = tiling_to_arcs(t, n)
        H = count_ham_paths(adj, n)
        is_gs = all(t[p] == t[rho_map[p]] for p in tl)
        Hvals.append(H)
        gridsym_flags.append(is_gs)
        if H > gmax:
            gmax = H
    n_max = sum(1 for h in Hvals if h == gmax)
    n_max_gs = sum(1 for h, g in zip(Hvals, gridsym_flags) if h == gmax and g)
    n_max_ngs = n_max - n_max_gs
    return gmax, m, 2 ** m, n_max, n_max_gs, n_max_ngs

# ---------------------------------------------------------------------------
# self-converse cross-check: rho-fixed tiling <=> phi-self-converse tournament
# ---------------------------------------------------------------------------
def selfconv_crosscheck(n):
    tl = tiles(n)
    m = len(tl)
    rho_map = {p: rho(*p, n) for p in tl}
    mism = 0
    cnt_fixed = 0
    cnt_sc = 0
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        is_fixed = all(t[p] == t[rho_map[p]] for p in tl)
        adj = tiling_to_arcs(t, n)
        cr = converse_relabel_arcs(adj, n)
        is_sc = arcs_equal(adj, cr, n)
        if is_fixed:
            cnt_fixed += 1
        if is_sc:
            cnt_sc += 1
        if is_fixed != is_sc:
            mism += 1
    return cnt_fixed, cnt_sc, mism

# ---------------------------------------------------------------------------
def main():
    print("=" * 72)
    print("ADVERSARIAL re-check of HYP-2688 (FRESH script, no reuse)")
    print("=" * 72)

    canon_global = {3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661, 9: 3357}

    # --- sanity: Hamiltonian counter on a known case ---
    # transitive tournament on n vertices has exactly 1 Hamiltonian path.
    for n in range(3, 8):
        tl = tiles(n)
        t = {p: 0 for p in tl}   # all forward -> transitive n->..->1 plus a->b all forward = transitive
        adj = tiling_to_arcs(t, n)
        H = count_ham_paths(adj, n)
        print(f"  [sanity] n={n} all-forward (transitive) H = {H}  (expect 1)  "
              f"{'OK' if H == 1 else 'FAIL'}")

    # --- self-converse cross-check (validates rho == converse_relabel fixed pts) ---
    print("\n[X] rho-fixed tiling  <=>  phi-self-converse tournament  (independent)")
    for n in range(3, 8):
        cf, cs, mism = selfconv_crosscheck(n)
        print(f"  n={n}: #rho-fixed={cf}  #phi-self-converse={cs}  mismatches={mism}  "
              f"2^half={2**half_n(n)}  {'OK' if (cf==cs==2**half_n(n) and mism==0) else 'FAIL'}")

    # --- Q1/Q2/Q4: grid-sym maxima ---
    print("\n[GS] grid-symmetric maxima (exhaust 2^half(n) grid-sym tilings)")
    print(f"{'n':>3} {'half':>5} {'2^half':>8} {'GS_maxH':>8} {'canon_max':>10} {'ATTAINS':>8}")
    gs_results = {}
    for n in [5, 6, 7, 8, 9]:
        hmax, spec, ngs = grid_sym_max(n)
        gs_results[n] = (hmax, spec, ngs)
        cm = canon_global[n]
        print(f"{n:>3} {half_n(n):>5} {ngs:>8} {hmax:>8} {cm:>10} "
              f"{str(hmax == cm):>8}")
        top = dict(sorted(spec.items(), reverse=True)[:5])
        print(f"        top grid-sym H spectrum (H:count) = {top}")

    # --- Q3: full-cube cross-check for n=5,6,7 ---
    print("\n[FC] FULL-CUBE exhaustion (all 2^m tilings) n=5,6,7")
    print(f"{'n':>3} {'m':>3} {'2^m':>7} {'globalH':>8} {'#max':>5} {'#max_GS':>8} "
          f"{'#max_nonGS':>11} {'canon':>6} {'match':>6}")
    for n in [5, 6, 7]:
        gmax, m, twom, nmax, nmaxgs, nmaxngs = full_cube_analysis(n)
        cm = canon_global[n]
        match = (gmax == cm) and (gmax == gs_results[n][0])
        print(f"{n:>3} {m:>3} {twom:>7} {gmax:>8} {nmax:>5} {nmaxgs:>8} "
              f"{nmaxngs:>11} {cm:>6} {str(match):>6}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
