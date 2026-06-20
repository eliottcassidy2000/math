#!/usr/bin/env python3
"""
half_tiling_gnomon_kpswf.py
================================================================
THREAD 3 -- HYP-2687: gnomon shells and the even-odd dichotomy.
kind-pasteur (wf) 2026-06-20.

Builds on 04-computation/half_tiling_framework_kps.py (REUSED helpers).

The rho-involution rho(a,b) = (n+1-b, n+1-a) acts on the staircase delta_{n-2}.
Its orbits are: d = floor((n-1)/2) FIXED cells (on anti-diagonal a+b = n+1) and
(m-d)/2 transposed PAIRS.  #orbits = half(n) = floor((n-1)^2/4).

For ODD n = 2k+1:   half = k^2 = 1+3+5+...+(2k-1)   (k odd GNOMON shells)
For EVEN n = 2k:    half = k(k-1)                     (pronic, NO central fixed cell)

This script:
 (a) Identifies the k gnomon shells explicitly as sets of rho-orbits, in (a,b)
     and pin-grid (r,c) coords; characterizes span/strip/distance-from-apex.
 (b) Lists the fixed diagonal cells {a+b=n+1} per n; identifies the tournament
     feature each controls (apex (n,1) = source->sink arc, etc.).
 (c) EVEN-ODD: computes c3 (# directed 3-cycles) parity across ALL tilings for
     n=4..7, split by grid-sym vs non-grid-sym; relates to THM-016/THM-017.
 (d) Mode-B shell test: does the outer gnomon (n -> n+2) correspond to one nested
     SC sub-tournament shell on the inner n-2 vertices?
"""

import itertools
from math import comb, floor
from collections import Counter, defaultdict

# ----------------------------------------------------------------------
# Reused machinery (identical conventions to half_tiling_framework_kps.py)
# ----------------------------------------------------------------------

def tiles(n):
    """Tiles (a,b) with a>=b+2, explorer order: for y=1..n-2: for x=n..y+2: (x,y)."""
    out = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            out.append((x, y))
    return out

def rc(a, b):
    """pin-grid coordinates: r = a-b-1 (the 'skip'/distance-1 strip index), c = b."""
    return (a-b-1, b)

def rho(a, b, n):
    """rho(a,b) = (n+1-b, n+1-a). Upper vertex first since b<a => n+1-b>n+1-a."""
    return (n+1-b, n+1-a)

def tiling_to_tournament(t, tile_list, n):
    T = {}
    for k in range(n, 1, -1):
        T[(k, k-1)] = 1
        T[(k-1, k)] = 0
    for (a, b) in tile_list:
        bit = t[(a, b)]
        if bit == 0:
            T[(a, b)] = 1; T[(b, a)] = 0
        else:
            T[(a, b)] = 0; T[(b, a)] = 1
    return T

def count_c3(T, n):
    """Number of directed 3-cycles (cyclic triangles) in tournament T on {1..n}."""
    c = 0
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            for k in range(j+1, n+1):
                # a directed 3-cycle on {i,j,k}: either i->j->k->i or i->k->j->i
                if T[(i, j)] and T[(j, k)] and T[(k, i)]:
                    c += 1
                elif T[(i, k)] and T[(k, j)] and T[(j, i)]:
                    c += 1
    return c

def involution_structure(n):
    tl = tiles(n)
    tset = set(tl)
    m = len(tl)
    fixed, pairs, seen = [], [], set()
    for (a, b) in tl:
        img = rho(a, b, n)
        assert img in tset, f"rho image {img} of {(a,b)} not a tile (n={n})"
        if img == (a, b):
            fixed.append((a, b))
        else:
            key = frozenset([(a, b), img])
            if key not in seen:
                seen.add(key)
                pairs.append(((a, b), img))
    d = len(fixed)
    half = d + len(pairs)
    return tl, m, d, fixed, pairs, half

def half_formula(n):
    return floor((n-1)**2 / 4) if n >= 1 else 0

# ----------------------------------------------------------------------
# (a) GNOMON SHELLS for odd n
# ----------------------------------------------------------------------
#
# The half-tiling fundamental domain is one "triangle half" of the staircase:
# the rho-orbit reps.  We choose canonical reps = the orbit element on the kept
# side a+b <= n+1 (anti-diagonal fixed cells included once).
#
# For odd n=2k+1 the domain has k^2 orbits.  We index orbits by their pin-grid
# (r,c) of the canonical rep.  The natural Cartesian gnomon coordinate is
#    (R, C) where we re-base the half-domain to a k x k square.
#
# Geometric realization of the gnomon: shell i (i=1..k) = the L-shaped layer
# {orbits with max(R,C) = i}.  Shell i has 2i-1 orbits.
#
# We MAP the half-domain to the k x k square as follows, then read off shells.

def canonical_orbit_reps(n):
    """Return list of (rep_tile, orbit) where orbit is a frozenset of 1 or 2 tiles.
       Canonical rep = the tile in the orbit with a+b <= n+1 (fixed cells: a+b=n+1)."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    reps = []
    for (a, b) in fixed:
        reps.append(((a, b), frozenset([(a, b)]), 'fixed'))
    for (p, q) in pairs:
        # pick the one with smaller a+b as the kept rep
        sp, sq = p[0]+p[1], q[0]+q[1]
        rep = p if sp <= sq else q
        reps.append((rep, frozenset([p, q]), 'pair'))
    return reps, fixed, pairs

def square_coords(n):
    """Map each rho-orbit of odd n=2k+1 into a k x k integer grid (R,C), 1<=R,C<=k.
       Construction: the kept half is {(a,b): a+b<=n+1}.  In pin coords (r,c) with
       r=a-b-1>=1, c=b>=1, r+c<=n-1.  The kept-half condition a+b<=n+1 is r+2c<=n-1.
       We build an explicit bijection to {1..k}^2 by a diagonal sweep so gnomon
       shells = L-layers."""
    assert n % 2 == 1
    k = (n-1)//2
    reps, fixed, pairs = canonical_orbit_reps(n)
    # Sort orbits by (a+b, then b) to get a stable canonical ordering.
    reps_sorted = sorted(reps, key=lambda x: (x[0][0]+x[0][1], x[0][1], x[0][0]))
    # We assign square coords by the *gnomon-distance* of each orbit.
    # gnomon distance g(orbit) = number of vertices "between the legs", made precise
    # below in shells_by_apex_distance; here we just return reps_sorted + k.
    return reps_sorted, k, fixed, pairs

def shells_by_apex_distance(n):
    """Partition the rho-orbits into k gnomon shells for odd n=2k+1.

    DEFINITION (apex-nested gnomons):  The apex tile is (n,1) (pin (n-2,1),
    the source->sink arc, on the fixed diagonal).  We define for each orbit a
    'shell index' s = the Mode-B nesting level it first appears at.

    Concretely: vertex set {1,...,n}.  Strip away the OUTER pair of vertices
    {1, n} (Mode-B step n -> n-2).  A tile (a,b) 'belongs to the inner n-2
    sub-tournament on {2,...,n-1}' iff 2 <= b and a <= n-1, i.e. it survives the
    peel.  Shell index of an orbit = how many peels until BOTH tiles of the orbit
    are gone = ceil based on min over orbit of the peel-survival count.

    Peel t (t=0,1,2,...) removes vertices {1,...,t} and {n-t+1,...,n}.
    Tile (a,b) survives peel t iff b > t and a <= n-t, i.e. t < b and t <= n-a,
    so t < min(b, n-a+1) => survives for t = 0,...,min(b, n+1-a)-1.
    survival(a,b) = min(b, n+1-a).  An orbit {(a,b),(n+1-b,n+1-a)} has the SAME
    survival for both tiles (rho preserves min(b,n+1-a)!) -- verified below.
    Shell index = (max survival over half-domain) - survival + 1, so the apex
    (survival=1 ... actually apex (n,1): min(1, n+1-n)=min(1,1)=1) is the OUTERMOST
    shell (index k), and the deepest-nested orbit is shell 1 (innermost).

    Returns: dict shell_index -> list of (rep_tile, orbit, survival).
    """
    reps, fixed, pairs = canonical_orbit_reps(n)
    def survival(a, b):
        return min(b, n+1-a)
    # verify rho-invariance of survival on every orbit
    for rep, orbit, kind in reps:
        svals = set()
        for (a, b) in orbit:
            svals.add(survival(a, b))
        assert len(svals) == 1, f"survival not rho-invariant on orbit {orbit} (n={n}): {svals}"
    smax = max(survival(*rep) for rep, _, _ in reps)
    shells = defaultdict(list)
    for rep, orbit, kind in reps:
        s = survival(*rep)
        shell_idx = smax - s + 1   # innermost (largest survival) = shell 1
        shells[shell_idx].append((rep, orbit, s, kind))
    return dict(shells), smax

# ----------------------------------------------------------------------
# (b) fixed diagonal cells and tournament features
# ----------------------------------------------------------------------

def fixed_cells_features(n):
    """List the fixed cells (a+b=n+1) with pin coords and the tournament feature
       each controls."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rows = []
    for (a, b) in sorted(fixed, key=lambda ab: -(ab[0]-ab[1])):
        r, c = rc(a, b)
        span = a - b
        # The arc spans vertices a..b; on the anti-diagonal a+b=n+1.
        # apex (n,1): the single arc between the global source-side (n) and sink (1).
        if (a, b) == (n, 1):
            feat = "APEX: arc between extreme vertices n and 1 (source<->sink); span n-1"
        else:
            feat = f"mid-diagonal arc span={span}; symmetric pivot at vertices ({a},{b})"
        rows.append((a, b, r, c, span, feat))
    return rows, d

# ----------------------------------------------------------------------
# (c) c3 parity across tilings
# ----------------------------------------------------------------------

def c3_parity_analysis(n):
    """Exhaustive over all 2^m tilings: distribution of c3, parity of c3, split by
       grid-sym vs non-grid-sym."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rho_map = {p: rho(*p, n) for p in tl}
    c3_dist = Counter()
    parity_all = Counter()           # 0/1 -> count
    parity_gs = Counter()
    parity_ngs = Counter()
    c3_gs = Counter()
    c3_ngs = Counter()
    total = 0
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        c = count_c3(T, n)
        c3_dist[c] += 1
        par = c & 1
        parity_all[par] += 1
        is_gs = all(t[p] == t[rho_map[p]] for p in tl)
        if is_gs:
            parity_gs[par] += 1
            c3_gs[c] += 1
        else:
            parity_ngs[par] += 1
            c3_ngs[c] += 1
        total += 1
    return {
        'n': n, 'm': m, 'total': total,
        'c3_dist': dict(sorted(c3_dist.items())),
        'parity_all': dict(parity_all),
        'parity_gs': dict(parity_gs),
        'parity_ngs': dict(parity_ngs),
        'c3_gs': dict(sorted(c3_gs.items())),
        'c3_ngs': dict(sorted(c3_ngs.items())),
        'both_parities_all': (parity_all.get(0,0)>0 and parity_all.get(1,0)>0),
        'both_parities_gs': (parity_gs.get(0,0)>0 and parity_gs.get(1,0)>0),
    }

# ----------------------------------------------------------------------
# (d) Mode-B nested SC shell test
# ----------------------------------------------------------------------

def converse_relabel(T, n):
    R = {}
    for (u, v), val in T.items():
        pu, pv = n+1-u, n+1-v
        if val == 1:
            R[(pv, pu)] = 1
            R[(pu, pv)] = 0
    return R

def is_phi_self_converse(t, tl, n):
    rho_map = {p: rho(*p, n) for p in tl}
    return all(t[p] == t[rho_map[p]] for p in tl)

def mode_b_shell_test(n):
    """Take a grid-sym (phi-self-converse) tiling on n vertices.  Peel the outer
       vertices {1,n} (Mode-B).  Question: is the restriction to {2,...,n-1} a
       phi'-self-converse tournament on those n-2 vertices, with phi'(i)=n+1-i
       (the SAME involution, now mapping {2..n-1} onto itself)?

       i.e. does each outer gnomon shell add exactly one nested SC layer, leaving
       the inner n-2 vertices carrying an SC sub-tournament?

       We test: for EVERY phi-self-converse tiling on n vertices, the inner
       sub-tournament on {2..n-1} is phi'-self-converse on those n-2 labels."""
    tl, m, d, fixed, pairs, half = involution_structure(n)
    rho_map = {p: rho(*p, n) for p in tl}
    inner_tiles = [(a, b) for (a, b) in tl if b >= 2 and a <= n-1]  # tiles of {2..n-1}
    # phi' on inner labels {2..n-1} is still i -> n+1-i (maps the set to itself).
    # On a tile (a,b) of the inner sub-tournament, the inner-rho is (n+1-b, n+1-a)
    # -- the SAME formula -- because phi'(i)=n+1-i on {2..n-1}.
    ok = 0
    tested = 0
    fails = []
    for bits in itertools.product([0, 1], repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        if not is_phi_self_converse(t, tl, n):
            continue
        tested += 1
        # inner is automatically phi'-self-converse because the inner tiles form
        # rho-orbits among themselves (rho maps inner tiles to inner tiles).
        inner_ok = all(t[(a, b)] == t[rho_map[(a, b)]] for (a, b) in inner_tiles)
        # also check that rho maps inner tiles to inner tiles (closure):
        closure = all(rho_map[(a, b)] in set(inner_tiles) for (a, b) in inner_tiles)
        if inner_ok and closure:
            ok += 1
        else:
            fails.append((bits, inner_ok, closure))
    return {'n': n, 'tested_gridsym': tested, 'inner_selfconv': ok,
            'n_inner': n-2, 'n_inner_tiles': len(inner_tiles),
            'all_ok': (ok == tested), 'num_fail': len(fails)}

def mode_b_shell_count(n):
    """Count the orbits added when going n-2 -> n (one Mode-B step = outer gnomon).
       For odd n=2k+1: outer shell = shell index k has 2k-1 orbits.
       Verify half(n) - half(n-2) = (outer gnomon width)."""
    return half_formula(n) - half_formula(n-2)

# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    print("="*74)
    print("HYP-2687  GNOMON SHELLS & EVEN-ODD DICHOTOMY  -- kind-pasteur(wf) 2026-06-20")
    print("="*74)

    # ---------- size sanity ----------
    print("\n[0] half(n) = floor((n-1)^2/4):  odd k^2 (gnomon square) / even k(k-1) (pronic)")
    print(f"{'n':>3} {'m':>4} {'d=fix':>6} {'pairs':>6} {'half':>5} {'shape':>22} {'half(n)-half(n-2)':>18}")
    for n in range(3, 14):
        tl, m, d, fixed, pairs, half = involution_structure(n)
        if n % 2 == 1:
            k = (n-1)//2; shape = f"odd k={k}: k^2={k*k}"
        else:
            k = n//2; shape = f"even k={k}: k(k-1)={k*(k-1)}"
        modeb = half - half_formula(n-2) if n >= 5 else half
        print(f"{n:>3} {m:>4} {d:>6} {len(pairs):>6} {half:>5} {shape:>22} {modeb:>18}")

    # ---------- (a) gnomon shells (odd n) ----------
    print("\n" + "="*74)
    print("[a] GNOMON SHELLS for ODD n  (shell i = nesting level; widths 1,3,5,...,2k-1)")
    print("    survival(a,b)=min(b, n+1-a) is rho-invariant; shell = smax - survival + 1")
    print("    shell 1 = innermost (deepest Mode-B core), shell k = outermost (apex layer)")
    for n in [3, 5, 7, 9]:
        shells, smax = shells_by_apex_distance(n)
        k = (n-1)//2
        print(f"\n  --- n={n} (k={k}, half={k*k}, smax={smax}) ---")
        widths = []
        for s_idx in sorted(shells):
            orbs = shells[s_idx]
            widths.append(len(orbs))
            descr = []
            for rep, orbit, surv, kind in sorted(orbs, key=lambda z: z[0]):
                r, c = rc(*rep)
                tag = 'FIX' if kind == 'fixed' else 'pair'
                if kind == 'pair':
                    other = [x for x in orbit if x != rep][0]
                    descr.append(f"{rep}~{other}[{tag},rc={r,c},surv={surv}]")
                else:
                    descr.append(f"{rep}[{tag},rc={r,c},surv={surv}]")
            print(f"    shell {s_idx} (width {len(orbs)}): " + "  ".join(descr))
        gnomon = [2*i-1 for i in range(1, k+1)]
        print(f"    -> shell widths {widths}  vs gnomon {gnomon}  "
              f"{'OK (1,3,5,... gnomon)' if widths==gnomon else 'MISMATCH'}")

    # ---------- (b) fixed diagonal cells ----------
    print("\n" + "="*74)
    print("[b] FIXED DIAGONAL CELLS {a+b=n+1}  (count d=floor((n-1)/2); apex=(n,1))")
    for n in range(3, 11):
        rows, d = fixed_cells_features(n)
        print(f"\n  --- n={n}: d={d} fixed cells (apex (n,1) always present) ---")
        for (a, b, r, c, span, feat) in rows:
            print(f"    (a,b)=({a},{b})  pin(r,c)=({r},{c})  span={span}  | {feat}")

    # ---------- (c) c3 parity ----------
    print("\n" + "="*74)
    print("[c] c3 = #directed 3-cycles: PARITY across ALL tilings, n=4..7")
    print("    (THM-016/017 'even-odd split' is the ALTERNATING-SUM identity; here we")
    print("     test whether VERTEX-count parity of n controls c3-PARITY spread.)")
    results_c = {}
    for n in range(3, 8):
        res = c3_parity_analysis(n)
        results_c[n] = res
        print(f"\n  --- n={n} (m={res['m']}, total tilings={res['total']}) ---")
        print(f"    c3 distribution        : {res['c3_dist']}")
        print(f"    parity (all)   even:odd : {res['parity_all']}  "
              f"both parities present={res['both_parities_all']}")
        print(f"    parity (grid-sym)       : {res['parity_gs']}  "
              f"both parities present={res['both_parities_gs']}")
        print(f"    parity (non-grid-sym)   : {res['parity_ngs']}")
        print(f"    c3 values grid-sym only : {sorted(res['c3_gs'].keys())}")
    print("\n  EVEN-ODD CONJECTURE TEST (c (a) variant): does odd n force c3 to take")
    print("  BOTH parities while even n constrains it to one?")
    for n in range(3, 8):
        r = results_c[n]
        par = 'EVEN n' if n % 2 == 0 else 'ODD  n'
        c3vals = list(r['c3_dist'].keys())
        c3par = set(v & 1 for v in c3vals)
        print(f"    n={n} [{par}]: c3 values={c3vals}  c3-parities present={sorted(c3par)}  "
              f"grid-sym c3-parities={sorted(set(v&1 for v in r['c3_gs']))}")

    # ---------- (d) Mode-B nested SC shell ----------
    print("\n" + "="*74)
    print("[d] MODE-B SHELL TEST: outer gnomon (n -> n+2) = one nested SC sub-tournament?")
    print("    For every phi-self-converse tiling on n vertices, is the inner")
    print("    sub-tournament on {2..n-1} phi'-self-converse on those n-2 labels?")
    for n in range(4, 9):
        res = mode_b_shell_test(n)
        widthadd = mode_b_shell_count(n)
        print(f"  n={n}: grid-sym tilings tested={res['tested_gridsym']:>4}  "
              f"inner({res['n_inner']})-self-converse={res['inner_selfconv']:>4}  "
              f"all_ok={res['all_ok']}  | half(n)-half(n-2)={widthadd} "
              f"({'=2k-1 outer gnomon' if n%2==1 else '=2(k-1) outer pronic ring'})")
    print("\n  Interpretation: the inner tiles {(a,b): 2<=b, a<=n-1} are closed under rho")
    print("  (rho maps {2..n-1}-tiles to {2..n-1}-tiles), so an SC tiling restricts to an")
    print("  SC tiling on n-2 vertices. The outer gnomon = the rho-orbits touching vertex")
    print("  1 or n. This is exactly the Mode-B descent G_n -> G_{n-2}.")

    # ---------- (c-proof) the c3-parity mechanism ----------
    print("\n" + "="*74)
    print("[c-PROOF] MECHANISM for the grid-sym c3-parity dichotomy")
    print("  phi(i)=n+1-i is an ANTI-automorphism of any phi-self-converse (grid-sym) T.")
    print("  phi maps each directed 3-cycle to a directed 3-cycle. A directed 3-cycle on")
    print("  vertex-set S is self-paired iff S is phi-INVARIANT; otherwise cycles pair off.")
    print("  => c3 parity = parity of #(directed 3-cycles on phi-invariant 3-sets).")
    print("  phi has a fixed vertex (n+1)/2 only at ODD n; a 3-set (odd size) is")
    print("  phi-invariant iff it contains that fixed vertex + one phi-pair.")
    print(f"  {'n':>3} {'parity(n)':>9} {'d=fixcells':>11} {'#phi-inv 3-sets':>16} {'c3 forced even?':>16}")
    for n in range(3, 14):
        d = floor((n-1)/2)
        def phi(i): return n+1-i
        inv3 = sum(1 for S in itertools.combinations(range(1, n+1), 3)
                   if set(phi(x) for x in S) == set(S))
        par = 'EVEN' if n % 2 == 0 else 'ODD'
        forced = (inv3 == 0)
        print(f"  {n:>3} {par:>9} {d:>11} {inv3:>16} {str(forced):>16}")
    print("  EVEN n: 0 phi-invariant 3-sets (no fixed vertex) => grid-sym c3 ALWAYS EVEN.")
    print("  ODD  n: exactly d=(n-1)/2 phi-invariant 3-sets {(n+1)/2, x, n+1-x} => c3 can be ODD.")
    print("  This is the geometric source of the even-odd split, localized at the fixed")
    print("  diagonal: the EXTRA fixed cell gained per Mode-B step (n->n+2) at odd n is what")
    print("  furnishes a phi-invariant 3-set and unlocks odd c3-parity.")

    print("\nDONE.")

if __name__ == "__main__":
    main()
