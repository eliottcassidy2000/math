#!/usr/bin/env python3
"""
edge_thickness_s20cs.py — Self-loop + edge thickness decomposition
kind-pasteur-2026-03-23-S20cs (overnight session)

THE KEY FORMULA:
  T_n = SL_n + sum_{edges e} thickness(e)

where:
  T_n = total transition orbits (Burnside-computable)
  SL_n = self-loop count (arc flips that stay in same merged class)
  thickness(e) = number of (tournament, arc) transitions generating edge e

If we can compute SL_n via Burnside, then:
  E(G_n) = number of edges with thickness > 0
  And since thickness -> 2 uniformly, E(G_n) ~ (T_n - SL_n) / 2

This script computes the EXACT decomposition for n=3..6 (and n=7 if feasible).
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations, combinations
from collections import defaultdict, Counter
from functools import reduce
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EDGE THICKNESS DECOMPOSITION: T_n = SL_n + sum(thickness)")
print("  kind-pasteur-2026-03-23-S20cs")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================
def tadj(n, bits):
    a = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): a[i][j] = 1
            else: a[j][i] = 1
            idx += 1
    return a

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def comp(a, n):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def ssq(a, n):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

# ============================================================================
# BURNSIDE T_n COMPUTATION
# ============================================================================

def compute_T_burnside(n):
    """Compute T_n via Burnside: count (tournament, arc-flip) pairs fixed by each sigma."""
    m = comb(n, 2)
    # For each permutation sigma, count Fix(sigma) * m_sigma
    # where m_sigma = number of arc positions fixed by sigma
    # Actually T_n = sum over all tournaments T of (number of distinct classes reachable)
    # This is different from Burnside...

    # Simpler: T_n = (1/n!) * sum_sigma Fix(sigma) * f(sigma)
    # where f(sigma) = number of arc orbits under sigma's action on arcs
    # No wait, that's not right either.

    # T_n (transition orbits) = number of orbits of S_n on (tournament, arc) pairs
    # By Burnside: T_n = (1/n!) * sum_sigma |{(T, (i,j)) : sigma fixes (T, (i,j))}|
    # sigma fixes (T, (i,j)) iff sigma(T) = T AND sigma({i,j}) = {i,j}
    # i.e., sigma is an automorphism of T that fixes the arc {i,j} as a set

    # This equals: sum_sigma (number of T fixed by sigma) * (number of arcs fixed by sigma)

    total = 0
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    # For efficiency, group by cycle type
    from collections import Counter as Ctr
    cycle_type_count = defaultdict(int)
    for sigma in permutations(range(n)):
        # Compute cycle type
        visited = [False]*n; cycles = []
        for i in range(n):
            if visited[i]: continue
            c = 0; j = i
            while not visited[j]:
                visited[j] = True; j = sigma[j]; c += 1
            cycles.append(c)
        ct = tuple(sorted(cycles, reverse=True))
        cycle_type_count[ct] += 1

    for ct, count in cycle_type_count.items():
        # For this cycle type, compute Fix(sigma) and arc_fix(sigma)
        # Use one representative sigma with this cycle type

        # Build representative sigma
        sigma = list(range(n)); pos = 0
        for c in ct:
            # cycle of length c starting at pos
            for i in range(c-1):
                sigma[pos+i] = pos+i+1
            sigma[pos+c-1] = pos
            pos += c

        # Fix(sigma) = number of tournaments fixed by sigma
        # = 2^(number of arc orbits with all-odd cycles)... actually:
        # Arc (i,j) maps to (sigma(i), sigma(j)). Group arcs into orbits.
        arc_orbits = []
        visited_arcs = set()
        for i in range(n):
            for j in range(i+1, n):
                if (i,j) in visited_arcs: continue
                orbit = set()
                a, b = i, j
                while (min(a,b), max(a,b)) not in orbit:
                    orbit.add((min(a,b), max(a,b)))
                    a, b = sigma[a], sigma[b]
                for x in orbit: visited_arcs.add(x)
                arc_orbits.append(orbit)

        # Check if all cycles in sigma are odd (for Fix(sigma) > 0)
        all_odd = all(c % 2 == 1 for c in ct)
        if all_odd:
            # Fix(sigma) = 2^(number of arc orbits that are "free")
            # Actually for tournaments: an arc orbit of size k contributes
            # 1 free bit if the orbit maps arc to arc (not to reverse)
            # For all-odd-cycle permutations: all arc orbits are "consistent"
            # Fix(sigma) = 2^(number of arc orbits)

            # Need to check: does the orbit reverse arcs or not?
            n_free = 0
            for orbit in arc_orbits:
                # Check if orbit is self-reversing
                (a0, b0) = list(orbit)[0]
                # Trace: (a0,b0) -> (sigma(a0), sigma(b0)) -> ...
                # Check if at some point we get (b,a) instead of (a,b)
                a, b = a0, b0
                reverses = False
                for _ in range(len(orbit)):
                    a, b = sigma[a], sigma[b]
                    if (min(a,b), max(a,b)) in orbit:
                        # Check direction
                        if a > b:  # reversed direction
                            reverses = True
                    if (min(a,b), max(a,b)) == (min(a0,b0), max(a0,b0)):
                        break
                # For tournaments: if orbit doesn't reverse, 2 choices (free bit)
                # If orbit reverses, the tournament value is forced (1 choice)
                # Actually for all-odd sigma: no orbit reverses
                n_free += 1

            fix_sigma = 2 ** n_free
        else:
            fix_sigma = 0

        # Number of arcs fixed by sigma (as sets): {i,j} fixed iff sigma({i,j}) = {i,j}
        # i.e., either sigma(i)=i and sigma(j)=j, or sigma(i)=j and sigma(j)=i
        arcs_fixed = 0
        for i in range(n):
            for j in range(i+1, n):
                if (sigma[i] == i and sigma[j] == j) or (sigma[i] == j and sigma[j] == i):
                    arcs_fixed += 1

        total += count * fix_sigma * arcs_fixed

    T_n = total // factorial(n)
    return T_n


# ============================================================================
# EXACT DECOMPOSITION
# ============================================================================

def exact_decomposition(n):
    m = comb(n,2); total = 1<<m; t0 = time.time()
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
    iso = defaultdict(list)
    for b in range(total):
        a = tadj(n, b); iso[canon(a, n)].append(b)

    cl = []; c2c = {}
    for idx, (cn, mem) in enumerate(sorted(iso.items())):
        a = tadj(n, mem[0])
        cl.append({'cid':idx,'canon':cn,'mem':mem,'adj':a,
                   'H':Hdp(a,n),'sc':canon(comp(a,n),n)==cn,
                   'aut':factorial(n)//len(mem),'size':len(mem)})
        c2c[cn] = idx
    for d in cl:
        cc = canon(comp(d['adj'],n),n)
        d['comp_cid'] = c2c.get(cc,-1)

    # Merge map
    mi = {}; mid = 0
    for d in cl:
        ci = d['cid']
        if ci in mi: continue
        cp = d['comp_cid']; mi[ci] = mid
        if cp != ci and cp >= 0: mi[cp] = mid
        mid += 1
    Vm = mid

    mc = {}
    for d in cl:
        mv = mi[d['cid']]
        if mv not in mc:
            mc[mv] = {'H':d['H'], 'sc':d['sc'], 'aut':d['aut'], 'size':d['size']}

    # Count ALL transitions: for each tournament, flip each arc, classify result
    # Track: self-loops, edge thicknesses (per merged edge)
    self_loops = 0  # total (tournament, arc) pairs that stay in same merged class
    edge_thickness = Counter()  # (min_merged, max_merged) -> count of transitions
    edge_color = {}  # edge -> 'blue' or 'black'

    # Also track per-class self-loop count
    class_self_loops = Counter()

    for d in cl:
        ci = d['cid']
        src = mi[ci]
        for b in d['mem']:
            for k in range(m):
                fb = b ^ (1 << k)
                fa = tadj(n, fb)
                fc = canon(fa, n)
                nb_cid = c2c.get(fc)
                if nb_cid is None: continue
                tgt = mi[nb_cid]

                if tgt == src:
                    self_loops += 1
                    class_self_loops[src] += 1
                else:
                    e = (min(src, tgt), max(src, tgt))
                    edge_thickness[e] += 1
                    if e not in edge_color:
                        edge_color[e] = 'blue' if mc[src]['sc'] == mc[tgt]['sc'] else 'black'

    # Total transitions
    T_total = self_loops + sum(edge_thickness.values())

    # T_n should equal T_total / n! ... wait, T_total counts raw (tournament, arc) pairs
    # T_n (transition orbits) = T_total / n! ... no, that's not right either
    # Each tournament appears |orbit| = n!/|Aut(T)| times
    # And each has m arc flips
    # So T_total = sum over classes of |orbit| * m = 2^m * m
    # T_n (orbits) = T_total / n! ... no.
    # Actually T_n as defined by Burnside = number of (iso_class, arc_position) orbits
    # under diagonal S_n action.

    # Let me just report the raw counts
    print(f"\n  n={n}: V_merged={Vm}, {len(cl)} iso classes [{time.time()-t0:.1f}s]")
    print(f"  Total (tournament, arc) transitions: {total * m}")
    print(f"    = {total} tournaments x {m} arcs")
    print(f"  Self-loops (stay in same merged class): {self_loops}")
    print(f"  Cross-transitions (change merged class): {sum(edge_thickness.values())}")
    print(f"  Edges (distinct pairs): {len(edge_thickness)}")

    # Thickness distribution
    thick_dist = Counter(edge_thickness.values())
    print(f"\n  EDGE THICKNESS DISTRIBUTION:")
    for t in sorted(thick_dist.keys()):
        E_at_t = thick_dist[t]
        is_blue = sum(1 for e, th in edge_thickness.items() if th == t and edge_color.get(e) == 'blue')
        is_black = E_at_t - is_blue
        print(f"    thickness={t}: {E_at_t} edges ({is_blue} blue, {is_black} black)")

    # Verify: self_loops + sum(thicknesses) = total transitions
    total_cross = sum(edge_thickness.values())
    print(f"\n  VERIFICATION:")
    print(f"    self_loops + sum(thickness) = {self_loops} + {total_cross} = {self_loops + total_cross}")
    print(f"    total (T,arc) pairs = {total * m}")
    match = (self_loops + total_cross == total * m)
    print(f"    Match? {'YES' if match else 'NO'}")

    # Per-class analysis
    print(f"\n  PER-CLASS SELF-LOOP ANALYSIS:")
    print(f"    {'mid':>4} {'H':>4} {'SC':>3} {'|Aut|':>5} {'size':>6} {'SL':>6} {'SL/size':>8} {'SL/(size*m)':>12}")
    for mv in sorted(mc.keys()):
        d = mc[mv]
        sl = class_self_loops.get(mv, 0)
        sl_per = sl / d['size'] if d['size'] > 0 else 0
        sl_frac = sl / (d['size'] * m) if d['size'] > 0 else 0
        if Vm <= 15 or d['sc']:  # show all at small n, just SC at larger n
            print(f"    {mv:4d} {d['H']:4d} {'Y' if d['sc'] else 'N':>3} {d['aut']:5d} "
                  f"{d['size']:6d} {sl:6d} {sl_per:8.1f} {sl_frac:12.4f}")

    # Key ratios
    E = len(edge_thickness)
    avg_thickness = sum(edge_thickness.values()) / E if E > 0 else 0
    print(f"\n  KEY RATIOS:")
    print(f"    E (edges) = {E}")
    print(f"    SL (self-loops) = {self_loops}")
    print(f"    Total transitions = {total * m}")
    print(f"    Avg thickness = {avg_thickness:.4f}")
    print(f"    SL / total = {self_loops / (total * m):.4f}")
    print(f"    (total - SL) / (2*E) = {(total*m - self_loops) / (2*E):.4f}" if E > 0 else "")

    # Blue vs black thickness
    blue_thick = [edge_thickness[e] for e in edge_thickness if edge_color.get(e) == 'blue']
    black_thick = [edge_thickness[e] for e in edge_thickness if edge_color.get(e) == 'black']
    if blue_thick:
        print(f"\n    Blue edges: {len(blue_thick)}, avg thickness = {np.mean(blue_thick):.2f}, "
              f"min={min(blue_thick)}, max={max(blue_thick)}")
    if black_thick:
        print(f"    Black edges: {len(black_thick)}, avg thickness = {np.mean(black_thick):.2f}, "
              f"min={min(black_thick)}, max={max(black_thick)}")

    # Self-loop fraction by SC/NS type
    sc_sl = sum(class_self_loops.get(mv, 0) for mv in mc if mc[mv]['sc'])
    ns_sl = sum(class_self_loops.get(mv, 0) for mv in mc if not mc[mv]['sc'])
    sc_total = sum(mc[mv]['size'] * m for mv in mc if mc[mv]['sc'])
    ns_total = sum(mc[mv]['size'] * m for mv in mc if not mc[mv]['sc'])
    print(f"\n    SC self-loop fraction: {sc_sl}/{sc_total} = {sc_sl/sc_total:.4f}" if sc_total > 0 else "")
    print(f"    NS self-loop fraction: {ns_sl}/{ns_total} = {ns_sl/ns_total:.4f}" if ns_total > 0 else "")

    return {
        'n': n, 'V': Vm, 'E': E, 'SL': self_loops,
        'total_trans': total * m, 'avg_thick': avg_thickness,
        'thick_dist': thick_dist, 'blue_thick': blue_thick, 'black_thick': black_thick,
        'class_sl': dict(class_self_loops), 'mc': mc
    }


# ============================================================================
# MAIN
# ============================================================================

all_results = {}

for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")
    all_results[n] = exact_decomposition(n)


# ============================================================================
# CROSS-n SYNTHESIS
# ============================================================================

print(f"\n\n{'='*80}")
print("  CROSS-n SYNTHESIS")
print(f"{'='*80}")

print(f"\n  {'n':>3} {'V_m':>6} {'E':>6} {'SL':>8} {'total':>10} {'SL/tot':>8} "
      f"{'avg_th':>8} {'(tot-SL)/2E':>12}")
print(f"  {'-'*70}")
for n in [3,4,5,6]:
    r = all_results[n]
    ratio = (r['total_trans'] - r['SL']) / (2*r['E']) if r['E'] > 0 else 0
    print(f"  {n:3d} {r['V']:6d} {r['E']:6d} {r['SL']:8d} {r['total_trans']:10d} "
          f"{r['SL']/r['total_trans']:8.4f} {r['avg_thick']:8.2f} {ratio:12.4f}")

# Self-loop sequence
print(f"\n  Self-loop sequence SL_n: {', '.join(str(all_results[n]['SL']) for n in [3,4,5,6])}")
print(f"  Edge sequence E_n: {', '.join(str(all_results[n]['E']) for n in [3,4,5,6])}")

# The formula: E = (total - SL) / avg_thickness
# Can we predict SL from Burnside?
sl_ratios = []
for n in [3,4,5,6]:
    m = comb(n,2)
    sl_ratios.append(f"{all_results[n]['SL']/(2**m * m):.6f}")
print(f"\n  SL_n / (2^m * m): {', '.join(sl_ratios)}")

# SL per tournament (averaged)
sl_per = []
for n in [3,4,5,6]:
    sl_per.append(f"{all_results[n]['SL']/2**comb(n,2):.2f}")
print(f"  SL per tournament: {', '.join(sl_per)}")
print(f"  (should be ~neutral arcs = n-1 for transitive, varies for others)")

# Thickness distribution summary
print(f"\n  THICKNESS DISTRIBUTION SUMMARY:")
for n in [3,4,5,6]:
    td = all_results[n]['thick_dist']
    print(f"    n={n}: {dict(sorted(td.items()))}")

# Key question: is avg_thickness EXACTLY 2 * (something)?
at2 = [f"{all_results[n]['avg_thick']/2:.4f}" for n in [3,4,5,6]]
print(f"\n  avg_thickness / 2: {', '.join(at2)}")

# The magic formula check: E = (total - SL) / avg_thickness
print(f"\n  CHECK: E = (total - SL) / avg_thickness ?")
for n in [3,4,5,6]:
    r = all_results[n]
    predicted_E = (r['total_trans'] - r['SL']) / r['avg_thick'] if r['avg_thick'] > 0 else 0
    print(f"    n={n}: predicted E = {predicted_E:.1f}, actual E = {r['E']}  "
          f"{'EXACT' if abs(predicted_E - r['E']) < 0.01 else 'APPROX'}")


print(f"\n  DONE.")
print("=" * 80)
