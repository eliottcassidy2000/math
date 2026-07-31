#!/usr/bin/env python3
"""Lane E independent audit of MISTAKE-329 / MSG-2937 (mixed Lorenz double-2 repair).

Independent reconstruction, NOT importing repo referee logic for the core
computation (a separate repo cross-check import is done at the end).

Verifies for the two claimed missed rows:
  row1: F=(1,4,5,7,9,11), L=D=194040, S=55392, ds=(2,2,194040), 27696>13860
  row2: F=(1,5,7,8,9,11), L=D=388080, S=109044, ds=(2,2,388080), 54522>27720
  1. L = lcm(14 v) and |S_D| at D=L, from the raw danger-comb definition.
  2. Reflection x -> L-1-x invariance of the safe set and exact parity balance.
  3. fibre_cap(D,d3,2) = (D//lcm(d3,2)) * ceil(ceil(d3/7)/gcd(d3,2)).
  4. Exact parity kill fires; coarse |S|>2*ambient_capacity does not.
  5. The rows pass common-u and the full phase-free Lorenz screen for
     shape (2,2,D) over every proper nontrivial divisor quotient q|D
     (so they genuinely reach the double-2 branch as survivors).
Then a global superset scan over all 3003 bodies (arc arithmetic, no bitsets)
finds every support-hard (row, shape=(2,2,d3)) occurrence on which the coarse
and exact double-2 tests disagree, to confirm the flip set is exactly these 2.
"""

from math import gcd, lcm
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import sys

CUTOFF = Q(26, 31)          # SUPPORT_CUTOFF = (36/91)/(558/1183) = 26/31


def build_safe_bool(F):
    """Direct definition: danger = union of teeth around k*L/v, width L/(14v)."""
    L = lcm(*(14 * v for v in F))
    danger = bytearray(L)
    for v in F:
        half = L // (14 * v)
        period = L // v
        for k in range(v + 1):
            c = k * period
            lo, hi = max(0, c - half), min(L, c + half)
            for x in range(lo, hi):
                danger[x] = 1
    safe = [x for x in range(L) if not danger[x]]
    return L, safe


def ceil_div(a, b):
    return -(-a // b)


def fibre_cap(D, d, q):
    """Sharp q-fibre cap of the enlarged denominator-d needle: max single
    q-fiber load = H*ceil(ell/g) with g=gcd(d,q), H=D//lcm(d,q), ell=ceil(d/7)."""
    ell = ceil_div(d, 7)
    g = gcd(d, q)
    H = D // lcm(d, q)
    return H * ceil_div(ell, g)


def capacity(D, d):
    return (D // d) * ceil_div(d, 7)


def needle_lorenz_parameters(D, d, q):
    ell = ceil_div(d, 7)
    g = gcd(d, q)
    H = D // lcm(d, q)
    A, r = divmod(ell, g)
    high = r * (q // g)
    assert H * (A * q + high) == (D // d) * ell
    return (H * A, H, high)


def loads_mod(points, q):
    """Load of each residue class mod q, from an explicit point list."""
    loads = [0] * q
    for x in points:
        loads[x % q] += 1
    return loads


def lorenz_value(sorted_desc, s):
    return sum(sorted_desc[:s])


def full_lorenz_check(D, points, ds):
    """Return list of (q, s, left, right) positive violations for shape ds."""
    violations = []
    divs = [q for q in range(2, D) if D % q == 0]
    for q in divs:
        loads = sorted(loads_mod(points, q), reverse=True)
        profiles = [needle_lorenz_parameters(D, d, q) for d in ds]
        base = sum(p[0] for p in profiles)
        # test every s (exhaustive; stronger than breakpoint set)
        cum = 0
        for s in range(1, q + 1):
            cum += loads[s - 1]
            right = base * s + sum(h * min(s, hc) for _b, h, hc in profiles)
            if cum > right:
                violations.append((q, s, cum, right))
                break  # one violation per q suffices
    return violations


def check_row(F, D_exp, S_exp, d3, cap_exp):
    L, safe = build_safe_bool(F)
    assert L == D_exp, (L, D_exp)
    S = len(safe)
    assert S == S_exp, (S, S_exp)
    # projection at D=L is identity
    # reflection invariance and parity balance
    safeset = set(safe)
    assert all((L - 1 - x) in safeset for x in safe), "not reflection invariant"
    even = sum(1 for x in safe if x % 2 == 0)
    odd = S - even
    assert even == odd == S // 2, (even, odd)
    # caps
    pcap = fibre_cap(D_exp, d3, 2)
    assert pcap == cap_exp, (pcap, cap_exp)
    amb = capacity(D_exp, d3)
    exact_kill = (S // 2) > pcap
    coarse_kill = S > 2 * amb
    # support-hard?
    hard = Q(S, D_exp) <= CUTOFF
    # common-u for shape (2,2,d3)
    caps = [capacity(D_exp, 2), capacity(D_exp, 2), capacity(D_exp, d3)]
    commonu = any(
        d in (2, 3) and S > (sum(caps) - caps[i])
        for i, d in enumerate((2, 2, d3))
    )
    # baseline all-top test: top ceil(d/7) class loads per mask
    def top_load(d):
        ell = ceil_div(d, 7)
        loads = sorted(loads_mod(safe, d), reverse=True)
        return sum(loads[:ell])
    baseline = S <= top_load(2) + top_load(2) + top_load(d3)
    # full Lorenz screen for shape (2,2,d3)
    viol = full_lorenz_check(D_exp, safe, (2, 2, d3))
    print(f"F={F} L=D={L} S={S} d3={d3}")
    print(f"  parity balance: even={even} odd={odd} (S/2={S//2})")
    print(f"  reflection-invariant: True")
    print(f"  ambient capacity(D,d3)={amb}  q=2 fibre cap={pcap}")
    print(f"  support-hard (S/D<=26/31): {hard}  baseline all-top: {baseline}")
    print(f"  common-u fires: {commonu}")
    print(f"  Lorenz violations for (2,2,{d3}): {viol[:3]}{' NONE -> survivor' if not viol else ''}")
    print(f"  EXACT parity test: S/2={S//2} > {pcap} -> kill={exact_kill}")
    print(f"  COARSE test: S={S} > 2*amb={2*amb} -> kill={coarse_kill}")
    ok = (not viol) and (not commonu) and hard and baseline and exact_kill and not coarse_kill
    print(f"  row verdict: {'CONFIRMED missed kill' if ok else 'MISMATCH'}")
    return ok


def merged_arc_length_mod(ranges, D):
    """|image mod D| of half-open integer ranges, by cyclic arc merge (no bitset)."""
    pieces = []
    for left, right in ranges:
        length = right - left
        if length >= D:
            return D
        r = left % D
        e = r + length
        if e <= D:
            pieces.append((r, e))
        else:
            pieces.append((r, D))
            pieces.append((0, e - D))
    pieces.sort()
    total = 0
    cl, cr = pieces[0]
    for left, right in pieces[1:]:
        if left <= cr:
            cr = max(cr, right)
        else:
            total += cr - cl
            cl, cr = left, right
    total += cr - cl
    return total


def safe_cell_ranges(F):
    """Merged safe ranges (independent reimplementation of the definition)."""
    L = lcm(*(14 * v for v in F))
    danger = []
    for v in F:
        half = L // (14 * v)
        period = L // v
        for k in range(v + 1):
            c = k * period
            danger.append((max(0, c - half), min(L, c + half)))
    danger.sort()
    merged = []
    for left, right in danger:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    safe = []
    cur = 0
    for left, right in merged:
        if cur < left:
            safe.append((cur, left))
        cur = max(cur, right)
    if cur < L:
        safe.append((cur, L))
    return L, safe


def global_flip_scan():
    """All support-hard rows x shapes (2,2,d3), lcm=D, where coarse and exact
    double-2 tests disagree (either direction).  Arc arithmetic only."""
    flips = []
    n_bodies = 0
    n_rows = 0
    n_hard = 0
    for F in combinations(range(1, 15), 6):
        n_bodies += 1
        L, ranges = safe_cell_ranges(F)
        d = 1
        divs = []
        while d * d <= L:
            if L % d == 0:
                divs.append(d)
                if d != L // d:
                    divs.append(L // d)
            d += 1
        for D in divs:
            n_rows += 1
            if D % 2:
                continue
            S = merged_arc_length_mod(ranges, D)
            if Q(S, D) > CUTOFF:
                continue
            n_hard += 1
            cands = [D]
            if (D // 2) % 2 == 1 and D // 2 > 1:
                cands.append(D // 2)  # lcm(2, D/2)=D when D/2 odd
            for d3 in cands:
                if lcm(2, 2, d3) != D:
                    continue
                pcap = fibre_cap(D, d3, 2)
                amb = capacity(D, d3)
                exact_kill = (S % 2 == 0) and (S // 2 > pcap)
                coarse_kill = S > 2 * amb
                if exact_kill != coarse_kill:
                    flips.append((F, D, S, d3, S // 2, pcap, 2 * amb))
    print(f"bodies={n_bodies} divisor_rows={n_rows} support_hard_even_D={n_hard}")
    print(f"flip occurrences (coarse vs exact disagree): {len(flips)}")
    for f in flips:
        print(f"  F={f[0]} D={f[1]} S={f[2]} d3={f[3]} S/2={f[4]} pcap={f[5]} 2amb={f[6]}")
    return flips


def repo_cross_check():
    sys.path.insert(0, "/home/claude/math/04-computation")
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location(
        "sup",
        "/home/claude/math/04-computation/"
        "lrc14_two_drift_body_projection_support_thm2928.py",
    )
    sup = module_from_spec(spec)
    spec.loader.exec_module(sup)
    for F, D, S in (
        ((1, 4, 5, 7, 9, 11), 194040, 55392),
        ((1, 5, 7, 8, 9, 11), 388080, 109044),
    ):
        L, ranges = sup.safe_cell_ranges(F)
        assert L == D
        got = sup.support_size_bitset(D, ranges)
        assert got == S, (F, got, S)
        print(f"repo cross-check F={F}: L={L} S={got} OK")


def main():
    ok1 = check_row((1, 4, 5, 7, 9, 11), 194040, 55392, 194040, 13860)
    print()
    ok2 = check_row((1, 5, 7, 8, 9, 11), 388080, 109044, 388080, 27720)
    print()
    flips = global_flip_scan()
    print()
    repo_cross_check()
    print()
    claimed = {((1, 4, 5, 7, 9, 11), 194040), ((1, 5, 7, 8, 9, 11), 388080)}
    found = {(f[0], f[1]) for f in flips}
    print(f"rows confirmed: {ok1 and ok2}")
    print(f"flip set == claimed two rows: {found == claimed}")
    print("LANE-E VERDICT:",
          "PASS" if (ok1 and ok2 and found == claimed) else "FAIL")


if __name__ == "__main__":
    main()
