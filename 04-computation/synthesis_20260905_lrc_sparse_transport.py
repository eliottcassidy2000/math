#!/usr/bin/env python3
"""Exact audits of sparse-comb transport; no floating-point decisions.

The analytic theorem is in the matching report.  The selected certificate is
computed without overlap lengths; physical intersection is a separate control.
Run from the repository root with python -B (also under -O).
"""
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import gcd, lcm
from random import Random
from hashlib import sha256
import json
import sys

CHECKS = 0
PAIRS = ((0, 1), (0, 2), (1, 2))


def check(value, detail):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(detail)


def sheet(a, j, scale, radius=Q(1, 14)):
    # scale is a common multiple of 3*radius.denominator*a.
    rr = scale * radius.numerator // (radius.denominator * a)
    out = []
    for k in range(a):
        center = (scale * k // a - scale * j // 3) % scale
        lo, hi = center - rr, center + rr
        if lo < 0:
            out.extend(((0, hi), (lo + scale, scale)))
        elif hi > scale:
            out.extend(((lo, scale), (0, hi - scale)))
        else:
            out.append((lo, hi))
    return sorted(out)


def meet(aa, bb):
    i = j = 0
    out = []
    while i < len(aa) and j < len(bb):
        lo, hi = max(aa[i][0], bb[j][0]), min(aa[i][1], bb[j][1])
        if lo < hi:
            out.append((lo, hi))
        aend, bend = aa[i][1], bb[j][1]
        i += aend <= bend
        j += bend <= aend
    return out


def network(aa, bb):
    # Records contact and endpoint-free capacities only.
    i = j = 0
    loads_a = [0] * len(aa)
    loads_b = [0] * len(bb)
    edges = []
    total = 0
    while i < len(aa) and j < len(bb):
        if max(aa[i][0], bb[j][0]) < min(aa[i][1], bb[j][1]):
            cap = min(aa[i][1] - aa[i][0], bb[j][1] - bb[j][0])
            loads_a[i] += cap
            loads_b[j] += cap
            total += cap
            edges.append((i, j))
        aend, bend = aa[i][1], bb[j][1]
        i += aend <= bend
        j += bend <= aend
    feasible = (all(v <= b-a for v, (a, b) in zip(loads_a, aa)) and
                all(v <= b-a for v, (a, b) in zip(loads_b, bb)))
    return total, feasible, edges


def row(w, radius=Q(1, 14), audit_reference=False):
    scale = 3 * radius.denominator * lcm(*w)
    sheets = {(i, j): sheet(a, j, scale, radius) for i, a in enumerate(w)
              for j in range(3)}
    uppers = [0, 0, 0]
    controls = [0, 0, 0]
    branches = 0
    for perm in permutations((0, 1, 2)):
        for pi, (i, j) in enumerate(PAIRS):
            k = 3-i-j
            aa = meet(sheets[i, perm[i]], sheets[j, perm[j]])
            bb = sheets[k, perm[k]]
            cap, feasible, edges = network(aa, bb)
            check(feasible, ("edge-min capacity violation", w, perm, (i, j)))
            delta_a = 2*scale*radius.numerator//(radius.denominator*max(w[i], w[j]))
            check(all(y[0]-x[1] >= 2*delta_a for x, y in zip(aa, aa[1:])),
                  ("pair inherited separation", w, perm, (i, j)))
            deg_a = [0]*len(aa)
            deg_b = [0]*len(bb)
            for a, b in edges:
                deg_a[a] += 1
                deg_b[b] += 1
            check(all(deg_a[a] == 1 or deg_b[b] == 1 for a, b in edges),
                  ("contact graph is a star forest", w, perm, (i, j)))
            branches += sum(d >= 2 for d in deg_a + deg_b)
            uppers[pi] += cap
            # Physical control follows selection of this contact certificate.
            controls[pi] += sum(b-a for a, b in meet(aa, bb))
            if audit_reference:
                from lrc14_third_sheet_component_network_thm4409 import bipartite_capacity
                ref = bipartite_capacity(aa, bb, edges)
                check(ref == cap, ("independent augmenting-flow", w, perm, (i, j)))
    check(len(set(controls)) == 1, ("pair-independent physical mass", w))
    check(all(c >= controls[0] for c in uppers), ("certificate upper", w))
    upper = tuple(Q(v, scale) for v in uppers)
    mass = Q(controls[0], scale)
    return upper, mass, branches


def raw_roof(w):
    """Independent raw-lattice formula, retaining exact positive support.

    Deliberately no physical interval code or contact construction is called.
    """
    r = Q(3, 14)
    bounds = [int(r*(w[j]+w[k])) for j, k in ((1, 2), (0, 2), (0, 1))]
    roof = [2*r/a for a in w]
    out = [Q(0)]*3
    mass = Q(0)
    carriers = []
    for x in range(-bounds[0], bounds[0]+1):
        for y in range(-bounds[1], bounds[1]+1):
            znum = -w[0]*x-w[1]*y
            if znum % w[2]:
                continue
            c = (x, y, znum//w[2])
            if abs(c[2]) > bounds[2] or any(v % 3 == 0 for v in c):
                continue
            pairs = [r/w[i]+r/w[j]-Q(abs(c[3-i-j]), w[i]*w[j]) for i, j in PAIRS]
            length = min(roof+pairs)
            if length <= 0:
                continue
            kernels = [min(roof+[p]) for p in pairs]
            check(min(kernels) == length, ("pointwise pair selection", w, c))
            mass += length
            for i, value in enumerate(kernels):
                out[i] += value
            carriers.append((c, str(length), tuple(map(str, kernels))))
    return tuple(out), mass, carriers


def main():
    sys.stdout.reconfigure(newline="\n")
    # Universe A independently replays all THM-4409 triples with integer grids.
    values = [v for v in range(1, 80, 2) if v % 3]
    base = [w for w in combinations(values, 3) if gcd(*w) == 1]
    check(len(base) == 2910, "base universe cardinality")
    # Universe B is deterministic hostile probing, not a complete height census.
    rng = Random(20260905)
    samples = set()
    vals = [v for v in range(1, 1002, 2) if v % 3]
    while len(samples) < 400:
        w = tuple(sorted(rng.sample(vals, 3)))
        if gcd(*w) == 1:
            samples.add(w)
    # Near equality resonance, very unbalanced scales, and shifted fixed gaps.
    for t in (1, 5, 11, 29, 71, 149, 311):
        for offsets in ((0, 2, -2), (2, 0, 2), (0, -2, 2), (2, 4, -2)):
            w = tuple(sorted(c*t+d for c, d in zip((1, 5, 11), offsets)))
            if len(set(w)) == 3 and all(v > 0 and v % 2 and v % 3 for v in w) and gcd(*w) == 1:
                samples.add(w)
    for t in (79, 101, 211, 503, 1001, 2003):
        for a, b in ((1, 5), (1, 11), (5, 11), (11, 13)):
            w = (a, b, t)
            if len(set(w)) == 3 and t % 3 and gcd(*w) == 1:
                samples.add(w)
    records = []
    base_exact = equality = branches = 0
    strict_min = None
    for w in base:
        caps, mass, bran = row(w)
        chosen = min(caps)
        check(chosen <= Q(6, 77), ("base sharp certificate", w, caps))
        base_exact += chosen == mass
        equality += chosen == Q(6, 77)
        branches += bran
        records.append((w, tuple(map(str, caps)), str(mass)))
    check((base_exact, equality) == (1747, 1), "inherited exactness/equality controls")
    sample_failures = []
    sample_max = None
    for w in sorted(samples):
        caps, mass, bran = row(w)
        branches += bran
        chosen = min(caps)
        if chosen > Q(6, 77):
            sample_failures.append((w, str(chosen), str(mass)))
        candidate = (chosen, w, mass)
        sample_max = max(sample_max, candidate) if sample_max else candidate
        records.append((w, tuple(map(str, caps)), str(mass)))
    # Actual max-flow replay on equality, strict loss, and branching controls.
    for w in ((1, 5, 11), (1, 19, 79), (1, 5, 2003), (11, 13, 1001)):
        row(w, audit_reference=True)
    # The new raw-roof formula is audited independently over a complete low
    # universe plus the canonical strict-loss and equality controls.
    raw_values = [v for v in range(1, 44, 2) if v % 3]
    raw_universe = [w for w in combinations(raw_values, 3) if gcd(*w) == 1]
    raw_universe.extend(((1, 19, 79), (1, 5, 101), (11, 13, 101)))
    for w in raw_universe:
        caps, mass, _ = row(w)
        raw_caps, raw_mass, _ = raw_roof(w)
        check(raw_caps == caps and raw_mass == mass, ("raw-roof/interval dictionary", w))
    hostile_caps, hostile_mass, hostile_carriers = raw_roof((1, 19, 79))
    # Sharpness of the generic gap coefficient two.
    aa = [(Q(9, 10), Q(13, 5))]
    bb = [(Q(0), Q(1)), (Q(5, 2), Q(7, 2))]
    cap, feasible, edges = network(aa, bb)
    check(cap == 2 and not feasible and aa[0][1]-aa[0][0] == Q(17, 10),
          "gap=3/2 hostile must exceed source capacity")
    # Nonconstant weights can break edgewise feasibility despite gap factor six.
    # I=[.9,7.1], J1=[0,1], J2=[7,8]; phi=1 on outside caps,
    # phi=1/100 on I.  Each leaf has >1/10 mass, center only 31/500.
    center_mass = Q(31, 500)
    leaf_mass = Q(9, 10) + Q(1, 1000)
    weighted_edge_sum = 2*min(center_mass, leaf_mass)
    check(weighted_edge_sum > center_mass, "arbitrary weight hostile")
    print("status=PROVED_ANALYTICALLY report; FINITE_EXACT independent computational checks")
    print("theorem=gap >= 2*maximum component length on each shore makes edge minima feasible")
    print("comb_application=every positive integer speed, arbitrary phases, radius<=1/6, every finite intersection partition")
    print("base_universe=2910 primitive distinct odd ternary-unit triples through 79")
    print("base_min_capacity_exact="+str(base_exact)+"; target_equalities="+str(equality))
    print("sample_universe=deterministic stress probes; count="+str(len(samples))+"; max_height="+str(max(w[-1] for w in samples)))
    print("sample_target_failures="+str(sample_failures))
    print("sample_max_capacity="+str(sample_max))
    print("branching_vertices_checked="+str(branches))
    print("maxflow_independent_controls=4 triples, all 18 sheet/pair choices")
    print("raw_roof_independent_controls="+str(len(raw_universe))+" triples")
    print("strict_loss_hostile=(1,19,79); capacities="+str(hostile_caps)+"; mass="+str(hostile_mass))
    print("strict_loss_raw_carriers="+str(hostile_carriers))
    print("generic_gap_hostile=gap 3/2, edge-min sum 2, maxflow 17/10")
    print("weighted_gap_six_hostile=center 31/500, edge-min sum 31/250")
    print("checks="+str(CHECKS))
    print("rows_sha256="+sha256(json.dumps(records, sort_keys=True).encode()).hexdigest())
    print("LRC14=OPEN; all-height sharp 6/77 bound=OPEN")


if __name__ == "__main__":
    main()
