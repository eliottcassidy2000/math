#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S93 -- HYP-3840 part 3: the tight-set slope landscape across moduli.
For q = 5..16: enumerate tight (q-1)-sets (perm-lift, perm-2lift, dup+drop, dup+drop+relift;
THM-593 filters), report the slope floor vs C_AP(q) and whether rigidity holds.
Answers: at which moduli is "the AP is the slope minimizer" true?  (q=8 known beater.)
"""
from fractions import Fraction as F
from math import gcd
import itertools

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def sweep_q(Q, jmax_single=6, jmax_multi=3):
    UNITS = [a for a in range(1, Q) if gcd(a, Q) == 1]
    NONUNITS = [a for a in range(1, Q) if gcd(a, Q) > 1]
    RQ = F(1, Q)
    BASE = list(range(1, Q))

    def unit_witness_filter(S):
        for a in UNITS:
            t = F(a, Q)
            if min(dist(v * t) for v in S) != RQ:
                return False
        return True

    def is_tight(S):
        Sl = sorted(set(S)); dens = set()
        for i, v in enumerate(Sl):
            for w in Sl[i:]:
                dens.add(v + w)
                if w > v: dens.add(w - v)
        for den in sorted(dens):
            for m in range(1, den):
                t = F(m, den)
                if all(dist(v * t) > RQ for v in Sl):
                    return False
        return True

    found = {}
    def consider(S, tag):
        S = tuple(sorted(S))
        if len(S) != Q - 1 or len(set(S)) != Q - 1 or S in found or S == tuple(BASE):
            return
        if not unit_witness_filter(S): return
        if not is_tight(S): return
        found[S] = tag

    for x in BASE:
        for j in range(1, jmax_single + 1):
            consider([e for e in BASE if e != x] + [x + Q * j], f"perm-lift {x}+{Q*j}")
    for x, y in itertools.combinations(BASE, 2):
        for j1 in range(1, jmax_multi + 1):
            for j2 in range(1, jmax_multi + 1):
                consider([e for e in BASE if e not in (x, y)] + [x + Q * j1, y + Q * j2],
                         f"perm-2lift {x},{y}")
    for v in NONUNITS:
        for s in range(1, Q):
            if s == v: continue
            for j in range(1, jmax_single + 1):
                consider([e for e in BASE if e != v] + [s + Q * j],
                         f"drop{v} dup{s}@{s+Q*j}")
    for v in NONUNITS:
        keep = [e for e in BASE if e != v]
        for s in range(1, Q):
            if s == v: continue
            for j1 in range(1, jmax_multi + 1):
                for x in keep:
                    for j2 in range(1, jmax_multi + 1):
                        consider([e for e in keep if e != x] + [x + Q * j2, s + Q * j1],
                                 f"drop{v} dup{s}@{s+Q*j1} relift{x}->{x+Q*j2}")
    return found

import functools
from lonely_profile import profile  # direct exact slope (run from 04-computation/ or with path)

def slope_direct(S, Q):
    return profile(sorted(S), F(1, Q)).critical_slope(Q)

def slope_formula(S, Q):
    """(2/q) sum 1/v_max(u): exact iff argmax f_S = unit fractions (clean tight sets);
    a lower bound otherwise (extra witnesses only add components)."""
    UNITS = [a for a in range(1, Q) if gcd(a, Q) == 1]
    vm = {}
    for v in S:
        u = v % Q
        if gcd(u, Q) == 1:
            vm[u] = max(vm.get(u, 0), v)
    return F(2, Q) * sum(F(1, vm[u]) for u in UNITS)

print("PRIMITIVE (gcd=1) tight sets only; slopes verified DIRECTLY (formula cross-checked)")
print(f"{'q':>3} {'C_AP(q)':>12} {'#tight(non-AP)':>15} {'slope floor':>14} {'rigid?':>7}  beaters")
for Q in range(5, 17):
    c_ap = F(2, Q) * sum(F(1, a) for a in range(1, Q) if gcd(a, Q) == 1)
    found = {S: tag for S, tag in sweep_q(Q).items()
             if functools.reduce(gcd, S) == 1}
    slopes = {S: slope_direct(list(S), Q) for S in found}
    floor = min([c_ap] + list(slopes.values()))
    beaters = [(list(S), str(c)) for S, c in slopes.items() if c < c_ap]
    rigid = "YES" if all(c == c_ap for c in slopes.values()) else "no"
    print(f"{Q:>3} {str(c_ap):>12} {len(found):>15} {str(floor):>14} {rigid:>7}  "
          f"{beaters if beaters else ''}", flush=True)
    for S, tag in sorted(found.items()):
        cf = slope_formula(list(S), Q)
        note = "formula=direct" if cf == slopes[S] else f"formula {cf} < direct (extra witnesses)"
        print(f"      {tag:45s} slope={slopes[S]}  [{note}]  S={list(S)}", flush=True)
