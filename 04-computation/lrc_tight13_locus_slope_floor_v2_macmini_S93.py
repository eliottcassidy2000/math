#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S93 -- HYP-3840 part 2 (v2, fast): tight-13-set locus mod 14 and the
exact slope floor at N=14.  Uses THM-593 as a filter:
  tight => every unit a/q is a witness with min exactly 1/q  (Cor A2)  [6-point filter]
  tight => dropped residues are non-units                    (Cor A3)
Families (residue classification HYP-3750 + element relifts, q=14):
  (A) perm-lift: residues={1..13}, one element lifted x -> x+14j, j<=6
  (A2) perm-2lift: two elements lifted, j<=3 each
  (B) dup+drop, no relift: drop v (non-unit), add s+14j, j<=6
  (B2) dup+drop + ONE relift: drop v (non-unit), dup s via s+14j1, relift x -> x+14j2
       (this is the q=8-beater pattern), j1,j2 <= 3
Tightness test: unit-witness filter (necessary), then early-exit scan of the binding-pair
grid for any min > 1/14 (if none, M = 1/14 exactly since 1/14 attained).
Output flushed section by section (lesson from the buffered-timeout loss).
"""
from fractions import Fraction as F
from math import gcd
import sys, itertools

Q = 14
UNITS = [a for a in range(1, Q) if gcd(a, Q) == 1]
NONUNITS = [a for a in range(1, Q) if gcd(a, Q) > 1]
RQ = F(1, Q)

def dist(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def unit_witness_filter(S):
    """Necessary for tight: min at every unit a/14 is exactly 1/14."""
    for a in UNITS:
        t = F(a, Q)
        if min(dist(v * t) for v in S) != RQ:
            return False
    return True

def is_tight(S):
    """Assuming unit_witness_filter passed (so M >= 1/14 attained), verify no t has
    min > 1/14 by scanning binding-pair candidates with early exit."""
    Sl = sorted(set(S)); dens = set()
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            dens.add(v + w)
            if w > v: dens.add(w - v)
    for den in sorted(dens):
        for m in range(1, den):
            t = F(m, den)
            mn = RQ
            ok_gt = True
            for v in Sl:
                d = dist(v * t)
                if d <= RQ:
                    ok_gt = False
                    break
            if ok_gt:
                return False  # found min > 1/14 => M > 1/14
    return True

def v_max_map(S):
    vm = {}
    for v in S:
        u = v % Q
        if gcd(u, Q) == 1:
            vm[u] = max(vm.get(u, 0), v)
    return vm

def slope_thm593(S):
    """c_S = (2/q) sum_{u unit} 1/v_max(u)  (THM-593 Part B)."""
    vm = v_max_map(S)
    if set(vm) != set(UNITS):
        return None
    return F(2, Q) * sum(F(1, vm[u]) for u in UNITS)

def slope_direct(S):
    """Independent check: measure in last linear cell below 1/14."""
    pts = set()
    Sl = sorted(set(S))
    for i, v in enumerate(Sl):
        for w in Sl[i:]:
            for den in ({v + w} | ({w - v} if w > v else set())):
                d = 1
                while F(d, den) < RQ:
                    pts.add(F(d, den)); d += 1
    last = max(pts) if pts else F(0)
    r1 = last + (RQ - last) / 3
    iv = []
    for v in Sl:
        half = r1 / v
        for k in range(v):
            c = F(k, v); a, b = c - half, c + half
            if a < 0: iv.append((a + 1, F(1))); iv.append((F(0), b))
            elif b > 1: iv.append((a, F(1))); iv.append((F(0), b - 1))
            else: iv.append((a, b))
    iv.sort(); tot = F(0); ca = cb = None
    for a, b in iv:
        if cb is None: ca, cb = a, b
        elif a <= cb: cb = max(cb, b)
        else: tot += cb - ca; ca, cb = a, b
    tot += cb - ca
    return (1 - tot) / (1 - Q * r1)

BASE = list(range(1, Q))
found = {}

def consider(S, tag):
    S = tuple(sorted(S))
    if len(S) != 13 or S in found or S == tuple(BASE):
        return
    if len(set(S)) != 13: return
    if not unit_witness_filter(S): return
    if not is_tight(S): return
    found[S] = tag
    print(f"  TIGHT [{tag}]  S={list(S)}", flush=True)

print("=" * 78, flush=True)
print("(A) perm-lift x->x+14j, j<=6", flush=True)
for x in BASE:
    for j in range(1, 7):
        consider([e for e in BASE if e != x] + [x + Q * j], f"perm-lift {x}->{x+Q*j}")

print("(A2) perm-2lift, j<=3 each", flush=True)
for x, y in itertools.combinations(BASE, 2):
    for j1 in range(1, 4):
        for j2 in range(1, 4):
            consider([e for e in BASE if e not in (x, y)] + [x + Q * j1, y + Q * j2],
                     f"perm-2lift {x}->{x+Q*j1},{y}->{y+Q*j2}")

print("(B) dup+drop (drop non-unit v, dup s via s+14j), j<=6", flush=True)
for v in NONUNITS:
    for s in range(1, Q):
        if s == v: continue
        for j in range(1, 7):
            consider([e for e in BASE if e != v] + [s + Q * j],
                     f"drop{v} dup{s} via {s+Q*j}")

print("(B2) dup+drop + one relift, j<=3 each  (the q=8-beater pattern)", flush=True)
for v in NONUNITS:
    keep = [e for e in BASE if e != v]
    for s in range(1, Q):
        if s == v: continue
        for j1 in range(1, 4):
            dup_el = s + Q * j1
            for x in keep:
                for j2 in range(1, 4):
                    S = [e for e in keep if e != x] + [x + Q * j2, dup_el]
                    consider(S, f"drop{v} dup{s} via {dup_el} relift {x}->{x+Q*j2}")

print("(sanity) drops of UNIT residues (lemma predicts none tight), sample j<=2", flush=True)
lemma_violations = 0
for v in UNITS:
    for s in range(1, Q):
        if s == v: continue
        for j in (1, 2):
            S = tuple(sorted([e for e in BASE if e != v] + [s + Q * j]))
            if len(set(S)) == 13 and unit_witness_filter(S) and is_tight(S):
                lemma_violations += 1
                print(f"  !! LEMMA VIOLATION drop{v} dup{s}", flush=True)
print(f"unit-drop tight sets found: {lemma_violations} (THM-593 predicts 0)", flush=True)

print("=" * 78, flush=True)
c_ap = F(2, Q) * sum(F(1, a) for a in UNITS)
print(f"C_AP(14) = {c_ap} = {float(c_ap):.6f}", flush=True)
print(f"non-AP tight 13-sets found: {len(found)}", flush=True)
floor = c_ap
for S, tag in sorted(found.items()):
    c593 = slope_thm593(list(S))
    cdir = slope_direct(list(S))
    match = "OK" if c593 == cdir else f"MISMATCH direct={cdir}"
    rel = "= C_AP (rigid)" if cdir == c_ap else (f"< C_AP  <<<< BEATER" if cdir < c_ap else "> C_AP")
    floor = min(floor, cdir)
    print(f"  slope {cdir} = {float(cdir):.6f}  {rel}  [thm593 {match}]  {tag}", flush=True)
    print(f"     S={list(S)}", flush=True)
print(f"\nSLOPE FLOOR at N=14 over enumerated tight locus: {floor} = {float(floor):.6f}", flush=True)
print("(families: perm-lift j<=6, perm-2lift j<=3, dup+drop j<=6, dup+drop+relift j<=3)", flush=True)
