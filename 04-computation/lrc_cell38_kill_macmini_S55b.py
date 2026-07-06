#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S55b -- HYP-4222: THE 3/38 CELL KILL (with kps-4197's fold).

An attainer W (primitive, M(W) = 3/38 exactly) must satisfy, ALL AT ONCE:
  (A)  38-grid tight cover: for every odd unit a mod 38, some runner within
       dist_38 <= 3 (else margin 4/38 = 2/19 > 2/25 at that a, contradiction
       with M < 2/25 -- wait: M = 3/38 means NO t beats 3/38, so any a with
       all runners dist >= 4 gives margin >= 4/38 = 2/19 > 3/38: kills
       attainment).  In ban language (evens ban {0,+-1}-halves = 1 class;
       odds ban +-w^{-1} AND +-3w^{-1} = 2 classes): the level-4 bans COVER
       all 9 dilation classes mod 19.
  (B)  the level-3 witness: SOME class a0 escapes the level-3 bans (evens:
       same class as level 4!; odds: only +-w^{-1}), with the binder: an odd
       runner at exactly +-3 (its 3-ban covers a0).
  (C)  39-grid tightness: for every unit b mod 39, some runner within
       dist_39 <= 3 (39*3/38 = 3.07 => <= 3).  39 = 3*13: the band folds to
       (mod 3, mod 13) = {(0,0), (1,1), (2,2), (0,3)} up to sign.
  (D)  the gap profile's 13-pinning (pair-covering mod 13 or a 13-multiple)
       and mod-19 pinning (= (A)'s shadow), covering 2..12, spread, 24-comp.

THIS SCRIPT: enumerate the JOINT template space (mod 38, mod 39) at the level
of the class data the constraints see, and test satisfiability of
(A) + (B) + (C) + the mod-13 pair pinning + parity/3-part consistency.
The space: per element, (parity, r19, r3, r13) -- we enumerate the MULTISET
structure smartly: mod-19 pair-multiplicities (9 pairs, sum 12, >= 1 each by
(A)-shadow... not forced >= 1: (A) is about the level-4 ban cover which
includes odds' 3-bans), so we enumerate directly over per-element class
choices with symmetry reduction and constraint propagation.

PRAGMATIC REDUCTION (exact, no loss): all constraints are invariant under
the JOINT dilation action (a, b) -> (ua, vb)?? NO: the element data couple
r19/r13/r3/parity through the INTEGER w.  The honest finite test: enumerate
w-values as residues mod 2*3*13*19 = 1482 (the CRT class determines ALL
constraint participation), i.e. choose 12 residues in Z/1482 (with the
no-0-mod-38 and covering-shadow conditions), and check (A)-(D).  The raw
space C(1482,12) is out of reach -- but the SEARCH is for SATISFIABILITY:
we run (i) structured hill-climbing / randomized completion over the
constraint graph to hunt ANY satisfying template (if found: the fold does
NOT kill -- report the surviving template for the integral-lever phase);
(ii) an EXACT feasibility count on the dominant sub-symmetry: fix the
witness class a0 = [1] (dilation normalization) and the binder, propagate.

VERDICT semantics:
  - SAT template found  => fold-level survivor; print it (integral phase needed).
  - extensive UNSAT     => strong evidence the 38x39 fold alone kills; then
                           attempt the exact pigeonhole count for the proof.
"""
from math import gcd
import random, sys, time
from itertools import product

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(3838)

P19, P13 = 19, 13
def d19(x): x %= 19; return min(x, 19 - x)
def d38(x): x %= 38; return min(x, 38 - x)
def d39(x): x %= 39; return min(x, 39 - x)
def d13(x): x %= 13; return min(x, 13 - x)

UNITS38 = [a for a in range(1, 38) if gcd(a, 38) == 1]
UNITS39 = [b for b in range(1, 39) if gcd(b, 39) == 1]

def constraints(ws):
    """ws: 12 integer representatives (use residues mod 1482 lifted small).
    Returns dict of constraint truth."""
    # (A) 38-grid: every odd unit a: some w with d38(w*a) <= 3
    A = all(any(d38(w * a) <= 3 for w in ws) for a in UNITS38)
    # (B) witness: some odd unit a0: all d38(w*a0) >= 3, with equality (binder)
    B = False
    for a in UNITS38:
        ds = [d38(w * a) for w in ws]
        if min(ds) == 3:
            B = True
            break
    # (C) 39-grid: every unit b: some w with d39(w*b) <= 3
    C = all(any(d39(w * b) <= 3 for w in ws) for b in UNITS39)
    # (D13) mod-13 pair pinning (or 13-mult)
    if any(w % 13 == 0 for w in ws):
        D13 = True
    else:
        D13 = all(any(w % 13 in (b, 13 - b) for w in ws) for b in range(1, 7))
    # no runner 0 mod 38 (else no witness at 38) and 0 mod 39 shadow fine
    Z = all(w % 38 != 0 for w in ws)
    return A, B, C, D13, Z

def dq(x, q):
    x %= q
    return min(x, q - x)

def sat_score(ws):
    A, B, C, D13, Z = constraints(ws)
    s = 0
    s += sum(1 for a in UNITS38 if not any(d38(w * a) <= 3 for w in ws))
    s += 0 if B else 5
    s += sum(1 for b in UNITS39 if not any(d39(w * b) <= 3 for w in ws))
    s += 0 if Z else 50
    # FULL level-3/38 pinning: q <= 12: covering (some w == 0 mod q);
    # q in 13..25: near-unit (every unit a: some w*a in {0,+-1} mod q)
    for q in range(2, 13):
        if not any(w % q == 0 for w in ws):
            s += 3
    for q in range(13, 26):
        if any(w % q == 0 for w in ws):
            continue
        for b in range(1, q // 2 + 1):
            if gcd(b, q) != 1:
                continue
            if not any(w % q in (b, q - b) for w in ws):
                s += 1
    return s

log("PHASE 1b: hunt with the FULL level-3/38 pinning stack (covering 2..12 + near-unit 13..25)")
from math import lcm as _lcm
M = 1
for q in range(2, 40):
    M = _lcm(M, q)          # residues mod lcm(2..39): all constraints fold here
best = None
found = []
t0 = time.time()
restarts = 0
while time.time() - t0 < 900 and len(found) < 5:
    restarts += 1
    ws = [random.randrange(1, M) for _ in range(12)]
    ws = [w if w % 38 else w + 1 for w in ws]
    cur = sat_score(ws)
    for step in range(4000):
        if cur == 0:
            break
        i = random.randrange(12)
        old = ws[i]
        ws[i] = random.randrange(1, M)
        if ws[i] % 38 == 0:
            ws[i] += 1
        new = sat_score(ws)
        if new <= cur:
            cur = new
        else:
            ws[i] = old
    if best is None or cur < best[0]:
        best = (cur, list(ws))
    if cur == 0:
        found.append(list(ws))
        log(f"  SAT TEMPLATE: {sorted(w % M for w in ws)}")
log(f"restarts: {restarts}; best violation score: {best[0]}; SAT templates: {len(found)}")
if found:
    log("=> the 38x39+13 fold does NOT kill alone; integral levers needed on these:")
    for ws in found:
        A, B, C, D13, Z = constraints(ws)
        log(f"   {sorted(w % M for w in ws)}  A={A} B={B} C={C} D13={D13}")
else:
    log("=> NO satisfying template found (score floor > 0): the joint fold LOOKS unsatisfiable;")
    log("   next: the exact pigeonhole (phase 2, separate run) for the proof.")
log(f"[t = {time.time()-T0:.0f}s]")
