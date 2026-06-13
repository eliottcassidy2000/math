#!/usr/bin/env python3
"""claudebox-2026-06-11-S7 part B — the FIBERED-SHELL lattice closes the joint failures.

Fibered-shell witness criterion (to prove + verify): for q = d*m with d | 14, t = a/q,
gcd(a,q)=1: t is a strict witness (M > 1/14) iff every v in S has
(v*a mod q) not in +-{0,1,...,floor(q/14)}.
For a d-core runner v = d*w this is (w*a mod m) not in +-{0..floor(m/14)}: for m <= 13
the band is {0} — a pure divisor condition on the core/d. So the d-core enters its own
shell-m problem with the band rescaled by d: the THM-421 clock (m=1), the S643 window
(m -> infinity), and the 3-tower 27/9/3 (d=1) are all faces of ONE lattice
Q = {d*m : d in {1,2,7,14}, m <= 27}.

Test: re-generate the part-3 joint failures (same seeds/families), and for each find
the MINIMAL modulus q admitting a twisted witness, and whether q in Q. Also adversarial
push: try to construct configs with NO witness in Q (blocking the whole lattice)."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

TARGET = F(1, 14)

def band(q):
    b = q // 14
    if 14 * b == q: pass  # j = q/14 gives exactly 1/14, not strict: it IS in the band
    return b

def witness_at_q(S, q):
    """Return twist a (gcd(a,q)=1) with all (v*a mod q) outside +-{0..floor(q/14)}, or None."""
    B = q // 14
    prof = [v % q for v in S]
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        ok = True
        for p in prof:
            r = (p * a) % q
            r = min(r, q - r)
            if r <= B: ok = False; break
        if ok: return a
    return None

def min_witness_modulus(S, qmax=400):
    for q in range(2, qmax + 1):
        a = witness_at_q(S, q)
        if a is not None: return q, a
    return None, None

def normF(x):
    fr = x - int(x)
    if fr < 0: fr += 1
    return min(fr, 1 - fr)

# 0. PROVE-CHECK the criterion itself: t=a/q strict witness <-> band condition (random check)
random.seed(7)
bad = 0
for _ in range(3000):
    q = random.randint(2, 120); a = random.randint(1, q - 1)
    if gcd(a, q) != 1: continue
    v = random.randint(1, 400)
    lhs = normF(F(v * a, q)) > TARGET
    r = (v * a) % q; r = min(r, q - r)
    rhs = r > q // 14
    if lhs != rhs: bad += 1
print(f"criterion check (||va/q|| > 1/14 <-> (va mod q) outside +-floor(q/14)): mismatches = {bad}")

Q_LATTICE = sorted({d * m for d in (1, 2, 7, 14) for m in range(1, 28)} - {1})
print(f"fibered lattice Q (d|14, m<=27): {len(Q_LATTICE)} moduli, max {max(Q_LATTICE)}")

# 1. regenerate the joint failures from part 3 (same construction, fast head tests inline)
def head_B_fails(S):
    if any(v % 27 == 0 for v in S): pass
    else:
        units = [v % 27 for v in S if v % 3]
        if any(witness_at_q([u], 27) for u in [0]):  # placeholder, real test below
            pass
    # exact reuse: head B succeeds iff shell 3 or 9 free of multiples, or shell-27 unit twist
    if all(v % 3 for v in S): return False
    if all(v % 9 for v in S): return False
    if any(v % 27 == 0 for v in S): return True
    units = [v % 27 for v in S if v % 3]
    for a in range(1, 27):
        if gcd(a, 27) != 1: continue
        if all((a * u) % 27 not in (1, 26) for u in units): return False
    return True

def safe_window(S_other, t0):
    L, R = t0 - 1, t0 + 1
    for v in S_other:
        x = v * t0; kf = int(x)
        cl, cr = [], []
        for k in (kf - 1, kf, kf + 1):
            lo, hi = k - TARGET, k + TARGET
            if hi <= x: cl.append(hi)
            if lo >= x: cr.append(lo)
            if lo < x < hi: return None
        L = max(L, max(cl) / v if cl else t0 - 1)
        R = min(R, min(cr) / v if cr else t0 + 1)
        if L >= R: return None
    return (L, R)

def core_dodge_in(core, L, R):
    arcs = []
    for v in core:
        for k in range(int(v * L) - 1, int(v * R) + 2):
            a, b = F(k - TARGET) / v, F(k + TARGET) / v
            if b > L and a < R: arcs.append((a, b))
    arcs.sort(); cur = L
    for a, b in arcs:
        if a > cur: return (cur + min(a, R)) / 2
        cur = max(cur, b)
        if cur >= R: return None
    return (cur + R) / 2 if cur < R else None

def head_A_fails(S):
    for d in (14, 7, 2):
        core = [v for v in S if v % d == 0]
        rest = [v for v in S if v % d]
        if not core: return False
        for b in range(1, d):
            if gcd(b, d) != 1: continue
            win = safe_window(rest, F(b, d))
            if win and core_dodge_in(core, *win) is not None: return False
    return True

def joint_fail(S):
    return head_A_fails(S) and head_B_fails(S)

jf = []
for r in range(1, 1200):
    if r % 7 == 0: continue
    S = sorted([7 * k for k in range(1, 13)] + [r])
    if len(set(S)) != 13 or reduce(gcd, S) != 1: continue
    if joint_fail(S): jf.append(S)
for r in range(1, 1200):
    if r % 7 == 0 or r == 54: continue
    S = sorted([14 * k for k in range(1, 7)] + [7, 21, 35, 49, 77, 54] + [r])
    if len(set(S)) != 13 or reduce(gcd, S) != 1: continue
    if joint_fail(S): jf.append(S)
random.seed(14)
for trial in range(4000):
    c = random.randint(6, 12)
    ks = random.sample(range(1, 19), c)
    if not any((7 * k) % 14 == 0 for k in ks): ks[0] = 2 * random.randint(1, 9)
    S = [7 * k for k in ks] + random.sample([v for v in range(1, 120) if v % 7], 13 - c)
    S = sorted(S)
    if len(set(S)) != 13 or reduce(gcd, S) != 1: continue
    if not any(v % 14 == 0 for v in S): continue
    if joint_fail(S): jf.append(S)
print(f"joint failures regenerated: {len(jf)}")

# 2. minimal witness modulus for each; coverage by the fibered lattice Q
from collections import Counter
minq = Counter(); inQ = 0; outQ_examples = []
covered_by_Q = 0
for S in jf:
    q, a = min_witness_modulus(S, 400)
    minq[q] += 1
    if q in Q_LATTICE: inQ += 1
    # independent: does ANY Q-modulus work?
    okQ = any(witness_at_q(S, qq) for qq in Q_LATTICE)
    if okQ: covered_by_Q += 1
    elif len(outQ_examples) < 5: outQ_examples.append((S, q, a))
print(f"minimal witness moduli histogram: {dict(sorted(minq.items()))}")
print(f"minimal modulus already in Q: {inQ}/{len(jf)}")
print(f"covered by SOME Q-modulus:    {covered_by_Q}/{len(jf)}")
for S, q, a in outQ_examples:
    print(f"  NOT covered by Q: {S}  (min witness q={q}, a={a})")

# 3. adversarial: block the whole lattice Q? A config needs, for EVERY q in Q, either a
# runner killing all twists (multiple of q for the relevant stratum) or a +-band cover.
# Try greedy: start from a rich blocking skeleton and tune the free runners mod the
# prime shells in Q (13*7=91 needs mult of 91 or cover; band +-6 per runner => 12 bad
# twist-classes per runner mod 91 — coverable? phi(91)=72, 13 runners * 12 = 156 >= 72).
# Random adversarial search over residue-tuned configs:
random.seed(91)
def blocks_q(S, q):
    return witness_at_q(S, q) is None
best = None
for trial in range(2000):
    # skeleton: multiples chosen to kill the small shells + divisor clocks
    S = [42, 54, 56, 80, 88, 117, 50, 119, 11 * 23, 13 * 29]  # 42=2*3*7(14? no)...
    # ensure multiple of 14 and 13 distinct, then 3 free tuned runners
    S = [14 * 3, 54, 56, 80, 88, 117, 50, 119, 253, 377]
    free = random.sample(range(1, 1500), 3)
    T = sorted(set(S + free))
    if len(T) != 13 or reduce(gcd, T) != 1: continue
    nb = sum(1 for q in Q_LATTICE if blocks_q(T, q))
    if best is None or nb > best[0]: best = (nb, T)
print(f"adversarial lattice blocking: best blocked {best[0]}/{len(Q_LATTICE)} Q-moduli")
print(f"  config: {best[1]}")
unblocked = [q for q in Q_LATTICE if not blocks_q(best[1], q)][:10]
print(f"  first unblocked Q-moduli: {unblocked}")
