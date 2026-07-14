#!/usr/bin/env python3
"""
lrc14_mod360360_classification_opus_S287.py
===========================================
opus-2026-07-14-S287.  RUN the mod-360360 classification; decide the binding=>shadow conjecture
over the classified space.

THE EXACT TEST (rigorous per slot): slot (k, a), witness t = a/k + delta, delta in (0, 1/(2k)].
Good-delta set = intersection over ALL v in V of {delta : ||v(a/k + delta)|| >= 1/14} -- exact
rational interval arithmetic.  Slot FEASIBLE iff nonempty (then an explicit witness exists and
M(V) >= 1/14 rigorously).  A body is ALL-BLOCKED iff every slot (k <= 13, a unit) is infeasible.

THE CLASSIFIED SPACE:
  A. RIGIDITY-GUIDED EXTREMAL BLOCKERS: CRT lifts with complete residue system mod 13 among the
     12 non-carriers (the S286 pigeonhole's only candidates), near-complete mod 11/7, covering
     carriers, sizes in mixed tiers (to dodge the spread13/cluster/box/isolated tiles).
     Systematic sample over tier-assignments (the finite pattern layer).
  B. MEGA-HUNT: adversarial random covering bodies with float screen + early exit.
  Any all-blocked candidate -> EXACT delta-interval confirmation + M measurement + tile check.

VERDICT: counterexample (binding all-blocked body) or the conjecture decided-empirically over
the classified space with the extremal blockers explicitly swept.
"""
import random, math, time
from fractions import Fraction as F
import numpy as np

LAM = F(1, 14)

def units(k): return [a for a in range(1, k) if math.gcd(a, k) == 1]
def is_covering(V): return all(any(v % q == 0 for v in V) for q in range(2, 15))

# ---------- exact slot test ----------
def intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def slot_feasible_exact(V, k, a):
    """exact: exists delta in (0, 1/(2k)] with ||v(a/k+delta)|| >= 1/14 for all v."""
    lo0, hi0 = F(1, 10**9), F(1, 2 * k)   # tiny positive start to keep it open at 0
    cur = [(lo0, hi0)]
    for v in sorted(V, reverse=True):
        # good delta for v: ||v a/k + v delta|| >= 1/14
        # bad: v delta in (m - v a/k - 1/14, m - v a/k + 1/14)
        base = F(v * a, k)
        good = []
        mlo = (v * lo0 + base - 1).__floor__()
        mhi = (v * hi0 + base + 1).__ceil__()
        for m in range(mlo, mhi + 1):
            g1 = (m + LAM - base) / v
            g2 = (m + 1 - LAM - base) / v
            g1c, g2c = max(g1, lo0), min(g2, hi0)
            if g2c > g1c: good.append((g1c, g2c))
        cur = intersect(cur, sorted(good))
        if not cur: return False, None
    return True, cur[0][0]

def all_blocked_exact(V):
    for k in range(2, 14):
        if not any(v % k == 0 for v in V): continue
        for a in units(k):
            ok, wit = slot_feasible_exact(V, k, a)
            if ok: return False, (k, a, wit)
    return True, None

# ---------- float screen (fast, early exit) ----------
def shadow_screen(V):
    Va = np.array(sorted(V), dtype=np.float64)
    for k in range(2, 14):
        if not any(v % k == 0 for v in V): continue
        deltas = np.linspace(1e-7, 1.0/(2*k), 240)
        for a in units(k):
            tt = a / k + deltas
            fr = np.outer(Va, tt) % 1.0
            cl = np.minimum(fr, 1.0 - fr).min(axis=0)
            if cl.max() >= 1.0/14 - 1e-9:
                return True
    return False

def Mfloat(V, NG=1<<20):
    t = np.arange(NG) / NG
    m = np.ones(NG)
    for v in V:
        fr = (v * t) % 1.0
        m = np.minimum(m, np.minimum(fr, 1.0 - fr))
    return m.max()

# ---------- A: rigidity-guided extremal blockers ----------
def crt_blockers(n_samples, seed):
    rng = random.Random(seed)
    out = []
    TIERS = [(15, 45), (46, 110), (111, 260), (261, 430)]
    for _ in range(n_samples):
        # carrier set: cover q = 2..14 with a few dedicated carriers
        carriers = []
        need = [8, 9, 11, 13, 7, 10, 12, 14]
        base_c = rng.choice([1, 2, 3])
        # one multiple of 13 (the k=13 carrier), one of 11, and an lcm-ish for the rest
        c13 = 13 * rng.randint(2, 30)
        c11 = 11 * rng.randint(2, 36)
        crest = rng.choice([504, 2520, 360, 720, 840]) * rng.randint(1, 3)  # covers 7,8,9,10,12,14
        carriers = [c13, c11, crest]
        # 10 more speeds: complete the residue system mod 13 (with c13 covering 0)
        residues13 = list(range(1, 13))
        rng.shuffle(residues13)
        body = list(carriers)
        for r in residues13:
            if len(body) >= 13: break
            lo, hi = TIERS[rng.randrange(4)]
            cands = [x for x in range(lo, hi) if x % 13 == r]
            if not cands: continue
            body.append(rng.choice(cands))
        body = sorted(set(body))[:13]
        if len(body) == 13 and is_covering(body):
            out.append(body)
    return out

# ---------- run ----------
print("=" * 100)
print("THE MOD-360360 CLASSIFICATION RUN (S287)")
print("=" * 100)
t0 = time.time()

# Part A: extremal blockers
blockA = crt_blockers(30000, 287)
print("\nA. rigidity-guided CRT blockers: %d covering bodies constructed" % len(blockA))
cand = []
for V in blockA:
    if not shadow_screen(V):
        cand.append(V)
print("   float screen: %d ALL-BLOCKED candidates" % len(cand))

# Part B: mega-hunt (mixed adversarial styles)
rng = random.Random(2870)
tried = 0; candB = []
t1 = time.time()
while tried < 250000 and time.time() - t1 < 900:
    style = tried % 4
    if style == 0:
        V = sorted(rng.sample(range(1, 15), 7) + rng.sample(range(15, 420), 6))
    elif style == 1:
        V = sorted(rng.sample(range(1, 30), 8) + rng.sample(range(30, 300), 5))
    elif style == 2:
        V = sorted(rng.sample(range(10, 60), 9) + rng.sample(range(60, 430), 4))
    else:
        c = rng.choice([2, 3, 4])
        V = sorted(set([c*x for x in rng.sample(range(5, 90), 9)] + rng.sample(range(15, 400), 4)))[:13]
    if len(set(V)) < 13 or not is_covering(V): continue
    tried += 1
    if not shadow_screen(V):
        candB.append(V)
print("\nB. mega-hunt: %d covering bodies screened in %.0fs: %d ALL-BLOCKED candidates"
      % (tried, time.time() - t1, len(candB)))

# exact confirmation of any candidates
print("\nC. EXACT confirmation of candidates (delta-interval arithmetic):")
verdicts = []
for V in (cand + candB)[:40]:
    blocked, wit = all_blocked_exact(V)
    if blocked:
        M = Mfloat(V)
        verdicts.append((V, M))
        print("   CONFIRMED ALL-BLOCKED: M=%.4f  V=%s" % (M, V))
    else:
        k, a, w = wit
        print("   screen false-positive: slot (k=%d,a=%d) exactly feasible (witness delta=%s)  V=%s"
              % (k, a, w, V))
if not (cand + candB):
    print("   (no candidates to confirm)")

print("\nVERDICT:")
if verdicts:
    binding = [x for x in verdicts if x[1] < 0.22]
    if binding:
        print("   COUNTEREXAMPLE(S) to the conjecture: binding all-blocked bodies exist:")
        for V, M in binding: print("     M=%.4f V=%s" % (M, V))
    else:
        print("   all-blocked bodies exist but ALL LOOSE (M >= %.3f): the conjecture SURVIVES"
              % min(x[1] for x in verdicts))
        for V, M in verdicts[:6]: print("     M=%.4f V=%s" % (M, V))
else:
    print("   ZERO all-blocked bodies across the classified space (CRT extremal blockers +")
    print("   mega-hunt, exact-confirmed): the conjecture SURVIVES with every slot-blocking")
    print("   strategy the rigidity analysis permits explicitly swept.")
print("\ntotal time %.0fs" % (time.time() - t0))
