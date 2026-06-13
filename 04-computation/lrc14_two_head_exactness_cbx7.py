#!/usr/bin/env python3
"""claudebox-2026-06-11-S7 — HYP-2436 Test 3: does C'(14) reduce EXACTLY to
C'(5)-on-the-3-core  UNION  the THM-421 mod-2/mod-7 fiber?

Parts:
  0. Primitivity hygiene: C'(14) as literally stated (no gcd condition) is false:
     2*{1..13} contains 14 and is tight. All work below uses gcd(S)=1.
  1. DEGENERATION LEMMA (exact, exhaustive over residues): at any unit twist a of
     shell 27, every runner with 3|v, 27 nmid v has ||v a/27|| >= 1/9 > 1/14.
     Hence at level 1/14 the 3-core carries only the divisor condition (27|v),
     i.e. the descended shell-9 problem has forbidden band {0}; C'(5) proper
     (level 1/5 on shell 9) has band {0,+-1}. The descent DEGENERATES: the
     3-core head is divisor-theoretic, NOT C'(5).
  2. NON-NECESSITY witness: explicit primitive multiple-of-14 config whose 3-core
     is the TIGHT n=5 AP {1,2,3,4} (so 'C'(5)-on-the-core' is blocked at the
     n=5 band) yet the config is loose with a shell-27 witness.
  3. NON-SUFFICIENCY hunt: exact implementations of
       HEAD A (fiber): d-clock windows, d in {2,7,14}, all b: exact R-safe window
                       around b/d, exact core-arc coverage check inside it.
       HEAD B (3-tower): shell-27 unit twists (band {0,+-1}) using Lemma 1,
                       shell-9 (band {0}), shell-3 (band {0}).
     Search structured + random primitive multiple-of-14 configs failing BOTH
     heads; verify each failure is still loose (exact M) and locate its witness.
Normalization: canon (THM-398): S = 13 distinct positive ints, M(S)=max_t min ||vt||,
LRC(14): M >= 1/14; C'(14): primitive S containing a multiple of 14 => M > 1/14.
Exact arithmetic throughout (integers / Fraction); no floats in any verdict.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools, sys

TARGET = F(1, 14)

def normF(x):  # ||x|| for Fraction x
    fr = x - int(x)
    if fr < 0: fr += 1
    return min(fr, 1 - fr)

def M_exact(S):
    """Exact M(S) = max_t min_v ||vt||. Optimum is at a crossing ||vt||=||wt||
    (denominator v+w or |v-w|) or a single-runner peak (denominator 2v)."""
    dens = set()
    for v in S:
        dens.add(2 * v)
        for w in S:
            if w > v:
                dens.add(v + w); dens.add(w - v)
    dens.discard(0)
    best_num, best_den, best_t = 0, 1, F(0)
    for d in sorted(dens):
        prof = [v % d for v in S]
        for k in range(1, d):
            m = d  # min over v of min(vk mod d, d - vk mod d), times 1/d
            for p in prof:
                r = (p * k) % d
                r = min(r, d - r)
                if r < m:
                    m = r
                    if m * best_den <= best_num * d: break
            else:
                pass
            if m * best_den > best_num * d:
                best_num, best_den, best_t = m, d, F(k, d)
    return F(best_num, best_den), best_t

# ---------- HEAD B: the 3-adic shell tower 27 -> 9 -> 3 ----------
def head_B(S):
    """Twisted witness on the 3-tower. Returns (True, t) or (False, reason)."""
    # shell 3 and shell 9: band {0} at level 1/14 (nonzero residue => >= 1/9 > 1/14)
    if all(v % 3 for v in S):
        return True, F(1, 3)
    if all(v % 9 for v in S):
        return True, F(1, 9)
    # shell 27, unit twists: danger = {27|v} u {units v with a*v = +-1 (mod 27)}
    if any(v % 27 == 0 for v in S):
        return False, "mult-of-27 blocks all shell-27 twists; 9,3 blocked too"
    units = [v % 27 for v in S if v % 3 != 0]
    classes = {min(u, 27 - u) for u in units}
    for a in range(1, 27):
        if gcd(a, 27) != 1: continue
        if all((a * u) % 27 not in (1, 26) for u in units):
            return True, F(a, 27)
    return False, f"unit classes {sorted(classes)} +-cover all 9 classes"

# ---------- HEAD A: divisor-clock windows (THM-421 / S643), d in {2,7,14} ----------
def safe_window(S_other, t0):
    """Maximal open interval (L,R) around t0 on which every v in S_other has
    ||vt|| > 1/14, computed exactly. Danger arcs of v: [(k-1/14)/v, (k+1/14)/v]."""
    L, R = t0 - 1, t0 + 1
    for v in S_other:
        x = v * t0  # Fraction
        # nearest danger interval edges around x in Z +- 1/14
        kf = int(x)  # floor for positive x
        cands_left, cands_right = [], []
        for k in (kf - 1, kf, kf + 1):
            lo, hi = k - TARGET, k + TARGET
            if hi <= x: cands_left.append(hi)
            if lo >= x: cands_right.append(lo)
            if lo < x < hi:  # t0 strictly inside v's danger: window empty
                return None
        lv = max(cands_left) / v if cands_left else t0 - 1
        rv = min(cands_right) / v if cands_right else t0 + 1
        L, R = max(L, lv), min(R, rv)
        if L >= R: return None
    return (L, R)

def core_dodge_in(core, L, R):
    """Exact: is there t in open (L,R) with ||vt|| > 1/14 for all v in core?
    core danger arcs are closed [(k-1/14)/v, (k+1/14)/v]; check union covers (L,R)."""
    if L >= R: return None
    arcs = []
    for v in core:
        k0 = int(v * L) - 1
        k1 = int(v * R) + 2
        for k in range(k0, k1):
            a, b = F(k - TARGET, 1) / v if False else (k - TARGET) / F(v), (k + TARGET) / F(v)
            a, b = F(k - TARGET) / v, F(k + TARGET) / v
            if b > L and a < R:
                arcs.append((a, b))
    arcs.sort()
    cur = L
    for a, b in arcs:
        if a > cur:  # gap (cur, a) is core-safe; any t inside works
            return (cur + min(a, R)) / 2
        cur = max(cur, b)
        if cur >= R: return None
    return (cur + R) / 2 if cur < R else None

def head_A(S):
    """Fiber head: for d in {2,7,14}, all b coprime d: window around b/d safe for
    non-multiples-of-d, dodge the d-core inside it. Returns (True,t) or (False,info)."""
    for d in (14, 7, 2):
        core = [v for v in S if v % d == 0]
        rest = [v for v in S if v % d]
        if not core:
            return True, F(1, d)  # d-clock itself is a clean witness
        for b in range(1, d):
            if gcd(b, d) != 1: continue
            t0 = F(b, d)
            win = safe_window(rest, t0)
            if win is None: continue
            t = core_dodge_in(core, win[0], win[1])
            if t is not None:
                # certify
                assert all(normF(v * t) > TARGET for v in S), (S, d, b, t)
                return True, t
    return False, "all d-windows (d=2,7,14) covered by the d-core arcs"

def is_primitive_mult14(S):
    return reduce(gcd, S) == 1 and any(v % 14 == 0 for v in S) and len(set(S)) == 13

def witness_shell(t):
    return t.denominator

# =================== PART 0: primitivity ===================
print("=" * 72)
print("PART 0 - primitivity is essential in C'(14) as stated")
S0 = [2 * k for k in range(1, 14)]
M0, t0 = M_exact(S0)
print(f"S = 2*{{1..13}} = {S0}")
print(f"  contains 14: {14 in S0};  gcd = 2;  M(S) = {M0} (= 1/14: {M0 == F(1,14)}) at t = {t0}")
print("  => 'multiple of 14 => loose' is FALSE without gcd(S)=1; C' must read PRIMITIVE.")

# =================== PART 1: degeneration lemma ===================
print("=" * 72)
print("PART 1 - the 3-core degeneration lemma (exhaustive over residues)")
viol = 0
for v in range(1, 27):           # residues mod 27 with 3|v, 27 nmid v
    if v % 3 != 0: continue
    for a in range(1, 27):
        if gcd(a, 27) != 1: continue
        if normF(F(v * a, 27)) < F(1, 9): viol += 1
print(f"  ||v*a/27|| >= 1/9 for all 3|v, 27 nmid v, all units a: violations = {viol}")
band_14 = [j for j in range(0, 5) if F(j, 9) <= F(1, 14)]
band_5  = [j for j in range(0, 5) if F(j, 9) <= F(1, 5)]
print(f"  descended shell-9 forbidden band at level 1/14: +-{band_14}  (divisor condition only)")
print(f"  C'(5) shell-9 forbidden band at level 1/5:      +-{band_5}   (genuine +-1 structure)")
print("  => the descent DEGENERATES at 1/14: the 3-core head is the divisor tower {3,9,27},")
print("     NOT C'(5). Bands agree only for 1/n in [1/9, 2/9), i.e. n in 5..9. n=14 is past it.")

# =================== PART 2: non-necessity ===================
print("=" * 72)
print("PART 2 - C'(5)-on-the-core is NOT necessary: core tight as n=5, config loose")
core = [3, 6, 9, 12]                      # /3 = {1,2,3,4} = the TIGHT n=5 AP
Mc, tc = M_exact([1, 2, 3, 4])
S2 = core + [1, 2, 4, 5, 7, 8, 10, 11, 28]   # 28 = the multiple of 14; units leave +-13 free
assert is_primitive_mult14(S2), S2
M2, t2 = M_exact(S2)
okB, tB = head_B(S2)
print(f"  core/3 = [1,2,3,4]: M = {Mc} (tight as n=5 problem: {Mc == F(1,5)})")
print(f"  S = {sorted(S2)}")
print(f"  head B (3-tower): {okB}, twist t = {tB}")
if okB:
    mn = min(normF(v * tB) for v in S2)
    print(f"    certified: min_v ||v t|| = {mn} > 1/14: {mn > TARGET}")
print(f"  exact M(S) = {M2} > 1/14: {M2 > TARGET}  (witness t* = {t2}, shell {witness_shell(t2)})")
print("  => the 3-core can be C'(5)-TIGHT while C'(14) looseness holds: not necessary.")

# =================== PART 3: non-sufficiency hunt ===================
print("=" * 72)
print("PART 3 - hunting primitive mult-of-14 configs failing BOTH heads")
random.seed(14)

def test(S, label, verbose=True):
    S = sorted(S)
    if not is_primitive_mult14(S): return None
    okA, infoA = head_A(S)
    okB, infoB = head_B(S)
    if okA or okB: return ("covered", okA, okB)
    M, t = M_exact(S)
    if verbose:
        print(f"  JOINT FAILURE: {S}")
        print(f"    head A: {infoA}")
        print(f"    head B: {infoB}")
        print(f"    exact M = {M} (> 1/14: {M > TARGET}); witness t* = {t}, shell {t.denominator}")
    return ("joint-fail", M, t)

# family 1: 7*{1..12} + r  (rich 7-core, single stranger r)
fails = []
for r in range(1, 1200):
    if r % 7 == 0: continue
    S = [7 * k for k in range(1, 13)] + [r]
    if len(set(S)) != 13: continue
    res = test(S, f"7*{{1..12}}+{r}", verbose=False)
    if res and res[0] == "joint-fail":
        fails.append(("fam1", sorted(S), res[1], res[2]))
print(f"family 1 (7*{{1..12}} + r, r<1200): joint failures = {len(fails)}")
for tag, S, M, t in fails[:6]:
    print(f"  {S}  M={M}  t*={t} (shell {t.denominator})  loose: {M > TARGET}")

# family 2: kill head B with 54=2*27, rich 14-core
n_f2 = 0
for r in range(1, 1200):
    if r % 7 == 0 or r == 54: continue
    S = [14 * k for k in range(1, 7)] + [7, 21, 35, 49, 77, 54] + [r]
    if len(set(S)) != 13: continue
    res = test(S, f"fam2+{r}", verbose=False)
    if res and res[0] == "joint-fail":
        n_f2 += 1
        fails.append(("fam2", sorted(S), res[1], res[2]))
print(f"family 2 (14*{{1..6}} u odd-7-mults u {{54}} u r): joint failures = {n_f2}")
for tag, S, M, t in [f for f in fails if f[0] == 'fam2'][:6]:
    print(f"  {S}  M={M}  t*={t} (shell {t.denominator})  loose: {M > TARGET}")

# family 3: random structured: dense 7-divisible cores of varying size + tuned units
n_f3 = 0; tried = 0
for trial in range(4000):
    c = random.randint(6, 12)                       # size of the 7-core
    ks = random.sample(range(1, 19), c)
    if not any((7 * k) % 14 == 0 for k in ks):      # ensure a multiple of 14
        ks[0] = 2 * random.randint(1, 9)
    corev = [7 * k for k in ks]
    restv = random.sample([v for v in range(1, 120) if v % 7], 13 - c)
    S = corev + restv
    if len(set(S)) != 13 or reduce(gcd, S) != 1: continue
    tried += 1
    res = test(S, "fam3", verbose=False)
    if res and res[0] == "joint-fail":
        n_f3 += 1
        if n_f3 <= 8: fails.append(("fam3", sorted(S), res[1], res[2]))
print(f"family 3 (random, {tried} primitive tries): joint failures = {n_f3}")
for tag, S, M, t in [f for f in fails if f[0] == 'fam3'][:8]:
    print(f"  {S}  M={M}  t*={t} (shell {t.denominator})  loose: {M > TARGET}")

# verdict + witness landscape
print("=" * 72)
jf = [f for f in fails]
print(f"TOTAL joint failures found: {len(jf)}")
if jf:
    all_loose = all(M > TARGET for _, _, M, _ in jf)
    shells = sorted({t.denominator for _, _, _, t in jf})
    print(f"  all loose (C' holds on them): {all_loose}")
    print(f"  witness shells (denominators of t*): {shells}")
    print("  => the union [C'(5)-3-core head] u [THM-421 fiber head] is NOT sufficient:")
    print("     these configs are loose only via witnesses OUTSIDE both heads.")
else:
    print("  none found in these families - heads may be jointly sufficient on them;")
    print("  the exactness question stays open on this evidence.")
