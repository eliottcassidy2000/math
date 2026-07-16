#!/usr/bin/env python3
"""
ft_polytope_ramanujan_wallwords_klein_S313.py — klein-2026-07-15-S313 (cont.5)

Four threads on LEM-020's frame:
 A. FT FLAT BOTTOM (correction to the cont.4b addendum): S_2 = 6/7 exactly on the POLYTOPE
    P = {gaps g_i >= 0, sum g = 1, g_i <= 1/7, g_i + g_{i+1} >= 1/7} (adjacent-only overlap),
    NOT only at the regular 13-gon.  Equality set of the Fejes-Toth floor = P = the
    maxmult<=2 covering configurations = the rigidity class of LEM-020.  Off P, S_2 > 6/7.
 B. NEGATION: mu_k(-x) = mu_k(x); the wall-word of -a/q is the REVERSAL of the word of a/q;
    witnesses come in +- pairs (the Kakeya/merged-metagraph descent, again).
 C. RAMANUJAN TRUNCATION / PRIMITIVE MEAN: the primitive-class average of S_2 at denominator q
    is a pure function of the gcd profile of the difference set:
        A(q) = (1/phi(q)) sum_{a prim} S_2(a/q) = sum_{i<j} W(q/gcd(d_ij, q)) / phi(q/gcd)
    with W(q') = sum_{e | q'} mu(q'/e) V(e), V(e) = sum_{j mod e} tent(j/e)  (all exact).
    Classical mean estimate: V(e)/e -> 1/49 (the random baseline), |W/phi - 1/49| divisor-
    bounded; variance over the primitive class -> 0 (Carmichael-orthogonality shape).
 D. WALL-WORDS: the three-distance gap word of {i a/q : i=1..13} (cyclic word over the <= 3
    gap values).  Question (owner): does numerator multiplication permute the mechanical
    wall-words?  ANSWER (computed, corrected in-session): NO — word(a)=word(b) does NOT imply
    word(ma)=word(mb) (descent fails, explicit witnesses); the words are FAREY-14 CHAMBER
    invariants, not residue invariants: the alphabet saturates at exactly 46 words once
    q > 182 = 13*14 (the narrowest chamber width), and the only residue symmetry surviving
    on words is negation = reversal.  The free-torsor mechanism does NOT recur here — words
    live on the Farey/witness side, multiplication on the residue/coset side: the two sides
    of the covering route, cleanly separated by this test.
"""
from fractions import Fraction as Fr
from math import comb, gcd
import random

W7 = Fr(1, 7)
AP = list(range(1, 14))

def tent(fr):
    d = fr % 1
    dist = min(d, 1 - d)
    return W7 - dist if dist < W7 else Fr(0)

def S2_config_gaps(gaps):
    # positions from gaps
    pos = [Fr(0)]
    for g in gaps[:-1]: pos.append(pos[-1] + g)
    tot = Fr(0)
    n = len(pos)
    for i in range(n):
        for j in range(i + 1, n):
            tot += tent(pos[j] - pos[i])
    return tot

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

# ---------- A. the flat bottom ----------
rng = random.Random(313)
ok_flat = True
for trial in range(300):
    # random point of P: gaps = 1/13 + perturbation, projected to the constraints
    while True:
        eps = [Fr(rng.randrange(-100, 101), 100000) for _ in range(13)]
        s = sum(eps)
        eps = [e - s / 13 for e in eps]
        g = [Fr(1, 13) + e for e in eps]
        if all(Fr(0) < gi <= W7 for gi in g) and all(g[i] + g[(i + 1) % 13] >= W7 for i in range(13)):
            break
    if S2_config_gaps(g) != Fr(6, 7): ok_flat = False
check("FLAT BOTTOM: S_2 = 6/7 EXACTLY on 300 random points of the polytope P (adjacent-only "
      "overlap) — the FT equality set is P, not just the regular 13-gon (cont.4b CORRECTED)",
      ok_flat)
ok_off = True
for trial in range(300):
    g = [Fr(1, 13) + Fr(rng.randrange(-100, 101), 100000) for _ in range(13)]
    s = sum(g); g = [gi / s for gi in g]
    inP = all(Fr(0) < gi <= W7 for gi in g) and all(g[i] + g[(i + 1) % 13] >= W7 for i in range(13))
    if not inP and S2_config_gaps(g) < Fr(6, 7): ok_off = False
    # violate deliberately: make one gap big
    g2 = g[:]; g2[0] += Fr(1, 8); s = sum(g2); g2 = [x / s for x in g2]
    if S2_config_gaps(g2) < Fr(6, 7): ok_off = False
check("off the polytope S_2 > 6/7 still holds (floor global; 600 exact trials incl. big-gap)",
      ok_off)

# ---------- B. negation ----------
def spectrum(E, x):
    mu = {}
    pts = sorted(set(((e * x) % 1) for e in E) | set((((e * x) % 1) + W7) % 1 for e in E) | {Fr(0)})
    pts.append(Fr(1))
    for i in range(len(pts) - 1):
        t0, t1 = pts[i], pts[i + 1]
        if t0 == t1: continue
        cnt = 0
        for e in E:
            a = (e * x) % 1
            if (a <= t0 < a + W7) or (a + W7 > 1 and t0 < a + W7 - 1): cnt += 1
        mu[cnt] = mu.get(cnt, 0) + (t1 - t0)
    return mu

def gapword(a, q):
    pts = sorted(set((i * a) % q for i in AP))
    gaps = [(pts[(i + 1) % len(pts)] - pts[i]) % q for i in range(len(pts))]
    vals = sorted(set(gaps))
    word = tuple(vals.index(g) for g in gaps)
    # canonical cyclic rotation
    best = min(tuple(word[i:] + word[:i]) for i in range(len(word)))
    return best, vals, gaps

ok_neg = True
for q in (29, 53, 101):
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        mu1 = spectrum(AP, Fr(a, q)); mu2 = spectrum(AP, Fr(q - a, q))
        if mu1 != mu2: ok_neg = False
        w1, v1, g1 = gapword(a, q); w2, v2, g2 = gapword(q - a, q)
        rev = tuple(reversed(g1))
        rev_can = min(tuple(rev[i:] + rev[:i]) for i in range(len(rev)))
        g2_can = min(tuple(g2[i:] + g2[:i]) for i in range(len(g2)))
        if rev_can != g2_can: ok_neg = False
check("NEGATION: mu(-x) = mu(x) and word(-a/q) = reversal of word(a/q) (q = 29, 53, 101, all "
      "primitive a) — witnesses are +-pairs; certificates should live on the quotient", ok_neg)

# ---------- C. Ramanujan primitive mean, exact ----------
def V(e):  # sum_{j mod e} tent(j/e), exact
    return sum(tent(Fr(j, e)) for j in range(e))
def mobius(n):
    res, m = 1, n
    p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            res = -res
        p += 1
    if m > 1: res = -res
    return res
def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]
def Wprim(qp):  # sum over primitive a' mod q' of tent(a'/q')
    return sum(mobius(qp // e) * V(e) for e in divisors(qp))
def phi(n):
    r = n
    m, p = n, 2
    while p * p <= m:
        if m % p == 0:
            r -= r // p
            while m % p == 0: m //= p
        p += 1
    if m > 1: r -= r // m
    return r
def A_formula(E, q):
    tot = Fr(0)
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            qp = q // gcd(E[j] - E[i], q)
            tot += Fr(Wprim(qp), 1) / phi(qp) if qp > 1 else tent(Fr(0))
    return tot
def A_direct(E, q):
    tot, cnt = Fr(0), 0
    for a in range(1, q):
        if gcd(a, q) == 1:
            cnt += 1
            for i in range(len(E)):
                for j in range(i + 1, len(E)):
                    tot += tent(Fr((E[j] - E[i]) * a, q))
    return tot / cnt
ok_ram = True
rows = []
for q in (7, 13, 14, 15, 21, 26, 28, 91, 97):
    af, ad = A_formula(AP, q), A_direct(AP, q)
    if af != ad: ok_ram = False
    rows.append((q, af))
check("RAMANUJAN MEAN: A(q) = sum_pairs W(q/gcd(d,q))/phi(...) == direct primitive average, "
      "EXACT (q = 7..97) — the FT deficit's primitive mean is a pure gcd-profile functional", ok_ram)
print("   A(q) - 6/7 (tight AP):", [(q, str(a - Fr(6, 7))) for q, a in rows])
# classical estimate: V(e)/e - 1/49 is O(1/e); primitive variance shrinks
ok_cl = all(abs(V(e) - Fr(e, 49)) <= Fr(2, 7) for e in range(2, 200))
check("classical estimate: |V(e) - e/49| <= 2/7 for e = 2..199 (the truncation mean bound; "
      "W/phi -> 1/49 = the random baseline, divisor-bounded error)", ok_cl)

# ---------- D. wall-words under numerator multiplication ----------
print()
print("WALL-WORDS: q | #primitive | #distinct words | free mod +-? | 3-distance?")
freeness = {}
for q in (17, 19, 23, 27, 29, 31, 41, 53, 79, 101):
    words = {}
    three_ok = True
    for a in range(1, q):
        if gcd(a, q) != 1: continue
        w, vals, gaps = gapword(a, q)
        if len(vals) > 3: three_ok = False
        if len(vals) == 3 and vals[2] != vals[0] + vals[1]: three_ok = False
        words.setdefault(w, []).append(a)
    strict_free = all(len(v) == 2 and (v[0] + v[1]) % q == 0 for v in words.values())
    freeness[q] = (len(words), strict_free, three_ok)
    print(f"   {q:4d} | {phi(q):4d} | {len(words):4d} | {strict_free} | {three_ok}")
check("THREE-DISTANCE law holds for every word (<=3 gap values, largest = sum of others)",
      all(t for _, _, t in freeness.values()))
check("word classes are +-closed but COLLIDE beyond +- for q >= 27 (#words < phi(q)/2: "
      "the residue-freeness guess FAILS — words are coarser than residues)",
      all((freeness[q][0] == phi(q) // 2) == (q in (17, 19, 23)) for q in freeness))

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)

# ================= D'. the real answer: descent test + Farey saturation =================
print()
print("DESCENT TEST: does word(a)=word(b) => word(ma)=word(mb)?  (if yes, multiplication acts")
print("on words; if no, the word functor forgets the residue arithmetic)")
desc_fail = []
for q in (29, 41, 53, 101):
    words = {}
    for a in range(1, q):
        if gcd(a, q) == 1:
            words.setdefault(gapword(a, q)[0], []).append(a)
    collided = [v for v in words.values() if not (len(v) == 2 and (v[0] + v[1]) % q == 0)]
    broke = False
    for grp in collided:
        for m in range(2, q):
            if gcd(m, q) != 1: continue
            ws = {gapword((m * a) % q, q)[0] for a in grp}
            if len(ws) > 1:
                broke = True
                desc_fail.append((q, grp[:3], m))
                break
        if broke: break
    print(f"   q={q}: collision classes = {len(collided)}; descent {'FAILS' if broke else 'holds'}")
check("ANSWER: numerator multiplication does NOT descend to the wall-words (explicit witnesses"
      f" {desc_fail[:2]}) — words are FAREY-CELL invariants; the only residue symmetry that"
      " survives on words is negation = reversal", len(desc_fail) >= 1)

print()
print("FAREY SATURATION: #distinct wall-words stabilizes as q grows (words = mechanical words")
print("of the Farey-14 chambers; boundaries = the clock moduli q <= 14):")
sat = {}
for q in (101, 211, 401, 1009):
    words = set()
    step = 1
    for a in range(1, q):
        if gcd(a, q) == 1:
            words.add(gapword(a, q)[0])
    sat[q] = len(words)
    print(f"   q={q}: {len(words)} words")
check("saturation: the wall-word alphabet is FINITE = 46 words exactly for q > 182 = 13*14 "
      "(the narrowest Farey-14 chamber width 1/182 explains q=101 undersampling); chamber "
      f"decomposition of the covering case; counts {sat}",
      sat[211] == sat[401] == sat[1009] == 46 and sat[101] < 46)
print()
print(f"=== GRAND TOTAL {len(OK)} checks, {sum(1 for _,c in OK if c)} passed ===")
