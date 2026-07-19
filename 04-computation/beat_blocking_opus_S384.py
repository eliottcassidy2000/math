# opus-2026-07-17-S384 -- HYP-7690: DOES THE BEAT CONSTRAINT BREAK THE BLOCKING?
#
# THM-1110: 13 free speeds can block ANY single modulus q, because each kills
# k_q = #{w in W_q\{0} : gcd(w,q)=1} numerators and 13*k_q >= phi(q) always.
# THM-1170: the optimum sits at a BEAT frequency q = v_i +- v_j.
#
# THE CLAIMED BITE.  W_q is symmetric (r in W_q <=> -r in W_q).  At q = v_i + v_j
# we have v_i = -v_j (mod q), so for every p,  v_j p = -v_i p (mod q), and
#     v_i p in W_q  <=>  v_j p in W_q.
# The two speeds kill IDENTICAL numerator sets -- they are redundant, so 13
# speeds act as at most 12 distinct blockers at their own beat frequency.
# THM-1110's counting theorem fires when (#distinct blockers)*k_q < phi(q), and
# max over q of phi(q)/k_q is 12 (at q = 90).  So the drop 13 -> 12 lands EXACTLY
# on the boundary.  This run measures what actually happens.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
from collections import Counter

def W(q): return set([0] + [r % q for j in range(1, (q-1)//14 + 1) for r in (j, -j)])
def killset(v, q, Wq):
    return frozenset(p for p in range(1, q) if gcd(p, q) == 1 and (v*p) % q in Wq)
def phi(n): return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)
def k_q(q, Wq): return sum(1 for w in Wq if w != 0 and gcd(w, q) == 1)

print("(1) THE REDUNDANCY CLAIM: at q = v_i + v_j, do v_i and v_j kill the same set?")
random.seed(384)
bad = 0; n = 0
for _ in range(200):
    a, b = random.sample(range(1, 200), 2)
    q = a + b
    if q <= 14: continue
    Wq = W(q)
    n += 1
    if killset(a, q, Wq) != killset(b, q, Wq): bad += 1
print(f"    {n} beat pairs tested; kill-sets differ in {bad} cases")
print("    (0 confirms the pairing is exact -- W_q symmetric and v_j = -v_i mod q)")

print()
print("(2) DOES THE SAME HOLD AT DIFFERENCE BEATS q = |v_i - v_j|?")
bad2 = 0; n2 = 0
for _ in range(200):
    a, b = random.sample(range(1, 200), 2)
    q = abs(a - b)
    if q <= 14: continue
    Wq = W(q)
    n2 += 1
    if killset(a, q, Wq) != killset(b, q, Wq): bad2 += 1
print(f"    {n2} difference beats tested; kill-sets differ in {bad2} cases")
print("    (at q = |v_i - v_j| we have v_i = v_j mod q -- identical, not antipodal)")

print()
print("(3) THE DECISIVE COUNT: over a family's OWN beat frequencies, how close does")
print("    (#distinct blockers) * k_q  come to  phi(q)?   Theorem fires if <.")
fires = 0; fams = 0; best_margins = []
for _ in range(30):
    V = sorted(random.sample(range(1, 90), 13))
    fams += 1
    beats = set()
    for a, b in combinations(V, 2):
        beats.add(a+b)
        if abs(a-b) > 14: beats.add(abs(a-b))
    for v in V: beats.add(2*v)
    best = None
    for q in sorted(beats):
        if q <= 14: continue
        if any(v % q == 0 for v in V): continue      # cheaply blocked, skip
        Wq = W(q); kq = k_q(q, Wq)
        if kq == 0: continue
        distinct = len({killset(v, q, Wq) for v in V})
        margin = phi(q) - distinct * kq              # > 0 means theorem FIRES
        if best is None or margin > best[0]: best = (margin, q, distinct, phi(q), kq)
    if best:
        best_margins.append(best)
        if best[0] > 0: fires += 1
print(f"    families where the counting theorem FIRES at some own-beat q: {fires}/{fams}")
if best_margins:
    bm = max(best_margins, key=lambda t: t[0])
    print(f"    best case seen: q={bm[1]}  phi={bm[3]}  k_q={bm[4]}  distinct blockers={bm[2]}"
          f"  margin={bm[0]}")
    md = sorted(t[0] for t in best_margins)[len(best_margins)//2]
    print(f"    median best margin over families: {md}   (needs > 0 to certify)")
