"""
lrc14_gap_vs_longestAP_opus_S181.py  (opus-2026-07-09-S181)

DEEP-DIVE TEST (Freiman neighborhood): the good-period dichotomy keys on longest-AP, but Freiman theory
says the true 'near-tight' parameter is the DOUBLING constant K=|S+S|/|S| (equivalently additive energy).
These DIVERGE on multi-dimensional GENERALIZED APs (GAPs): a 2-D GAP {a+b*P} has SMALL doubling (very
structured) but SHORT 1-D longest-AP (only the side length). QUESTION: is a 2-D GAP near-tight (small
lonely measure L) DESPITE a short longest-AP? If yes, longest-AP is the WRONG dichotomy parameter and
doubling |S+S| is right -- and GAP-structured sets are a blind spot for a longest-AP-keyed dichotomy.
"""
import numpy as np
from collections import Counter

NG = 120011
TAU = (np.arange(NG) + 0.5) / NG
H = 1.0 / 14
MAIN = (6.0 / 7.0) ** 13

def E(S):
    r = Counter(a + b for a in S for b in S); return sum(c * c for c in r.values())
def sumset(S): return len(set(a + b for a in S for b in S))
def longest_ap(S):
    Sset = set(S); S = sorted(S); best = 1
    for i, a in enumerate(S):
        for b in S[i+1:]:
            d = b - a
            if a - d in Sset: continue
            L = 2; x = b + d
            while x in Sset: L += 1; x += d
            best = max(best, L)
    return best
def lonely(S):
    M = np.zeros(NG)
    for v in S:
        d = np.abs(((v * TAU + 0.5) % 1.0) - 0.5); M += (d <= H)
    return float((M == 0).mean())

def gap2d(A, B, P):
    # shift by +1 so all speeds are NONZERO (speed 0 would trivially cover the origin => fake L=0)
    S = sorted(set(1 + a + b*P for a in range(A) for b in range(B)))
    return S[:13] if len(S) >= 13 else None

def Mval(S, NG2=600011):
    """M(S) = max_tau min_i ||v_i tau||, fine 1-D search (the actual loneliness)."""
    t = (np.arange(1, NG2) ) / NG2
    mn = np.full(NG2 - 1, 1.0)
    for v in S:
        d = np.abs(((v * t + 0.5) % 1.0) - 0.5)
        mn = np.minimum(mn, d)
    return float(mn.max())

print("="*104)
print(f"GAP vs longest-AP: is a small-doubling 2-D GAP near-tight (small L) despite a SHORT longest-AP?")
print(f"(6/7)^13={MAIN:.4f}  1/14=0.0714   [near-tight = small L; loose = L~0.13]")
print("="*104)
print(f"  {'set':>34} {'|S+S|':>6} {'K=|S+S|/13':>10} {'E(S)':>6} {'longestAP':>9} {'L lonely':>9} {'M=maxmin':>9} {'14*M':>6}")
tests = []
tests.append(("AP {1..13} (1-D, tight)", list(range(1,14))))
tests.append(("Sidon (dissoc)", [1,2,5,11,22,33,45,60,78,95,110,130,150]))
# 2-D GAPs: small doubling, short 1-D AP
tests.append(("2-D GAP 4x4 P=5", gap2d(4,4,5)))
tests.append(("2-D GAP 4x4 P=6", gap2d(4,4,6)))
tests.append(("2-D GAP 3x5 P=7", gap2d(3,5,7)))
tests.append(("2-D GAP 5x3 P=20", gap2d(5,3,20)))
tests.append(("2-D GAP 4x4 P=20", gap2d(4,4,20)))
tests.append(("2-D GAP 7x2 P=15", gap2d(7,2,15)))
tests.append(("2-D GAP 3x5 P=40 (wide)", gap2d(3,5,40)))
for name, S in tests:
    if S is None or len(set(S)) != 13: 
        print(f"  {name:>34}  (skip, size {len(set(S)) if S else 0})"); continue
    S = sorted(set(S))
    M = Mval(S)
    print(f"  {name:>34} {sumset(S):>6} {sumset(S)/13:>10.3f} {E(S):>6} {longest_ap(S):>9} {lonely(S):>9.4f} {M:>9.5f} {14*M:>6.3f}")
print("="*104)
print("READING: if a 2-D GAP row shows SHORT longestAP (<=5) but SMALL K & SMALL L (near-tight), then the")
print("dichotomy must key on the DOUBLING |S+S|, not longest-AP -- GAPs are a longest-AP blind spot.")
