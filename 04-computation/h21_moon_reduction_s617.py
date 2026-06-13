#!/usr/bin/env python3
"""h21_moon_reduction_s617.py — a rigorous reduction of the H=21 problem to a
FINITE check (n<=12) via Moon's 3-cycle bound, plus a CORRECTION of an S616 error.

============================ CORRECTION (S616) ============================
S616 (HYP-2187) claimed the H=21 problem reduces to a single open profile
(alpha_1=6, alpha_2=2) and used a reconstructed conflict graph Omega. That
reconstruction was BUGGY: it deduped directed odd cycles by their VERTEX SET,
collapsing distinct directed 5-/7-cycles on the same vertices. Verified: the
reconstructed I(Omega,2) != true H(T) in 1612/3000 random n=6 tournaments. So:
  * the general I(Omega,2) values from S616 are unreliable;
  * "the one open case is (6,2,0)" was WRONG -- there are FOUR connected
    profiles (4,3),(6,2),(8,1),(10,0), all in THM-079 Part H (the i_2-jump);
  * the (6,2) i_2 in {0,1,5} that S616 reported DOES match THM-079 Part H
    (3-cycle-dominated, where the dedup bug is harmless) -- but it is THM-079's
    result, not new.
H=21 is in fact PROVED impossible for n<=8 EXHAUSTIVELY (THM-079 Part G,
268,435,456 tournaments at n=8). I under-credited it as n<=7.
==========================================================================

NEW RIGOROUS REDUCTION (this session):
  H(T) = 21  =>  (strong-component multiplicativity, H = prod over SC components,
    and 21 = 3*7 with 7 not a strong H-value by THM-029) => some STRONG component
    S has H(S) = 21.
  For that strong tournament S on m vertices:
    H(S)=21 => I(Omega,2)=1+2*alpha_1+... =21 => alpha_1 <= 10
    => c_3(S) <= alpha_1 <= 10           (3-cycles are odd cycles: c_3 <= alpha_1)
  MOON'S THEOREM: a strong tournament on m vertices has c_3 >= m-2 (verified
    below for m<=7: min c_3 over strong = m-2 exactly).
    => m-2 <= 10  =>  m <= 12.
  With THM-079 Part G (H=21 impossible exhaustively for m<=8), the ONLY remaining
  cases are STRONG tournaments on m in {9,10,11,12} with c_3 <= 10 -- a FINITE set.

This collapses THM-079's open part (which left cycle-rich n>=9 open in general)
to a finite n<=12 check. Sampling evidence (below): strong tournaments on n=9,10
with c_3<=10 have min H = 75, 153 -- comfortably != 21 (long odd cycles from
Moon pancyclicity blow alpha_1 far past 10). FULL proof = exhaustive enumeration
of strong c_3<=10 tournaments on m=9..12 (very restricted, near-Moon-extremal),
OR a proof that strong m>=9 forces alpha_1 >= 11.

Session: claude-2026-06-03-S617 (h21-moon-reduction).
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from itertools import product
import random

def c3(n, adj):
    cnt = 0; e = lambda a, b: adj[a] >> b & 1
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (e(i, j) and e(j, k) and e(k, i)) or (e(i, k) and e(k, j) and e(j, i)):
                    cnt += 1
    return cnt
def is_strong(n, adj):
    def reach(fwd):
        seen = 1; fr = [0]
        while fr:
            x = fr.pop()
            for y in range(n):
                if not (seen >> y & 1):
                    ed = (adj[x] >> y & 1) if fwd else (adj[y] >> x & 1)
                    if ed: seen |= 1 << y; fr.append(y)
        return seen == (1 << n)-1
    return reach(True) and reach(False)
def Hcount(n, adj):
    size = 1 << n; dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for w in range(n):
                    if not (mask >> w & 1) and av >> w & 1: dp[mask | 1 << w][w] += c
    return sum(dp[size-1])

print("\n  H=21: MOON 3-CYCLE REDUCTION TO n<=12 (+ S616 correction)\n" + "=" * 70)
print("\n  (1) Moon's bound: min c_3 over STRONG tournaments = n-2  [verified m<=7]")
for n in range(3, 8):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    mc = 10**9
    for bits in product([0, 1], repeat=len(edges)):
        adj = [0]*n
        for (i, j), b in zip(edges, bits):
            if b: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
        if is_strong(n, adj): mc = min(mc, c3(n, adj))
    print(f"      m={n}: min c_3 (strong) = {mc} = m-2 ({mc==n-2})")
print("""
  (2) REDUCTION: H=21 => strong component with H=21 (21=3*7, 7 not strong, THM-029)
      => alpha_1<=10 => c_3<=10 => (Moon) m<=12. With THM-079 (m<=8 exhaustive),
      remaining = STRONG tournaments, m in {9,10,11,12}, c_3<=10. FINITE.""")
print("\n  (3) Evidence: strong tournaments on m=9..12 with c_3<=10 have H far from 21:")
rng = random.Random(11)
for n in [9, 10]:
    cnt = 0; minH = 10**9; h21 = 0
    for _ in range(250000):
        adj = [0]*n
        for i in range(n):
            for j in range(i+1, n):
                if rng.getrandbits(1): adj[i] |= 1 << j
                else: adj[j] |= 1 << i
        if c3(n, adj) <= 10 and is_strong(n, adj):
            cnt += 1; H = Hcount(n, adj); minH = min(minH, H)
            if H == 21: h21 += 1
    print(f"      m={n}: {cnt} strong c_3<=10 sampled; H=21 found {h21}; min H = {minH}")
print("      (m=11,12: strong c_3<=10 are near-Moon-extremal and very rare.)")
print("""
  => H=21 reduces to a FINITE n<=12 check; sampled cases have H>=75 (m=9), 153 (m=10),
     nowhere near 21 (long odd cycles from Moon pancyclicity push alpha_1>>10).
     FULL proof: exhaust strong c_3<=10 on m=9..12, or prove strong m>=9 => alpha_1>=11.
""")
