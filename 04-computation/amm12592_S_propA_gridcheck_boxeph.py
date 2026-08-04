"""Lane E2 -- machine check of Proposition A's arithmetic (static
uniform-envelope closure) on a hostile integer grid.

Hypotheses: 3 <= m, 6(m+2) <= d, 6 m B <= d, B >= 1; uniform interior bound
Abar = 2B (the proved cell-1 transient ceiling).  Checked implications
(exact Fractions):

  (P1) delta=1 interior pump, any 2 <= t <= m:
         2*Abar*t/d + Abar*t*(t-1)/(d*(d-t+1)) <= 4/5        [decay >= 1/5]
  (P2) delta=0 interior pump, any 2 <= t <= m:
         Abar*t/(d-t) <= 4/5                                  [decay >= 1/5]
  (P3) freeze, delta=1 mate:  2*Abar + Abar*m/(d-m) <= d/(m+1)
  (P4) freeze, delta=0:       Abar <= (d-1-m)/(m+1)
  (P5) freeze, delta=1 extreme: Abar <= d*(d-m-1)/((m+1)*(m+2))
  (P6) cell-1 ceiling: B*(1 + (B+2)/(4*d)) <= 2*B   (i.e. B + 2 <= 4d)
  (P7) deadlines vs budget (>= d - 4, parallel not sequential):
         max( 10*B,                       [cells 2..m extinct]
              B//2 + 2*isqrt(B*d) + 2,    [cell 1: transient + drain integral]
              (d + B)//2 + 1 )            [debt clock, |j0| <= d + B]
         <=  d - 4

Grid: d over 7 decades incl. the real feed-end degrees at R = 2^7..2^20,
m over {3, 5, 12, d//12, d//6 - 2}, B = d // (6m) (the extremal allowed B)
and B = 1 (degenerate end).
"""
import sys
from fractions import Fraction as F
from math import isqrt

fails = []
def chk(name, cond, ctx):
    if not cond:
        fails.append((name, ctx))

DS = [60, 99, 196, 393, 786, 1613, 3226, 6451, 12903, 25806, 103224,
      825792, 10**6]
for d in DS:
    ms = sorted({3, 5, 12, max(3, d // 12), max(3, d // 6 - 2)})
    for m in ms:
        if not (3 <= m and 6 * (m + 2) <= d):
            continue
        for B in {1, d // (6 * m)}:
            if B < 1 or 6 * m * B > d:
                continue
            Ab = 2 * B
            for t in range(2, m + 1):
                chk("P1", F(2 * Ab * t, d) + F(Ab * t * (t - 1),
                    d * (d - t + 1)) <= F(4, 5), (d, m, B, t))
                chk("P2", F(Ab * t, d - t) <= F(4, 5), (d, m, B, t))
            chk("P3", F(2 * Ab) + F(Ab * m, d - m) <= F(d, m + 1), (d, m, B))
            chk("P4", Ab <= F(d - 1 - m, m + 1), (d, m, B))
            chk("P5", Ab <= F(d * (d - m - 1), (m + 1) * (m + 2)), (d, m, B))
            chk("P6", B + 2 <= 4 * d, (d, m, B))
            chk("P7", max(10 * B, B // 2 + 2 * isqrt(B * d) + 2,
                          (d + B) // 2 + 1) <= d - 4, (d, m, B))

if fails:
    print("FAIL %d cases; first 10:" % len(fails))
    for f in fails[:10]:
        print(" ", f)
    sys.exit(1)
print("Proposition A grid check: ALL PASS (%d degree values, extremal B)"
      % len(DS))
