from fractions import Fraction as F
from math import gcd
from itertools import combinations
from lrc import measure_B, measure_intersection, arc_endpoints

# ---------- PART 1: FIRST MOMENT ----------
print("="*70)
print("PART 1: FIRST MOMENT  E[N] = sum_i measure(B_i) = 2m/n = 2(n-1)/n")
print("="*70)

speed_sets = [
    [1,2,3,4],        # classic, m=4 n=5
    [1,3,4,7],
    [2,3,5,7],
    [1,2,3,4,5],      # m=5 n=6
    [1,3,4,5,9],
    [3,5,7,11,13],
]
for speeds in speed_sets:
    m = len(speeds); n = m+1
    EN = sum(measure_B(v,n) for v in speeds)
    print(f"  speeds={speeds} n={n}: E[N]={EN}={float(EN):.4f}  (2m/n={F(2*m,n)})")
    assert EN == F(2*m, n)
print("  All E[N] == 2m/n confirmed.  -> ~2 for large n.")

# ---------- PART 2: SECOND MOMENT / PAIRWISE OVERLAP ----------
print()
print("="*70)
print("PART 2: PAIRWISE OVERLAP  measure(B_i ∩ B_j), correlation ratio")
print("="*70)
print("  correlation ratio  r = measure(B_i∩B_j) / (2/n)^2  = measure * n^2/4")
print("  r=1 iff independent; r>1 over-correlated (resonant); r<1 under.")
print()

def overlap(vi, vj, n):
    return measure_intersection([vi,vj], n)

def ratio(vi, vj, n):
    o = overlap(vi,vj,n)
    return o * F(n*n,4), o

# Scan many pairs at a fixed n to see the pattern.
def scan(n, vmax):
    rows = []
    for vi in range(1, vmax+1):
        for vj in range(vi+1, vmax+1):
            r, o = ratio(vi, vj, n)
            g = gcd(vi,vj)
            rows.append((float(r), vi, vj, g, o))
    return rows

for n in [5, 7]:
    print(f"--- n={n}, scanning pairs (vi,vj), vmax=20 ---")
    rows = scan(n, 20)
    rows.sort()
    print(f"  total pairs: {len(rows)}; ratio min={rows[0][0]:.4f} max={rows[-1][0]:.4f}")
    over = [x for x in rows if x[0] > 1.0000001]
    under = [x for x in rows if x[0] < 0.9999999]
    eq = [x for x in rows if abs(x[0]-1)<1e-9]
    print(f"  #over-correlated(r>1): {len(over)}   #under(r<1): {len(under)}   #exactly 1: {len(eq)}")
    print("  Most over-correlated (top 8):")
    for r,vi,vj,g,o in rows[-8:][::-1]:
        print(f"     v=({vi},{vj}) gcd={g} v_i+v_j={vi+vj} |v_i-v_j|={abs(vi-vj)}  ratio={r:.4f} meas={o}")
    print("  Most under-correlated (bottom 5):")
    for r,vi,vj,g,o in rows[:5]:
        print(f"     v=({vi},{vj}) gcd={g} ratio={r:.4f} meas={o}")
    print()
