from fractions import Fraction as F
from math import exp
from itertools import combinations
from lrc import measure_B, measure_intersection, lonely_measure

# ---------- PART 3: LONELY MEASURE mu vs independence baseline ----------
print("="*70)
print("PART 3: LONELY MEASURE mu = measure{N=0}  vs  (1-2/n)^m")
print("="*70)

def baseline(m,n):
    return (F(1) - F(2,n))**m

def moments(speeds,n):
    """E[N], E[N(N-1)], Var(N)."""
    m=len(speeds)
    EN=sum(measure_B(v,n) for v in speeds)
    # E[N(N-1)] = sum_{i!=j} measure(B_i∩B_j) = 2 * sum_{i<j} overlap
    S2=F(0)
    for i,j in combinations(range(m),2):
        S2 += measure_intersection([speeds[i],speeds[j]],n)
    ENN1=2*S2
    EN2=ENN1+EN
    Var=EN2-EN*EN
    return EN,ENN1,Var

# Examine known/structured speed sets.
print("\nm=4 (n=5):  baseline (1-2/5)^4 = (3/5)^4 =", baseline(4,5), "=", float(baseline(4,5)))
print(" e^-2 =", round(exp(-2),5))
test4=[
 [1,2,3,4],[1,2,3,5],[1,3,5,7],[1,2,5,7],[1,4,5,6],[2,3,4,5],
 [1,2,4,8],[1,2,4,7],[1,5,7,11],[3,4,5,7],[1,2,3,7],[1,6,7,11],
]
print(f"\n{'speeds':22}{'mu':>10}{'mu(float)':>11}{'baseline':>10}{'mu/base':>9}{'Var':>9}")
res4=[]
for s in test4:
    n=5; mu=lonely_measure(s,n); EN,ENN1,Var=moments(s,n)
    b=baseline(4,n); ratio=float(mu/b) if b else 0
    res4.append((ratio,s,mu,Var))
    print(f"{str(s):22}{str(mu):>10}{float(mu):>11.5f}{float(b):>10.5f}{ratio:>9.4f}{float(Var):>9.4f}")

# brute-force search smallest mu/baseline for m=4, n=5 over speed sets up to 12
print("\nSearch m=4,n=5 over all 4-subsets of {1..14}: smallest mu/baseline (hardest sets):")
n=5; b=baseline(4,n); records=[]
for s in combinations(range(1,15),4):
    mu=lonely_measure(list(s),n)
    records.append((float(mu/b),float(mu),s))
records.sort()
zero=[r for r in records if r[1]==0]
print(f"  total sets={len(records)}; sets with mu=0 (NO lonely time, would refute LRC): {len(zero)}")
print("  smallest mu/baseline:")
for r,mu,s in records[:12]:
    print(f"     {s}  mu={mu:.5f}  mu/base={r:.4f}")
print("  largest mu/baseline:")
for r,mu,s in records[-4:]:
    print(f"     {s}  mu={mu:.5f}  mu/base={r:.4f}")
