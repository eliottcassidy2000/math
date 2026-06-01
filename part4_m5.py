from fractions import Fraction as F
from math import exp
from itertools import combinations
from lrc import measure_B, measure_intersection, lonely_measure

def baseline(m,n): return (F(1)-F(2,n))**m
def moments(speeds,n):
    m=len(speeds)
    EN=sum(measure_B(v,n) for v in speeds)
    S2=sum(measure_intersection([speeds[i],speeds[j]],n) for i,j in combinations(range(m),2))
    ENN1=2*S2; EN2=ENN1+EN; Var=EN2-EN*EN
    return EN,ENN1,Var

# ---- m=5, n=6 ----
print("="*70)
print("PART 3 (m=5, n=6): baseline (1-2/6)^5 = (2/3)^5 =", baseline(5,6),"=",float(baseline(5,6)))
print("="*70)
n=6; b=baseline(5,n); records=[]
for s in combinations(range(1,16),5):
    mu=lonely_measure(list(s),n)
    EN,ENN1,Var=moments(list(s),n)
    records.append((float(mu/b),float(mu),float(Var),s))
records.sort()
zero=[r for r in records if r[1]==0]
print(f"  total sets={len(records)}; mu=0 sets (LRC-tight): {len(zero)}")
print("  the mu=0 sets:", [r[3] for r in zero][:20])
print("  smallest positive mu/baseline:")
shown=0
for r,mu,Var,s in records:
    if mu>0:
        print(f"     {s}  mu={mu:.5f}  mu/base={r:.4f}  Var={Var:.4f}")
        shown+=1
        if shown>=10: break
print("  largest mu/baseline:")
for r,mu,Var,s in records[-4:]:
    print(f"     {s}  mu={mu:.5f}  mu/base={r:.4f}  Var={Var:.4f}")

# ---- PART 4: relation between Var/correlation and mu ----
print()
print("="*70)
print("PART 4: Does high pairwise correlation (resonance) <=> small mu?")
print("="*70)
# Define total excess correlation C = sum_{i<j} (ratio_ij - 1) = (E[N(N-1)] - 2*C(m,2)*(2/n)^2)/(2/n)^2 *...
# Simpler: 'excess' = E[N(N-1)] - m(m-1)*(2/n)^2  (positive => net over-correlation).
n=6
def excess(s):
    m=len(s)
    ENN1=2*sum(measure_intersection([s[i],s[j]],n) for i,j in combinations(range(m),2))
    indep=m*(m-1)*(F(2,n))**2
    return float(ENN1-indep)
rows=[]
for s in combinations(range(1,14),5):
    mu=float(lonely_measure(list(s),n))
    rows.append((excess(s),mu,s))
# correlation between excess and mu
import statistics
xs=[r[0] for r in rows]; ys=[r[1] for r in rows]
mx=statistics.mean(xs); my=statistics.mean(ys)
cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))/len(xs)
sx=statistics.pstdev(xs); sy=statistics.pstdev(ys)
print(f"  Pearson corr(excess-correlation, mu) over {len(rows)} sets (n=6,m=5): {cov/(sx*sy):.4f}")
print(f"  (negative => more over-correlation/resonance tends to SMALLER lonely measure)")
# among mu=0 sets, are they all over-correlated?
z=[r for r in rows if r[1]==0]
print(f"  mu=0 sets: {len(z)}; their excess values: min={min(r[0] for r in z):.4f} max={max(r[0] for r in z):.4f} mean={statistics.mean([r[0] for r in z]):.4f}")
pos=[r for r in rows if r[1]>0]
print(f"  mu>0 sets: excess mean={statistics.mean([r[0] for r in pos]):.4f}")
