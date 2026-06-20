import sys, itertools
from fractions import Fraction
if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
from lrc14_verify_glaisher_indep import meas_S7, maxrun, richC2, odd_part

def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    sx=sum((x-mx)**2 for x in xs)**.5; sy=sum((y-my)**2 for y in ys)**.5
    return cov/(sx*sy) if sx*sy else 0

# n_distinct_odd_parts
def ndistinct_odd(E): return len({odd_part(e) for e in E})

print("="*70)
print("PART 4: headline correlations at k=8 (verify claim's numbers)")
print("="*70)
vals=[(meas_S7(E),E) for E in itertools.combinations(range(1,14),8)]
p0s=[float(v) for v,E in vals]
runs=[maxrun(E) for v,E in vals]
c2=[richC2(E) for v,E in vals]
nod=[ndistinct_odd(E) for v,E in vals]
print(f"  corr(p0, contiguity/maxrun) = {pearson(p0s,runs):+.4f}  (claimed +0.268)")
print(f"  corr(p0, richness C2)       = {pearson(p0s,c2):+.4f}  (claimed +0.156)")
print(f"  corr(p0, #distinct-odd)     = {pearson(p0s,nod):+.4f}  (claimed -0.190)")

import collections
print("\n  mean p0 by max-run:")
byrun=collections.defaultdict(list)
for v,E in vals: byrun[maxrun(E)].append(float(v))
for r in sorted(byrun): print(f"     run={r}: mean_p0={sum(byrun[r])/len(byrun[r]):.4f} (n={len(byrun[r])})")

print()
print("="*70)
print("PART 5: cap comparison — is consec p0 below the EXACT caps?")
print("  cap_8=2243/5880, cap_9=1979/4004, cap_10=55/91")
print("="*70)
caps={8:Fraction(2243,5880),9:Fraction(1979,4004),10:Fraction(55,91)}
for k in (8,9,10):
    pc=meas_S7(tuple(range(1,k+1)))
    cap=caps[k]
    print(f"  k={k}: consec p0={float(pc):.5f}={pc}  cap={float(cap):.5f}={cap}  "
          f"p0 {'<=' if pc<=cap else '>'} cap  margin={float(cap-pc):+.5f}")
print("  NOTE: if consec p0 > cap, the sector route does NOT close for that k;")
print("        the binding question is whether MAX over E of p0 <= cap.")

print()
print("="*70)
print("PART 6: SUMMARY — does dyadic richness EXPLAIN consec-maximality?")
print("="*70)
# The decisive test: among the 40 single-swaps from consec{1..8}, rank the
# drop in p0 against (a) whether the removed element was in the pow2 chain,
# (b) the position of the removed element (interior hole vs boundary).
base=list(range(1,9)); p0b=meas_S7(base)
pow2={1,2,4,8}
rows=[]
for i in range(8):
    rem=base[i]
    for new in range(9,14):
        E=base[:i]+[new]+base[i+1:]
        drop=float(p0b-meas_S7(E))
        rows.append((drop,rem,rem in pow2,new))
# average drop when removing a pow2 element vs non-pow2
dp_pow2=[d for d,r,ip,n in rows if ip]
dp_non=[d for d,r,ip,n in rows if not ip]
print(f"  avg p0-drop removing a POW2 element {{1,2,4,8}}: {sum(dp_pow2)/len(dp_pow2):.4f}")
print(f"  avg p0-drop removing a NON-pow2 element       : {sum(dp_non)/len(dp_non):.4f}")
print(f"  => removing chain elements hurts LESS than non-chain: "
      f"{sum(dp_pow2)/len(dp_pow2) < sum(dp_non)/len(dp_non)}")
print("  This INVERTS the dyadic-chain mechanism prediction.")
