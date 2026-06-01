from fractions import Fraction as F
from math import gcd
from lrc import measure_intersection

# Derive exact formula for measure(B_i ∩ B_j).
# Claim: by a change of variable, the overlap depends only on g=gcd(vi,vj),
# a=vi/g, b=vj/g (coprime), and n. In fact since substituting t->t/g maps
# the structure, and B is periodic, the measure(B_a ∩ B_b) with gcd(a,b)=1
# determines everything: measure(B_{ga} ∩ B_{gb}) = measure(B_a ∩ B_b).
# Verify that scaling invariance.

print("Scaling invariance: measure(B_{g a} ∩ B_{g b}) == measure(B_a ∩ B_b) ?")
for n in [5,7,11]:
    for (a,b) in [(1,2),(2,3),(1,3),(3,5),(1,5)]:
        base = measure_intersection([a,b], n)
        for g in [2,3,4,5]:
            assert measure_intersection([g*a,g*b], n) == base, (n,a,b,g)
    print(f"  n={n}: scaling invariance holds.")

# So overlap is a function of coprime (a,b) and n only. The arcs of B_a (period 1/a)
# vs B_b. The relevant object: on the torus, ||a t||<1/n and ||b t||<1/n.
# Let u = a t mod 1, then b t = (b/a) u ... messy. Instead use the 2D lattice:
# The fraction of t with ||a t||<1/n & ||b t||<1/n. Since gcd(a,b)=1, the map
# t -> (a t mod 1, b t mod 1) is NOT equidistributed on the 2-torus (it's a line),
# but over t in [0,1) the pair (frac(a t), frac(b t)) traces the closed geodesic
# of slope b/a, hitting the box [-1/n,1/n]^2 (centered, mod 1).
# Exact count: measure = (1/(ab)) * #{(k,l): the line segment passes through cell}...
# Simplest: just compute exactly and tabulate vs (a mod n, b mod n).

print()
print("Hypothesis: for coprime (a,b), overlap depends on (a mod n, b mod n).")
print("Tabulate correlation ratio r = measure * n^2/4 for coprime pairs, n=7:")
n=7
seen={}
import itertools
def coprime(a,b): return gcd(a,b)==1
rows=[]
for a in range(1,40):
    for b in range(a+1,40):
        if not coprime(a,b): continue
        r = measure_intersection([a,b],n)*F(n*n,4)
        key=(a%n, b%n)
        rows.append((a,b,a%n,b%n,float(r)))
# group by residue pair (unordered, and also ± since ||.|| symmetric)
from collections import defaultdict
groups=defaultdict(set)
for a,b,ar,br,r in rows:
    # canonical residue: ||.|| symmetric under negation, so map residue x to min(x, n-x)
    def fold(x): return min(x%n, (-x)%n)
    key=tuple(sorted((fold(ar),fold(br))))
    groups[key].add(round(r,6))
print("  residue-pair (folded mod n, |.|) -> set of ratios observed:")
for key in sorted(groups):
    print(f"     {key}: {sorted(groups[key])}")
