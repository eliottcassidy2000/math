from s3_residue_colorings import *
F = MQ([-3, -11])
w1 = F.elt({0: Fr(1, 2), 1: Fr(1, 2)}); w3 = F.elt({0: Fr(5, 6), 2: Fr(1, 6)}); w3b = F.conj(w3)
gens = []; z = F.one()
for j in range(6):
    gens += [z, F.mul(z, w3), F.mul(z, w3b)]; z = F.mul(z, w1)
pts = build_patch(F, gens, 3, 2500)
edges, cand, diffs = exact_edges(F, pts)
G = set(gens) | set(tuple(-c for c in g) for g in gens)
allu = set()
for d in diffs:
    allu.add(d); allu.add(tuple(-c for c in d))
extra = [u for u in allu if u not in G]
print("unit vectors found:", len(allu), " generators:", len(G), " extra:", len(extra))
import cmath
for u in sorted(extra, key=lambda u: cmath.phase(F.numval(u)))[:12]:
    print("  ", [str(c) for c in u], " angle(deg)=%.4f" % math.degrees(cmath.phase(F.numval(u))))
# Z-span check: are the extras in Z-span of the 18? (denominators)
print("max denominator among extras:", max(max(c.denominator for c in u) for u in extra))
