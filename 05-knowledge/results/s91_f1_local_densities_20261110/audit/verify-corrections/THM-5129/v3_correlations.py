#!/usr/bin/env python3
"""Verifier (THM-5129 (D)): Pearson correlations of the covered share of [h, 3h/2] with the number of strong
classes, recomputed from the two rows printed in the session's source_supply.out (sections (5) and (6))."""
import io, re, sys
P = ("C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/source_supply.out")
txt = io.open(P, encoding="utf-8").read()
sec5 = txt.split("(5) h:")[1].split("h with the whole window")[0]
sec6 = txt.split("(6) h:")[1].split("Pearson")[0]
share = {int(a): int(b) / int(c) for a, b, c in re.findall(r"(\d+):(\d+)/(\d+)", sec5)}
m = {int(a): int(b) for a, b in re.findall(r"(\d+):(\d+)\b(?!/)", sec6.split("\n", 1)[1])}
def pearson(hs):
    xs = [share[h] for h in hs]; ys = [m[h] for h in hs]; n = len(hs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs); syy = sum((y - my) ** 2 for y in ys)
    return sxy / (sxx * syy) ** 0.5
ok = True
for lo, want in ((23, 0.68), (15, 0.66), (3, 0.41)):
    hs = [h for h in sorted(share) if lo <= h <= 73]
    r = pearson(hs)
    print("%d <= h <= 73: %d attained offsets, Pearson correlation %.4f (quoted %.2f)" % (lo, len(hs), r, want))
    ok &= abs(r - want) < 0.005
two = [h for h in sorted(m) if m[h] == 2]
print("offsets <= 73 with exactly two strong classes:", two)
ok &= two == [9, 17, 19, 35, 39, 47, 63]
print("m(67) =", m[67], "; m(27) =", m[27], "; m(55) =", m[55])
print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
