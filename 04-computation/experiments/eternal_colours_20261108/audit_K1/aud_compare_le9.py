#!/usr/bin/env python3
# audit K1: compare, order by order (n <= 9), the self-converse value sets and strong flags of
#   (a) THM-4922's census (re-verified class by class by aud_census.c), and
#   (b) the auditor's anti-invariant enumeration aud_sc.c (modes full and red).
import sys, os
W = sys.argv[1]
ok = True
for n in range(1, 10):
    cen = {}
    for line in open(os.path.join(W, f"census_vals_{n}.txt")):
        _, h, s, ss, a = map(int, line.split())
        if s: cen[h] = int(ss > 0)
    for mode in ("full", "red"):
        an = {}
        for line in open(os.path.join(W, f"anti_{n}_{mode}.txt")):
            h, c, sc = map(int, line.split()); an[h] = int(sc > 0)
        same = an == cen
        ok &= same
        print(f"order {n} {mode:4s}: census SC values {len(cen)} (strong {sum(cen.values())}); anti-invariant {len(an)} "
              f"(strong {sum(an.values())}); identical incl. strong flags: {same}")
        if not same:
            print("   only census:", sorted(set(cen) - set(an))[:20], " only anti:", sorted(set(an) - set(cen))[:20],
                  " flag diffs:", [h for h in cen if h in an and cen[h] != an[h]][:20])
print("ALL ORDERS <= 9 AGREE" if ok else "DISAGREEMENT")
