"""f_section4.py - section 4 table of PROOF.md from the auditor's own counts (d_image.c, window [2^22, 2^23)).
union miss = fraction of odd N in the window hit by no member; product = sum_c n_c prod_s (1 - v_sc/n_c) / sum_c n_c."""
import json, sys, math
from collections import defaultdict
res = {}
for f in ('bankA_all', 'bankA_nobal', 'bankA_generic34', 'bankA_r2d16'):
    res[f] = {}
    for M in (3, 945, 12288):
        n = {}; prod = defaultdict(lambda: 1.0); uh = {}; seeds = set()
        for l in open(f'out_w22_{f}_M{M}.txt'):
            t = l.split()
            if t[0] == 'S':
                h, R, e, d, MM, c, nc, p, v = map(int, t[1:]); n[c] = nc; prod[c] *= (1 - v / nc); seeds.add((h, R, e, d))
            else:
                MM, c, nc, u = map(int, t[1:]); uh[c] = u
        N = sum(n.values())
        res[f][M] = {'seeds': len(seeds), 'union_miss': round(sum(n[c] - uh[c] for c in n) / N, 5),
                     'product_of_member_miss': round(sum(n[c] * prod[c] for c in n) / N, 5)}
        if M == 3:
            res[f][M]['union_miss_by_class'] = [round((n[c] - uh[c]) / n[c], 4) for c in sorted(n)]
print(json.dumps(res, indent=1))
json.dump(res, open('out_section4_w22.json', 'w'), indent=1)
