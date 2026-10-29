"""Merge the n=7 k=3 census slices, compare every class with the producer's census, and summarise glue losses."""
import re
from compare_census import parse_prod, parse_mine
lines = []
for f in ('out/census_n7_k3.out', 'out/census_n7_k3_s1.out', 'out/census_n7_k3_s2.out', 'out/census_n7_k3_s3.out'):
    lines += [l for l in open(f) if re.match(r'^[01]+ ', l)]
open('out/census_n7_k3_merged.out', 'w').writelines(lines)
P = parse_prod('../out/census_n7_k3.out')
M = parse_mine('out/census_n7_k3_merged.out')
assert set(P) == set(M), (len(P), len(M))
bad = [b for b in P if (P[b][0], P[b][1], P[b][2]) != M[b][:3] or P[b][3] != M[b][3]]
print('n=7 k=3: %d classes; disagreements with the producer (Smith M, Smith N, spur2/3, kills): %d' % (len(M), len(bad)))
loss = sum('gluediff=True' in l for l in lines)
spur = sum(not re.search(r'spur2=0 spur3=0', l) for l in lines)
hidden = sum('gluediff=True' in l and re.search(r'spur2=0 spur3=0', l) is not None for l in lines)
kills = sum('kills=True' in l for l in lines)
print('Smith form of N differs from M: %d classes; spur2+spur3 > 0: %d; differs with spur = 0: %d; boundary killers: %d'
      % (loss, spur, hidden, kills))
