# Same rank test as identify_new.py for newforms of non-2-power levels (candidate_forms_odd.out).
import re
exec(open('identify_new.py').read().split('results = []')[0].replace('print(f"{len(qs)} primes', 'pass  # print(f"{len(qs)} primes'))
oforms = {}
for l in open('candidate_forms_odd.out'):
    lab, p, tr, nr = l.split(); oforms.setdefault(lab, {})[int(p)] = int(tr)
print(f"{len(qs)} primes; {len(oforms)} orbits; new dimension {rBN - rB}")
hits = []
for lab, d in oforms.items():
    k = lab_k(lab)
    cols = [[q ** j * d[q] for q in qs] for j in range(0, (6 - k) // 2 + 1)]
    if not any(any(c) for c in cols): continue
    rBC = rank(base + cols); dd = rank(base + cols + newcols) - rBC
    if dd < rBN - rB: hits.append((dd, lab))
print("orbits that reduce the new dimension:", sorted(hits) if hits else 'none')
