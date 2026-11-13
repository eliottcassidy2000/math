"""Single ternary clock in the nested-funnel (automatically ready) language:
   core prefix -> inner funnel L1 b^j R1 -> fixed suffix words (the last one is the outer funnel) -> ready z.
   Ready count R(j) = f . c1(j) - 1 with f = (1,0,0,0,1,0) M_suffix.  Want R(j)+1 = c 3^j for all j>=1."""
import sys, itertools
import numpy as np
sys.path.insert(0, 'scripts')
from funnel import step6, SINGLETON
from atlas import cores, words_upto
from twoclock_search import contexts, mid_matrix
maxcore, inctx, sufmax, sufword = map(int, sys.argv[1:5])
C = cores(maxcore)
ICTX = contexts(inctx)
ws = words_upto(sufword)
sufs = []
def gen(prefix, budget):
    if prefix: sufs.append(tuple(prefix))
    for w in ws:
        if len(w)+1 <= budget:
            gen(prefix+[w], budget-len(w)-1)
gen([], sufmax)
print('cores', len(C), 'inner', len(ICTX), 'suffixes', len(sufs), flush=True)
F = []
for s in sufs:
    M = mid_matrix(s)
    f = np.array([1,0,0,0,1,0], dtype=object).dot(M)
    F.append((s, f))
FA = np.array([[int(x) for x in f] for s, f in F], dtype=np.int64)
hits = []
for seq, (c0, nv) in C.items():
    for (L1, b1, R1) in ICTX:
        cs = [np.array(step6(c0, L1 + b1*j + R1), dtype=object) for j in range(1, 10)]
        D = np.array([[int(x) for x in (cs[j+1] - 3*cs[j])] for j in range(8)], dtype=np.int64).T
        ok = np.nonzero(np.abs(FA.dot(D)).max(axis=1) == 0)[0]
        for i in ok:
            s, f = F[i]
            vals = []
            for j in range(1, 7):
                st = c0
                for w in [L1 + b1*j + R1] + list(s): st = step6(st, w)
                vals.append(st[0] + st[4] - 1)
            if all(vals[t+1] == 3*vals[t] + 2 for t in range(5)):
                nvert = nv + len(L1) + len(R1) + sum(len(w)+1 for w in s) + 2
                hits.append((vals[0] + 1, nvert, seq, (L1, b1, R1), s))
hits.sort(key=lambda h: (h[0], h[1]))
print('hits', len(hits))
seen = set()
for h in hits:
    if h[0] in seen: continue
    seen.add(h[0])
    print('R(j)+1 = %d*3^(j-1) ;  order(j=1)=%d  core=%s inner=%s suffix=%s' % (h[0], h[1], h[2], h[3], h[4]))
