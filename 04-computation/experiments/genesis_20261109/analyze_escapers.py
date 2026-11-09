"""Classes, converse pairing, decomposability defect and deletion labels of order-12 escapers (genesis lane)."""
import subprocess, sys
from collections import Counter, defaultdict
from grammar import profile, classify
from defect import defsub, twin
from verify_escapers import adj_from, lam

def canon_many(codes, n=12):
    inp = ''.join('%d %s\n' % (n, c) for c in codes)
    r = subprocess.run(['./canon'], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    return [ln.split()[1] for ln in r if ln.strip()]

codes = []
for fn in sys.argv[1:]:
    for line in open(fn):
        p = line.split()
        if not p: continue
        if p[0] == 'E': codes.append(p[1])
        elif p[0] == '12': codes.append(p[1])
can = canon_many(codes)
classes = sorted(set(can))
conv = dict(zip(classes, canon_many([''.join('1' if b == '0' else '0' for b in c) for c in classes])))
sc = [c for c in classes if conv[c] == c]
print('escaper records', len(codes), 'classes', len(classes), 'self-converse', len(sc),
      'classes whose converse is also found', sum(1 for c in classes if conv[c] in set(classes)))
prof = profile([(12, c) for c in classes])
decs = Counter(); hv = defaultdict(list); dellab = Counter(); allok = True
for c, p in zip(classes, prof):
    labels, st = classify(p)
    allok &= (st == 'ESCAPER')
    out = adj_from(c, 12)
    d, M = defsub(out, 12)
    dec = min(p['lam'], d)
    decs[(p['lam'], d, twin(out, 12))] += 1
    hv[p['H']].append(c)
    # labels of the 12 deletions (each is an order-11 non-creator): SUM / SUB / SUM2 / FCV(at order 11) / none
    sub = []
    for v in range(12):
        vs = [x for x in range(12) if x != v]
        code11 = ''.join('1' if (out[vs[i]] >> vs[j]) & 1 else '0' for i in range(11) for j in range(i + 1, 11))
        sub.append((11, code11))
    labs = set()
    for q in profile(sub):
        l, _ = classify(q)
        labs.add(','.join(l) if l else 'ESC11')
    dellab[tuple(sorted(labs))] += 1
print('all re-classified as escapers:', allok)
print('H range', min(hv), max(hv), ' distinct H', len(hv), ' classes per H:', Counter(len(v) for v in hv.values()))
print('(lambda, module defect, twin) -> classes:', dict(decs))
print('label sets of the 12 deletions (order 11):')
for k, v in dellab.most_common(): print('  ', v, k)
