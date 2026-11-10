# a09b: compare the audit's brute-force gadget matrices with the author's gadget tables (as matrix sets),
# count distinct matrices per maximal order, write the distinct-matrix files used by the audit closures,
# and test the structural facts the exactness argument needs: E_u = M[0][2] >= 1, s_u = -M[0][1] <= h_G = M[0][0].
import sys, hashlib
SCR, AUTH = sys.argv[1], sys.argv[2]
def load(fn):
    rows = []
    for l in open(fn):
        p = l.split()
        if len(p) < 13: continue
        rows.append((int(p[0]), tuple(int(x) for x in p[4:13])))
    return rows
mine5 = load(SCR + '/agad5_labelled.txt'); mine7 = load(SCR + '/agad7_classes.txt')
auth5 = load(AUTH + '/gadgets5.txt'); auth7 = load(AUTH + '/gadgets7c.txt')
for name, rows in (('audit labelled <=5', mine5), ('audit classes <=7', mine7), ('author gadgets5', auth5), ('author gadgets7c', auth7)):
    print('%-20s rows %6d' % (name, len(rows)) + ''.join('  distinct(order<=%d) %d' % (G, len(set(m for g, m in rows if g <= G))) for G in range(1, 8) if any(g == G for g, m in rows)))
for G in range(1, 6):
    a = set(m for g, m in mine5 if g <= G); b = set(m for g, m in mine7 if g <= G); c = set(m for g, m in auth5 if g <= G)
    print('order <= %d: audit labelled == audit classes: %s ; audit == author gadgets5: %s' % (G, a == b, a == c))
for G in range(1, 8):
    a = set(m for g, m in mine7 if g <= G); c = set(m for g, m in auth7 if g <= G)
    print('order <= %d: audit classes == author gadgets7c: %s (%d vs %d)' % (G, a == c, len(a), len(c)))
allm = set(m for g, m in mine7) | set(m for g, m in mine5)
print('structural facts over all %d distinct matrices: min E_u = %d ; violations of 0 <= s_u <= h_G: %d ; min h_G = %d'
      % (len(allm), min(m[2] for m in allm), sum(1 for m in allm if not (0 <= -m[1] <= m[0])), min(m[0] for m in allm)))
for G in (5, 6, 7):
    ms = sorted(set(m for g, m in mine7 if g <= G))
    with open(SCR + '/amat_le%d.txt' % G, 'w') as f:
        for m in ms: f.write(' '.join(map(str, m)) + '\n')
    print('wrote amat_le%d.txt: %d distinct matrices, sha256 %s' % (G, len(ms), hashlib.sha256(open(SCR + '/amat_le%d.txt' % G, 'rb').read()).hexdigest()[:16]))
