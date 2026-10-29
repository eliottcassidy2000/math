"""Conjecture check: coker M_k(T) is torsion-free (every nonzero invariant factor 1) iff every strong component of T
has at most one Hamiltonian cycle (a strong component of order 1 counts as having none)."""
import re, itertools, sys, glob
from runs import tour_from_bits

def ham_cycles_sub(A, verts):
    verts = list(verts)
    m = len(verts)
    if m < 3:
        return 0
    cnt = 0
    for w in itertools.permutations(verts[1:]):
        cyc = (verts[0],) + w
        if all(A[cyc[i]][cyc[(i + 1) % m]] for i in range(m)):
            cnt += 1
    return cnt

def strong_components(A):
    n = len(A)
    reach = [{i} for i in range(n)]
    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(n):
                if A[i][j] and not reach[j] <= reach[i]:
                    reach[i] |= reach[j]; changed = True
    comps = {}
    for i in range(n):
        key = frozenset(j for j in reach[i] if i in reach[j])
        comps[key] = 1
    return list(comps)

tally = {}
viol = []
for f in sorted(glob.glob('out/census_n*_k*.out')):
    n, k = map(int, re.findall(r'n(\d+)_k(\d+)', f)[0])
    if k >= n:
        continue
    for l in open(f):
        if not re.match(r'^[01]{3,}\s', l):
            continue
        m = re.search(r' M=(.*?)( unglued=| P=|$)', l)
        if not m:
            continue
        bits = l.split()[0]
        A = tour_from_bits(n, bits)
        comps = strong_components(A)
        hcs = [ham_cycles_sub(A, c) for c in comps]
        pred = all(h <= 1 for h in hcs)
        tf = all(t.startswith('1^') or t.startswith('0^') for t in m.group(1).split())
        key = (n, k)
        tally.setdefault(key, [0, 0])
        tally[key][0] += 1
        if pred != tf:
            tally[key][1] += 1
            viol.append((n, k, bits, hcs, m.group(1)[:80]))
for key in sorted(tally):
    print('n=%d k=%d: classes=%d, violations of [torsion-free <=> every strong component has <= 1 HC]: %d' % (key + tuple(tally[key])))
for v in viol[:40]:
    print('  VIOLATION', v)
