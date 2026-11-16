#!/usr/bin/env python3
"""maximisers.py: classify the maximising pairs (W, l) listed by `hratio4 ... print` up to isomorphism.
usage: maximisers.py print_file [print_file ...]
Also tests the block structure: B = the 7-vertex tournament of the word 0 0 7 6 0 followed by the letter 6; for each
length n the set of canonical forms of the maximisers' T_(W l) is compared with the set of canonical forms of T_(W') => B
(every arc from W' into B) over all literal words W' of length n - 6 (n - 5 vertices; the empty tournament at n = 5).
Only the pairs attaining the largest listed value per (kind, length) are kept, so that the eight restricted listings
(first letter fixed; at lengths 5-6 some first letters attain only a smaller local maximum) can be classified together.
For every listed pair the literal span-3 tournament T_W (first vertex, then one vertex per letter: the new vertex beats
the older vertices at the gaps in the letter, loses to all others) and T_{W l} are built; nauty's labelg gives canonical
digraph6 strings.  Per length and per kind (maxH / maxZ) the script reports the number of listed pairs, the number of
distinct canonical forms of T_{W l}, whether T_W is transitive, and which ranks of the transitive order of T_W the new
vertex beats (ranks 1 = source ... n+1 = sink)."""
import sys, subprocess
from collections import defaultdict
def build(word):
    n = len(word) + 1; adj = [[0] * n for _ in range(n)]
    for i, l in enumerate(word, start=1):
        for g in (1, 2, 3):
            j = i - g
            if j >= 0:
                if l >> (g - 1) & 1: adj[i][j] = 1
                else: adj[j][i] = 1
        for j in range(0, i - 3): adj[j][i] = 1
    return adj
def digraph6(adj):
    n = len(adj); bits = ''.join(str(adj[i][j]) for i in range(n) for j in range(n))
    bits += '0' * (-len(bits) % 6)
    return '&' + chr(63 + n) + ''.join(chr(63 + int(bits[k:k + 6], 2)) for k in range(0, len(bits), 6))
def canon(strings):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(strings) + '\n', capture_output=True, text=True, check=True).stdout.split()
    assert len(out) == len(strings); return out
def transitive_order(adj):
    n = len(adj); sc = [sum(adj[i]) for i in range(n)]
    if sorted(sc) != list(range(n)): return None
    return sorted(range(n), key=lambda i: -sc[i])     # rank 1 = source
def block_join(wprime):
    """T_(W') => B: vertices of W' first, then the block B = T(0 0 7 6 0 + letter 6); all arcs from W' to B."""
    A = build(wprime) if wprime is not None else []; k = len(A)
    Bm = build([0, 0, 7, 6, 0, 6]); n = k + len(Bm)
    adj = [[0] * n for _ in range(n)]
    for i in range(k):
        for j in range(k): adj[i][j] = A[i][j]
        for j in range(k, n): adj[i][j] = 1
    for i in range(len(Bm)):
        for j in range(len(Bm)): adj[k + i][k + j] = Bm[i][j]
    return adj
import itertools
groups = defaultdict(list)
for fn in sys.argv[1:]:
    for line in open(fn):
        line = line.strip()
        if not (line.startswith('maxH') or line.startswith('maxZ')): continue
        kind, _, rest = line.partition(' length ')
        n = int(rest.split(':')[0]); body = rest.split(':', 1)[1]
        letter = int(body.split('letter')[1].split(',')[0])
        wpart = body.split('word')[1].split('(')[0].split()
        word = list(map(int, wpart)); frac = line.split('(')[-1].rstrip(')')
        groups[(kind, n)].append((word, letter, frac))
from fractions import Fraction
for (kind, n) in sorted(groups):
    allitems = groups[(kind, n)]
    top = max(Fraction(f) for _, _, f in allitems)
    items = [it for it in allitems if Fraction(it[2]) == top]
    dropped = len(allitems) - len(items)
    full = [digraph6(build(w + [l])) for w, l, _ in items]
    base = [digraph6(build(w)) for w, l, _ in items]
    cf = canon(full); cb = canon(base)
    ranks = set(); trans = True
    for (w, l, _) in items:
        adj = build(w); order = transitive_order(adj)
        if order is None: trans = False; continue
        m = len(w) + 1
        beaten = [m - g for g in (1, 2, 3) if l >> (g - 1) & 1 and m - g >= 0]   # the new vertex beats older vertex m-g
        ranks.add(tuple(sorted(order.index(v) + 1 for v in beaten)))
    fr = sorted(set(f for _, _, f in items))
    block_note = ""
    if kind == 'maxH' and n >= 5:
        if n == 5: preds = [digraph6(block_join(None))]
        else: preds = [digraph6(block_join(list(wp))) for wp in itertools.product(range(8), repeat=n - 6)]
        cp = set(canon(preds)); block_note = f"; block structure T_(W') => B over all W' of length {n - 6}: {len(cp)} classes, equal to the maximiser classes: {cp == set(cf)}"
    print(f"{kind} length {n}: {len(items)} pairs at the maximum {fr}" + (f" ({dropped} listed pairs below it dropped)" if dropped else "") + f", distinct T_(W l) up to isomorphism: {len(set(cf))}, distinct T_W: {len(set(cb))}, T_W transitive for all: {trans}, ranks beaten by the new vertex: {sorted(ranks)}" + block_note)
