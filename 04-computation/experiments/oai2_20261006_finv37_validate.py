#!/usr/bin/env python3
"""Independent validation of the certified UNSAT clause set for the fully gap-determined Erdos-592 game Q_gap(3,7)
(THM-4560 / THM-470 C; mac-mini-2026-10-06-oaimath2, written independently of the reader's finv_unsat_audit.py).

Variables 1..1098 = the lex-positive gap vectors of {-6..6}^3 in lexicographic order (an edge x < y of [7]^3 is blue iff
its gap y - x is in E). Every clause must be a genuine constraint:
  * negative clause {-a, -b, -c}: a triangle, i.e. d_a + d_b = d_c with x, x + d_a, x + d_a + d_b all in [7]^3 for some x
    (per coordinate max(|a_i|, |b_i|, |a_i + b_i|) <= 6);
  * positive clause: contains every gap of a binary subgrid of [7]^3 (8 leaves stored in the .leaves file; tree check).
If all clauses are genuine, an UNSAT verdict for the clause set (any solver) is an UNSAT verdict for Q_gap(3,7).
Run: python3 oai2_20261006_finv37_validate.py <cnf> <leaves>   (gunzip the files in oai2_20261006_readers/hindman_ramsey/)
"""
import sys, itertools
from collections import Counter
CNF, LEAVES = sys.argv[1], sys.argv[2]
t, n = 7, 3
gaps = [d for d in itertools.product(range(-(t-1), t), repeat=n) if next((x for x in d if x), 0) > 0]
gaps.sort()
assert len(gaps) == 1098
idx = {d: i + 1 for i, d in enumerate(gaps)}
vec = {i + 1: d for i, d in enumerate(gaps)}
# sanity of the inferred numbering against the reader's first clauses: (0,0,1) -> 1, (0,1,-6) -> 7
assert idx[(0, 0, 1)] == 1 and idx[(0, 0, 2)] == 2 and idx[(0, 1, -6)] == 7 and idx[(1, -1, -1)] == 155
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def realizable(a, b):
    return all(max(abs(x), abs(y), abs(x + y)) <= t - 1 for x, y in zip(a, b))
neg_ok = neg_bad = 0
pos_cnf = Counter()
with open(CNF) as f:
    header = f.readline().split()
    assert header[:2] == ['p', 'cnf'] and int(header[2]) == 1098
    m = int(header[3])
    count = 0
    for line in f:
        lits = [int(x) for x in line.split()]
        assert lits[-1] == 0
        lits = lits[:-1]
        count += 1
        if all(l < 0 for l in lits):
            vs = [vec[-l] for l in lits]
            ok = False
            if len(vs) == 3:
                for a, b, c in itertools.permutations(vs):
                    if add(a, b) == c and realizable(a, b):
                        ok = True
                        break
            neg_ok += ok
            neg_bad += (not ok)
        elif all(l > 0 for l in lits):
            pos_cnf[tuple(sorted(set(lits)))] += 1
        else:
            raise SystemExit("mixed clause: %s" % lits)
    assert count == m, (count, m)
def is_binary_subgrid(pts):
    # pts: 8 distinct points of [t]^3; check the tree structure: 2 first coordinates, under each 2 second coordinates,
    # under each (first, second) 2 third coordinates
    if len(set(pts)) != 8 or any(not (0 <= x < t) for p in pts for x in p): return False
    A = {p[0] for p in pts}
    if len(A) != 2: return False
    for a in A:
        B = {p[1] for p in pts if p[0] == a}
        if len(B) != 2: return False
        for b in B:
            C = {p[2] for p in pts if p[0] == a and p[1] == b}
            if len(C) != 2: return False
    return True
pos_ok = pos_bad = 0
pos_leaves = Counter()
with open(LEAVES) as f:
    for line in f:
        lit_part, leaf_part = line.split('|')
        lits = tuple(sorted(set(int(x) for x in lit_part.split())))
        pts = [tuple(int(v) for v in s.split(',')) for s in leaf_part.split()]
        ok = is_binary_subgrid(pts)
        if ok:
            gs = set()
            for p, q in itertools.combinations(sorted(pts), 2):
                gs.add(idx[tuple(y - x for x, y in zip(p, q))])
            ok = gs.issubset(set(lits))
        pos_ok += ok
        pos_bad += (not ok)
        pos_leaves[lits] += 1
missing = sum((pos_cnf - pos_leaves).values())
print(f"clauses {m}: negative {neg_ok} valid / {neg_bad} invalid; positive (from .leaves) {pos_ok} valid / {pos_bad} invalid; "
      f"positive clauses in CNF without a certified subgrid: {missing}")
print("VALID" if neg_bad == 0 and pos_bad == 0 and missing == 0 else "INVALID")
