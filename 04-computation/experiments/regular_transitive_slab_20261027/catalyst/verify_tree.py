"""Independent check of a merge-tree certificate (THM-4820): reads 'MERGE a b c = x y z + x2 y2 z2' lines (count vectors of
unit charges 1, omega, omega^2), checks that they form one binary tree whose leaves are single unit charges and whose
root is the given count vector, and recomputes the cost sum of max(0, Q(x+y) - Q(x) - Q(y)) with
Q(a + b omega) = a^2 - ab + b^2 evaluated on the Eisenstein VALUES (not on the count formula used by the DP)."""
import sys
from collections import Counter
def value(u):   # count vector (p, q, r) -> Eisenstein integer p + q w + r w^2 = (p - r) + (q - r) w
    p, q, r = u; return (p - r, q - r)
def Q(z): a, b = z; return a * a - a * b + b * b
merges = []
for line in open(sys.argv[1]):
    if not line.startswith('MERGE'): continue
    t = list(map(int, line.replace('=', ' ').replace('+', ' ').split()[1:]))
    merges.append((tuple(t[0:3]), tuple(t[3:6]), tuple(t[6:9])))
root = tuple(map(int, sys.argv[2].split(',')))
need = Counter([root]); cost = 0; leaves = Counter()
parents = Counter(m[0] for m in merges)
pool = Counter([root]); used = Counter()
# walk: every parent occurrence must be produced exactly once; children are either merged further or are unit leaves
avail = Counter(m[0] for m in merges)
todo = [root]; mi = Counter(); bymerge = {}
for mg in merges: bymerge.setdefault(mg[0], []).append(mg)
while todo:
    u = todo.pop()
    if sum(u) == 1: leaves[u] += 1; continue
    lst = bymerge.get(u)
    assert lst, f"no merge produces {u}"
    mg = lst.pop(); _, l, r = mg
    assert tuple(a + b for a, b in zip(l, r)) == u and sum(l) > 0 and sum(r) > 0
    x, y = value(l), value(r); s = (x[0] + y[0], x[1] + y[1])
    cost += max(0, Q(s) - Q(x) - Q(y)); todo += [l, r]
assert all(len(v) == 0 for v in bymerge.values()), "unused merges"
assert leaves == Counter({(1, 0, 0): root[0], (0, 1, 0): root[1], (0, 0, 1): root[2]}), leaves
print(f"tree for {root}: {len(merges)} merges, {sum(leaves.values())} unit leaves, root value {value(root)}, total positive interaction cost {cost}")
