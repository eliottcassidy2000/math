"""Independent verifier for merge-tree certificates in the 'MERGE a b c = x y z + x2 y2 z2' format (THM-4820 audit).

Checks, written without reference to verify_tree.py:
  * every non-blank line parses as a MERGE line; parent = left + right; both children nonzero and nonnegative;
  * multiset accounting: the non-unit children are exactly the parents minus the root (each produced once), so the
    lines form a single binary tree (sizes strictly decrease along child links, hence no cycles);
  * an explicit top-down assignment of lines to nodes uses every line exactly once;
  * the unit-leaf multiset equals the root count vector;
  * cost recomputed in the sqrt(-3) model: value(p,q,r) = (X + Y sqrt(-3))/2, X = 2p - q - r, Y = q - r,
    norm = (X^2 + 3Y^2)/4 (exact integer), beta = norm(parent) - norm(left) - norm(right);
    cross-checked against the count formula 3 u.v - |u||v|;
  * potential bookkeeping: Phi(root) = floor(norm/3) and the slack of every merge; prints the merges with slack > 0.
Usage: python3 audit_tree.py FILE a,b,c
"""
import re
import sys
from collections import Counter, deque

fn, root = sys.argv[1], tuple(map(int, sys.argv[2].split(",")))
pat = re.compile(r"^MERGE\s+(\d+)\s+(\d+)\s+(\d+)\s*=\s*(\d+)\s+(\d+)\s+(\d+)\s*\+\s*(\d+)\s+(\d+)\s+(\d+)\s*$")
lines = []
for raw in open(fn):
    if not raw.strip():
        continue
    m = pat.match(raw.strip())
    assert m, f"unparsable line: {raw!r}"
    t = tuple(map(int, m.groups()))
    lines.append((t[0:3], t[3:6], t[6:9]))


def norm4(u):  # 4 * norm in the sqrt(-3) model
    p, q, r = u
    X, Y = 2 * p - q - r, q - r
    return X * X + 3 * Y * Y


def norm(u):
    n4 = norm4(u)
    assert n4 % 4 == 0
    return n4 // 4


ok = True
for P, L, Rr in lines:
    ok &= all(a + b == c for a, b, c in zip(L, Rr, P))
    ok &= sum(L) > 0 and sum(Rr) > 0 and min(L + Rr) >= 0
assert ok, "a line is not a valid merge"
parents = Counter(P for P, _, _ in lines)
children = Counter(ch for _, L, Rr in lines for ch in (L, Rr) if sum(ch) > 1)
units = Counter(ch for _, L, Rr in lines for ch in (L, Rr) if sum(ch) == 1)
assert parents - children == Counter({root: 1}) and not (children - parents), "parent/child multisets do not form one tree"
assert units == Counter({(1, 0, 0): root[0], (0, 1, 0): root[1], (0, 0, 1): root[2]}), f"leaf multiset {units}"
# explicit top-down assignment
avail = {}
for i, (P, L, Rr) in enumerate(lines):
    avail.setdefault(P, []).append(i)
used = [False] * len(lines)
dq = deque([root])
nleaves = 0
while dq:
    u = dq.popleft()
    if sum(u) == 1:
        nleaves += 1
        continue
    i = avail[u].pop()
    used[i] = True
    dq.extend(lines[i][1:])
assert all(used), "unused lines"
cost = 0
slack_tot = 0
pos = Counter()
slack_merges = []
for P, L, Rr in lines:
    beta = norm(P) - norm(L) - norm(Rr)
    bc = 3 * sum(a * b for a, b in zip(L, Rr)) - sum(L) * sum(Rr)
    assert beta == bc, "value beta != count-formula beta"
    bp = max(beta, 0)
    cost += bp
    pos[beta] += 1
    sl = bp - (norm(P) // 3 - norm(L) // 3 - norm(Rr) // 3)
    assert sl >= 0
    slack_tot += sl
    if sl:
        slack_merges.append((P, L, Rr, beta, sl))
X, Y = 2 * root[0] - root[1] - root[2], root[1] - root[2]
print(f"{fn}: {len(lines)} merges, {nleaves} unit leaves, root count {root}, root value (X+Y sqrt(-3))/2 with "
      f"X={X}, Y={Y} (= {X // 2 if Y == 0 else 'non-rational'}), norm {norm(root)}")
print(f"  total cost sum beta+ = {cost};  Phi(root) = {norm(root) // 3};  total slack = {slack_tot};  "
      f"Phi + slack = {norm(root) // 3 + slack_tot}")
print(f"  beta histogram: {dict(sorted(pos.items()))};  sum beta = {sum(b * c for b, c in pos.items())} "
      f"(= Q(root) - #leaves = {norm(root) - nleaves})")
print(f"  merges with positive slack: {slack_merges}")
