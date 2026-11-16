"""(1) Atoms of the multiplicative monoid of counts (odd numbers except 7, 21; true up to 20,405 by the census).
(2) Order-6 tournaments: for each count H, the isomorphism classes and how converse pairs them."""
import itertools, sys
sys.path.insert(0, '.')
from w45_hk import H_of
N = 20405
S = set(x for x in range(1, N + 1, 2) if x not in (7, 21))
def is_prime(n):
    if n < 2: return False
    d = 2
    while d * d <= n:
        if n % d == 0: return False
        d += 1
    return True
atoms = []
for x in sorted(S):
    if x == 1: continue
    red = False
    d = 3
    while d * d <= x:
        if x % d == 0 and d in S and (x // d) in S: red = True; break
        d += 2
    if not red: atoms.append(x)
pred = set(p for p in range(3, N + 1, 2) if is_prime(p) and p != 7) | set(7 * p for p in range(5, N // 7 + 1) if is_prime(p)) | {63, 343}
pred = set(a for a in pred if a <= N)
print("atoms of the count monoid up to", N, "match {primes != 7} U {7p : p >= 5} U {63, 343}:", set(atoms) == pred)
print("   atoms divisible by 7 (smallest):", sorted(a for a in atoms if a % 7 == 0)[:8])
# (2) order-6 classes
n = 6; pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
perms = list(itertools.permutations(range(n)))
def canon(adj):
    best = None
    for p in perms:
        key = tuple(adj[p[i]][p[j]] for i, j in pairs)
        if best is None or key < best: best = key
    return best
classes = {}
for bits in range(1 << len(pairs)):
    out = [0] * n; adj = [[0] * n for _ in range(n)]
    for t, (i, j) in enumerate(pairs):
        if bits >> t & 1: out[i] |= 1 << j; adj[i][j] = 1
        else: out[j] |= 1 << i; adj[j][i] = 1
    h = H_of(n, out)
    if h not in (33, 45, 37, 41, 25):      # restrict canonicalisation to a few counts of interest
        continue
    c = canon(adj)
    if c in classes: continue
    conv = [[adj[j][i] if i != j else 0 for j in range(n)] for i in range(n)]
    classes[c] = (h, canon(conv))
by_h = {}
for c, (h, cc) in classes.items():
    by_h.setdefault(h, []).append((c, cc))
for h in sorted(by_h):
    cl = by_h[h]; selfc = sum(1 for c, cc in cl if c == cc); pairs_ = (len(cl) - selfc) // 2
    print(f"order 6, H = {h}: {len(cl)} classes = {selfc} self-converse + {pairs_} converse pairs")
