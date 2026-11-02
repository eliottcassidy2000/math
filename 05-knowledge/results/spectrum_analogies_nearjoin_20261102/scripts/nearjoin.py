"""Near-join family: T' = T -> M_k -> z (transitive chain M_k of k >= 1 vertices after T, then z) with the single arc z -> t0
reversed.  Claim: H(T') = h(T) + 2^k * B(T, t0), B(T,t0) = sum over S subset of V(T)-t0 of h(T[S]) * hs(T[S^c u t0], t0),
hs(U, t) = Hamiltonian paths of U starting at t.  (1) verify the formula by brute force on all tournaments of order <= 5,
k = 1, 2; (2) census of (h, B) over all tournaments of order <= 7 and all t0; (3) coverage of odd numbers <= LIMIT by
values h + 2^k B (k >= 1) and by the multiplicative closure with the strong values below 625 (3,5,9,11,13,15,17,19,23..)."""
import subprocess, itertools, sys
LIMIT = 100000
def decode(code, n):
    adj = [[0] * n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[idx] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj
def hpaths(adj, verts, start=None):
    """Hamiltonian paths of the induced subtournament on verts (list), optionally starting at start."""
    n = len(verts); pos = {v: i for i, v in enumerate(verts)}
    if n == 0: return 1 if start is None else 0
    f = [[0] * n for _ in range(1 << n)]
    for i, v in enumerate(verts):
        if start is None or v == start: f[1 << i][i] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for i in range(n):
            c = row[i]
            if not c: continue
            for j in range(n):
                if adj[verts[i]][verts[j]] and not (S >> j & 1): f[S | 1 << j][j] += c
    return sum(f[(1 << n) - 1])
def B(adj, n, t0):
    others = [v for v in range(n) if v != t0]; tot = 0
    for r in range(len(others) + 1):
        for S in itertools.combinations(others, r):
            Sc = [v for v in others if v not in S]
            tot += hpaths(adj, list(S)) * hpaths(adj, Sc + [t0], start=t0)
    return tot
def build_Tprime(adj, n, t0, k):
    N = n + k + 1; A = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n): A[i][j] = adj[i][j]
    for i in range(n):
        for j in range(n, N): A[i][j] = 1          # T -> M, T -> z
    for i in range(n, N):
        for j in range(i + 1, N): A[i][j] = 1      # chain M_1 -> ... -> M_k -> z
    A[t0][N - 1] = 0; A[N - 1][t0] = 1             # reversed arc z -> t0
    return A
gen = lambda n: subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
# (1) verification
bad = 0; checks = 0
for n in range(1, 6):
    codes = gen(n) if n >= 2 else ['']
    for code in codes:
        adj = decode(code, n) if n >= 2 else [[0]]
        h = hpaths(adj, list(range(n)))
        for t0 in range(n):
            b = B(adj, n, t0)
            for k in (1, 2):
                Tp = build_Tprime(adj, n, t0, k); N = n + k + 1
                H = hpaths(Tp, list(range(N))); checks += 1
                if H != h + 2 ** k * b: bad += 1
print('(1) formula H(T\') = h + 2^k B checked on %d cases (orders <= 5, k = 1, 2): mismatches %d' % (checks, bad))
# (2) census of (h, B)
pairs = {}
for n in range(2, 8):
    for code in gen(n):
        adj = decode(code, n); h = hpaths(adj, list(range(n)))
        for t0 in range(n):
            b = B(adj, n, t0); pairs.setdefault((h, b), (n, code, t0))
print('(2) distinct (h, B) pairs over orders 2..7: %d; smallest B per h for h <= 45:' % len(pairs),
      {h: min(b for (hh, b) in pairs if hh == h) for h in sorted({hh for hh, _ in pairs}) if h <= 45})
# (3) coverage of odd numbers
covered = {}
for (h, b), wit in pairs.items():
    k = 1
    while h + (b << k) <= LIMIT:
        v = h + (b << k)
        if v % 2 == 1 and v not in covered: covered[v] = (h, b, k, wit)
        k += 1
odd_missing = [v for v in range(1, LIMIT + 1, 2) if v not in covered]
print('(3) odd numbers <= %d NOT of the form h + 2^k B (k >= 1) with (h, B) from orders <= 7: %d; first 40: %s' % (LIMIT, len(odd_missing), odd_missing[:40]))
# multiplicative closure with small attained values (join of components): use the census values of orders <= 7 as atoms
atoms = sorted({v for v in covered if v <= 2000} | {3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45})
reach = set(covered)
for a in atoms:
    for v in list(reach):
        if a * v <= LIMIT: reach.add(a * v)
odd_missing2 = [v for v in range(1, LIMIT + 1, 2) if v not in reach]
print('    after multiplying by attained atoms: %d odd numbers <= %d still uncovered; first 40: %s' % (len(odd_missing2), LIMIT, odd_missing2[:40]))
print('    residues mod 8 of the uncovered (counts):', {r: sum(1 for v in odd_missing2 if v % 8 == r) for r in (1, 3, 5, 7)})
