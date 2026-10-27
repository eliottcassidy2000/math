"""checks_misc.py -- small independent checks for the THM-4823 / PROBLEM_CLASSIFICATION section 5 / sigma-four audit.
  1. Witnesses at n = 5 that #S4 is not a function of the score sequence (from gentourng classes).
  2. H, c3 and directed 4-cycles of the four 4-vertex types; H = 1 + 2 c3 on all 64 labelled 4-tournaments.
  3. Nuclear norm of 2(E_ij - E_ji); rank change of one arc reversal (rank is 2-Lipschitz per reversal).
  4. Paley skew matrices for primes q = 3 mod 4: circulant, singular values sqrt(q) (mult q-1) and 0, nuclear norm
     (q-1) sqrt(q) = the odd-order bound; scale comparison of the three edit lower bounds (Lipschitz, Frobenius, duality).
  5. Regular class counts at orders 3..11 (gentourng -d -D).
Usage: python3 checks_misc.py
"""
import subprocess, itertools, math
import numpy as np

def classes(n, extra=()):
    out = subprocess.run(["gentourng", "-q", *extra, str(n)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = np.zeros((n, n), dtype=int); k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': A[i, j] = 1
                else: A[j, i] = 1
                k += 1
        yield line, A

def s4count(A):
    n = len(A); c = 0
    for q in itertools.combinations(range(n), 4):
        sc = sorted(int(A[u, list(q)].sum()) for u in q)
        if sc == [1, 1, 2, 2]: c += 1
    return c

print("1. n = 5: score sequences with more than one #S4 value")
by = {}
for line, A in classes(5):
    s = tuple(int(x) for x in sorted(A.sum(1)))
    by.setdefault(s, []).append((s4count(A), line))
for s, L in sorted(by.items()):
    vals = sorted(set(v for v, _ in L))
    if len(vals) > 1:
        print("   scores", s, "#S4 values", vals, "classes", L)

print("2. the four 4-vertex types")
def ham_paths(A):
    n = len(A); return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i], p[i + 1]] for i in range(n - 1)))
def c3(A):
    n = len(A); return sum(1 for a, b, c in itertools.combinations(range(n), 3) if (A[a, b] and A[b, c] and A[c, a]) or (A[a, c] and A[c, b] and A[b, a]))
def c4(A):
    n = len(A); cnt = 0
    for p in itertools.permutations(range(n)):
        if p[0] == 0 and all(A[p[i], p[(i + 1) % 4]] for i in range(4)): cnt += 1
    return cnt
names = {(0, 1, 2, 3): 'TT4', (1, 1, 1, 3): 'W4', (0, 2, 2, 2): 'L4', (1, 1, 2, 2): 'S4'}
for line, A in classes(4):
    s = tuple(int(x) for x in sorted(A.sum(1)))
    print(f"   {names[s]}: scores {s}  H = {ham_paths(A)}  c3 = {c3(A)}  directed 4-cycles = {c4(A)}  tr(A^4)/4 = {np.trace(np.linalg.matrix_power(A, 4)) // 4}")
bad = 0
E = list(itertools.combinations(range(4), 2))
for mask in range(64):
    A = np.zeros((4, 4), dtype=int)
    for b, (i, j) in enumerate(E):
        if mask >> b & 1: A[i, j] = 1
        else: A[j, i] = 1
    if ham_paths(A) != 1 + 2 * c3(A): bad += 1
print(f"   H = 1 + 2 c3 on all 64 labelled 4-tournaments: failures {bad}")

print("3. one arc reversal")
n = 6; D = np.zeros((n, n)); D[1, 4] = 2; D[4, 1] = -2
print(f"   singular values of 2(E_ij - E_ji): {np.round(np.linalg.svd(D, compute_uv=False)[:3], 6)}; nuclear norm {np.linalg.svd(D, compute_uv=False).sum():.6f}; rank {np.linalg.matrix_rank(D)}")

print("4. Paley skew matrices")
def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
for q in [7, 11, 19, 23, 31, 43, 47, 59]:
    S = np.array([[legendre(j - i, q) for j in range(q)] for i in range(q)], dtype=float)
    circ = all(S[i, j] == S[(i + 1) % q, (j + 1) % q] for i in range(q) for j in range(q))
    skew = np.allclose(S, -S.T)
    sv = np.linalg.svd(S, compute_uv=False)
    nz = sv[sv > 1e-8]
    print(f"   q={q}: circulant {circ}, skew {skew}, rank {len(nz)}, nonzero singular values all sqrt(q): {np.allclose(nz, math.sqrt(q))}, "
          f"nuclear {sv.sum():.4f} = (q-1)sqrt(q) {(q - 1) * math.sqrt(q):.4f}; S S^T = qI - J: {np.allclose(S @ S.T, q * np.eye(q) - np.ones((q, q)))}")
print("   edit lower bounds against parses with nuclear norm U = (13/12) q log2 q (Session75 b = 3 constant):")
for q in [103, 1031, 10007, 100003]:
    nu = (q - 1) * math.sqrt(q); U = 13 / 12 * q * math.log2(q)
    lip = (nu - U) / 4                      # one reversal moves the nuclear norm by <= 4
    frob = (nu - U) ** 2 / (8 * (q - 1))    # ||A-B||_F^2 = 8m and rank <= q-1
    dual = (q * (q - 1) - math.sqrt(q) * U) / 4   # <A,B> = n(n-1) - 4m <= ||A||_op ||B||_*
    print(f"   q={q}: Lipschitz {lip:.3g} (~q^1.5/4 = {q ** 1.5 / 4:.3g}); Frobenius {frob:.3g} (~q^2/8 = {q * q / 8:.3g}); duality {dual:.3g} (~q^2/4 = {q * q / 4:.3g})")

print("5. regular classes (gentourng -d m -D m on 2m+1 vertices)")
for m in range(1, 6):
    nn = 2 * m + 1
    r = subprocess.run(["gentourng", "-u", f"-d{m}", f"-D{m}", str(nn)], capture_output=True, text=True)
    print(f"   n={nn}: {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else r.stdout.strip()}")

print("6. the TT4-free classes at orders 5, 6, 7 (is the order-7 one doubly regular, i.e. P7; is the order-6 one P7 - v?)")
def tt4free(A):
    n = len(A)
    for q in itertools.combinations(range(n), 4):
        if sorted(int(A[u, list(q)].sum()) for u in q) == [0, 1, 2, 3]: return False
    return True
def canon(A):
    n = len(A); best = None
    for p in itertools.permutations(range(n)):
        key = tuple(int(A[p[i], p[j]]) for i in range(n) for j in range(n))
        if best is None or key < best: best = key
    return best
P7 = np.array([[1 if (j - i) % 7 in (1, 2, 4) else 0 for j in range(7)] for i in range(7)])
P7m = P7[1:, 1:]
for nn in (5, 6, 7):
    found = [A for _, A in classes(nn) if tt4free(A)]
    for A in found:
        sc = sorted(int(x) for x in A.sum(1))
        msg = f"   n={nn}: TT4-free class, scores {sc}"
        if nn == 7:
            co = {int((A[u] & A[v]).sum()) for u in range(7) for v in range(7) if u != v}
            msg += f", common out-neighbours of distinct pairs {co} (doubly regular iff {{1}}), isomorphic to P7: {canon(A) == canon(P7)}"
        if nn == 6:
            msg += f", isomorphic to P7 - v: {canon(A) == canon(P7m)}"
        print(msg)

print("7. skew-Hadamard doubling D = [[H, H], [-H^T, H^T]] (sigma-four section 2(b) extension)")
H = np.array([[1, 1], [-1, 1]])
for k in range(2, 8):
    H = np.block([[H, H], [-H.T, H.T]]); N = len(H); S = H - np.eye(N)
    nu = np.linalg.svd(S, compute_uv=False).sum()
    print(f"   order {N}: H + H^T = 2I {np.allclose(H + H.T, 2 * np.eye(N))}, H H^T = N I {np.allclose(H @ H.T, N * np.eye(N))}, "
          f"S = H - I is a tournament skew matrix {set(np.unique(S[~np.eye(N, dtype=bool)])) == {-1.0, 1.0}}, nuclear norm {nu:.4f} = N sqrt(N-1) {N * math.sqrt(N - 1):.4f}")

print("8. Goldbach row: every even N is a sum of two units modulo m (all m <= 300)")
bad = []
for m in range(2, 301):
    U = [u for u in range(m) if math.gcd(u, m) == 1]
    sums = {(u + v) % m for u in U for v in U}
    evens = {(2 * t) % m for t in range(m)}
    if not evens <= sums: bad.append(m)
print(f"   failures: {bad}")
print("9. Sloane row: 2 is a primitive root mod 3^k (k <= 12), so 2^n mod 3^k is equidistributed over the units; units whose k low ternary digits avoid 0 number 2^k, a share (2/3)^(k-1)")
for k in range(1, 13):
    M = 3 ** k; o = 1; x = 2
    while x != 1: x = x * 2 % M; o += 1
    nz = sum(1 for r in range(M) if r % 3 and all((r // 3 ** i) % 3 for i in range(k)))
    print(f"   k={k}: ord_M(2) = {o} = phi = {2 * 3 ** (k - 1)}: {o == 2 * 3 ** (k - 1)}; zero-free units {nz} = 2^k {nz == 2 ** k}" if k <= 9 else f"   k={k}: ord_M(2) = {o} = phi: {o == 2 * 3 ** (k - 1)}")
print("10. labelled W4-free counts from Theorem C: W_lab(n) = n! + sum_{k=3..n} C(n,k)(n-k)!(2^(k-1)(k-1)! - k!)")
print("   ", [math.factorial(n) + sum(math.comb(n, k) * math.factorial(n - k) * (2 ** (k - 1) * math.factorial(k - 1) - math.factorial(k)) for k in range(3, n + 1)) for n in range(1, 11)])
