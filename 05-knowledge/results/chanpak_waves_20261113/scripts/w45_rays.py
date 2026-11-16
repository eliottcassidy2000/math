"""Search sandwich modules (|G_u|+|G_d| <= 4) for rays and digit pairs.
Ray: Hankel rank <= 2 of n -> h M^n w on the reachable space (n >= 1), i.e. H_{n+1} = lam H_n + charge.
Digit pair: modules M0, M1 whose joint Hankel matrix (words u, v with |v| >= 1) has rank 2."""
import random, itertools, sys
from fractions import Fraction
from w45_sandwich import regs, module_matrix, apply, HIDX, randtour, tours

def rank(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(9)) for j in range(9)] for i in range(9)]

def hrow(M):           # the functional h M as a row vector
    return [sum(M[i][j] for i in HIDX) for j in range(9)]

rng = random.Random(11)
samples = [regs(n, randtour(n, rng)) for n in (2, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7)]
# basis of the physical span
basis = []
for w in samples:
    if rank(basis + [w]) > len(basis): basis.append(w)
print("physical span dimension:", len(basis))

modules = []
for a in range(0, 5):
    for b in range(0, 5 - a):
        if a + b == 0: continue
        for ou in tours(a):
            for od in tours(b):
                for xu in range(1 << a):
                    for xd in range(1 << b):
                        modules.append(((a, tuple(ou), xu, b, tuple(od), xd), module_matrix(a, ou, xu, b, od, xd)))
print("modules:", len(modules))

def seq(M, w, L):
    out = []; v = w
    for _ in range(L):
        v = apply(M, v); out.append(sum(v[i] for i in HIDX))
    return out

rays = {}
for key, M in modules:
    rows = []
    for w in basis:
        s = seq(M, w, 8)
        rows.append(s)
    # Hankel rank of n -> f(n) for each start: test recurrence of order 2 common to all starts
    # columns: windows; simpler: rank of matrix [f_{n}, f_{n+1}, f_{n+2}] over n and starts
    trip = [[s[n], s[n + 1], s[n + 2]] for s in rows for n in range(0, 6)]
    if rank(trip) <= 2:
        # find lam: f_{n+2} = (lam+1) f_{n+1} - lam f_n
        lam = None
        for s in rows:
            for n in range(0, 6):
                d = s[n + 1] - s[n]
                if d != 0:
                    lam = Fraction(s[n + 2] - s[n + 1], d); break
            if lam is not None: break
        rays.setdefault(lam, []).append((key, M))
print("rays by lambda:", {str(k): len(v) for k, v in rays.items()})
# ---- digit pairs among exact rays (joint Hankel rank 2) ----
def f(Ms, word, w):
    v = w
    for ch in reversed(word):
        v = apply(Ms[ch], v)
    return sum(v[i] for i in HIDX)
def charges(M, lam):
    ds = []
    for w in basis:
        s = seq(M, w, 8)
        c = {s[n + 1] - lam * s[n] for n in range(0, 7)}
        if len(c) != 1: return None
        ds.append(c.pop())
    return ds
lefts = [''.join(p) for L in range(0, 4) for p in itertools.product('01', repeat=L)]
rights = [''.join(p) for L in range(1, 4) for p in itertools.product('01', repeat=L)]
for lam in (Fraction(2), Fraction(3)):
    good = [(key, M, ds) for key, M in rays.get(lam, []) for ds in [charges(M, lam)] if ds is not None]
    nz = [g for g in good if any(d != 0 for d in g[2])]
    found = 0
    for (k0, M0, d0), (k1, M1, d1) in itertools.combinations(good, 2):
        if d0 == d1: continue
        Ms = {'0': M0, '1': M1}
        rows = [[f(Ms, u + v, w) for v in rights for w in basis] for u in lefts]
        if rank(rows) <= 2: found += 1
    print(f"lambda={lam}: {len(good)} exact rays H -> lam H + charge, {len(nz)} with nonzero charge; digit pairs (joint Hankel rank 2): {found}")
