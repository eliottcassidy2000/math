"""
audit_S25_hn.py -- independent re-derivation for the audit of THM-4628 (session opus-2026-10-13-S25). Own code.

(1) Decomposition groups in (Z/55)^x: for every residue l coprime to 55, the order of l and whether -1 = 54 lies in <l>
    (conjugation-fixed primes); the fixed residues with their orders; the fixed rational primes <= 2500 with residue
    degree f and q = l^(f/2); checks: 2 is NOT fixed (order 20, 2^10 = 34 mod 55), no fixed prime has order 4 or 20,
    the ramified primes 5 and 11 are not fixed either (their decomposition groups do not contain (-1, -1)).
(2) F_(q^2) = F_q[sqrt d] for q in QS: norm-one torus N1 (size q+1, symmetric, no zero), the eigenvalues
    lambda_a = sum_{z in N1} psi(Tr(a z)) for EVERY a in F_(q^2), the check that lambda_a depends only on N(a), the identity
    lambda_a = -K_q(N(a)) with an independent Kloosterman sum K_q(c) = sum_{x in F_q^x} e((x + c/x)/q), lambda_min,
    Hoffman 1 + (q+1)/|lambda_min|, Weil floor 1 + (q+1)/(2 sqrt q); multiplicity of the eigenvalue q+1 (= 1 iff connected).
(3) Small q (7, 11, 23): the full adjacency spectrum of Cay(F_(q^2), N1) by numpy.linalg.eigvalsh versus the character sums
    (sign convention: least adjacency eigenvalue = -max_c K_q(c)).
(4) Weil floors at q = 19^5, 29^5, 79^5 and the least q with 1 + (q+1)/(2 sqrt q) > 6.
Run: python3 audit_S25_hn.py
"""
import math
import numpy as np

# ---------- (1) decomposition groups ----------
def order_mod(l, m):
    k, x = 1, l % m
    while x != 1:
        x = (x * l) % m
        k += 1
    return k

print('# (1) (Z/55)^x: residues l whose decomposition group <l> contains conjugation -1 = 54')
fixed_res = {}
for l in range(1, 55):
    if math.gcd(l, 55) != 1:
        continue
    o = order_mod(l, 55)
    has = any(pow(l, j, 55) == 54 for j in range(o))
    if has:
        fixed_res[l] = o
print('   fixed residues (l mod 55: order):', fixed_res)
print('   orders occurring among fixed residues:', sorted(set(fixed_res.values())))
print('   fixed residues of order 10 are exactly l = 4 mod 5 with l a primitive root mod 11:',
      all((l % 5 == 4 and order_mod(l, 11) == 10) for l, o in fixed_res.items() if o == 10),
      '; residues:', sorted(l for l, o in fixed_res.items() if o == 10))
print('   2 mod 55: order', order_mod(2, 55), ', 2^10 mod 55 =', pow(2, 10, 55), ', fixed:', 2 in fixed_res)
print('   any residue of order 20 fixed:', any(o == 20 for o in fixed_res.values()),
      '; any of order 4 fixed:', any(o == 4 for o in fixed_res.values()))
# ramified primes: D_5 = (Z/5)^x x <5 mod 11>, D_11 = <11 mod 5> x (Z/11)^x
sub5 = {pow(5, j, 11) for j in range(10)}
sub11 = {pow(11, j, 5) for j in range(4)}
print(f'   ramified: D_5 contains (-1,-1) iff 10 in <5 mod 11> = {sorted(sub5)}: {10 in sub5};  D_11 contains (-1,-1) iff 4 in <11 mod 5> = {sorted(sub11)}: {4 in sub11}')
primes = [p for p in range(2, 2501) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
fixed_primes = [(p, fixed_res[p % 55]) for p in primes if p % 55 in fixed_res]
print('   fixed primes <= 2500 (p, residue degree f = order):')
print('   ', [(p, f, f'q={p}^{f//2}') for p, f in fixed_primes])
deg2 = [p for p, f in fixed_primes if f == 2]
print('   residue-degree-two fixed primes <= 2000:', [p for p in deg2 if p <= 2000])
print('   least q over all fixed primes:', min(p ** (f // 2) for p, f in fixed_primes), '(109 from the degree-two list; 19^5 =', 19 ** 5, ')')

# ---------- (2) torus spectra ----------
def nonsquare(q):
    return next(d for d in range(2, q) if pow(d, (q - 1) // 2, q) == q - 1)

def kloosterman(q):
    x = np.arange(1, q)
    inv = np.array([pow(int(v), q - 2, q) for v in x])
    K = np.zeros(q)
    for c in range(1, q):
        K[c] = np.cos(2 * np.pi * ((x + c * inv) % q) / q).sum()
    return K

def torus_spectrum(q):
    d = nonsquare(q)
    a = np.arange(q)
    A, B = np.meshgrid(a, a, indexing='ij')
    norm = (A * A - d * B * B) % q
    tx, ty = np.nonzero(norm == 1)
    assert len(tx) == q + 1, 'torus size'
    S = set(zip(tx.tolist(), ty.tolist()))
    assert (0, 0) not in S and all(((-x) % q, (-y) % q) in S for x, y in S), 'symmetric, loopless'
    # lambda_(u,v) = sum_{(x,y) in N1} cos(2 pi * 2(u x + d v y)/q), for every (u, v)
    lam = np.zeros((q, q))
    for u in range(q):
        # shape (q, q+1): v index, torus index
        arg = (2 * (u * tx[None, :] + d * a[:, None] * ty[None, :])) % q
        lam[u] = np.cos(2 * np.pi * arg / q).sum(axis=1)
    return d, norm, lam

for q in (109, 439):
    d, norm, lam = torus_spectrum(q)
    K = kloosterman(q)
    spread = max(np.ptp(lam[norm == c]) for c in range(1, q))
    err = max(abs(lam[norm == c][0] + K[c]) for c in range(1, q))
    lmin = lam[norm != 0].min()
    mult_top = int(np.sum(np.abs(lam - (q + 1)) < 1e-6))
    hoff = 1 + (q + 1) / abs(lmin)
    weil = 1 + (q + 1) / (2 * math.sqrt(q))
    print(f'\n# (2) q = {q}, d = {d}: |N1| = {q+1}; max spread of lambda within a norm class {spread:.2e}; max |lambda + K(N(a))| = {err:.2e}')
    print(f'   lambda_min = {lmin:.4f}  (= -max_c K(c) = {-K[1:].max():.4f});  max |K| = {np.abs(K[1:]).max():.4f} <= 2 sqrt q = {2*math.sqrt(q):.4f}')
    print(f'   multiplicity of eigenvalue q+1 = {mult_top} (connected iff 1);  Hoffman bound = {hoff:.4f} -> chi >= {math.ceil(hoff - 1e-9)};  Weil floor = {weil:.4f}')
    print(f'   the trivial-norm class: lambda at a = 0 is {lam[0,0]:.1f} = degree q+1')

# ---------- (3) adjacency spectrum for small q ----------
print('\n# (3) adjacency spectrum of Cay(F_(q^2), N1) versus character sums, small q')
for q in (7, 11, 23):
    d, norm, lam = torus_spectrum(q)
    K = kloosterman(q)
    a = np.arange(q)
    A, B = np.meshgrid(a, a, indexing='ij')
    tx, ty = np.nonzero(norm == 1)
    n = q * q
    adj = np.zeros((n, n))
    idx = lambda x, y: x * q + y
    for x in range(q):
        for y in range(q):
            for sx, sy in zip(tx, ty):
                adj[idx(x, y), idx((x + sx) % q, (y + sy) % q)] = 1
    assert np.array_equal(adj, adj.T) and np.all(adj.sum(axis=1) == q + 1)
    ev = np.sort(np.linalg.eigvalsh(adj))
    cs = np.sort(lam.ravel())
    print(f'   q = {q}: regular of degree {q+1}; max |sorted adjacency eigenvalues - sorted character sums| = {np.abs(ev - cs).max():.2e}; '
          f'least adjacency eigenvalue {ev[0]:.6f} = -max K(c) = {-K[1:].max():.6f}; Hoffman {1 + (q+1)/abs(ev[0]):.4f}')

# ---------- (4) Weil floors ----------
print('\n# (4) Weil floors 1 + (q+1)/(2 sqrt q)')
for l in (19, 29, 79):
    q = l ** 5
    print(f'   q = {l}^5 = {q}: {1 + (q + 1) / (2 * math.sqrt(q)):.1f}')
least = next(q for q in range(2, 200) if 1 + (q + 1) / (2 * math.sqrt(q)) > 6)
print(f'   least integer q with floor > 6: {least} (floor at 97: {1 + 98/(2*math.sqrt(97)):.4f}; at 101: {1 + 102/(2*math.sqrt(101)):.4f}; at 109: {1 + 110/(2*math.sqrt(109)):.4f})')
print('done')
