#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S706
Generalizing THM-412 (density quantization of the unit-distance layer) from
2D lattices to ALL CM orders, and characterizing the density-quantum SPECTRUM
as the Euler-totient ladder M(N)=max{m: phi(m)|N}.

THM-412 (S703): for a 2D lattice with squared-norm form Q and w proper automorphs
(= #roots of unity in the imaginary quadratic order, w in {2,4,6}), the order-w
rotation group acts FREELY on each norm-D shell (D>0), so w | r_Q(D); the
unit-distance density r_Q(D)/2 is quantized in steps of w/2, capped at 3 in 2D.

CLAIM (THM-416, this session): the SAME free-action mechanism holds in EVERY CM
field. If K is a CM field of degree N=2d with w = #roots of unity mu(K), embed
the ring of integers O_K (or any CM order) in R^{2d} by the canonical
(Minkowski) embedding. Multiplication by a root of unity zeta is an ISOMETRY of
that embedding (|sigma(zeta)|=1 for every archimedean place) and fixes only 0
(zeta*g=g => (zeta-1)g=0 => g=0). Hence the cyclic group mu(K) of order w acts
freely on every positive norm-shell, so  w | r(D)  for all D>0, and the
unit-distance density quantum on the CM lattice is w/2.

Therefore the achievable density-quantum spectrum in dimension/degree N=2d is
   { w/2 : w = #mu(K), K a CM field with [K:Q] | N }  =  { m/2 : phi(m) | N, m even }.
The MAXIMAL quantum at degree N is  M(N)/2  where  M(N) = max{m : phi(m) | N}.
2D (N=2) caps at M(2)=6, quantum 3 (= triangular/Eisenstein, the "everything is
the triangle" constant). Higher quanta require a DIMENSION (rank) jump -> the
totient ladder is the precise "spectrum between triangular and CM".

VERIFICATION DONE HERE:
 (A) w | r(D) re-confirmed for 2D orders (disc -3,-4,-15,-23,-7,-8,-11).
 (B) w | r(D) verified for degree-4 cyclotomic lattices Z[zeta_5] (w=10),
     Z[zeta_8] (w=8), Z[zeta_12] (w=12), and degree-6 Z[zeta_7] (w=14),
     Z[zeta_9] (w=18) -- the genuinely NEW (>2D) cases.
     Gram matrix = Ramanujan-sum matrix c_m((i-j) mod m); norms are integers.
     Only norm-shells that are STABLE between box radius B and B+1 (hence
     fully captured = closed under the group action) are trusted.
 (C) The totient ladder M(N)/2 tabulated for N=2..24.
"""
import numpy as np
from math import gcd, isqrt
from itertools import product

def euler_phi(m):
    return sum(1 for a in range(1, m+1) if gcd(a, m) == 1)

def units_mod(m):
    return [a for a in range(1, m+1) if gcd(a, m) == 1]

def ramanujan_sum(k, m):
    """c_m(k) = sum_{a in (Z/m)^*} cos(2 pi a k / m), an integer."""
    s = sum(np.cos(2*np.pi*a*(k % m)/m) for a in units_mod(m))
    si = int(round(s))
    assert abs(s - si) < 1e-7, (k, m, s)
    return si

def cyclotomic_gram(m):
    """Integral Gram matrix of Z[zeta_m] in the canonical embedding.
    Basis 1, zeta, ..., zeta^{phi(m)-1}. G[i][j] = c_m((i-j) mod m)."""
    d = euler_phi(m)
    G = np.zeros((d, d), dtype=np.int64)
    for i in range(d):
        for j in range(d):
            G[i][j] = ramanujan_sum(i - j, m)
    return G

def quad_gram(a, b, c):
    """2D form a x^2 + b xy + c y^2 -> Gram [[2a,b],[b,2c]] (so c^T G c /1 = 2*Q)."""
    # use Q directly as the quadratic form value to keep integer norms = Q
    return np.array([[a, b/2],[b/2, c]], dtype=float)

def shell_counts(G, B, integer_form=None):
    """Enumerate lattice vectors c in [-B,B]^dim, return dict norm->count.
    norm = c^T G c (rounded to int)."""
    dim = G.shape[0]
    counts = {}
    for c in product(range(-B, B+1), repeat=dim):
        v = np.array(c, dtype=float)
        val = v @ G @ v
        iv = int(round(val))
        assert abs(val - iv) < 1e-6, (c, val)
        if iv <= 0:
            continue
        counts[iv] = counts.get(iv, 0) + 1
    return counts

def stable_shells(G, B):
    """Return shells PROVABLY fully captured by the box [-B,B]^dim.
    Rigorous capture bound: c^T G c = D >= lambda_min(G)*|c|^2, so any vector of
    squared-norm D has every coordinate |c_j| <= |c| <= sqrt(D/lambda_min). Hence
    D <= lambda_min*B^2 guarantees ALL norm-D vectors lie in the box -> the shell
    count is exact. (No truncation artifact; the group action cannot move a
    captured vector outside a fully-captured ball.)"""
    lam_min = float(np.linalg.eigvalsh(G).min())
    D_safe = lam_min * (B**2)
    cB = shell_counts(G, B)
    return {D: n for D, n in cB.items() if D <= D_safe + 1e-9}

def w_of_cyclotomic(m):
    """#roots of unity in Q(zeta_m) = m if m even else 2m."""
    return m if m % 2 == 0 else 2*m

print("="*78)
print("THM-416 verification: density quantization w | r(D) for CM orders")
print("="*78)

fail = 0

# ---- (A) 2D orders: re-confirm THM-412 ----
print("\n[A] 2D lattices (re-confirm THM-412):")
twoD = [
    ("disc -3 (triangular, Eisenstein)", (1,1,1), 6),
    ("disc -4 (square, Gaussian)",       (1,0,1), 4),
    ("disc -7",                          (1,1,2), 2),
    ("disc -8",                          (1,0,2), 2),
    ("disc -11",                         (1,1,3), 2),
    ("disc -15 (h=2, the density-5 rung)",(1,1,4), 2),
    ("disc -23 (h=3)",                   (1,1,6), 2),
]
for name,(a,b,c),w in twoD:
    G = quad_gram(a,b,c)
    st = stable_shells(G, 12)
    bad = [(D,n) for D,n in sorted(st.items()) if n % w != 0]
    status = "OK" if not bad else f"FAIL {bad[:3]}"
    if bad: fail += 1
    first = sorted(st.items())[:6]
    print(f"  {name:38s} w={w}  shells(first6)={first}  -> {status}")

# ---- (B) higher CM: cyclotomic lattices (THE NEW CASES) ----
print("\n[B] CM cyclotomic lattices Z[zeta_m] (degree phi(m) > 2 -- NEW):")
cyclo = [5, 8, 12, 7, 9]   # degrees 4,4,4,6,6
for m in cyclo:
    d = euler_phi(m)
    w = w_of_cyclotomic(m)
    G = cyclotomic_gram(m)
    B = 5 if d <= 4 else 3
    st = stable_shells(G, B)
    bad = [(D,n) for D,n in sorted(st.items()) if n % w != 0]
    status = "OK" if not bad else f"FAIL {bad[:5]}"
    if bad: fail += 1
    first = sorted(st.items())[:6]
    print(f"  Z[zeta_{m:<2}] deg={d}  w(#roots of unity)={w:<3} "
          f"quantum w/2={w//2:<2} stable-shells={len(st)} first6={first} -> {status}")

# ---- (C) the totient ladder M(N) = max{m: phi(m) | N} ----
print("\n[C] Density-quantum ladder  M(N)=max{m: phi(m)|N},  quantum=M(N)/2:")
print("    N=2d (dim/degree) | M(N)=max #roots of unity | quantum M(N)/2 | witness field")
def M_of(N, cap=10000):
    best = 1
    for m in range(1, cap):
        if N % euler_phi(m) == 0:
            best = m
    return best
ladder = []
for N in range(2, 26, 2):
    M = M_of(N)
    ladder.append((N, M, M//2))
    # witness: Q(zeta_M) has degree phi(M) | N
    print(f"      N={N:2d}  M(N)={M:3d}  quantum={M//2:3d}   Q(zeta_{M})  [phi={euler_phi(M)}]")

print("\n  Maximal-quantum sequence by degree 2,4,6,8,10,12,...:",
      [q for _,_,q in ladder])
print("  Maximal-w (roots of unity) sequence:               ",
      [M for _,M,_ in ladder])

# resonance flag (honest): does quantum 21 (=forbidden H) appear, and at which degree?
for N,M,q in ladder:
    if q == 21:
        print(f"\n  [RESONANCE, honest] quantum 21 first appears at degree N={N} "
              f"(w={M}=2*21, Q(zeta_{M})). 21 = forbidden tournament H-value (THM-079).")
    if q == 3:
        pass

# ---- (D) LRC corollary: quantum(Z[zeta_{2n-1}]) = 2n-1 = worry-set rosette order ----
print("\n[D] LRC corollary: quantum of the worry-set field Z[zeta_C], C=2n-1:")
print("    n  | C=2n-1 | deg=phi(C) | w=2C | quantum=w/2 | == C ?")
for n in range(3, 15):
    C = 2*n - 1
    d = euler_phi(C)
    w = w_of_cyclotomic(C)        # C odd -> 2C
    q = w // 2
    print(f"   {n:2d}  |  {C:3d}   |    {d:3d}     | {w:3d}  |    {q:3d}      | {'YES' if q==C else 'NO'}"
          + ("   <- n=14, C=27=3^3 (V* split)" if n == 14 else ""))

print("\n" + "="*78)
print(f"TOTAL FAILURES: {fail}")
print("PASS" if fail == 0 else "SOME FAILED")
print("="*78)
