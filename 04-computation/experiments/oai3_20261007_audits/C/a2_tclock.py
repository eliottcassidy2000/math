# Independent audit implementation of the THM-4569 Terras-clock chain.
# State (j, c) with c an exact Fraction in Z[1/3]; transitions coded straight from the THM table:
#  (p,e)=(0,0): (j, c/2); (0,1): (j+1, (3c+1)/2); (1,0): (j, (3c+1-3^j)/2); (1,1): (j-1, (c-3^(j-1))/2)
# where e = c mod 2 computed 2-adically (c = u/3^k, 3^k odd => parity of u).
import random, sys
from fractions import Fraction as F
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def pw3(k): return F(3)**k
def par(c):  # c in Z[1/3]: c = u / 3^k  -> c mod 2 = u mod 2
    assert c.denominator & (c.denominator - 1) != 0 or c.denominator == 1
    d = c.denominator
    while d % 3 == 0: d //= 3
    assert d == 1, c
    return c.numerator & 1
def step(j, c, p):
    e = par(c)
    if p == 0 and e == 0: return j, c / 2
    if p == 0 and e == 1: return j + 1, (3*c + 1) / 2
    if p == 1 and e == 0: return j, (3*c + 1 - pw3(j)) / 2
    return j - 1, (c - pw3(j - 1)) / 2

# ---- (A) transitions vs direct 2-adic iteration of x and y (mod 2^K)
def T(z): return z >> 1 if z % 2 == 0 else (3*z + 1) >> 1
rng = random.Random(20261007)
K = 500; steps_ok = 0; bad = 0; up = down = flat = 0
for trial in range(250):
    j0 = rng.randint(-5, 5); c0 = F(rng.randint(-10**5, 10**5), 3**rng.randint(0, 4))
    M = 1 << K
    y = rng.getrandbits(K)
    inv3 = pow(3, -1, M)
    def emb(q, m):  # Fraction in Z[1/3] -> Z/2^m
        return (q.numerator * pow(q.denominator, -1, 1 << m)) % (1 << m)
    x = (emb(pw3(j0), K) * y + emb(c0, K)) % M
    j, c = j0, c0
    for s in range(400):
        m = K - s
        if (x - emb(pw3(j), m) * y - emb(c, m)) % (1 << m) != 0: bad += 1; break
        if (x & 1) != ((y & 1) ^ par(c)): bad += 1; break
        jn, cn = step(j, c, y & 1)
        if jn > j: up += 1
        elif jn < j: down += 1
        else: flat += 1
        j, c = jn, cn
        x, y = T(x) % (1 << (m-1)), T(y) % (1 << (m-1))
        steps_ok += 1
print(f"(A) relation x_s = 3^j_s y_s + c_s and parity rule checked on {steps_ok} steps; failures {bad}; j up {up} down {down} flat {flat}")
# merge precursors and forced prefix of S19 lag-1 (q in 5+16Z_2 => parities 1,0,0,0)
st = (1, F(2))
for p in (1, 0, 0, 0): st = step(*st, p)
print("(A) S19 lag-1 start (1,2) after forced prefix 1000 ->", st)
q = 5 + 16*12345; ps = []
for _ in range(4): ps.append(q & 1); q = T(q)
print("    parities of q=5+16k for 4 steps:", ps)
# x = y+1: first step
print("    (0,1) after p=0:", step(0, F(1), 0), " after p=1:", step(0, F(1), 1))

# ---- (B) value iteration on the box B(J,C) = {|j|<=J, |c 3^max(0,-j)| <= C 3^|j|}
def box_states(J, C):
    S = []
    for j in range(-J, J + 1):
        d = max(0, -j); R = C * 3**abs(j)
        for N in range(-R, R + 1):
            S.append((j, F(N, 3**d)))
    return S
def inbox(j, c, J, C):
    if abs(j) > J: return False
    N = c * pw3(max(0, -j))
    assert N.denominator == 1
    return abs(N.numerator) <= C * 3**abs(j)
def vi(J, C, sweeps, report):
    S = box_states(J, C); idx = {s: i for i, s in enumerate(S)}; n = len(S)
    OUT = n
    nxt = np.full((n, 2), OUT, dtype=np.int64)
    for i, (j, c) in enumerate(S):
        for p in (0, 1):
            jn, cn = step(j, c, p)
            if inbox(jn, cn, J, C): nxt[i, p] = idx[(jn, cn)]
    z = idx[(0, F(0))]
    V = np.zeros(n + 1)
    for it in range(sweeps):
        Vn = 0.5 * (V[nxt[:, 0]] + V[nxt[:, 1]]); Vn[z] = 1.0
        V[:n] = Vn
    # exact absorption-within-box probabilities via sparse linear solve (limit of VI)
    rows, cols, vals = [], [], []
    b = np.zeros(n)
    for i in range(n):
        rows.append(i); cols.append(i); vals.append(1.0)
        if i == z: b[i] = 1.0; continue
        for p in (0, 1):
            k = nxt[i, p]
            if k == OUT: continue
            rows.append(i); cols.append(k); vals.append(-0.5)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    h = spla.spsolve(A.tocsc(), b)
    small = [idx[(j, F(N, 3**max(0, -j)))] for j in (-1, 0, 1) for N in range(-10*3**abs(j), 10*3**abs(j) + 1) if (j, N) != (0, 0)]
    out = {k: (V[idx[s]], h[idx[s]]) for k, s in report.items()}
    out['min B(1,10)'] = (V[small].min(), h[small].min())
    print(f"(B) box B({J},{C}): {n} states, {sweeps} sweeps:  " + ";  ".join(f"{k}: VI {a:.6f} / solve {b_:.6f}" for k, (a, b_) in out.items()))
rep = {'(0,1)': (0, F(1)), '(2,1)': (2, F(1)), '(1,1)': (1, F(1))}
vi(3, 30, 4000, rep)
vi(4, 100, 4000, rep)
