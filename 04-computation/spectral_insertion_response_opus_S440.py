#!/usr/bin/env python3
"""
THE SPECTRAL INSERTION-RESPONSE (opus-2026-07-20-S440) -- how the combinatorial a = vertex-insertion
(THM-1900) acts on the SKEW characteristic polynomial char_S, for GENERAL tournaments. Complements
kps THM-1880 (which does only the transitive tournament's a/b Chebyshev-Pell frame) by giving the
a-action on ALL tournaments' spectra, and bridges combinatorial-a (insert a vertex) to algebraic-a.

Skew matrix S: S_ij = +1 if i->j, -1 if j->i, 0 diag. Insert vertex u beating exactly P (sign
vector s: s_j=+1 iff j in P). New skew S' borders S with s. Claim (bordered/Schur):

    char_S(T + u_P)(x) = x * char_S(T)(x) + s^T adj(xI - S) s .

We verify this exactly, characterize the correction B_P(x)=s^T adj(xI-S) s, show source=sink
(=> the complement b preserves char_S), confirm eigenvalue INTERLACING (Cauchy on iS), and recover
kps's transitive recursion E_n = x*E_{n-1} + B as the special case.
"""
import itertools, sympy as sp

x = sp.symbols('x')

def skew(adj, n):
    S = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i, j] = 1 if adj[i][j] else -1
    return S

def edges_iter(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(pairs)):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: adj[i][j] = 1
            else:               adj[j][i] = 1
        yield adj

def char_skew(S, n):
    return sp.expand((x*sp.eye(n) - S).det())

def border(S, m, P):
    """add vertex m beating exactly P (subset of range(m))."""
    Sp = sp.zeros(m+1, m+1)
    Sp[:m, :m] = S
    for j in range(m):
        sj = 1 if j in P else -1
        Sp[m, j] = sj
        Sp[j, m] = -sj
    return Sp

print("SPECTRAL INSERTION-RESPONSE: char_S(T+u_P) = x*char_S(T) + s^T adj(xI-S) s")
print("="*74)
bordered_ok = True
source_eq_sink = True
deriv_hits = 0; deriv_tot = 0
examples = []
for m in range(2, 5):                       # T on m vertices; T+u on m+1
    adjs = list(edges_iter(m))
    for Ti, adj in enumerate(adjs):
        S = skew(adj, m)
        cS = char_skew(S, m)
        adjc = (x*sp.eye(m) - S).adjugate()
        cSprime = sp.diff(cS, x)             # char_S'(x) = sum_i char_{S-i}
        for P in (frozenset(p) for r in range(m+1) for p in itertools.combinations(range(m), r)):
            s = sp.Matrix([1 if j in P else -1 for j in range(m)])
            B = sp.expand((s.T * adjc * s)[0, 0])
            lhs = char_skew(border(S, m, P), m+1)
            if sp.expand(lhs - (x*cS + B)) != 0:
                bordered_ok = False
            # source (P=all) vs sink (P=empty)
        # source/sink check
        s_src = sp.Matrix([1]*m); s_snk = sp.Matrix([-1]*m)
        B_src = sp.expand((s_src.T*adjc*s_src)[0,0]); B_snk = sp.expand((s_snk.T*adjc*s_snk)[0,0])
        if B_src != B_snk: source_eq_sink = False
        # is the DIAGONAL part of the correction char_S'(x)? (s_i^2=1 always => diag sum = char_S')
        diag = sp.expand(sum(adjc[i, i] for i in range(m)))
        deriv_tot += 1
        if sp.expand(diag - cSprime) == 0: deriv_hits += 1
        if m == 3 and Ti in (0, 6):
            examples.append((m, Ti, sp.factor(cS), sp.factor(B_src)))
    print(f"  m={m}: {len(adjs)} tournaments checked")

print(f"\n bordered recursion char_S(T+u_P)=x*char_S+s^T adj s holds for ALL (T,P): {bordered_ok}")
print(f" source-insertion char_S == sink-insertion char_S (=> complement b preserves char_S): {source_eq_sink}")
print(f" diagonal cofactor sum == char_S'(x) (so B_P = char_S' + signed off-diagonal): {deriv_hits}/{deriv_tot}")
for m, Ti, cS, B in examples:
    print(f"   e.g. m={m} T#{Ti}: char_S={cS},  source correction B={B}")

# ---- recover kps's transitive recursion E_n = x*E_{n-1} + B_source ----
print("\n" + "="*74)
print("TRANSITIVE TOWER (kps THM-1880 recovered as the source-insertion special case)")
def transitive_skew(n):
    adj = [[1 if i < j else 0 for j in range(n)] for i in range(n)]  # i->j iff i<j
    return skew(adj, n)
prev = None
for n in range(1, 7):
    S = transitive_skew(n); E = char_skew(S, n)
    Eexpect = sp.expand(((x+1)**n + (x-1)**n)/2)
    tag = "= ((x+1)^n+(x-1)^n)/2  OK" if sp.expand(E - Eexpect) == 0 else "MISMATCH"
    # recursion from level n-1
    rec = ""
    if prev is not None:
        Sm = transitive_skew(n-1); adjc = (x*sp.eye(n-1)-Sm).adjugate()
        B = sp.expand((sp.Matrix([1]*(n-1)).T*adjc*sp.Matrix([1]*(n-1)))[0,0])
        rec = "  x*E_{n-1}+B = E_n OK" if sp.expand(x*prev + B - E) == 0 else "  REC-MISMATCH"
    print(f"  n={n}: E_n={sp.expand(E)}  {tag}{rec}")
    prev = E

# ---- interlacing (spectral insertion-response) ----
print("\n" + "="*74)
print("INTERLACING: eigenvalues of S' interlace those of S (Cauchy on iS) -- numeric check")
import numpy as np
def numeigs(adj, n):
    S = np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    return sorted((np.linalg.eigvals(S)/1j).real)   # imaginary parts (skew => purely imaginary)
random_ok = True
import random; random.seed(7)
for _ in range(2000):
    m = 5
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if random.random() < .5: adj[i][j]=1
            else: adj[j][i]=1
    P = set(k for k in range(m) if random.random()<.5)
    big = [row[:]+[0] for row in adj]+[[0]*(m+1)]
    for j in range(m):
        if j in P: big[m][j]=1
        else: big[j][m]=1
    mu = numeigs(adj, m); nu = numeigs(big, m+1)
    # nu (m+1 vals) must interlace mu (m vals): nu_k <= mu_k <= nu_{k+1}
    for k in range(m):
        if not (nu[k] - 1e-9 <= mu[k] <= nu[k+1] + 1e-9): random_ok = False
print(f"  eigenvalues of T+u_P interlace those of T on the imaginary axis (2000 random n=5): {random_ok}")
print("\n READING: a = vertex-insertion acts on char_S by x*(.)+B_P; B_P=char_S'+signed-offdiag;")
print(" source=sink => complement b is char_S-invariant; the new spectrum INTERLACES the old.")
print(" Transitive tower recovers kps THM-1880. This is the a-operation on ALL tournaments' spectra.")
