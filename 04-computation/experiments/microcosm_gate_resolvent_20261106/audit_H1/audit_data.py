# audit_H1: FINITE-EXACT spot checks: phi minimum over strong classes n = 5..9 (gentourng -c), HC = 7 / HC = 21 class counts and gates, gated spectra H(Gamma T) for cores n <= 7, note section 4 face formula and its determination by L P_T.
import itertools
import math
import random
import subprocess
import sys
from fractions import Fraction
import numpy as np
from audit_tournament_lib import (GENTOURNG, from_upper, classes, rand_tour, transitive, scores, gate, phi, cover_poly)

HP = sys.argv[1] if len(sys.argv) > 1 else "./hp_count"
rng = random.Random(211)
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def upper(A):
    n = len(A)
    return "".join("1" if A[i][j] else "0" for i in range(n) for j in range(i + 1, n))


def run_counts(lines, cyc=False):
    inp = "".join(lines)
    args = [HP] + (["-c"] if cyc else [])
    out = subprocess.run(args, input=inp, capture_output=True, text=True, check=True).stdout.split()
    return [int(x) for x in out]


def has_gate(A):
    n = len(A)
    sc = scores(A)
    indeg = [n - 1 - x for x in sc]
    return any(A[a][b] and sc[a] == 1 and indeg[b] == 1 for a in range(n) for b in range(n))


# ---------- (a) phi minimum over strong tournaments ----------
print("(a) min / max of phi over strong classes")
for n in range(5, 10):
    strs = subprocess.run([GENTOURNG, "-q", "-c", str(n)], capture_output=True, text=True, check=True).stdout.split()
    m = len(strs)
    iu = [(i, j) for i in range(n) for j in range(i + 1, n)]
    bits = np.array([[1 if ch == '1' else -1 for ch in st] for st in strs], dtype=float)   # M[i][j] for i<j
    Ms = np.zeros((m, n, n))
    for k, (i, j) in enumerate(iu):
        Ms[:, i, j] = bits[:, k]
        Ms[:, j, i] = -bits[:, k]
    sys_ = np.eye(n)[None, :, :] - Ms
    x = np.linalg.solve(sys_, np.ones((m, n, 1)))[:, :, 0]
    ph = x.sum(axis=1)
    lo = ph.min()
    cand = np.where(ph < lo + 1e-7)[0]
    exact = [phi(from_upper(n, strs[c])) for c in cand]
    mn = min(exact)
    mins = [c for c, e in zip(cand, exact) if e == mn]
    sc = [sorted(scores(from_upper(n, strs[c]))) for c in mins]
    print(f"  n = {n}: {m} strong classes; min phi = {mn} attained by {len(mins)} class(es), scores {sc}; "
          f"max phi = {ph.max():.6f}")

# ---------- (b) HC = 7 and HC = 21 ----------
print("(b) Hamiltonian-cycle counts 7 and 21 among strong classes")
for n in range(3, 9):
    strs = subprocess.run([GENTOURNG, "-q", "-c", str(n)], capture_output=True, text=True, check=True).stdout.split()
    hc = run_counts([f"{n} {st}\n" for st in strs], cyc=True)
    c7 = [st for st, h in zip(strs, hc) if h == 7]
    c21 = [st for st, h in zip(strs, hc) if h == 21]
    g7 = sum(has_gate(from_upper(n, st)) for st in c7)
    g21 = sum(has_gate(from_upper(n, st)) for st in c21)
    # sanity: HC(Gamma T) = H(T) on gated classes
    gated = [st for st in strs if has_gate(from_upper(n, st))]
    print(f"  n = {n}: {len(strs)} strong classes; HC = 7: {len(c7)} classes ({g7} gated); "
          f"HC = 21: {len(c21)} classes ({g21} gated); gated strong classes: {len(gated)}")

# ---------- (c) gated spectra ----------
print("(c) gated spectra H(Gamma T) over all core classes")
for n in range(1, 8):
    cl = classes(n)
    vals = sorted(set(run_counts([f"{n + 2} {upper(gate(A))}\n" for A in cl])))
    second = vals[1] if len(vals) > 1 else None
    check(vals[0] == 2 ** n + 1, f"min gated n={n}")
    if n >= 3:
        check(second == 3 + 3 * 2 ** (n - 1), f"second gated n={n}")
    print(f"  n = {n}: {len(vals)} distinct; first values {vals[:6]}; 2^n+1 = {2**n+1}; 3+3*2^(n-1) = {3 + 3 * 2 ** (n - 1) if n >= 1 else None}")
small = sorted(v for n in range(1, 8) for v in set(run_counts([f"{n + 2} {upper(gate(A))}\n" for A in classes(n)])) if v <= 40)
print("  gated values <= 40 over core orders 1..7:", sorted(set(small)))


# ---------- (d) face formula ----------
def desc_poly(A, verts):
    """backward-adjacency (descent) distribution over all orderings of verts: list of counts by k."""
    k = len(verts)
    if k == 0:
        return [1]
    out = [0] * k
    for p in itertools.permutations(verts):
        b = sum(1 for i in range(k - 1) if A[p[i + 1]][p[i]])
        out[b] += 1
    return out


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def face_formula(A):
    n = len(A)
    tot = [0]
    cache = {}
    def dp(S):
        if S not in cache:
            cache[S] = desc_poly(A, [v for v in range(n) if (S >> v) & 1])
        return cache[S]
    for lab in itertools.product(range(3), repeat=n):
        SA = sum(1 << v for v in range(n) if lab[v] == 0)
        SB = sum(1 << v for v in range(n) if lab[v] == 1)
        SC = sum(1 << v for v in range(n) if lab[v] == 2)
        base = pmul(pmul(dp(SA), dp(SB)), dp(SC))
        e1 = (SA != 0) + (SB == 0) + (SC != 0)
        e2 = 2 * (SB != 0)
        jw = [0] * 4
        jw[e1] += 1
        jw[e2] += 1
        tot = padd(tot, pmul(base, jw))
    return tot


def desc_from_cover(co, n):
    """invert P(q) = sum_f m_f C(q+f, n) (Session210 (7)) for m_f, return descent polynomial a_k = m_(n-1-k)."""
    P = lambda q: sum(c * q ** e for e, c in enumerate(co))
    # evaluations at q = 1..n are triangular in reversed f order
    m = [0] * n
    for q in range(1, n + 1):
        f = n - q   # C(q+f, n) = 1 at q + f = n; larger f already known
        val = P(q) - sum(m[g] * math.comb(q + g, n) for g in range(f + 1, n))
        m[f] = val
    return [m[n - 1 - k] for k in range(n)]


def L_poly(co):
    out = [0] * (len(co) + 2)
    for e, c in enumerate(co):
        if c:
            out[e + 2] += c
            for i in range(1, e + 1, 2):
                out[e - i + 1] += 2 * math.comb(e, i) * c
    return out


nface = 0
for n in range(1, 7):
    for trial in range(3 if n < 6 else 2):
        A = rand_tour(n, rng)
        G = gate(A)
        brute = desc_poly(G, list(range(n + 2)))
        ff = face_formula(A)
        while len(ff) < len(brute):
            ff.append(0)
        ff = ff[:len(brute)] if all(x == 0 for x in ff[len(brute):]) else ff
        check(ff == brute, f"face formula n={n}")
        a = desc_poly(A, list(range(n)))
        a1 = a[1] if len(a) > 1 else 0
        check(brute[0] == (n + 2) * a[0] + a1, "H(Gamma T) = (n+2) a0 + a1")
        # determination by Session210's operator: A_(Gamma T) from L P_T alone
        co = cover_poly(A)
        check(desc_from_cover(L_poly(co), n + 2) == brute, f"A_Gamma from L P_T n={n}")
        nface += 1
print(f"(d) face formula A_Gamma(T)(x) = sum_(A,B,C) A_T[A]A_T[B]A_T[C](x^([A]+[B=0]+[C]) + x^(2[B])) vs brute force, "
      f"H(Gamma T) = (n+2)a_0 + a_1, and A_Gamma(T) recovered from L P_T by the binomial transform: {nface} random cores n = 1..6")
print("FAILURES:", len(fails))
