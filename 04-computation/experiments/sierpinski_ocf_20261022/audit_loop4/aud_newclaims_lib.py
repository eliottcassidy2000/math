# AUDIT of the claims added to HYP-9371 during the audit (doubled-Paley DRT(23) class; (C1) at n = 31):
#  (a) rebuild the doubled-Paley DRT(23) independently (bordered P_11 -> skew-Hadamard K of order 12 ->
#      [[K, K], [-K^T, K^T]] -> normalise at every row, delete it), certify DRT and nauty class = drt23_doubled.json;
#  (b) own exact forest numbers (aud_fwdpoly, 2 primes): H, Delta N_8 vs 7168 Delta TT5, per-L shares, relative H;
#      the annealed TT5 law's out-of-sample error; own 5-census (true iso classes) direction;
#  (c) n = 31: twice-doubled Paley skew-Hadamard matrix of order 32 from P_7; own nauty classes, TT5 and 5-census
#      direction against P_31.
import sys, json, subprocess, itertools
from fractions import Fraction
from math import factorial
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE + '/audit_loop4')
from aud_task1 import F_mod, lagrange_coeffs, crt2, forest_N, P1, P2, code_adj, tt

def paley_S(q):
    sq = {(t * t) % q for t in range(1, q)}
    return np.array([[0 if i == j else (1 if (j - i) % q in sq else -1) for j in range(q)] for i in range(q)])
def skew_hadamard_from_drt(S):
    n = S.shape[0]; B = np.zeros((n + 1, n + 1), dtype=int); B[1:, 1:] = S; B[0, 1:] = 1; B[1:, 0] = -1
    K = np.eye(n + 1, dtype=int) + B
    assert (K @ K.T == (n + 1) * np.eye(n + 1, dtype=int)).all()
    return K
def doubled(K):
    m = K.shape[0]; K2 = np.block([[K, K], [-K.T, K.T]])
    assert (K2 @ K2.T == 2 * m * np.eye(2 * m, dtype=int)).all() and (K2 + K2.T == 2 * np.eye(2 * m, dtype=int)).all()
    return K2
def drts_from(K):
    m = K.shape[0]; B = K - np.eye(m, dtype=int); out = []
    for r in range(m):
        d = np.array([B[r, j] if j != r else 1 for j in range(m)])
        Br = d[:, None] * B * d[None, :]
        keep = [j for j in range(m) if j != r]
        D = Br[np.ix_(keep, keep)]
        assert (D @ D.T == (m - 1) * np.eye(m - 1, dtype=int) - 1).all()
        out.append(D)
    return out
def d6(S):
    n = S.shape[0]; bits = [1 if S[i, j] == 1 else 0 for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        s += chr(63 + v)
    return s
def canon(Ss):
    return subprocess.run(['labelg', '-q'], input='\n'.join(d6(S) for S in Ss) + '\n', capture_output=True, text=True).stdout.split()
def S_to_code(S):
    n = S.shape[0]; return ''.join('1' if S[i, j] == 1 else '0' for i in range(n) for j in range(i + 1, n))
# own 5-census by isomorphism class (as in aud_census5.py)
PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
perms = list(itertools.permutations(range(5)))
def adj10(c):
    A = [[0] * 5 for _ in range(5)]
    for b, (i, j) in enumerate(PAIRS):
        if c >> b & 1: A[i][j] = 1
        else: A[j][i] = 1
    return A
def code10(A): return sum(1 << b for b, (i, j) in enumerate(PAIRS) if A[i][j])
CAN = {c: min(code10([[adj10(c)[p[i]][p[j]] for j in range(5)] for i in range(5)]) for p in perms) for c in range(1024)}
CLS = sorted(set(CAN.values()))
def census5(S):
    n = S.shape[0]; cnt = {c: 0 for c in CLS}
    for X in itertools.combinations(range(n), 5):
        c = 0
        for b, (i, j) in enumerate(PAIRS):
            if S[X[i], X[j]] == 1: c |= 1 << b
        cnt[CAN[c]] += 1
    return cnt
TTc = 0   # canonical code of the transitive 5-tournament
VEC = None

