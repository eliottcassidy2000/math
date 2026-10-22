# Audit2 (E): the doubly regular tournaments of drt_results.out, checked independently.
#  - S S^T = nI - J (equivalently S^2 = J - nI) for every listed code;
#  - nauty (labelg) canonical forms: which codes are isomorphic to the cyclotomic Paley tournament, to each other,
#    and to the converses; self-converse test;
#  - H by Held-Karp (n = 19) / inclusion-exclusion (n = 23) with hcount.c (passed as argv[1]); compared with the file;
#  - complete 4- and 5-vertex subtournament censuses (isomorphism type by brute-force canonical form), P vs D;
#  - spectral checks: x_k^2 = 1/n on 1^perp, rho = ((n+1)/(n-1))^((n-1)/2), two-loop y_pred identical.
import sys, itertools, subprocess
from collections import Counter
import numpy as np
sys.path.insert(0, '.')
from common import code_to_A, skew, paley_A
from e_predictor import two_loop_dense

HC = sys.argv[1]
H_IE = len(sys.argv) > 2 and sys.argv[2] == 'withIE'

def d6(A):
    n = A.shape[0]
    bits = [int(A[i, j]) for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        s += chr(63 + v)
    return s

def canon(As):
    inp = '\n'.join(d6(A) for A in As) + '\n'
    out = subprocess.run(['labelg', '-q'], input=inp, capture_output=True, text=True).stdout.split()
    assert len(out) == len(As)
    return out

def A_to_code(A):
    n = A.shape[0]
    return ''.join('1' if A[i, j] else '0' for i in range(n) for j in range(i + 1, n))

def sub_canon(A, X):
    k = len(X); best = None
    for p in itertools.permutations(X):
        c = tuple(int(A[p[i], p[j]]) for i in range(k) for j in range(k) if i != j)
        if best is None or c < best: best = c
    return best

def census(A, k):
    n = A.shape[0]
    return Counter(sub_canon(A, X) for X in itertools.combinations(range(n), k))

rows = []
for line in open('../drt_results.out'):
    if line.startswith('#') or not line.strip(): continue
    f = line.split()
    rows.append((int(f[0]), f[1], int(f[2]), f[-1]))

print('code checks:')
by_n = {}
for n, lab, Hfile, code in rows:
    A = code_to_A(code); S = skew(A)
    drt = (S @ S.T == n * np.eye(n, dtype=np.int64) - np.ones((n, n), dtype=np.int64)).all()
    reg = all(int(A[i].sum()) == (n - 1) // 2 for i in range(n))
    by_n.setdefault(n, []).append((lab, Hfile, A))
    print(f"  {lab}: n={n}, regular={reg}, S S^T = nI - J: {drt}")

for n, lst in sorted(by_n.items()):
    P = paley_A(n)
    As = [P] + [A for _, _, A in lst] + [A.T.copy() for _, _, A in lst]
    cf = canon(As)
    cP = cf[0]; m = len(lst)
    print(f"\norder {n}: nauty classes (index 0 = cyclotomic Paley; then the listed codes; then their converses)")
    for i, (lab, Hfile, A) in enumerate(lst):
        iso_P = cf[1 + i] == cP
        selfconv = cf[1 + i] == cf[1 + m + i]
        same = [lst[j][0] for j in range(m) if j != i and cf[1 + j] == cf[1 + i]]
        conv_of = [lst[j][0] for j in range(m) if cf[1 + m + j] == cf[1 + i] and j != i]
        print(f"  {lab}: isomorphic to Paley: {iso_P}; self-converse: {selfconv}; isomorphic to other listed: {same}; "
              f"isomorphic to the converse of: {conv_of}")
    # H
    for lab, Hfile, A in lst:
        if n <= 21:
            Hm = int(subprocess.run([HC, 'hk'], input=A_to_code(A) + '\n', capture_output=True, text=True).stdout)
        elif H_IE:
            Hm = int(subprocess.run([HC, 'ie'], input=A_to_code(A) + '\n', capture_output=True, text=True).stdout)
        else:
            Hm = None
        rel = (Hfile - lst[0][1]) / lst[0][1]
        print(f"  {lab}: H(file) = {Hfile}  H(recomputed) = {Hm}  {'MATCH' if Hm == Hfile else ('not recomputed' if Hm is None else 'MISMATCH')}  "
              f"(H - H(P))/H(P) = {rel:+.4e}, n^3 * that = {n ** 3 * rel:+.4f}", flush=True)
    # spectra and two-loop
    for lab, Hfile, A in lst:
        S = skew(A).astype(float)
        x2 = np.sort(np.abs(np.linalg.eigvals(S / n).imag) ** 2)
        t = two_loop_dense(skew(A))
        print(f"  {lab}: x_k^2 in [{x2[1]:.12f}, {x2[-1]:.12f}] (1/n = {1 / n:.12f}); rho = {t['rho']:.12f} "
              f"(formula {((n + 1) / (n - 1)) ** ((n - 1) / 2):.12f}); y_pred = {t['y_pred']:.12f}; n*dumbbell = {t['nDB']:+.1e}")
    # census
    for k in (4, 5):
        cs = [census(A, k) for _, _, A in lst]
        keys = sorted(set().union(*cs))
        print(f"  {k}-vertex census ({len(keys)} isomorphism types present):")
        for i, (lab, _, _) in enumerate(lst):
            print(f"     {lab}: {[cs[i].get(key, 0) for key in keys]}")
        print(f"     all equal: {all(c == cs[0] for c in cs)}", flush=True)
