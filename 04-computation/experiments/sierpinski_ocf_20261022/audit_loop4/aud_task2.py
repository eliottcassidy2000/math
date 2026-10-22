# AUDIT task 2: the excess-3 engine (loop_engine4.series) on regular tournaments NOT used in excess3_validate.py,
# against the exact blow-up scaling audit2/b_mscaling.mscaling (exact rational X_3 coefficients):
#   * gentourng regular 9-tournament class 8 (not vertex-transitive; excess3_validate used R_9 and classes 2, 14), L <= 10;
#   * a random regular 11-tournament from an own 3-cycle-reversal chain with seed 2026 (excess3_validate: seed 5), L <= 8.
# Non-isomorphism to the tournaments of excess3_validate.py is certified with nauty labelg.
# Also: the series is recomputed with a second contour (K = 48, r = 0.7) to show contour independence.
import sys, time, subprocess
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE); sys.path.insert(0, HERE + '/audit2')
import loop_engine4 as LE4
from b_mscaling import mscaling
from common import rot_A, paley_A, skew, regular_classes, code_to_A

def d6(A):
    n = A.shape[0]; bits = [int(A[i][j]) for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        s += chr(63 + v)
    return s

def canon(As):
    return subprocess.run(['labelg', '-q'], input='\n'.join(d6(A) for A in As) + '\n', capture_output=True, text=True).stdout.split()

cl9 = regular_classes(9)
A9 = code_to_A(cl9[8])
# validated n = 9 tournaments: R_9, classes 2 and 14
c9 = canon([A9, rot_A(9), code_to_A(cl9[2]), code_to_A(cl9[14])])
print(f"n = 9 class 8: canonical form differs from R_9, class 2, class 14: {c9[0] not in c9[1:]}")

# own random regular 11-tournament: start from R_11, 3-cycle reversals, seed 2026
rng = np.random.default_rng(2026); n = 11
A = rot_A(11).copy()
for _ in range(6000):
    u, v, w = rng.choice(n, 3, replace=False)
    if A[u, v] and A[v, w] and A[w, u]:
        A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
assert all(A.sum(axis=1) == 5)
# the seed-5 tournament of excess3_validate.py
rng5 = np.random.default_rng(5)
B = np.array([[1 if 0 < (j - i) % n <= 5 else 0 for j in range(n)] for i in range(n)])
for _ in range(4000):
    u, v, w = rng5.choice(n, 3, replace=False)
    if B[u, v] and B[v, w] and B[w, u]: B[u, v] = B[v, w] = B[w, u] = 0; B[v, u] = B[w, v] = B[u, w] = 1
c11 = canon([A, B, paley_A(11)])
S11 = skew(A)
S4 = np.linalg.matrix_power(S11, 4)
print(f"n = 11 seed 2026: differs from the seed-5 tournament and P_11: {c11[0] not in c11[1:]}; diag(S^4) values {sorted(set(np.diag(S4)))} (non-constant => not vertex-transitive)")

worst = 0.0
for name, S, Lmax in (('regular 9-tournament, gentourng class 8', skew(A9), 10), ('random regular 11-tournament, seed 2026', S11, 8)):
    n = S.shape[0]; Sh = S / n
    t0 = time.time(); X = mscaling(S, Lmax, verbose=False); tm = time.time() - t0
    t0 = time.time(); co = LE4.series(Sh, 3, Lmax, K=32, r=0.5); te = time.time() - t0
    co2 = LE4.series(Sh, 3, Lmax, K=48, r=0.7)
    print(f"{name}: m-scaling {tm:.0f}s, engine {te:.0f}s")
    for L in range(2, Lmax + 1, 2):
        ex = X[L].get(3, 0)
        d = co[L] - float(ex); worst = max(worst, abs(d))
        print(f"   L={L:2d}: exact X_3^(L) = {str(ex):>40s} = {float(ex):+.13e}; engine (K=32,r=0.5) {co[L]:+.13e} diff {d:+.1e}; engine (K=48,r=0.7) diff {co2[L] - float(ex):+.1e}")
    # also the other excesses from the exact m-scaling at the top order, for orientation
    print(f"   exact excess profile at L = {Lmax}: " + ', '.join(f"X_{e}: {float(v):+.6e}" for e, v in sorted(X[Lmax].items())))
print(f"max |engine - exact| (K=32, r=0.5): {worst:.1e}")
