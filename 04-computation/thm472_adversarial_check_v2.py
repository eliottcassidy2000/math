#!/usr/bin/env python3
"""Adversarial verification of THM-472 (tournament Hadamard ceiling), v2.

Fresh code, written independently of hadamard_det_census_macmini_s2.py
AND of thm472_adversarial_check.py (a parallel checker found on disk).

Checks, for n = 3..7 (exhaustive over all 2^C(n,2) labeled tournaments):
  - max det(I+S) and comparison with the claimed bound
    (n+1)^((n-1)/2) for odd n, n^(n/2) for even n
  - number of iso classes attaining the max
  - n=7: the 512-attainers are exactly the iso classes of the switching
    class of QR_7 (and there are exactly 6 of them)
  - n=7: every attainer, after the kernel-vector switching D = diag(w),
    has row sums 0, score (n-1)/2 = 3, common-out-neighbor count
    (n-3)/4 = 1 for every pair, and S'S'^T = nI - J
  - bordering: I + S'(QR_7) borders to an 8x8 skew-Hadamard matrix
  - n=4: equality det = 16 <=> SS^T = 3I (skew conference)
  - monotonicity of p -> (1 + c/p)^p (numeric spot check of the calculus)
"""
import numpy as np
from itertools import permutations, combinations
from fractions import Fraction

def pairs(n):
    return list(combinations(range(n), 2))

def all_tournaments_S(n, batch=1 << 15):
    """Yield batches of skew matrices S for all 2^C(n,2) tournaments."""
    P = pairs(n)
    m = len(P)
    total = 1 << m
    for start in range(0, total, batch):
        idx = np.arange(start, min(start + batch, total), dtype=np.uint64)
        K = len(idx)
        S = np.zeros((K, n, n), dtype=np.float64)
        for b, (i, j) in enumerate(P):
            bit = ((idx >> np.uint64(b)) & np.uint64(1)).astype(np.float64)
            s = 2.0 * bit - 1.0          # +1 means i->j
            S[:, i, j] = s
            S[:, j, i] = -s
        yield idx, S

def exact_det_int(M):
    """Exact integer determinant via Fraction Gaussian elimination."""
    n = M.shape[0]
    A = [[Fraction(int(M[i, j])) for j in range(n)] for i in range(n)]
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if f:
                A[r] = [A[r][k] - f * A[c][k] for k in range(n)]
    assert det.denominator == 1
    return int(det)

def mask_to_A(mask, n):
    P = pairs(n)
    A = np.zeros((n, n), dtype=np.int64)
    for b, (i, j) in enumerate(P):
        if (mask >> b) & 1:
            A[i, j] = 1
        else:
            A[j, i] = 1
    return A

def canon(A, perms_ix, pow2):
    """Canonical form: min over all vertex permutations of packed adjacency bits."""
    best = None
    for p in perms_ix:
        B = A[np.ix_(p, p)]
        v = int(B.reshape(-1) @ pow2)
        if best is None or v < best:
            best = v
    return best

def run_n(n, bound, expect_max, expect_classes=None):
    print(f"--- n = {n}:  claimed bound {bound}, expected max {expect_max}")
    I = np.eye(n)
    best = -1
    attainer_masks = []
    for idx, S in all_tournaments_S(n):
        dets = np.rint(np.linalg.det(I + S)).astype(np.int64)
        mx = int(dets.max())
        if mx > best:
            best = mx
            attainer_masks = []
        if mx == best:
            attainer_masks.extend(int(i) for i in idx[dets == best])
    print(f"    max det(I+S) over all 2^{n*(n-1)//2} tournaments: {best}"
          f"  (claim: {expect_max})  bound respected: {best <= bound}")
    assert best == expect_max, f"REFUTED: max at n={n} is {best}, not {expect_max}"
    # exact-arithmetic spot check on a few attainers
    for mask in attainer_masks[:3]:
        A = mask_to_A(mask, n)
        S = A - A.T
        d = exact_det_int(np.eye(n, dtype=np.int64) + S)
        assert d == best, f"float det wrong: exact {d} != {best}"
    # iso classes among attainers
    perms_ix = [np.array(p) for p in permutations(range(n))]
    pow2 = (1 << np.arange(n * n, dtype=np.int64)[::-1]).astype(object)
    canons = set()
    for mask in attainer_masks:
        canons.add(canon(mask_to_A(mask, n), perms_ix, pow2))
    print(f"    labeled attainers: {len(attainer_masks)}, iso classes: {len(canons)}")
    if expect_classes is not None:
        assert len(canons) == expect_classes, \
            f"REFUTED: {len(canons)} classes, claim {expect_classes}"
    return attainer_masks, canons, perms_ix, pow2

# ---------------- monotonicity of p log(1 + c/p) ----------------
print("=== check 0: monotonicity of p -> (1+c/p)^p ===")
ok = True
for c in [1.0, 3.0, 10.0, 21.0, 28.0]:
    vals = [p * np.log1p(c / p) for p in np.arange(0.5, 40, 0.25)]
    ok &= all(b > a for a, b in zip(vals, vals[1:]))
print(f"    strictly increasing for c in test set: {ok}")
assert ok

# ---------------- n = 3..7 census ----------------
print("=== check 1: exhaustive determinant census ===")
run_n(3, 4, 4, 2)
run_n(4, 16, 16, 2)
run_n(5, 36, 32, 8)
run_n(6, 216, 160, 4)
att7, canons7, perms7, pow27 = run_n(7, 512, 512, 6)

# ---------------- n=7: switching class of QR_7 ----------------
print("=== check 2: 512-attainers at n=7 == switching class of QR_7 ===")
n = 7
QR = {1, 2, 4}  # quadratic residues mod 7
A_qr = np.zeros((n, n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % 7 in QR:
            A_qr[i, j] = 1
S_qr = A_qr - A_qr.T
sw_canons = set()
for wmask in range(1 << n):
    w = np.array([1 if (wmask >> k) & 1 == 0 else -1 for k in range(n)])
    Ssw = S_qr * np.outer(w, w)
    Asw = (Ssw == 1).astype(np.int64)
    sw_canons.add(canon(Asw, perms7, pow27))
print(f"    iso classes in switching class of QR_7: {len(sw_canons)}")
print(f"    equals the 512-attainer class set: {sw_canons == canons7}")
assert sw_canons == canons7 and len(sw_canons) == 6

# ---------------- n=7: forced DRT structure of every attainer ----------------
print("=== check 3: kernel switching => row sums 0, scores 3, pair-dom 1, S'S'^T = 7I - J ===")
J = np.ones((7, 7))
checked = 0
for mask in att7:
    A = mask_to_A(mask, 7)
    S = (A - A.T).astype(np.float64)
    G = S @ S.T
    evals, evecs = np.linalg.eigh(G)
    assert abs(evals[0]) < 1e-8 and abs(evals[1] - 7) < 1e-8
    v = evecs[:, 0]
    assert np.allclose(v**2, 1.0 / 7, atol=1e-9), "v_i^2 != 1/n"
    w = np.sign(v)
    Sp = S * np.outer(w, w)
    assert np.allclose(Sp.sum(axis=1), 0), "row sums of S' not 0"
    assert np.allclose(Sp @ Sp.T, 7 * np.eye(7) - J), "S'S'^T != nI - J"
    Ap = (Sp > 0.5).astype(int)
    assert all(Ap[i].sum() == 3 for i in range(7)), "score != (n-1)/2"
    for i in range(7):
        for j in range(7):
            if i != j:
                assert int((Ap[i] & Ap[j]).sum()) == 1, "pair-domination != (n-3)/4"
    checked += 1
print(f"    all {checked} labeled attainers pass the forced-DRT checks")

# ---------------- bordering ----------------
print("=== check 4: explicit border of I+S'(QR_7) to skew-Hadamard order 8 ===")
Sp = S_qr.astype(np.int64)               # QR_7 already has row sums 0
assert (Sp.sum(axis=1) == 0).all()
H = np.zeros((8, 8), dtype=np.int64)
H[0, 0] = 1
H[0, 1:] = 1
H[1:, 0] = -1
H[1:, 1:] = np.eye(7, dtype=np.int64) + Sp
assert np.array_equal(H @ H.T, 8 * np.eye(8, dtype=np.int64)), "HH^T != 8I"
assert np.array_equal(H + H.T, 2 * np.eye(8, dtype=np.int64)), "not skew type"
assert set(np.unique(H)) == {-1, 1}
print("    H = [[1, 1^T], [-1, I+S']] is skew-Hadamard of order 8: PASS")

# ---------------- even n = 4: equality <=> SS^T = (n-1)I ----------------
print("=== check 5: n=4 equality forces SS^T = 3I (skew conference) ===")
for mask in range(1 << 6):
    A = mask_to_A(mask, 4)
    S = (A - A.T).astype(np.float64)
    d = int(round(np.linalg.det(np.eye(4) + S)))
    is_conf = np.allclose(S @ S.T, 3 * np.eye(4))
    assert (d == 16) == is_conf
print("    det = 16 <=> SS^T = 3I, for all 64 tournaments: PASS")

print()
print("ALL CHECKS PASSED")
