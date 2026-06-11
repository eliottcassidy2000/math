#!/usr/bin/env python3
"""
thm472_adversarial_check.py — adversarial verification of THM-472
(tournament Hadamard ceiling). FRESH code, written independently of
hadamard_det_census_macmini_s2.py.

Checks:
  A. Brute force over ALL labeled tournaments n=3..7:
     max det(I+S) = 4, 16, 32, 160, 512 (exact integer recheck of attainers).
  B. Iso-class counts of attainers (n=5: 8, n=6: 4, n=7: 6).
  C. n=7 attainers == canonical forms of the switching class of QR_7.
  D. Equality analysis: every n=7 attainer has SS^T = 7I - ww^T, w in {+-1}^7,
     and D S D is doubly regular (scores 3, pair-domination 1).
  E. QR_7 spectrum {0} u {+-i sqrt7 x3}; bordering to skew-Hadamard order 8.
  F. n=6 attainers: det(S) = 81 (|Pf| = 9).
  G. n=4 attainers: SS^T = 3I (skew conference); n=8 via gentourng:
     max 4096, 4 classes, all SS^T = 7I.
  H. Monotonicity sanity for p log(1+c/p).

Provenance: adversarial-check subagent 2026-06-11 (proof check of THM-472).
"""
import itertools, subprocess, sys
from fractions import Fraction
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

def pairs_of(n):
    return list(itertools.combinations(range(n), 2))

def S_from_mask(mask, n, PAIRS):
    S = np.zeros((n, n), dtype=np.int64)
    for k, (i, j) in enumerate(PAIRS):
        v = 1 if (mask >> k) & 1 else -1
        S[i, j] = v
        S[j, i] = -v
    return S

def mask_from_S(S, n, PAIRS):
    m = 0
    for k, (i, j) in enumerate(PAIRS):
        if S[i, j] == 1:
            m |= (1 << k)
    return m

def exact_det(M):
    """Bareiss fraction-free determinant, exact integers."""
    n = len(M)
    A = [[int(x) for x in row] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            piv = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def all_dets(n):
    """det(I+S) for all 2^C(n,2) tournaments, batched numpy (rounded)."""
    PAIRS = pairs_of(n)
    m = len(PAIRS)
    total = 1 << m
    best = -1
    attainers = []
    B = 1 << 17
    I = np.eye(n)
    for start in range(0, total, B):
        cnt = min(B, total - start)
        masks = np.arange(start, start + cnt, dtype=np.int64)
        S = np.zeros((cnt, n, n))
        for k, (i, j) in enumerate(PAIRS):
            bit = ((masks >> k) & 1) * 2 - 1
            S[:, i, j] = bit
            S[:, j, i] = -bit
        d = np.rint(np.linalg.det(I + S)).astype(np.int64)
        mx = int(d.max())
        if mx > best:
            best = mx
            attainers = [int(x) for x in masks[d == mx]]
        elif mx == best:
            attainers.extend(int(x) for x in masks[d == mx])
    return best, attainers

def canon(mask, n, PAIRS, perms):
    """Canonical form: min over relabelings of the pair-bitmask."""
    S = S_from_mask(mask, n, PAIRS)
    best = None
    for p in perms:
        Sp = S[np.ix_(p, p)]
        mm = mask_from_S(Sp, n, PAIRS)
        if best is None or mm < best:
            best = mm
    return best

def canon_set(masks, n):
    PAIRS = pairs_of(n)
    perms = [list(p) for p in itertools.permutations(range(n))]
    return set(canon(mk, n, PAIRS, perms) for mk in masks)

claimed_max = {3: 4, 4: 16, 5: 32, 6: 160, 7: 512}
claimed_nclass = {3: 2, 4: 2, 5: 8, 6: 4, 7: 6}
results = {}

print("== A/B: brute force max det(I+S), all labeled tournaments ==")
for n in range(3, 8):
    PAIRS = pairs_of(n)
    best, att = all_dets(n)
    # exact integer recheck of every attainer
    ok_exact = all(exact_det(np.eye(n, dtype=np.int64) + S_from_mask(mk, n, PAIRS)) == best
                   for mk in att[:200])  # cap exact recheck at 200 reps
    nb_odd = (n + 1) ** ((n - 1) // 2) if n % 2 == 1 else n ** (n // 2)
    cs = canon_set(att, n)
    results[n] = (best, att, cs)
    print(f"n={n}: max={best} (claimed {claimed_max[n]}, bound {nb_odd}), "
          f"labeled attainers={len(att)}, iso classes={len(cs)} "
          f"(claimed {claimed_nclass[n]}), exact-recheck={'OK' if ok_exact else 'FAIL'}")
    assert best == claimed_max[n], f"MAX MISMATCH at n={n}"
    assert len(cs) == claimed_nclass[n], f"CLASS COUNT MISMATCH at n={n}"
    assert best <= nb_odd, f"BOUND VIOLATED at n={n}"

print("\n== C: n=7 attainers vs switching class of QR_7 ==")
n = 7
PAIRS7 = pairs_of(7)
QR = {1, 2, 4}
Sqr = np.zeros((7, 7), dtype=np.int64)
for i in range(7):
    for j in range(7):
        if i != j:
            Sqr[i, j] = 1 if ((j - i) % 7) in QR else -1
assert (Sqr.T == -Sqr).all()
sw_masks = []
for bits in range(64):  # w_0 = +1 wlog (w and -w give same DSD)
    w = np.array([1] + [1 if (bits >> t) & 1 else -1 for t in range(6)])
    D = np.diag(w)
    sw_masks.append(mask_from_S(D @ Sqr @ D, 7, PAIRS7))
cs_qr = canon_set(sw_masks, 7)
cs_att = results[7][2]
print(f"switching class of QR_7: {len(set(sw_masks))} labeled switchings, "
      f"{len(cs_qr)} iso classes; equal to attainer set: {cs_qr == cs_att}")
assert cs_qr == cs_att, "ATTAINERS != SWITCHING CLASS OF QR_7"

print("\n== D: equality analysis on every n=7 attainer ==")
att7 = results[7][1]
ok = True
for mk in att7:
    S = S_from_mask(mk, 7, PAIRS7)
    G = S @ S.T
    # SS^T must equal 7I - ww^T: so ww^T = 7I - G; diag must be 1, rank-1 +-1
    W = 7 * np.eye(7, dtype=np.int64) - G
    if not (np.diag(W) == 1).all() or not (np.abs(W) == 1).all():
        ok = False; break
    w = W[0]  # first row of ww^T with w_0=+1 normalization
    if not (np.outer(w, w) == W).all():
        ok = False; break
    D = np.diag(w)
    S2 = D @ S @ D
    if not (S2 @ S2.T == 7 * np.eye(7, dtype=np.int64) - np.ones((7, 7), dtype=np.int64)).all():
        ok = False; break
    if not (S2.sum(axis=1) == 0).all():
        ok = False; break
    A = (S2 == 1).astype(int)
    if not (A.sum(axis=1) == 3).all():
        ok = False; break
    for i in range(7):
        for j in range(7):
            if i != j and int((A[i] & A[j]).sum()) != 1:
                ok = False
assert ok, "EQUALITY ANALYSIS FAILS ON SOME ATTAINER"
print(f"all {len(att7)} labeled attainers: SS^T = 7I - ww^T, switch is regular "
      f"(scores 3) and doubly regular (pair-domination 1): OK")

print("\n== E: QR_7 spectrum and bordering to skew-Hadamard of order 8 ==")
ev = np.linalg.eigvals(Sqr.astype(float))
print("eigenvalues of S(QR_7):", sorted(np.round(ev.imag, 6)))
assert np.allclose(sorted(ev.imag), [-np.sqrt(7)] * 3 + [0] + [np.sqrt(7)] * 3)
M = np.eye(7, dtype=np.int64) + Sqr
H = np.zeros((8, 8), dtype=np.int64)
H[0, 0] = 1
H[0, 1:] = 1
H[1:, 0] = -1
H[1:, 1:] = M
assert (np.abs(H) == 1).all()
assert (H + H.T == 2 * np.eye(8, dtype=np.int64)).all(), "not skew type"
assert (H @ H.T == 8 * np.eye(8, dtype=np.int64)).all(), "not Hadamard"
print("H = border(I+S(QR_7)) is a skew-Hadamard matrix of order 8: OK")

print("\n== F: n=6 attainers, Pfaffian check ==")
PAIRS6 = pairs_of(6)
pf_ok = True
for mk in results[6][1]:
    S = S_from_mask(mk, 6, PAIRS6)
    ds = exact_det(S)
    if ds != 81:
        pf_ok = False
print("all n=6 attainers have det(S) = 81 (|Pf| = 9):", pf_ok)
assert pf_ok

print("\n== G: n=4 skew conference + n=8 via gentourng ==")
PAIRS4 = pairs_of(4)
for mk in results[4][1]:
    S = S_from_mask(mk, 4, PAIRS4)
    assert (S @ S.T == 3 * np.eye(4, dtype=np.int64)).all()
print(f"all {len(results[4][1])} n=4 labeled attainers satisfy SS^T = 3I: OK")

def read_gentourng(n):
    m = n * (n - 1) // 2
    r = subprocess.run(["gentourng", str(n)], capture_output=True, text=True)
    out = []
    for line in (r.stdout or r.stderr).splitlines():
        line = line.strip()
        if len(line) == m and set(line) <= {"0", "1"}:
            out.append(line)
    return out

reps8 = read_gentourng(8)
print(f"gentourng 8: {len(reps8)} iso classes")
PAIRS8 = pairs_of(8)
best8, n8att = -1, 0
conf_ok = True
for rep in reps8:
    S = np.zeros((8, 8), dtype=np.int64)
    # gentourng emits upper triangle row by row: pairs (0,1),(0,2),...,(0,7),(1,2),...
    k = 0
    for i in range(8):
        for j in range(i + 1, 8):
            v = 1 if rep[k] == "1" else -1
            S[i, j] = v
            S[j, i] = -v
            k += 1
    d = exact_det(np.eye(8, dtype=np.int64) + S)
    if d > best8:
        best8, n8att = d, 1
        conf_ok = (S @ S.T == 7 * np.eye(8, dtype=np.int64)).all()
    elif d == best8:
        n8att += 1
        conf_ok = conf_ok and (S @ S.T == 7 * np.eye(8, dtype=np.int64)).all()
print(f"n=8: max det = {best8} (claimed 4096), classes attaining = {n8att} "
      f"(claimed 4), all skew-conference: {bool(conf_ok)}")
assert best8 == 4096 and n8att == 4 and conf_ok

print("\n== H: monotonicity sanity p log(1+c/p) ==")
import math
for n in (5, 7, 9, 11):
    c = n * (n - 1) / 2
    vals = [p * math.log(1 + c / p) for p in range(1, (n - 1) // 2 + 1)]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
print("p log(1+c/p) strictly increasing on p = 1..(n-1)/2 for n=5,7,9,11: OK")

print("\nALL CHECKS PASSED")
