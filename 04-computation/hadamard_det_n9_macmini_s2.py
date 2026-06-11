#!/usr/bin/env python3
"""
macmini-2026-06-10-S2

Hadamard-determinant census at n=9: first new case beyond the parent n<=8 run.

For a tournament T with 0/1 adjacency A (A[i][j]=1 iff i->j), let
S = A - A^T (skew +-1 matrix, zero diagonal).  det(I+S) = 2^(n-1) * d(T),
d an iso/switching/reversal invariant.

This script enumerates ALL 191536 iso classes at n=9 via nauty gentourng and
computes:

  1. d-spectrum (distribution of d = det(I+S)/256 per class).
  2. Ceiling: n=9 == 1 mod 4 so no DRT exists; the odd-n ceiling
     (n+1)^((n-1)/2) = 10^4 = 10000 is NOT attainable.  Record the actual max
     det ("tournament Barba number" B_t(9)) and ALL classes attaining it:
     score sequences, regularity, SS^T structure (proximity to nI-J = 9I-J),
     and H(T) of the maximizers.
  3. Floor: count d=1 classes; verify per-class equivalence with local
     transitivity (no 3-cycle in any in/out-neighborhood).
  4. Pearson(d, H) over all classes; d of the H-maximizer class.
  5. Monte Carlo labeled average of d at n=9 (200k samples) vs I(9) = 2620
     (involution number).

MISTAKE-028 discipline: primary determinant path is batched numpy float det
with integer rounding, VERIFIED by exact fraction-free Bareiss (Python ints)
on a 1000-class random sample and on ALL extremal classes (max-d and d=1).
Vectorized H and local-transitivity code is verified against trusted scalar
reference implementations on random samples.
"""

from __future__ import annotations

import itertools
import subprocess
import time
from collections import Counter

import numpy as np

N = 9
M = N * (N - 1) // 2  # 36
PAIRS = list(itertools.combinations(range(N), 2))
POW = 1 << (N - 1)  # 256
GENTOURNG = "/opt/homebrew/bin/gentourng"
RNG = np.random.default_rng(20260610)

t0 = time.time()


def stamp(msg: str) -> None:
    print(f"[{time.time() - t0:8.1f}s] {msg}", flush=True)


# ----------------------------------------------------------------------------
# 0. Enumerate iso classes via gentourng and build adjacency tensor
# ----------------------------------------------------------------------------

def read_gentourng_matrix() -> np.ndarray:
    """Return bit matrix (Nclass, 36) uint8; column k = pair PAIRS[k] (MSB-first)."""
    result = subprocess.run([GENTOURNG, str(N)], capture_output=True, check=True)
    lines = [ln for ln in result.stdout.split(b"\n") if len(ln) == M]
    raw = np.frombuffer(b"".join(lines), dtype=np.uint8).reshape(len(lines), M)
    return (raw == ord("1")).astype(np.uint8)


def adjacency_from_bits(bits: np.ndarray) -> np.ndarray:
    """(Nclass, 36) bits -> (Nclass, 9, 9) adjacency, A[:,i,j]=1 iff i->j."""
    nc = bits.shape[0]
    A = np.zeros((nc, N, N), dtype=np.uint8)
    for k, (i, j) in enumerate(PAIRS):
        A[:, i, j] = bits[:, k]
        A[:, j, i] = 1 - bits[:, k]
    return A


stamp("running gentourng 9 ...")
bits = read_gentourng_matrix()
NC = bits.shape[0]
stamp(f"gentourng gave {NC} iso-class representatives (expected A000568(9)=191536)")
assert NC == 191536, f"class count mismatch: {NC}"
A = adjacency_from_bits(bits)

# sanity: A + A^T = J - I for every class (check on full tensor)
assert np.all(A + A.transpose(0, 2, 1) == (1 - np.eye(N, dtype=np.uint8))), "adjacency broken"

S_all = A.astype(np.int8) - A.transpose(0, 2, 1).astype(np.int8)  # skew +-1


# ----------------------------------------------------------------------------
# 1. Determinants: numpy float primary, exact Bareiss verification
# ----------------------------------------------------------------------------

def bareiss_det(mat: list[list[int]]) -> int:
    """Exact integer determinant, fraction-free Bareiss. mat is modified copy."""
    a = [row[:] for row in mat]
    n = len(a)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    a[k], a[r] = a[r], a[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]


def det_I_plus_S_exact(idx: int) -> int:
    mat = (np.eye(N, dtype=np.int64) + S_all[idx]).tolist()
    return bareiss_det([[int(x) for x in row] for row in mat])


stamp("computing det(I+S) for all classes (batched numpy float64)...")
M_all = np.eye(N, dtype=np.float64)[None, :, :] + S_all.astype(np.float64)
dets_f = np.empty(NC, dtype=np.float64)
CH = 65536
for lo in range(0, NC, CH):
    dets_f[lo:lo + CH] = np.linalg.det(M_all[lo:lo + CH])
del M_all
dets = np.rint(dets_f).astype(np.int64)
max_float_err = float(np.max(np.abs(dets_f - dets)))
stamp(f"float->int rounding: max |residual| = {max_float_err:.3e} (must be << 0.5)")
assert max_float_err < 1e-3, "float determinant rounding unsafe"
assert np.all(dets > 0), "nonpositive det(I+S) found (impossible for tournaments)"
assert np.all(dets % POW == 0), "det not divisible by 2^(n-1)=256 somewhere"
d_all = dets // POW

# ----------------------------------------------------------------------------
# 2. H(T) for all classes: batched bitmask DP, verified against scalar DP
# ----------------------------------------------------------------------------

stamp("computing H(T) for all classes (batched bitmask DP)...")
FULL = (1 << N) - 1
mask_order = sorted(range(1, FULL + 1), key=lambda m: (bin(m).count("1"), m))
H_all = np.empty(NC, dtype=np.int64)
B = 2048
Af = A.astype(np.int64)
for lo in range(0, NC, B):
    hi = min(lo + B, NC)
    nb = hi - lo
    Ab = Af[lo:hi]  # (nb, 9, 9)
    dp = np.zeros((nb, FULL + 1, N), dtype=np.int64)
    for v in range(N):
        dp[:, 1 << v, v] = 1
    for mask in mask_order:
        pc = bin(mask).count("1")
        if pc < 2:
            continue
        mm = mask
        while mm:
            lsb = mm & -mm
            v = lsb.bit_length() - 1
            mm ^= lsb
            pm = mask ^ lsb
            # dp[:, mask, v] = sum_u dp[:, pm, u] * A[u -> v]
            dp[:, mask, v] = np.einsum("bu,bu->b", dp[:, pm, :], Ab[:, :, v])
    H_all[lo:hi] = dp[:, FULL, :].sum(axis=1)
    del dp
stamp("H DP done")


def hamiltonian_paths_scalar(idx: int) -> int:
    rows = [0] * N
    for i in range(N):
        for j in range(N):
            if A[idx, i, j]:
                rows[i] |= 1 << j
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1 << N):
        missing = FULL ^ mask
        for v in range(N):
            val = dp[mask][v]
            if not val:
                continue
            nxt = rows[v] & missing
            while nxt:
                lsb = nxt & -nxt
                w = lsb.bit_length() - 1
                dp[mask | lsb][w] += val
                nxt ^= lsb
    return sum(dp[FULL])


# ----------------------------------------------------------------------------
# 3. Local transitivity for all classes (vectorized), verified vs scalar
# ----------------------------------------------------------------------------

stamp("computing local transitivity for all classes (vectorized)...")
Ai = A.astype(np.int64)
loc_trans = np.ones(NC, dtype=bool)
for v in range(N):
    for direction in ("out", "in"):
        nbr = Ai[:, v, :] if direction == "out" else Ai[:, :, v]  # (NC, 9)
        sub = Ai * nbr[:, :, None] * nbr[:, None, :]
        deg = sub.sum(axis=2)                       # (NC, 9) out-deg inside nbhd
        k = nbr.sum(axis=1)                         # (NC,)
        triples = k * (k - 1) * (k - 2) // 6        # C(k,3)
        trans_triples = (deg * (deg - 1) // 2).sum(axis=1)  # sum C(deg,2)
        c3 = triples - trans_triples                # number of 3-cycles in nbhd
        assert np.all(c3 >= 0)
        loc_trans &= (c3 == 0)
        del sub


def locally_transitive_scalar(idx: int) -> bool:
    rows = A[idx]
    for v in range(N):
        for nb in (np.flatnonzero(rows[v]), np.flatnonzero(rows[:, v])):
            for (a, b, c) in itertools.combinations(nb, 3):
                # 3-cycle test on triple {a,b,c}
                fwd = rows[a, b] + rows[b, c] + rows[c, a]
                if fwd == 3 or fwd == 0:
                    return False
    return True


# ----------------------------------------------------------------------------
# 4. Verification block (MISTAKE-028 discipline)
# ----------------------------------------------------------------------------

stamp("verification: exact Bareiss on 1000 random classes...")
sample = RNG.choice(NC, size=1000, replace=False)
for idx in sample:
    assert det_I_plus_S_exact(int(idx)) == int(dets[idx]), f"Bareiss mismatch at class {idx}"
stamp("  1000/1000 random-class dets match exactly")

stamp("verification: scalar H on 50 random classes...")
for idx in sample[:50]:
    assert hamiltonian_paths_scalar(int(idx)) == int(H_all[idx]), f"H mismatch at class {idx}"
stamp("  50/50 H values match scalar DP")

stamp("verification: scalar local transitivity on 60 random classes...")
ver_idx = list(sample[:40]) + list(np.flatnonzero(loc_trans)[:20])
for idx in ver_idx:
    assert locally_transitive_scalar(int(idx)) == bool(loc_trans[idx]), f"LT mismatch at {idx}"
stamp(f"  {len(ver_idx)}/{len(ver_idx)} local-transitivity flags match scalar check")

# ----------------------------------------------------------------------------
# 5. Results
# ----------------------------------------------------------------------------

print()
print("=" * 78)
print("n=9 HADAMARD-DETERMINANT CENSUS  (all 191536 iso classes)")
print("=" * 78)

# --- 1. d-spectrum -----------------------------------------------------------
spec = Counter(d_all.tolist())
print()
print(f"[1] d-SPECTRUM at n=9   (d = det(I+S)/256;  {len(spec)} distinct values)")
print(f"    {'d':>5} | {'#classes':>9}")
for dv in sorted(spec):
    print(f"    {dv:>5} | {spec[dv]:>9}")
print(f"    total classes = {sum(spec.values())}")
print(f"    parity check: all d odd? {all(dv % 2 == 1 for dv in spec)}")

# --- 2. ceiling --------------------------------------------------------------
dmax = int(d_all.max())
det_max = dmax * POW
maximizers = np.flatnonzero(d_all == dmax)
print()
print(f"[2] CEILING at n=9 (n == 1 mod 4, no DRT exists)")
print(f"    DRT/odd-n bound (n+1)^((n-1)/2) = 10^4 = 10000  -> NOT attained")
print(f"    tournament Barba number B_t(9) = max det(I+S) = {det_max}  (d = {dmax})")
print(f"    ratio to bound: {det_max}/10000 = {det_max/10000:.4f}")
print(f"    number of classes attaining max: {len(maximizers)}")

target = 9 * np.eye(N, dtype=np.int64) - np.ones((N, N), dtype=np.int64)  # nI - J
for rank, idx in enumerate(maximizers, 1):
    idx = int(idx)
    # exact verification of every maximizer
    de = det_I_plus_S_exact(idx)
    assert de == det_max, f"maximizer {idx} exact det mismatch"
    Sm = S_all[idx].astype(np.int64)
    SST = Sm @ Sm.T
    off = SST[~np.eye(N, dtype=bool)]
    off_hist = Counter(off.tolist())
    frob2 = int(((SST - target) ** 2).sum())
    scores = tuple(sorted(int(x) for x in A[idx].sum(axis=1)))
    regular = scores == tuple([4] * 9)
    print(f"    --- maximizer #{rank}: class index {idx} (gentourng order) ---")
    print(f"        det(I+S) = {de} (exact Bareiss verified), d = {de // POW}")
    print(f"        score sequence = {scores}   regular = {regular}")
    print(f"        H(T) = {int(H_all[idx])}")
    print(f"        SS^T diag = {sorted(set(int(x) for x in np.diag(SST)))} (always n-1=8)")
    print(f"        SS^T off-diag histogram = {dict(sorted(off_hist.items()))}")
    print(f"        ||SS^T - (9I-J)||_F^2 = {frob2}"
          f"   (0 would mean DRT; impossible at n=9)")
    print(f"        SS^T =")
    for row in SST:
        print(f"          {' '.join(f'{int(x):>3}' for x in row)}")
    print(f"        upper-triangle arc bits (gentourng line) = "
          f"{''.join(str(int(b)) for b in bits[idx])}")

# --- 3. floor ----------------------------------------------------------------
d1 = np.flatnonzero(d_all == 1)
lt = np.flatnonzero(loc_trans)
print()
print(f"[3] FLOOR at n=9")
print(f"    classes with d = 1: {len(d1)}")
print(f"    locally transitive classes: {len(lt)}")
same_set = (len(d1) == len(lt)) and bool(np.all(d1 == lt))
print(f"    per-class set equality (d=1 <=> locally transitive): {same_set}")
if not same_set:
    only_d1 = set(d1.tolist()) - set(lt.tolist())
    only_lt = set(lt.tolist()) - set(d1.tolist())
    print(f"    !! d=1 but not loc-trans: {sorted(only_d1)[:20]}")
    print(f"    !! loc-trans but not d=1: {sorted(only_lt)[:20]}")
stamp("verification: exact Bareiss on ALL d=1 classes...")
for idx in d1:
    assert det_I_plus_S_exact(int(idx)) == POW, f"d=1 class {idx} exact det mismatch"
print(f"    exact Bareiss verified det(I+S)=256 on all {len(d1)} floor classes")
print(f"    score sequences of d=1 classes: "
      f"{sorted(Counter(tuple(sorted(int(x) for x in A[int(i)].sum(axis=1))) for i in d1).items())}")

# --- 4. d vs H ---------------------------------------------------------------
print()
print(f"[4] d vs H ACROSS ALL {NC} CLASSES")
r = float(np.corrcoef(d_all.astype(np.float64), H_all.astype(np.float64))[0, 1])
print(f"    Pearson r(d, H) = {r:.6f}")
hmax = int(H_all.max())
hmax_idx = np.flatnonzero(H_all == hmax)
print(f"    H range: [{int(H_all.min())}, {hmax}]")
print(f"    H-maximizer classes ({len(hmax_idx)}): "
      f"{[(int(i), int(d_all[i])) for i in hmax_idx]}  (class, d)")
print(f"    d of H-maximizer class(es): {sorted(set(int(d_all[i]) for i in hmax_idx))}")
print(f"    H of d-maximizer class(es): {sorted(set(int(H_all[i]) for i in maximizers))}")
hmin_idx = np.flatnonzero(H_all == int(H_all.min()))
print(f"    H-minimizer: H = {int(H_all.min())} on {len(hmin_idx)} class(es), "
      f"d values {sorted(set(int(d_all[i]) for i in hmin_idx))}")
print(f"    mean d per class = {float(d_all.mean()):.4f}, "
      f"mean H per class = {float(H_all.mean()):.4f}")
# d-spectrum conditioned on H parity sanity (H always odd by Redei)
assert np.all(H_all % 2 == 1), "Redei violated?!"
print(f"    Redei check: all H odd = True")

# --- 5. Monte Carlo labeled average ------------------------------------------
print()
print(f"[5] MONTE CARLO LABELED AVERAGE at n=9 (200000 samples)")
NS = 200_000
acc = []
iu, ju = np.triu_indices(N, k=1)
for lo in range(0, NS, 50_000):
    nb = min(50_000, NS - lo)
    signs = RNG.integers(0, 2, size=(nb, M), dtype=np.int8) * 2 - 1
    Sr = np.zeros((nb, N, N), dtype=np.float64)
    Sr[:, iu, ju] = signs
    Sr[:, ju, iu] = -signs
    Mr = Sr + np.eye(N)[None, :, :]
    dr = np.rint(np.linalg.det(Mr)).astype(np.int64)
    assert np.all(dr % POW == 0)
    acc.append(dr // POW)
d_samples = np.concatenate(acc)
mean_d = float(d_samples.mean())
sem = float(d_samples.std(ddof=1) / np.sqrt(NS))
# involution numbers I(n): I(0)=1, I(1)=1, I(n)=I(n-1)+(n-1)I(n-2)
inv = [1, 1]
for k in range(2, 10):
    inv.append(inv[-1] + (k - 1) * inv[-2])
print(f"    I(9) (involution number) = {inv[9]}")
print(f"    Monte Carlo mean d = {mean_d:.3f}  +- {sem:.3f} (SEM)")
print(f"    z-score vs I(9)={inv[9]}: {(mean_d - inv[9]) / sem:+.2f}")
print(f"    sample d range: [{int(d_samples.min())}, {int(d_samples.max())}]")

# --- labeled-weighted check of E[d] from the class census (exact, via Aut) ---
# (skipped: requires |Aut| per class; MC covers the labeled average instead)

print()
stamp("ALL DONE")
print()
print("Summary line: "
      f"B_t(9)={det_max} (d={dmax}, {len(maximizers)} classes), "
      f"floor d=1 on {len(d1)} classes = loc-trans count {len(lt)} "
      f"(set-equal: {same_set}), "
      f"Pearson(d,H)={r:.4f}, "
      f"MC mean d={mean_d:.1f} vs I(9)={inv[9]}, "
      f"distinct d values={len(spec)}")
