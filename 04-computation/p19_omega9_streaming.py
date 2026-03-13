#!/usr/bin/env python3
"""
Memory-efficient P_19 Ω_9 computation.

Instead of storing all 6M constraint columns in memory,
stream them through Gaussian elimination one at a time.

Key insight: we don't need to store all columns — just process
each column immediately and discard it after reduction.

Ω_8 = 331958 (already computed)
R_9 = 268340

opus-2026-03-13-S71b
"""

import time
import sys

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def stream_orbit_reps_d9(p, QR_list):
    """Generator yielding d=9 orbit reps one at a time (starting with 1)."""
    QR_set = set(QR_list)
    d = 9

    def backtrack(seq, ps_set, last_ps, depth):
        if depth == d:
            yield tuple(seq)
            return
        for s in QR_list:
            new_ps = (last_ps + s) % p
            if new_ps in ps_set or new_ps == 0:
                continue
            seq.append(s)
            ps_set.add(new_ps)
            yield from backtrack(seq, ps_set, new_ps, depth + 1)
            seq.pop()
            ps_set.discard(new_ps)

    yield from backtrack([1], {0, 1}, 1, 1)

def zm_orbit_class(seq, QR_list, p):
    best = seq
    for q in QR_list:
        if q == 1:
            continue
        scaled = tuple((q * s) % p for s in seq)
        if scaled < best:
            best = scaled
    return best

def compute_constraint_col(sigma, p, d, QR_list, QR_set, junk_map, q_mod):
    """Compute the constraint column for a single orbit rep."""
    col = {}
    for i in range(1, d):
        merged = (sigma[i-1] + sigma[i]) % p
        if merged not in QR_set and merged != 0:
            face = list(sigma)
            face[i-1] = merged
            del face[i]
            ft = zm_orbit_class(tuple(face), QR_list, p)
            if ft not in junk_map:
                junk_map[ft] = len(junk_map)
            row = junk_map[ft]
            col[row] = (col.get(row, 0) + (-1)**i) % q_mod
    return {k: v for k, v in col.items() if v != 0}

# Streaming Gaussian elimination
class StreamingGaussElim:
    def __init__(self, q):
        self.q = q
        self.pivot_col = {}  # row -> index in reduced
        self.reduced = []
        self.rank = 0
        self.n_processed = 0

    def process_column(self, col):
        """Process one column through Gaussian elimination."""
        c = dict(col)
        q = self.q

        while c:
            r = min(c.keys())
            v = c[r] % q
            if v == 0:
                del c[r]
                continue
            if r in self.pivot_col:
                pv = self.reduced[self.pivot_col[r]]
                pv_val = pv[r]
                inv_pv = pow(pv_val, q - 2, q)
                factor = (v * inv_pv) % q
                for pr, pval in pv.items():
                    c[pr] = (c.get(pr, 0) - factor * pval) % q
                c = {k: v % q for k, v in c.items() if v % q != 0}
            else:
                self.pivot_col[r] = len(self.reduced)
                self.reduced.append(c)
                self.rank += 1
                break

        self.n_processed += 1

# ============================================================
p = 19
m = 9
QR = get_QR(p)
QR_set = set(QR)
q_mod = 997
d = 9

print(f"P_{p} (m={m}), d={d}")
print(f"QR = {QR}")
print(f"Streaming Gaussian elimination mod {q_mod}\n")

# Known values
omega_8 = 331958
R_8 = 63618
R_9 = omega_8 - R_8  # = 268340
print(f"Ω_8 = {omega_8}, R_8 = {R_8}, R_9 = {R_9}")

junk_map = {}
ge = StreamingGaussElim(q_mod)

t0 = time.time()
n_reps = 0

for sigma in stream_orbit_reps_d9(p, QR):
    n_reps += 1
    col = compute_constraint_col(sigma, p, d, QR, QR_set, junk_map, q_mod)
    if col:
        ge.process_column(col)

    if n_reps % 100000 == 0:
        elapsed = time.time() - t0
        rate = n_reps / elapsed
        print(f"  {n_reps} reps, rank={ge.rank}, junk_orbits={len(junk_map)}, "
              f"{elapsed:.1f}s ({rate:.0f} reps/s)", flush=True)

t1 = time.time()
print(f"\nTotal: {n_reps} orbit reps, {len(junk_map)} junk orbits")
print(f"rank(C_9) = {ge.rank}")

omega_9 = n_reps - ge.rank
print(f"Ω_9^orb = {n_reps} - {ge.rank} = {omega_9}")

Budget = omega_9 - R_9
print(f"\nBudget_9 = Ω_9 - R_9 = {omega_9} - {R_9} = {Budget}")
print(f"Time: {t1-t0:.1f}s")

# β_9 = Budget - R_10. We need R_10 from top recursion.
# But we can check: is Budget consistent with β_9 = (m-3)/2 = 3?
# β_9 = Budget - R_10
# From palindromic-like considerations or direct computation of top half.
print(f"\nPredicted β_9^orb = (m-3)/2 = {(m-3)//2}")
print(f"Predicted β_9 = m(m-3)/2 = {m*(m-3)//2}")
print(f"If β_9^orb = 3, then R_10 = Budget - 3 = {Budget - 3}")
