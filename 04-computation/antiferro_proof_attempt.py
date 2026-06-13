"""
antiferro_proof_attempt.py — opus-2026-04-04-S3

Attempt to PROVE: H(ij)·H(0) < H(i)·H(j) for all tile pairs.

Equivalently: flipping two tiles simultaneously gives FEWER Hamiltonian paths
than the product of individual effects would predict (log-submodularity).

Strategy: Express H(ij)·H(0) - H(i)·H(j) in terms of the OCF decomposition
and show it's always negative.

H(0) = 1 (transitive tournament).
H(i) = 1 + c_i = 1 + 2^(s_i - 1)
H(j) = 1 + c_j = 1 + 2^(s_j - 1)
H(ij) = 1 + c_i + c_j + c_{ij}

So: H(ij)·H(0) - H(i)·H(j)
  = (1 + c_i + c_j + c_{ij}) · 1 - (1 + c_i)(1 + c_j)
  = 1 + c_i + c_j + c_{ij} - 1 - c_i - c_j - c_i·c_j
  = c_{ij} - c_i·c_j

So the antiferromagnetic property is: c_{ij} < c_i · c_j.

Since c_i = 2^(s_i-1) and c_j = 2^(s_j-1):
  c_i · c_j = 2^(s_i + s_j - 2)

And from the OCF (THM-287):
  c_{ij} = 2·A(i,j) + 4·B(i,j)

We need: 2·A(i,j) + 4·B(i,j) < 2^(s_i + s_j - 2)

For shared-vertex same-end: c_{ij} = -2^max(1,|Δs|-1) < 0 < 2^(s_i+s_j-2).
  Trivially true! ✓

For shared-vertex cross-end: c_{ij} = +2^(s_i+s_j-2) = c_i·c_j.
  We need STRICT inequality: c_{ij} < c_i·c_j.
  But c_{ij} = 2^(s_i+s_j-2) = c_i·c_j exactly!

  Wait, is the odds ratio exactly 1 for cross-end pairs?

Let me check: cross-end pair with x1 = y2.
H(0) = 1, H(i) = 1 + 2^(s_i-1), H(j) = 1 + 2^(s_j-1)
H(ij) = 1 + c_i + c_j + c_{ij}
     = 1 + 2^(s_i-1) + 2^(s_j-1) + 2^(s_i+s_j-2)
     = 1 + 2^(s_i-1) + 2^(s_j-1) + 2^(s_i-1)·2^(s_j-1)
     = (1 + 2^(s_i-1)) · (1 + 2^(s_j-1))
     = H(i) · H(j)

So H(ij) = H(i)·H(j)! The odds ratio is EXACTLY 1 for cross-end pairs!
This means: the cross-end interaction is ZERO in the log.

But wait, log(H(ij)·H(0)/(H(i)·H(j))) = log(1) = 0? But the data shows
cross-end pairs have λ₂ ≈ -1.099 (not zero)...

Let me re-examine. The Möbius inversion for log(H) includes higher-order corrections.
The pairwise log-coefficient λ_{ij} is NOT just log(H(ij)·H(0)/(H(i)·H(j))).
It's the FULL Möbius inversion of log(H) at the subset {i,j}, which equals:

λ_{ij} = log H(ij) - log H(i) - log H(j) + log H(0)
       = log(H(ij)·H(0) / (H(i)·H(j)))

But I claimed above that H(ij) = H(i)·H(j) for cross-end pairs.
If true, λ_{ij} = log(1) = 0. But data says λ₂ ≠ 0.

There must be an error in my claim. Let me verify numerically.
"""
import numpy as np
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles
from math import log

def verify_cross_end():
    for n in [5, 6]:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        print(f"\nn={n}: Cross-end pair verification")
        for i in range(m):
            for j in range(i+1, m):
                x1, y1 = tiles[i]
                x2, y2 = tiles[j]
                if x1 == y2 or y1 == x2:
                    s1, s2 = x1-y1, x2-y2

                    h00 = int(H[0])
                    h10 = int(H[1 << i])
                    h01 = int(H[1 << j])
                    h11 = int(H[(1 << i) | (1 << j)])

                    c_i = int(coeffs.get(1 << i, 0))
                    c_j = int(coeffs.get(1 << j, 0))
                    c_ij = int(coeffs.get((1 << i) | (1 << j), 0))

                    product = h10 * h01
                    actual = h11 * h00
                    diff = actual - product

                    odds = actual / product if product else 0
                    log_odds = log(odds) if odds > 0 else float('-inf')

                    print(f"  {tiles[i]}-{tiles[j]} s=({s1},{s2}): "
                          f"H(0)={h00}, H(i)={h10}, H(j)={h01}, H(ij)={h11} "
                          f"| c_i={c_i}, c_j={c_j}, c_ij={c_ij}, c_i*c_j={c_i*c_j} "
                          f"| H(ij)*H(0)={actual}, H(i)*H(j)={product}, diff={diff} "
                          f"| odds={odds:.6f}, log={log_odds:.6f}")

    # Also verify: the claim c_{ij} = c_i·c_j for cross-end
    print(f"\nDoes c_ij = c_i*c_j for cross-end pairs?")
    for n in [5, 6]:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)
        for i in range(m):
            for j in range(i+1, m):
                x1, y1 = tiles[i]
                x2, y2 = tiles[j]
                if x1 == y2 or y1 == x2:
                    c_i = int(coeffs.get(1 << i, 0))
                    c_j = int(coeffs.get(1 << j, 0))
                    c_ij = int(coeffs.get((1 << i) | (1 << j), 0))
                    match = c_ij == c_i * c_j
                    if not match:
                        print(f"  n={n}: {tiles[i]}-{tiles[j]}: c_ij={c_ij} vs c_i*c_j={c_i*c_j}: MISMATCH")

def prove_antiferro():
    """
    The core inequality: c_{ij} < c_i · c_j for ALL pairs.

    From OCF: c_i · c_j = 2^(s_i-1) · 2^(s_j-1) = 2^(s_i+s_j-2)
    And: c_{ij} = 2·A(i,j) + 4·B(i,j)

    Check: is c_{ij} always < c_i·c_j?
    Equivalent: is 2·A + 4·B < 2^(s_i+s_j-2)?
    """
    for n in [5, 6, 7]:
        H, tiles, m = compute_all_H(n)
        coeffs = mobius_inversion(H, m)

        print(f"\nn={n}: Antiferromagnetic inequality c_ij < c_i*c_j")
        violations = 0
        equalities = 0
        for i in range(m):
            for j in range(i+1, m):
                c_i = int(coeffs.get(1 << i, 0))
                c_j = int(coeffs.get(1 << j, 0))
                c_ij = int(coeffs.get((1 << i) | (1 << j), 0))
                prod = c_i * c_j

                if c_ij > prod:
                    violations += 1
                    print(f"  VIOLATION: {tiles[i]}-{tiles[j]}: c_ij={c_ij} > c_i*c_j={prod}")
                elif c_ij == prod:
                    equalities += 1

        total = m * (m-1) // 2
        strict = total - violations - equalities
        print(f"  Strict c_ij < c_i*c_j: {strict}/{total}")
        print(f"  Equality c_ij = c_i*c_j: {equalities}/{total}")
        print(f"  Violations: {violations}/{total}")

if __name__ == "__main__":
    verify_cross_end()
    prove_antiferro()
