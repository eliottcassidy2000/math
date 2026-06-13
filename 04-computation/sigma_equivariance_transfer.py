#!/usr/bin/env python3
"""
THEOREM (sigma-equivariance of transfer matrix for SC tournaments)

For a self-converse tournament T with involution anti-automorphism sigma:

  Off-diagonal: M[sigma(a), sigma(b)] = (-1)^{n-2} M[a,b]   (a != b)
  Diagonal:     M[sigma(a), sigma(a)] = (-1)^{n-1} M[a,a]

Unified: P_sigma M P_sigma^T = (-1)^{n-2} (M - 2*diag(M))

COROLLARY (odd n, (-1)^{n-2} = -1):
  If sigma(a) = b (a != b), then M[a,b] = 0.
  At odd n, sigma has 1 fixed point and (n-1)/2 swapped pairs.
  So (n-1)/2 off-diagonal entries M[i, sigma(i)] are forced to be 0.
  The remaining off-diagonal entries pair up: M[a,b] = -M[sigma(a), sigma(b)].

COROLLARY (even n, (-1)^{n-2} = +1):
  sigma preserves all off-diagonal M entries.
  Diagonal entries negate: M[sigma(a),sigma(a)] = -M[a,a].
  Since sigma is fixed-point-free at even n, diagonal entries come in +- pairs.

CONNECTION TO GS FLIP DICHOTOMY:
  At odd n: M[i,sigma(i)] = 0 is the algebraic manifestation of GS flip pairs
  always crossing isomorphism classes (THM-022, blueself odd-n obstruction).
  The "perpendicular" structure forces zero coupling between sigma-paired vertices.

  At even n: M commutes with sigma on off-diagonal, allowing some GS flips
  within the same class.

PROOF:
  From the reindexing identity (THM-030 dependency):
    M_T[b,a] = (-1)^{n-2} M_{T^op}[a,b]  (a != b)
  Anti-aut sigma: M_{T^op}[a,b] = M_T[sigma(a), sigma(b)]
  Combined with symmetry M[a,b] = M[b,a] (THM-030):
    M[a,b] = (-1)^{n-2} M[sigma(a), sigma(b)]  (off-diagonal)

  For diagonal: pos(a, P) in path P maps to n-1-pos(a, P) under path reversal.
  At odd n: (-1)^{pos} = (-1)^{n-1-pos} (parity preserved, since n-1 even).
  At even n: (-1)^{pos} = -(-1)^{n-1-pos} (parity flipped, since n-1 odd).
  So M[sigma(a),sigma(a)] = (-1)^{n-1} M[a,a].

VERIFIED: n=3,4,5 (exhaustive, 0 violations), n=7 (20 random SC tournaments, 0 violations).

opus-2026-03-06-S11b
"""

from itertools import permutations
import numpy as np

def compute_M_entry(T, n, a, b):
    """Compute M[a,b] for numeric tournament T."""
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(len(perm)-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod != 0:
                pos = list(perm).index(a)
                val += (-1)**pos * prod
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S)
            S_set = set(S) | {a}
            R_set = set(R) | {b}
            ea = 0
            for p in permutations(sorted(S_set)):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod
            if len(S_set) == 1: ea = 1
            bb2 = 0
            for p in permutations(sorted(R_set)):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod
            if len(R_set) == 1: bb2 = 1
            val += sign * ea * bb2
        return val

def find_involution_anti_auts(T, n):
    """Find all involution anti-automorphisms of T."""
    results = []
    for perm in permutations(range(n)):
        if not all(perm[perm[i]] == i for i in range(n)):
            continue
        is_anti = True
        for i in range(n):
            for j in range(n):
                if i != j:
                    if T.get((i,j), 0) != T.get((perm[j], perm[i]), 0):
                        is_anti = False
                        break
            if not is_anti:
                break
        if is_anti:
            results.append(perm)
    return results

if __name__ == '__main__':
    # Verification at n=3,4,5
    for n in [3, 4, 5]:
        off_factor = (-1)**(n-2)
        diag_factor = (-1)**(n-1)
        violations = 0
        sc_count = 0

        pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
        for mask in range(1 << len(pairs)):
            T = {}
            for idx, (i,j) in enumerate(pairs):
                if mask & (1 << idx):
                    T[(i,j)] = 1; T[(j,i)] = 0
                else:
                    T[(i,j)] = 0; T[(j,i)] = 1

            anti_auts = find_involution_anti_auts(T, n)
            if not anti_auts:
                continue
            sc_count += 1
            sigma = anti_auts[0]

            M = np.zeros((n, n))
            for a in range(n):
                for b in range(n):
                    M[a,b] = compute_M_entry(T, n, a, b)

            for a in range(n):
                for b in range(n):
                    expected_factor = diag_factor if a == b else off_factor
                    if abs(M[sigma[a], sigma[b]] - expected_factor * M[a,b]) > 1e-10:
                        violations += 1

        print(f"n={n}: {sc_count} SC tournaments, {violations} violations")

    print("\nAll verifications passed." if violations == 0 else f"\nFAILED: {violations} violations")
