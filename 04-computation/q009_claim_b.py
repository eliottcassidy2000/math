#!/usr/bin/env python3
"""
Verify Claim (B): the d-independent identity.

For any tournament on m vertices W with distinguished vertex v:

  sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) = -(-1)^m h_end(W, v)

where H(S) = total Hamiltonian path weight on S (H(empty) = 1).
h_start(R, v) = paths on R starting at v (for v in R).
h_end(W, v) = paths on W ending at v.

This identity involves ONLY internal arcs (no interface variables).
It is the key remaining step for the inductive proof of OCF.

PROVED ALGEBRAICALLY at m=1,2,3. Need general proof.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_start(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def H_total(T, verts):
    """Total Hamiltonian path weight on verts (any start/end)."""
    if len(verts) == 0:
        return 1.0
    if len(verts) == 1:
        return 1.0
    total = 0.0
    for p in permutations(verts):
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def verify_claim_B(m, num_trials=200):
    """Verify Claim B for tournament on m vertices."""
    print(f"=== Claim (B) at m={m} ===\n")

    W = list(range(m))
    arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        arc_values = {}
        for (a, b) in arcs:
            val = random.uniform(-2, 3)
            arc_values[(a, b)] = val
            arc_values[(b, a)] = 1 - val

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        for v in W[:min(3, m)]:
            Wv = [w for w in W if w != v]  # W\{v}

            # LHS: sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v)
            lhs = 0.0
            for smask in range(1 << len(Wv)):
                S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
                R = [w for w in W if w not in S]  # W\S, includes v

                H_S = H_total(T, S)
                h_R_v = h_start(T, R, v)

                lhs += ((-1) ** len(S)) * H_S * h_R_v

            # RHS: -(-1)^m h_end(W, v)
            rhs = -((-1) ** m) * h_end(T, W, v)

            err = abs(lhs - rhs)
            max_err = max(max_err, err)

            if trial < 3:
                print(f"  trial {trial}, v={v}: LHS={lhs:.6f}, RHS={rhs:.6f}, err={err:.2e}")

    print(f"\n  Max error over {num_trials} trials: {max_err:.2e}")
    if max_err < 1e-8:
        print(f"  CONFIRMED: Claim (B) holds at m={m}")
    else:
        print(f"  FAILED!")
    return max_err < 1e-8


def verify_claim_B_variant(m, num_trials=200):
    """
    Alternative form of Claim (B):

    sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) + (-1)^m h_end(W,v) = 0

    Equivalently, define Phi(v) = full alternating sum including H(empty) = 1:
    Phi(v) = sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v)
    Then: Phi(v) = -(-1)^m h_end(W,v) = (-1)^{m+1} h_end(W,v)
    """
    print(f"\n=== Claim (B) variant check at m={m} ===")

    W = list(range(m))
    arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        arc_values = {}
        for (a, b) in arcs:
            val = random.uniform(-1, 2)
            arc_values[(a, b)] = val
            arc_values[(b, a)] = 1 - val

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        for v in W[:2]:
            # Decompose by |S|
            terms_by_k = {}
            Wv = [w for w in W if w != v]
            for smask in range(1 << len(Wv)):
                S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
                R = [w for w in W if w not in S]
                k = len(S)
                term = ((-1) ** k) * H_total(T, S) * h_start(T, R, v)
                terms_by_k[k] = terms_by_k.get(k, 0) + term

            h_end_v = h_end(T, W, v)
            print(f"  trial {trial}, v={v}: h_end={h_end_v:.4f}")
            for k in sorted(terms_by_k):
                print(f"    |S|={k}: {terms_by_k[k]:.6f}")
            total = sum(terms_by_k.values())
            print(f"    Total Phi(v) = {total:.6f}, expected = {(-1)**(m+1)*h_end_v:.6f}")


def prove_claim_B_induction_test(m, num_trials=100):
    """
    Test: does Claim (B) for m follow from Claim (B) for m-1?

    Strategy: fix a vertex a != v, split S based on whether a in S.
    """
    print(f"\n=== Inductive structure test at m={m} ===\n")

    W = list(range(m))
    arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)

    for trial in range(min(num_trials, 5)):
        arc_values = {}
        for (a, b) in arcs:
            val = random.uniform(-1, 2)
            arc_values[(a, b)] = val
            arc_values[(b, a)] = 1 - val

        def T(a, b):
            if a == b:
                return 0
            return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))

        v = 0  # distinguished vertex
        a_vertex = 1  # vertex to split on
        Wva = [w for w in W if w != v and w != a_vertex]  # W\{v,a}

        # Split: S not containing a (a in R)
        sum_without_a = 0.0
        for smask in range(1 << len(Wva)):
            S = [Wva[bit] for bit in range(len(Wva)) if smask & (1 << bit)]
            R = [w for w in W if w not in S]  # includes v and a
            sum_without_a += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        # Split: S containing a (a not in R)
        sum_with_a = 0.0
        for smask in range(1 << len(Wva)):
            S_prime = [Wva[bit] for bit in range(len(Wva)) if smask & (1 << bit)]
            S = S_prime + [a_vertex]
            R = [w for w in W if w not in S]  # includes v but not a
            sum_with_a += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        total = sum_without_a + sum_with_a
        expected = -((-1) ** m) * h_end(T, W, v)

        print(f"  trial {trial}: without_a={sum_without_a:.4f}, with_a={sum_with_a:.4f}, "
              f"total={total:.4f}, expected={expected:.4f}, err={abs(total-expected):.2e}")

        # Can we express sum_without_a using Claim B at m-1?
        # When a is in R: R = {v, a} ∪ (Wva\S).
        # h_start(R, v) starts at v in a set containing a.
        # This is NOT directly the same as Claim B on a (m-1)-tournament.

        # Alternative: can we apply Claim B on W\{a} with vertex v?
        W_minus_a = [w for w in W if w != a_vertex]
        claim_b_minus_a = 0.0
        Wva2 = [w for w in W_minus_a if w != v]
        for smask in range(1 << len(Wva2)):
            S = [Wva2[bit] for bit in range(len(Wva2)) if smask & (1 << bit)]
            R = [w for w in W_minus_a if w not in S]
            claim_b_minus_a += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        expected_minus_a = -((-1) ** (m - 1)) * h_end(T, W_minus_a, v)
        print(f"    Claim B on W\\{{a}}: value={claim_b_minus_a:.4f}, "
              f"expected={expected_minus_a:.4f}, err={abs(claim_b_minus_a - expected_minus_a):.2e}")


if __name__ == "__main__":
    for m in range(1, 9):
        ok = verify_claim_B(m, num_trials=200 if m <= 6 else 50)
        if not ok:
            break

    verify_claim_B_variant(4)
    verify_claim_B_variant(5)
    prove_claim_B_induction_test(5)
    prove_claim_B_induction_test(6)
