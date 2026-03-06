#!/usr/bin/env python3
"""
PROOF of Claim (B) by induction on m.

Claim (B): For any tournament on m vertices W with distinguished vertex v:

  Phi(v) := sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)

PROOF:

Base case m=1: W = {v}. LHS = H(empty)*h_start({v},v) = 1*1 = 1. RHS = (-1)^2 * 1 = 1. Done.

Inductive step: Assume Claim B for all sizes < m. Let W' = W\{v} (size m-1).

Separate the S = W' term (where R = W\S = {v}):
  Phi(v) = (-1)^{m-1} H(W') * 1  +  sum_{S ⊊ W'} (-1)^|S| H(S) h_start(W\S, v)

For S ⊊ W', the set R = W\S has |R| >= 2 and contains v. Apply first-step decomposition:
  h_start(R, v) = sum_{u in R\{v}} T(v,u) h_start(R\{v}, u)

Note R\{v} = W'\S (since S doesn't contain v). Substituting:
  Phi(v) = (-1)^{m-1} H(W') + sum_{S ⊊ W'} (-1)^|S| H(S) sum_{u in W'\S} T(v,u) h_start(W'\S, u)

Exchange summation order (for each u in W', sum over S ⊆ W'\{u}):
  = (-1)^{m-1} H(W') + sum_{u in W'} T(v,u) [sum_{S ⊆ W'\{u}} (-1)^|S| H(S) h_start(W'\S, u)]

The inner sum [...] is EXACTLY Phi_{W'}(u) — Claim B for the (m-1)-vertex set W' with
distinguished vertex u! By the induction hypothesis:
  sum_{S ⊆ W'\{u}} (-1)^|S| H(S) h_start(W'\S, u) = (-1)^m h_end(W', u)

(Using |W'| = m-1, so (-1)^{(m-1)+1} = (-1)^m.)

Therefore:
  Phi(v) = (-1)^{m-1} H(W') + sum_{u in W'} T(v,u) * (-1)^m * h_end(W', u)
         = (-1)^{m-1} H(W') + (-1)^m * sum_{u in W'} T(v,u) h_end(W', u)

KEY STEP: Use T(v,u) = 1 - T(u,v) to evaluate the sum:
  sum_{u in W'} T(v,u) h_end(W', u) = sum_u (1-T(u,v)) h_end(W', u)
    = sum_u h_end(W', u) - sum_u T(u,v) h_end(W', u)
    = H(W') - h_end(W, v)

where the last equality uses:
  - sum_u h_end(W', u) = H(W')  (total Ham path weight = sum over ending vertices)
  - sum_u T(u,v) h_end(W', u) = h_end(W, v)  (any Ham path on W ending at v is a
    Ham path on W' ending at some u, followed by the step u→v)

Substituting back:
  Phi(v) = (-1)^{m-1} H(W') + (-1)^m (H(W') - h_end(W, v))
         = (-1)^{m-1} H(W') + (-1)^m H(W') + (-1)^{m+1} h_end(W, v)
         = [(-1)^{m-1} + (-1)^m] H(W') + (-1)^{m+1} h_end(W, v)
         = 0 + (-1)^{m+1} h_end(W, v)

since (-1)^{m-1} + (-1)^m = 0 always. QED.

This script verifies each step of the proof numerically.

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
    if len(verts) <= 1:
        return 1.0
    total = 0.0
    for p in permutations(verts):
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def verify_proof_step_by_step(m, num_trials=100):
    """Verify each step of the inductive proof numerically."""
    print(f"=== Proof verification at m={m} ===\n")

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
            return arc_values[(a, b)]

        v = W[0]
        W_prime = [w for w in W if w != v]

        # Step 1: Compute Phi(v) directly
        phi_direct = 0.0
        for smask in range(1 << len(W_prime)):
            S = [W_prime[bit] for bit in range(len(W_prime)) if smask & (1 << bit)]
            R = [w for w in W if w not in S]
            phi_direct += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        # Step 2: Separate boundary term
        H_W_prime = H_total(T, W_prime)
        boundary = ((-1) ** (m - 1)) * H_W_prime

        # Step 3: Apply induction hypothesis to inner sums
        induction_sum = 0.0
        for u in W_prime:
            # Inner sum = Phi_{W'}(u) = (-1)^m h_end(W', u) by induction
            phi_u_induction = ((-1) ** m) * h_end(T, W_prime, u)

            # Verify induction hypothesis holds for this sub-problem
            phi_u_direct = 0.0
            W_prime_minus_u = [w for w in W_prime if w != u]
            for smask in range(1 << len(W_prime_minus_u)):
                S = [W_prime_minus_u[bit] for bit in range(len(W_prime_minus_u))
                     if smask & (1 << bit)]
                R = [w for w in W_prime if w not in S]
                phi_u_direct += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, u)

            if trial == 0 and u == W_prime[0]:
                err_ih = abs(phi_u_direct - phi_u_induction)
                print(f"  Induction hypothesis check: Phi_{{W'}}({u}) direct={phi_u_direct:.6f}, "
                      f"induction={phi_u_induction:.6f}, err={err_ih:.2e}")

            induction_sum += T(v, u) * phi_u_induction

        # Step 4: Key identity sum_u T(v,u) h_end(W',u) = H(W') - h_end(W,v)
        sum_Tvu_hend = sum(T(v, u) * h_end(T, W_prime, u) for u in W_prime)
        h_end_W_v = h_end(T, W, v)
        key_identity_lhs = sum_Tvu_hend
        key_identity_rhs = H_W_prime - h_end_W_v

        if trial == 0:
            print(f"  Key identity: sum T(v,u)h_end(W',u) = {key_identity_lhs:.6f}, "
                  f"H(W')-h_end(W,v) = {key_identity_rhs:.6f}, "
                  f"err = {abs(key_identity_lhs - key_identity_rhs):.2e}")

        # Step 5: Combine
        phi_from_proof = boundary + induction_sum
        # = (-1)^{m-1} H(W') + (-1)^m (H(W') - h_end(W,v))
        # = (-1)^{m-1} H(W') + (-1)^m H(W') + (-1)^{m+1} h_end(W,v)
        # = 0 + (-1)^{m+1} h_end(W,v)
        phi_expected = ((-1) ** (m + 1)) * h_end_W_v

        err = abs(phi_direct - phi_expected)
        max_err = max(max_err, err)

        if trial < 3:
            print(f"  trial {trial}: Phi(v)={phi_direct:.6f}, "
                  f"(-1)^{{m+1}}h_end(W,v)={phi_expected:.6f}, err={err:.2e}")
            print(f"    boundary={boundary:.6f}, induction_sum={induction_sum:.6f}, "
                  f"proof_total={phi_from_proof:.6f}")

    print(f"\n  Max error: {max_err:.2e}")
    if max_err < 1e-8:
        print(f"  CONFIRMED at m={m}\n")
    return max_err < 1e-8


def verify_key_identity(m, num_trials=200):
    """
    Verify: sum_{u in W'} T(v,u) h_end(W', u) = H(W') - h_end(W, v)

    This follows from:
    1. sum_u h_end(W', u) = H(W')  (partition of Ham paths by ending vertex)
    2. sum_u T(u,v) h_end(W', u) = h_end(W, v)  (last-step decomposition)
    3. T(v,u) = 1 - T(u,v)  (tournament complement)
    """
    print(f"=== Key identity verification at m={m} ===")

    W = list(range(m))
    arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)
    max_err = [0.0, 0.0, 0.0]

    for trial in range(num_trials):
        arc_values = {}
        for (a, b) in arcs:
            val = random.uniform(-2, 3)
            arc_values[(a, b)] = val
            arc_values[(b, a)] = 1 - val

        def T(a, b):
            if a == b:
                return 0
            return arc_values[(a, b)]

        v = W[0]
        W_prime = [w for w in W if w != v]

        H_Wp = H_total(T, W_prime)
        h_end_W_v = h_end(T, W, v)

        # Identity 1: sum_u h_end(W', u) = H(W')
        sum_hend = sum(h_end(T, W_prime, u) for u in W_prime)
        max_err[0] = max(max_err[0], abs(sum_hend - H_Wp))

        # Identity 2: sum_u T(u,v) h_end(W', u) = h_end(W, v)
        sum_Tuv_hend = sum(T(u, v) * h_end(T, W_prime, u) for u in W_prime)
        max_err[1] = max(max_err[1], abs(sum_Tuv_hend - h_end_W_v))

        # Main identity: sum_u T(v,u) h_end(W', u) = H(W') - h_end(W, v)
        sum_Tvu_hend = sum(T(v, u) * h_end(T, W_prime, u) for u in W_prime)
        max_err[2] = max(max_err[2], abs(sum_Tvu_hend - (H_Wp - h_end_W_v)))

    print(f"  Identity 1 (sum h_end = H): max err = {max_err[0]:.2e}")
    print(f"  Identity 2 (last-step decomp): max err = {max_err[1]:.2e}")
    print(f"  Main identity: max err = {max_err[2]:.2e}\n")


def verify_cancellation(m, num_trials=200):
    """
    Verify the crucial cancellation: (-1)^{m-1} + (-1)^m = 0.

    This means the H(W') terms vanish:
    (-1)^{m-1} H(W') + (-1)^m H(W') = [(-1)^{m-1} + (-1)^m] H(W') = 0.
    """
    print(f"=== Cancellation check at m={m} ===")
    cancel = (-1) ** (m - 1) + (-1) ** m
    print(f"  (-1)^{{{m-1}}} + (-1)^{{{m}}} = {cancel}")
    assert cancel == 0
    print(f"  Confirmed: cancellation holds.\n")


if __name__ == "__main__":
    print("PROOF OF CLAIM (B) — NUMERICAL VERIFICATION\n")
    print("=" * 60)
    print("Each step of the inductive proof is verified numerically.")
    print("=" * 60 + "\n")

    # Verify the cancellation lemma
    for m in range(1, 10):
        verify_cancellation(m)

    # Verify key sub-identities
    for m in range(2, 8):
        verify_key_identity(m)

    # Verify full proof step by step
    for m in range(1, 8):
        verify_proof_step_by_step(m, num_trials=100 if m <= 5 else 30)
