# THM-026: Deletion-Sum Ratio Formula

**Status:** PROVED (algebraic identity)
**Date:** 2026-03-06
**Author:** kind-pasteur-2026-03-06-S18g

## Statement

For any tournament T on n vertices with H(T) > 0:

R(T) = sum_v H(T-v) / H(T) = n - E_weighted[|U(S)|]

where the weighted expectation is:

E_weighted[|U(S)|] = (sum_{S independent in Omega(T)} 2^{|S|} * |U(S)|) / (sum_{S independent} 2^{|S|})

and |U(S)| = |union of V(C) for C in S| = total tournament vertices covered by cycles in the independent set S.

Equivalently:

R(T) = n - U_sum / H(T)

where U_sum = 2 * sum_C |V(C)| * mu(C), the weighted cycle sum over all directed odd cycles C, with weights |V(C)| * mu(C) where mu(C) = H(T[V \ V(C)]).

## Proof

Step 1: By OCF, H(T-v) = I(Omega(T-v), 2). The cycles of T-v are exactly the cycles of T not passing through v. So Omega(T-v) = Omega(T) restricted to cycles avoiding v.

Step 2: An independent set S in Omega(T) contributes to H(T-v) iff all cycles in S avoid v, i.e., v not in U(S). This holds for exactly n - |U(S)| values of v.

Step 3: Therefore sum_v H(T-v) = sum_S 2^{|S|} * (n - |U(S)|) = n * H(T) - U_sum.

Step 4: Dividing by H(T): R(T) = n - U_sum / H(T) = n - E_weighted[|U(S)|].

Step 5: For the U_sum formula via Claim A: U_sum = sum_v [H(T) - H(T-v)] = sum_v delta_v. By Claim A, delta_v = 2 * sum_{C through v} mu(C). So U_sum = 2 * sum_C |V(C)| * mu(C) (each cycle C counted once per vertex it contains).

## Conjecture (OPEN-Q-017) — REFUTED AT n=7

The H-maximizer MINIMIZES R(T) among all tournaments with H(T) > 0.

Confirmed exhaustively: n=3,4,5,6. **FAILS at n=7**: tournaments with H=123 achieve R=1.585 < 5/3 = 1.667 (the maximizer's R). Many violations found. The R-minimizer at n=7 is NOT the H-maximizer.

Possible weaker versions to investigate:
- Does the maximizer minimize R within its score class?
- Is R monotonically related to H within score classes?

## Verification Scripts

- `04-computation/R_minimization_proof.py` — verifies R formula
- `04-computation/deletion_ratio_theorem.py` — exhaustive R-min check n=3-6
- `04-computation/R_min_n7_fast.py` — n=7 check (in progress)
- `04-computation/R_min_focused.py` — IP decomposition of U_sum
