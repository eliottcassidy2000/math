# THM-015: Swap Involution Polynomial Identity

**Type:** Theorem (proved at n<=6 by exhaustive symbolic verification)
**Certainty:** 5 -- PROVED at n<=6
**Status:** PROVED at n<=6
**Added by:** kind-pasteur-2026-03-05-S6
**Tags:** #ocf #arc-reversal #involution #proof #polynomial-identity

---

## Statement

Let T be a tournament on n vertices. Fix an arc i→j and let T' be the tournament obtained by flipping this arc to j→i.

Define the **swap involution** on Hamiltonian paths using arc i→j: swap the positions of i and j. A path is **unmatched** if the swapped path is invalid (blocked by a predecessor or successor).

Let U_T = #{unmatched T-paths} and U_{T'} = #{unmatched T'-paths}.

**Theorem:** For all n <= 5 (and conjecturally all n):

U_{T'} - U_T = delta_I

where delta_I = I(Omega(T'), 2) - I(Omega(T), 2) is the change in the independence polynomial.

Equivalently:

U_{T'} - U_T = 2 * sum_x s_x + 2*(D5-C5) + higher corrections

where s_x = 1 - T[x][i] - T[j][x] for x in V\{i,j}.

---

## Proof Method

**Polynomial identity over arc variables.** The key insight is that U_T and U_{T'} can be expressed as polynomials in the arc variables of T. Specifically:

1. Enumerate all "path shapes" — permutations with i→j (resp. j→i) consecutive
2. Each shape's validity is a monomial in arc variables
3. Each shape's blocking condition is a Boolean function of arc variables
4. U_T = sum over shapes of (validity AND blocking), a polynomial

The difference U_{T'} - U_T is obtained by the substitution i↔j applied to the arc variables, then subtracted.

**n=4 proof (by hand):**
With V\{i,j} = {a,b}, define p_x = T[x][i], q_x = T[j][x]:

U_T = p_a*q_b*(q_a+p_b-q_a*p_b) + p_b*q_a*(q_b+p_a-q_b*p_a) + p_a*q_a + p_b*q_b

After algebraic simplification, U_{T'} - U_T = 4 - 2*(p_a+p_b+q_a+q_b) = 2*(s_a+s_b).

This holds as a polynomial identity, not just over {0,1}.

**n=5 proof (symbolic verification):**
512 variable assignments checked exhaustively. U_{T'} - U_T = 2*sum(s_x) + 2*(D5-C5) for all assignments.

---

## Significance

This proof method converts the OCF identity (H(T) = I(Omega(T), 2)) into a finite polynomial identity that can be verified by exhaustive evaluation. Combined with the base case H(transitive) = 1 = I(empty, 2) and the arc-flip reachability of any tournament from the transitive tournament, this proves OCF for all n where the polynomial identity is verified.

**Proved so far:** OCF for n <= 6.
**Next target:** n = 7 (2^19 = 524288 cases, feasible but slower).

---

## Connection to THM-013 and THM-014

- THM-014 (swap involution): establishes U_{T'} - U_T = delta_H
- THM-013 (arc-flip formula): establishes delta_I = sum 2^k Delta(alpha_k)
- THM-015 (this theorem): proves delta_H = delta_I by showing U_{T'} - U_T = delta_I

Together these give: H(T) = I(Omega(T), 2) for all tournaments with n <= 5 (and n <= 6 pending verification).

---

## Verification Record

| n | Variable assignments | Path shapes | Result |
|---|---------------------|-------------|--------|
| 4 | 32 (2^5) | 6 | 32/32 PASS |
| 5 | 512 (2^9) | 24 | 512/512 PASS |
| 6 | 16384 (2^14) | 120 | 16384/16384 PASS |

---

## Open Question

Can this polynomial identity be proved for ALL n simultaneously? The identity
U_{T'} - U_T = delta_I holds as a multilinear polynomial identity in the arc variables.
A general proof would establish OCF (and hence Claim A) for all n.

Possible approaches:
1. Transfer matrix / permanent expansion of both sides
2. Inductive argument on n using the recursive structure of delta_I
3. Algebraic geometry (multilinear polynomials over {0,1} = F_2 affine variety)
