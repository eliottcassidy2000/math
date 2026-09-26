---
id: THM-4500
title: "Pairing Bellman contraction, complete finite policies, and attained natural-density optimum"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: crossroads-family-20260926, kind-pasteur
depends_on:
  - 01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md
  - 01-canon/theorems/THM-4496-pairing-kernel-convolution-and-log-density.md
related:
  - 05-knowledge/results/crossroads_family_20260926_flow.md
  - 01-canon/theorems/THM-4492-pairing-two-cutoff-density-separation.md
scripts:
  - 04-computation/experiments/crossroads_family_20260926_flow.py
  - 04-computation/experiments/crossroads_family_20260926_audit.py
outputs:
  - 05-knowledge/results/crossroads_family_20260926_flow.out
audit: "Root independently derived the local telescope and slow phase sweep. Geometry independently audited all proofs, including natural attainment; exact independent controls covered policies, child minima, L1 comparisons, telescopes, selector stability, and ancestor errors. Detailed scope and counts in the session audit."
---

# THM-4500 — Bellman contraction and an attained natural-density minimum

**PROVED.** Full proof: [flow note, sections 2–8](../../05-knowledge/results/crossroads_family_20260926_flow.md).
This theorem concerns the repository's **modified global two-step pairing
model**, not the ordinary Collatz map. No external novelty claim is made.

## Model and statement

A pairing is a binary sequence epsilon_i with constraints

    epsilon_(2k)+epsilon_(3k)=1,
    epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).

The harmless root bit is set to zero. These are precisely the global P2
constraints of THM-4491 (cited with the dependency path above). Let A_F(X)
count indices i<=X with epsilon_i=1. Let B_* be the attained minimum
logarithmic density / finite-kernel supremum identified by THM-4496.

With r=2/3 and probability Haar measure on Z_2 define

    (F d)(2k)   =1-r d(3k),
    (F d)(2k+1) =1+r min(d(3k+1),0)+r max(d(3k+2),0).

Then F contracts L1 by r and has integral 1 on every input. It has a
unique L1 fixed point d_*, essentially bounded by 3. For h>=1, if d_h=F^h(0),
e_h=||d_(h+1)-d_h||_1 and b_h=(1+integral min(d_h,0))/3, then

    B_*=(1+integral min(d_*,0))/3,
    b_h-e_h/2 <= B_* <= delta_h <= b_h+e_h/2,
    e_h <= (2/3)^h.                                  (1)

Here delta_h is the existing natural density of the explicit legal policy
that sets a free child bit to 1 iff d_h at that child's residue is negative.
Ties select zero. The period is at most 2^(h-1).

The exact finite-policy LP on M=2^s states has constraints

    0<=x_r<=1,
    3x_r+x_(2r/3)-x_((2r-1)/3)>=1,
    3x_r+x_(2r/3)-x_((2r+1)/3)<=2,                    (2)

where all arguments are residues modulo M. Every vertex is a pure coupled
free-choice policy. Its optimum V_s satisfies

    B_*<=V_(s+1)<=V_s<=B_*+(2/3)^(s+1).              (3)

There is **one global legal pairing** with

    A_F(X)=B_* X+O(X/(log X)^(1/3)).                   (4)

Consequently B_* is a minimum, not merely an infimum, of both existing
natural densities and upper natural densities among global P2 pairings.
The smaller minimum lower natural density alpha of THM-4491 is a different
quantity; THM-4492's obstruction to attaining alpha as a natural density
is not contradicted.

## Mechanisms and proof boundaries

For each fixed residue policy, an ancestor step has one reset branch and
two unresolved branches among three equally frequent CRT cylinders. Thus
unresolved mass after k steps is (2/3)^k. This proves existence of natural
density and, for period M and c=log_3 2,

    |A_F(X)-delta X| <=5 M^(1-c) X^c, X>=M.

For a bounded periodic potential d, telescope the local terms

    epsilon_i/i + sum_(child j) epsilon_j d_j/j -epsilon_i d_i/i.

Child reciprocal weights are (2/3)/i+O(i^-2). Local minimization has
parent-bit difference Fd-d and zero-parent mean
(integral d+integral min(d,0))/3. Summable errors and the outer fringe
contribute O(1), giving a universal harmonic lower certificate and a
matching sign-policy upper certificate. L1 contraction yields (1).
The contraction is **not** valid in the supremum norm; an exact hostile
is retained. Fixed-policy means, separately, do contract in that norm.

For (4), write t=log_(3/2)n and phi={t}. Interpolate S(h^3)=h linearly,
and at address n use the free bit of policy h=floor(S(t)+1-phi).
The actual ancestor satisfies p^j(n)=(2/3)^j n+O(1), so its log phase
stays nearly constant over k=floor((log_(3/2)X)^(1/3)) steps for n>=sqrtX.
Except for O(X/(log X)^(1/3)) indices near selector boundaries, those
ancestors use the same stationary policy. Resets erase earlier choices
on all but o(X/(log X)^(1/3)) further indices. The reference made by
patching stationary policies has the asserted density, using their
uniform discrepancy bounds on O(log X) intervals. The full proof counts
boundary bands in ordinary height, handles phase wraps, and gives
an exact integer-power implementation of the selector.

Abrupt height-band changes and pointwise choices of Haar equivalence-class
representatives do not justify this argument. Slow phase variation and
uniform reset estimates supply the missing control.

## Exact numerical certificate and scope

At h=24 the exact enclosure is

    750086585295922675/2369190669160808448 <= B_*
        <=750158422907862637/2369190669160808448.

Outward decimal bounds are 0.3166003458 < B_* < 0.3166306675.
The modulus-32 policy optimum is exactly 8209/25920, with a rational
primal/dual certificate. The deep certificate uses all 2^24 residual
states and guarded integer arithmetic. It is backed by the uniform
error theorem, not extrapolation from a numerical center.

No claim is made that B_* is rational, that a single finite-period policy
attains it, or that these modifications prove ordinary Collatz descent.
