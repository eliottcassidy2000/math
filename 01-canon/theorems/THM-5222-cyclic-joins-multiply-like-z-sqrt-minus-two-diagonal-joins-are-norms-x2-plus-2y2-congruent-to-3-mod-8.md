---
id: THM-5222
title: "Cyclic joins compose like Z[sqrt(-2)]: for A => B => c => A, H = (h_A + P_A)(h_B + P_B) + 2 h_A h_B = Re(beta_A conj(beta_B)) with beta = (h + P_2) + h sqrt(-2), where P_2 counts ordered covers by two disjoint nonempty paths; diagonal joins (B a copy of A) are norms x^2 + 2y^2, always congruent to 3 mod 8, so they realize a prime p only if p = 3 mod 8, and then only through the unique representation p = x^2 + 2y^2 with H(A) = y and P_2(A) = x - y"
status: >
  PROVED (path decomposition below). FINITE-EXACT: the formula matches on 200 of 200 random core pairs (orders 1-5); the
  diagonal primes found are all 3 mod 8. The general (non-diagonal) form is equivalent to the cyclic-join law
  3N = Re(alpha_A conj(alpha_B)), alpha = (f - 3h) + f sqrt(-2), of THM-5172 (C), with f = 2h + P_2. That identification
  of f, and the norm reading, are new here.
session: mac-mini-2026-11-15-norms (wave 34 of chanpak_waves_20261113)
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 34)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w34_cyclic_join_norm_law.py
related:
  - THM-5172 (C) (cyclic join in Z[sqrt(-2)]), THM-5133 (two-port kernels), HYP-9692 (band fullness)
---

# THM-5222 - cyclic joins and norms from Z[√−2]

## Proof

Let T be the cyclic join: A ⇒ B, B ⇒ c and c ⇒ A, with c a single vertex. Since c is visited once and A ⇒ B, every
Hamiltonian path has the form

    A₁ B₁ c A₂ B₂,

where (A₁, A₂) are disjoint paths covering A and (B₁, B₂) are disjoint paths covering B, each possibly empty.

- *Constraints.* c beats all of A, so a nonempty A₁ must be followed by a nonempty B₁. All of B beats c, so a nonempty
  B₂ must follow a nonempty A₂.
- *Counting the cases.*
  - A₁ = ∅ and A₂ a full path (h_A ways): B is arbitrary, giving 2h_B + P_B.
  - A₁ a full path and A₂ = ∅ (h_A ways): B₁ must be full, giving h_B.
  - Both A₁ and A₂ nonempty (P_A ways): B₁ ≠ ∅, giving h_B + P_B.
- *Total.* Summing, N = (h_A + P_A)(h_B + P_B) + 2h_A h_B. This is the real part of
  ((h_A + P_A) + h_A√−2)·conj((h_B + P_B) + h_B√−2).

## Diagonal joins

Take B a copy of A. Then N = (h + P₂)² + 2h², the norm of β = (h + P₂) + h√−2.
- *Residue.* By Rédei h is odd, so h² ≡ 1 (mod 8). Since N is odd, P₂ is even and (h + P₂)² ≡ 1 (mod 8). Hence
  N ≡ 3 (mod 8).
- *Primes.* A prime p is represented by x² + 2y² iff p ≡ 1 or 3 (mod 8), uniquely up to sign since Q(√−2) has class
  number 1. So a diagonal join realizes p only when p ≡ 3 (mod 8). Then x and y are odd, and the target becomes the
  **two-port equation** H(A) = y, P₂(A) = x − y.
- *Inert primes.* In any diagonal-join count, every prime q ≡ 5, 7 (mod 8) divides N to an even power.

## General joins: prime targets as linear two-port equations

Fix B, so that X_B = h_B + P_B and Y_B = h_B. Then N = X_B·(h_A + P_A) + 2Y_B·h_A is a linear form in the two ports
of A. A prime target p is reached as soon as the pair spectrum {(H(A), P₂(A))} meets the line
X_B·x + 2Y_B·y = p in a point with y = h_A and x = h_A + P_A.

So Chan–Pak for primes follows from a two-dimensional statement: the pair spectrum meets every such line for some
small fixed B.
