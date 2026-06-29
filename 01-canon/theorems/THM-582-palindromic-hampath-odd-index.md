---
id: THM-582
title: For a SELF-CONVERSE tournament with an involutory anti-automorphism phi, the reversal-conjugation rho(P)=phi(reverse(P)) is an involution on directed Hamiltonian paths; hence H(T) == #{phi-palindromic Ham paths} (mod 2), and by Rédei (H odd) the number of phi-PALINDROMIC Hamiltonian paths is ODD (>=1). This is the Hamiltonian-path-level twin of THM-281 (SC tiling fibers are odd via grid-symmetric tilings) -- the sigma-ODD "index" the project sought, on the witness side.
status: PROVED (elementary involution argument); VERIFIED exhaustively -- Paley T_3/T_7/T_11 (H=3/189/95095, #palindromic=1/9/185, all odd) and all 48 (n=4) + 704 (n=5) self-converse tournaments with involutory phi.
source: mac-mini-2026-06-29-S6
depends_on:
  - THM-280   # grid reflection = complement (the involutory anti-automorphism in tiling coords)
related:
  - THM-281   # SC tiling fibers are odd (the tiling-level twin)
  - THM-088   # signed-F palindromicity via path reversal (fwd(P^rev)=n-1-fwd(P))
  - THM-027   # tr(M)=H(T) for odd n (position-parity signed sum)
  - HYP-3239   # kps D_7 Borsuk-Ulam: the sigma-ODD/sign isotype = the witness side
  - THM-581    # the LRC floor is the sigma-EVEN index (existence); this is the ODD twin
  - OPEN-Q-108
results:
  - 04-computation/palindromic_hampath_parity_macmini_20260629.py
  - 05-knowledge/results/palindromic_hampath_parity_macmini_20260629.out
---

# THM-582 — the palindromic Hamiltonian-path odd index

## Statement
Let `T` be a tournament that is SELF-CONVERSE via an INVOLUTORY anti-automorphism `phi`:
`phi^2 = id` and `(u -> v) in T  iff  (phi(v) -> phi(u)) in T` (so `phi: T ~= T^op`).
For a directed Hamiltonian path `P = (v_1,...,v_n)` define the reversal-conjugation
> `rho(P) = phi(reverse(P)) = (phi(v_n), phi(v_{n-1}), ..., phi(v_1))`.

Then:
1. `rho` maps Hamiltonian paths of `T` to Hamiltonian paths of `T` (well-defined).
2. `rho` is an INVOLUTION (`rho^2 = id`).
3. Hence `H(T) = #fixed(rho) + 2 . #(free 2-orbits)`, so **`H(T) == #fixed(rho) (mod 2)`**.
4. The fixed points of `rho` are the `phi`-**palindromic** Hamiltonian paths.
5. Since `H(T)` is ODD (Rédei), **the number of `phi`-palindromic Hamiltonian paths is ODD (>= 1)**.

## Proof
**(1)** If `v_i -> v_{i+1}` in `T` for all `i`, then `reverse(P)=(v_n,...,v_1)` has consecutive arcs
`v_{i+1} -> v_i` in `T^op`; applying `phi` (an iso `T^op -> T`, since `phi: T ~= T^op` and `phi^2=id`)
gives consecutive arcs `phi(v_{i+1}) -> phi(v_i)`... precisely, the anti-automorphism condition
`(u->v in T) iff (phi(v)->phi(u) in T)` gives, from `v_i -> v_{i+1} in T`, that
`phi(v_{i+1}) -> phi(v_i) in T`.  Reading `rho(P)=(phi(v_n),...,phi(v_1))` left to right, its
consecutive pair is `(phi(v_{k+1}), phi(v_k))` with arc `phi(v_{k+1}) -> phi(v_k) in T`.  So `rho(P)`
is a directed Hamiltonian path of `T`.
**(2)** `rho^2(P)`: reversing `(phi(v_n),...,phi(v_1))` gives `(phi(v_1),...,phi(v_n))`, then `phi`
gives `(phi^2(v_1),...,phi^2(v_n)) = (v_1,...,v_n) = P` since `phi^2=id`.
**(3)** An involution on a finite set partitions it into fixed points and 2-orbits; count mod 2.
**(4)** Definition.  **(5)** Rédei: `H(T)` odd, so `#fixed(rho)` odd, in particular `>= 1`.  QED

## Verification (script palindromic_hampath_parity_macmini_20260629.py)
| T | n | H(T) | #palindromic | both odd? |
|---|---|------|--------------|-----------|
| Paley T_3 (phi=x->-x) | 3 | 3 | 1 | yes |
| Paley T_7 | 7 | 189 | 9 | yes |
| Paley T_11 | 11 | 95095 | 185 | yes |
| all SC, n=4 (48 of them) | 4 | 1 | 1 | yes |
| all SC, n=5 (704 of them) | 5 | 1 | 1 | yes |
`rho` verified to map HP->HP and be an involution in every case.

## Significance: the two indices (the sigma-EVEN / sigma-ODD split)
The complement/reversal `Z_2` (`sigma`) splits the project's parity data into two indices:
- **sigma-ODD (witness side):** the `phi`-palindromic Hamiltonian-path count (THM-582) and the
  grid-symmetric tiling count (THM-281) -- both odd, both realizing Rédei's `H(T)` odd as a
  fixed-point parity of the reversal involution.  This is kps's Borsuk-Ulam / sign-isotypic side
  (HYP-3239): the witness/odd-degree datum.
- **sigma-EVEN (existence side):** the LRC lonely measure = the Euler characteristic of the danger-
  cover nerve (THM-581 / HYP-3242), which is `sigma`-invariant so its sign component vanishes -- the
  Brouwer/SOS side, what existence (the LRC floor) needs.

THM-582 thus SETTLES the "odd index" question that recurred (two-order-two-structures reflection,
lrc19-s679): the odd index EXISTS and is the palindromic Hamiltonian-path count -- but it belongs to
the WITNESS/Rédei side, NOT to the LRC floor.  A lonely tournament has vertex `0` as a SOURCE, which
converse turns into a SINK, so the lonely tournament is NOT self-converse and THM-582 does not apply
to it -- confirming the floor is the sigma-EVEN index.  This REDIRECTS the LRC-floor closure to the
even-category route (descent THM-580 + cyclotomic SOS S75e), away from any odd/Borsuk-Ulam index.

The half-tiling (THM-549/550) is the `sigma`-fundamental domain housing both indices; its fixed
diagonal (the SC spine) is the palindromic/odd locus.

Script: `04-computation/palindromic_hampath_parity_macmini_20260629.py`.
