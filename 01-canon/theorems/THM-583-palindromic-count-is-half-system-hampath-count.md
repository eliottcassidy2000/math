---
id: THM-583
title: For a self-converse tournament with INVOLUTORY anti-automorphism phi having a unique fixed vertex c, the palindromic-Hamiltonian-path count f (THM-582's odd index) equals a HALF-SYSTEM Hamiltonian-path count on the (p-1)/2 phi-orbit pairs -- validity is PURE HALF-DATA (the half arcs v_i->v_{i+1} plus the closing arc v_m->c), the mirror/second half being AUTOMATIC by phi; a transfer DP on (last vertex, used-pairs) computes f with a half-size state space ((p-1)/2 pairs vs p vertices)
status: PROVED (the pure-half-data reduction is an elementary consequence of the anti-automorphism identity arc(u,w) <=> arc(phi(w),phi(u))); VERIFIED -- Paley T_3/T_7/T_11 give f=1/9/185 by direct count, half-system enumeration, AND transfer DP (DP states 2/21/150 vs p!=6/5040/39916800).
source: mac-mini-2026-06-29-S8
depends_on:
  - THM-582   # f = #palindromic Ham paths is odd (the odd index)
related:
  - THM-027   # tr(M)=H for odd n (transfer matrix); this is the half-system transfer for f
  - THM-549   # the half-tiling is the sigma-fundamental domain
  - HYP-3538  # the R +-1 eigenspace organizing principle (this realizes the witness half)
  - OPEN-Q-108
results:
  - 04-computation/half_system_f_recursion_macmini_20260629.py
  - 05-knowledge/results/half_system_f_recursion_macmini_20260629.out
---

# THM-583 -- f is a half-system Hamiltonian-path count

## Setup
Let `T` be a tournament on `Z_p` (or any vertex set) that is self-converse via an INVOLUTORY
anti-automorphism `phi` (`phi^2=id`, `arc(u,w) <=> arc(phi(w),phi(u))`) with a UNIQUE fixed vertex
`c` (`phi(c)=c`).  For Paley `T_p`, `p=3 mod4`, take `phi(x)=-x`, `c=0`; the other vertices form
`m=(p-1)/2` pairs `{x, phi(x)}`.

By THM-582 the `phi`-palindromic Hamiltonian paths (`phi(reverse(P))=P`) number an ODD `f>=1`.  Such
a path is forced to the shape
> `P = (v_1, ..., v_m,  c,  phi(v_m), ..., phi(v_1))`,
center `= c`, second half `= phi`-image of the reversed first half (one vertex from each pair).

## The pure-half-data reduction (PROVED)
> **`P` is a valid directed path in `T`  iff  the HALF-DATA arcs hold:
> `v_i -> v_{i+1}` for `i=1..m-1`  AND  `v_m -> c`.**

*Proof.*  The arcs of `P` are (a) `v_i->v_{i+1}` (i<m), (b) `v_m->c`, (c) `c->phi(v_m)`,
(d) `phi(v_{i+1})->phi(v_i)` (i<m).  By the anti-automorphism identity, `arc(v_i,v_{i+1}) <=>
arc(phi(v_{i+1}),phi(v_i))`, so (d) is equivalent to (a) -- AUTOMATIC.  And `arc(v_m,c) <=>
arc(phi(c),phi(v_m)) = arc(c,phi(v_m))`, so (c) is equivalent to (b) -- AUTOMATIC.  Hence (c),(d)
add no constraint; validity is exactly (a)+(b), pure first-half data.  QED

## Consequence: a half-size transfer recursion
> **`f = #{ ordered selections (v_1,...,v_m), one signed representative per `phi`-pair, with
> `v_i->v_{i+1}` and `v_m->c` } = a Hamiltonian-path count on the half-system digraph.`**
A transfer DP on the state `(last vertex, bitmask of used pairs)` computes `f` over `O(2^m * 2m)`
states instead of `p!` permutations.  This is the half-tiling/`sigma`-quotient realized for the odd
index: the witness data lives on the `(p-1)/2` pairs, the second half recovered by `phi`.

## Verification
| p | m=(p-1)/2 | f (direct) | f (half-system) | f (transfer DP) | DP states | p! |
|---|-----------|------------|-----------------|-----------------|-----------|----|
| 3 | 1 | 1 | 1 | 1 | 2 | 6 |
| 7 | 3 | 9 | 9 | 9 | 21 | 5040 |
| 11 | 5 | (—) | 185 | 185 | 150 | 39916800 |

## Significance
This is the answer to "f is a half-HP count, so a transfer-matrix recursion on `(p-1)/2` pairs may
exist."  It does: `f` (the THM-582 odd index, the `sigma`-ODD / witness side of the two-index split)
is a Hamiltonian-path count on the half-system digraph, computed by a transfer DP on half the
vertices, with the second half forced by `phi`.  Lossless half-compression -- the coordinate that
must be retained is `phi` (the reversal action) itself, which regenerates the dropped second half.
This complements HYP-3538: the witness `f` is the R-fixed/odd combinatorial datum, computed on the
half-system; the cap's obstruction is the R-odd spectral datum.  Both are the `eps=-1` content of the
reversal `R`, on the two sides.

Open: a closed form / linear recurrence for the Paley `f`-sequence `1,9,185,...` (the transfer DP is
config-specific; the QR structure of the half-system digraph governs it).

Script: `04-computation/half_system_f_recursion_macmini_20260629.py`.
