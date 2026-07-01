---
id: HYP-3804
title: A TOURNAMENT-INVARIANT SAFARI (creative new properties) -- three findings. (A) THE PACKING/COVERING ASYMMETRY: define the RAINBOW NUMBER R(n) = max dim of a subcube whose 2^k completions are ALL in DISTINCT iso classes (the packing DUAL of the flip-rank rho(n), HYP-3803, which is the covering min). VERIFIED (exhaustive n=3..6): R(n) = floor(log2|G_n|) EXACTLY (1,2,3,5), while rho(n)=1,2,4,7 EXCEEDS ceil(log2|G_n|) starting at n=6. Gap rho-R = 0,0,1,2 (grows) -- iso classes PACK at the information floor but COVER above the ceiling: covering is strictly harder than packing under the S_n action. (B) |Aut| max over classes = 3,3,5,9 for n=3..6. (C) THE SKEW-SPECTRUM IS A VERY WEAK INVARIANT: the eigenvalues of S=T-T^T (skew-symmetric +-1) give only 1,2,2,6 DISTINCT spectra among |G_n|=2,4,12,56 classes -- almost all tournaments are cospectral (clean NEGATIVE result: T-T^T spectrum is nearly useless for iso-separation). (D) THE QR-TRIANGLE DESIGN: the directed-triangle hypergraph of the quadratic-residue (doubly-regular) tournament on q=prime≡3 mod 4 vertices is a 2-(q,3,(q+1)/4) DESIGN (every vertex-pair in exactly (q+1)/4 directed triangles) -- VERIFIED q=3,7,11,19,23; the merely-rotational tournament {1,2,3} is NOT a design (pair-coverages 1,2,3) => the 2-design property CHARACTERIZES doubly-regularity. A clean tournament<->combinatorial-design bridge
status: MIXED (exhaustive small-n + established-adjacent designs + one negative). VERIFIED (canonicalization n=3..6): R(n)=1,2,3,5=floor(log2|G_n|) exhaustively; rho(n)=1,2,4,7 (HYP-3803); gap rho-R=0,0,1,2. Skew-spectrum distinct counts 1,2,2,6 (exact). QR-triangle 2-(q,3,(q+1)/4) design verified q in {3,7,11,19,23}; rotational-7 not a design. HONEST: R(n)=floor(log2) is a VERIFIED PATTERN (n<=6), conjectured general; the QR design is likely known (doubly-regular tournaments <-> skew-Hadamard/2-designs; cf. five_cycle_3design_s24e) -- freshly verified with the exact lambda=(q+1)/4; the skew-spectrum weakness is a documented negative. Exploratory, not LRC-tied (per the owner's invitation).
source: klein-2026-07-01-S72
depends_on:
  - HYP-3803   # S71: the flip-rank rho(n) (covering); R(n) is its packing dual
related:
  - THM-002    # OCF -- the directed triangles are the degree-3 OCF term; the QR design organizes them
  - HYP-3802   # the QR/Paley tournament on the atoms (doubly-regular, N(OCF)=80)
external: A000568; doubly-regular tournaments <-> skew-Hadamard matrices <-> 2-designs; cospectral graphs/tournaments; packing vs covering codes
results:
  - 04-computation/tournament_invariant_safari_klein.py
  - 05-knowledge/results/tournament_invariant_safari_klein.out
---

# HYP-3804 — a tournament-invariant safari (creative new properties)

Open exploration (owner's invitation: creative tournament properties, even if not LRC-tied), building on
the flip-rank cube model (HYP-3803). Three findings.

## (A) The packing/covering asymmetry: rainbow number vs flip-rank
Define the **rainbow number** `R(n)` = the max dimension of a subcube (fix arcs, free `k`) whose `2^k`
completions are **all in distinct iso classes** — the *packing* dual of the flip-rank `rho(n)` (the
*covering* minimum, HYP-3803). Bounds: `R(n) <= floor(log2|G_n|) <= ceil(log2|G_n|) <= rho(n)`.

| `n` | `\|G_n\|` | `R(n)` | `floor(log2)` | `rho(n)` | `ceil(log2)` | gap `rho-R` |
|---|---|---|---|---|---|---|
| 3 | 2 | 1 | 1 | 1 | 1 | 0 |
| 4 | 4 | 2 | 2 | 2 | 2 | 0 |
| 5 | 12 | 3 | 3 | 4 | 4 | 1 |
| 6 | 56 | 5 | 5 | 7 | 6 | 2 |

> **`R(n) = floor(log2|G_n|)` EXACTLY** (verified n=3..6): the iso classes always *pack* into a subcube at
> the information floor. But `rho(n)` *exceeds* `ceil(log2|G_n|)` starting at `n=6`, and the gap `rho-R`
> grows (0,0,1,2). **Covering all classes is strictly harder than packing distinct ones** under the `S_n`
> action: you can always find a perfect rainbow subcube, but you cannot always cover with the ceiling many
> bits. The symmetry group's folding obstructs covering (forces collisions you must absorb) but not packing
> (you have freedom to dodge collisions).

## (B) Symmetry: |Aut| distribution
Max `|Aut(T)|` over iso classes = `3, 3, 5, 9` for `n=3,4,5,6` (orbit size `= n!/|Aut|`). The most
symmetric tournaments are the circulant/rotational ones; at odd `n` the doubly-regular QR tournament leads.

## (C) The skew-spectrum is a very weak invariant (a clean negative)
For the skew-symmetric `+-1` matrix `S = T - T^T` (eigenvalues `+- i mu_k`), the multiset `{mu_k}` is an
iso-invariant. But it separates almost nothing: the number of **distinct** skew-spectra is `1, 2, 2, 6`
among `|G_n| = 2, 4, 12, 56` classes (`n=3..6`). Nearly all tournaments are **cospectral**. So the
`T-T^T` spectrum is nearly useless for distinguishing tournaments — worth recording as a dead-end invariant
(the skew-spectrum sees little beyond coarse cycle data).

## (D) The QR-triangle design (a tournament <-> design bridge)
The directed 3-cycles of a tournament form a 3-uniform hypergraph on the vertices. For the **quadratic-
residue (doubly-regular) tournament** on `q = prime ≡ 3 (mod 4)` vertices (`i -> j` iff `j-i` is a QR mod
`q`), this hypergraph is a **2-`(q, 3, (q+1)/4)` design**: every vertex-pair lies in exactly `(q+1)/4`
directed triangles. VERIFIED for `q = 3, 7, 11, 19, 23` (`lambda = 1, 2, 3, 5, 6`). The merely-rotational
tournament `{1,2,3}` on 7 vertices is **not** a design (pair-coverages `1,2,3`). So:
> **The directed-triangle hypergraph is a 2-design iff the tournament is doubly-regular** (QR), with
> `lambda = (q+1)/4`.
(Doubly-regular tournaments correspond to skew-Hadamard matrices; this is the design face of that, freshly
verified with the exact `lambda`. Cf. prior `five_cycle_3design_s24e` for the 5-cycle analogue.)

## Net
Three creative tournament properties: (A) a NEW packing/covering asymmetry (`R(n)=floor(log2)` vs
`rho(n)>ceil(log2)` from `n=6`), (B) the symmetry maxima, (C) the skew-spectrum's weakness (negative), (D)
the QR-triangle 2-design characterizing doubly-regularity. (A) is the most novel — the rainbow/flip-rank
duality quantifies how the `S_n` action makes covering harder than packing. Exploratory; not tied to the
LRC (per the owner's invitation), though (D)'s doubly-regular QR tournament is exactly the one on the LRC
atoms (HYP-3802).
