# Orlik–Solomon is a repo-wide pattern — and toric arrangements are the native LRC lens

> **SCOPE CORRECTION — MISTAKE-223 / THM-2047 (2026-07-21).** The braid and
> Shi point-count identities printed here are classical and the THM-1820
> Fourier relation-lattice formula is valid after Fejer regularization. Their
> advertised LRC identifications are not. `G_delta` is a thickened coordinate-
> slab complement pulled back to a one-parameter orbit, not a standard toric-
> arrangement complement; the Fourier kernel sum is not an arithmetic-Mobius
> layer sum; Shi walls do not encode the safe inequality; and the finite
> `N_R` triple census is not a Betti/Mobius invariant or a twelve-speed AP
> theorem. The exact arrangement carrier is the signed phase-height complex
> `v t +/- delta in Z`. Its horizontal slice preserves the pointwise predicate,
> and `chi(G_delta)>0` detects even isolated tight witnesses that volume misses.
> Read THM-2047 before using the exploratory text below.

*boxeph-2026-07-21-S209. Owner: explore creatively other areas where Orlik–Solomon could be leveraged;
think abstractly, look for similar structural patterns. Builds on boxeph-S208 (braid arrangement = NC2
Vandermonde; flat-localization → HYP-8775), THM-1820 (LRC relation-lattice pairing), THM-805 (staircase
Tutte / acyclic orientations), the LRC-permutohedron thread (s521/s525/s526). Four anchors verified in
`04-computation/orlik_solomon_across_the_repo_boxeph_S209.py`.*

## The abstract signature (what to look for)

Orlik–Solomon theory applies wherever you find this structure, independent of the objects' names:

1. a family of **conditions** (hyperplanes / equalities / resonances) with
2. an **intersection poset of flats** (or *layers*) that is a geometric (or arithmetic) lattice with a
   **Möbius function** `μ`, giving
3. a **count/invariant by inclusion–exclusion** over the poset — the **characteristic polynomial**
   `χ(t)=Σ_X μ(X) t^{dim X}` — computable by the **finite-field method** (`χ(q)=` #points off all
   conditions, `q≫0`), and
4. a **complement** whose topology is the **OS algebra** (`Poincaré = (−t)ⁿχ(−1/t)`, NBC basis), and
5. a **localization-at-a-flat** product factorization of everything near a degenerate flat.

Whenever a repo object is really "count / measure the configurations avoiding a lattice of coincidences,"
this machinery is available. It recurs in at least **four** guises — all verified.

## The repo map (four arrangement types)

| arrangement | repo object | what OS gives (verified) |
|---|---|---|
| **braid** `A_{n−1}` `{xᵢ=xⱼ}` | tournaments / NC2 Vandermonde (S208) | `χ=` falling factorial; regions `=n!=` transitive tournaments; **OS Poincaré `∏(1+kt)`**, Betti = Stirling-1st |
| **Shi / deformed braid** `{xᵢ−xⱼ∈{0,1}}` | LRC resonances `vt≈` integer | `χ=q(q−n)^{n−1}`; regions `=(n+1)^{n−1}` (parking functions), by **finite-field count** |
| **toric** (De Concini–Procesi) | LRC **relation lattice** `{k·v=0}` (THM-1820) | `\|G_δ\|=` complement volume `=` **arithmetic-Möbius sum over layers** `=` LRCMod mod-`q` |
| **generic** | figurate cutting (cake/bagel, S207) | regions `=ΣC(n,k)` (generic `χ`), deficit `=μ` of degenerate flats |

The point that unlocks the LRC side: **THM-1820's relation lattice `{k∈ℤⁿ:k·v=0}` is the poset of layers
of a *toric* arrangement** on the character torus, and De Concini–Procesi theory is the *toric analog of
Orlik–Solomon*. The repo already computes with it without the name — "LRCMod ladders = the relation lattice
tested mod `q`" (THM-1820) **is exactly the finite-field method on the torus**.

## Anchor 1 — tournament cohomology (a graded lens beyond char_A)

The complement of the braid arrangement (ordered configuration space) has OS Poincaré polynomial
`π(t)=∏_{k=1}^{n−1}(1+kt)`; its Betti numbers are the **unsigned Stirling numbers of the first kind**
`b_i=c(n,n−i)` (permutations with `n−i` cycles), top Betti `(n−1)!` (verified `n≤7`:
`n=4→[1,6,11,6]`, etc.). This is a **graded** invariant of ordering/tournament space, strictly finer than
the single datum char_A. Per-*tournament* refinement is available too: a tournament's comparability data
picks out a **graphic sub-arrangement** `A_G` with `χ_{A_G}=` chromatic polynomial and regions `=` acyclic
orientations — the deletion–restriction (Tutte) world THM-805 already lives in. So the tournament continuum
carries a cohomology, and char_A / the cycle spectrum are shadows of it.

## Anchor 2 — the finite-field method is the native LRC engine

Because LRC resonances are `{xᵢ−xⱼ = integer}` (a deformed braid arrangement) and the relation lattice is a
toric arrangement, the **characteristic (quasi-)polynomial is a Diophantine point-count** — precisely LRC's
number-theoretic flavor. Verified: braid `χ(q)=q(q−1)⋯(q−n+1)` and Shi `χ(q)=q(q−n)^{n−1}` reproduced by
raw `F_q` point-counts; Shi regions `(n+1)^{n−1}` are parking functions. This says the mod-`q` LRCMod
ladders are not an ad-hoc trick — they are computing an arrangement's characteristic quasi-polynomial, and
**Athanasiadis's finite-field method / Ehrhart quasi-polynomial theory is the systematic engine** for them.

## Anchor 3 — |G_δ| is a toric-complement volume (arithmetic Möbius)

Verified (matching THM-1820's bridge to `<3×10⁻³`): the good-set measure
`\|G_δ\| = ∫₀¹∏_j 1[‖v_jt‖≥δ]dt = Σ_{k·v=0} ∏_j ĝ(k_j)` — the left side is the **volume of the toric
arrangement's complement** (the loneliness set), the right side is the **arithmetic-Möbius sum over the
layers** `{k·v=0}`. Loneliness = the complement is nonempty = a chamber survives; covering = the complement
has measure/volume below threshold. This is the LRC restated as toric-arrangement region theory.

## Anchor 4 — tight = relation-richest = maximal Betti/Möbius mass

Verified: over primitive `n=3` speed triples, the **AP `(1,2,3)` uniquely maximizes relation richness**
`N_R` (matching THM-1820 B2 up to the ±sign convention — `N_R=4` counting `±k` once, `=8` with signs). In
arrangement language, `N_R` is the number of small **layers** — the toric arrangement's **Betti/Möbius mass**
— so *the LRC tight extremal is the configuration whose toric arrangement is richest in low-height layers*.
This is the same phenomenon as the braid side: the deepest, most-degenerate flat (all coordinates equal =
the AP/transitive vertex, the reify-ladder cold point) is where the Möbius mass concentrates.

## The leverage (what this buys, and what is still open)

- **(a) Tournament cohomology.** Use `π(t)=∏(1+kt)` and graphic-sub-arrangement deletion–restriction as a
  graded invariant of the continuum, finer than char_A; the cycle spectrum should be readable off the OS
  algebra.
- **(b) LRC characteristic quasi-polynomials.** Adopt the finite-field / Ehrhart engine explicitly for the
  LRCMod ladders — the toric arrangement's `χ(q)` *is* what those ladders compute; Wall A (AP-rigidity)
  becomes "the AP is the unique residue pattern maximizing the layer count," a point-count statement.
- **(c) The transferable trick — layer-localization.** The S208 win was that braid **flat-localization**
  factors the confluent Vandermonde into single-block L-P pieces (⇒ HYP-8775). The toric analog is
  **layer-localization** of the De Concini–Procesi arrangement: near a resonance layer, `\|G_δ\|` should
  factor into a lower-rank toric complement × a transverse braid factor. **This is proposed, not yet
  verified** — I confirmed the `\|G_δ\|=`complement-volume`=`arithmetic-Möbius identity, not the
  near-resonance factorization. If it holds, it is a concrete Wall-A tool: reduce the AP-extremality to a
  product of smaller resonance-block problems, exactly as the braid side reduced real-rootedness.

Honest scope: anchors 1–4 are verified identities/point-counts; the toric layer-localization (the actual
new lever for LRC) is a well-motivated conjecture by analogy with S208, and the general De Concini–Procesi
machinery over these specific lattices needs the arithmetic-Möbius computation carried out. The value here
is the **unification** — braid, Shi, toric, and generic arrangements are one OS pattern the repo keeps
re-encountering — and the identification of **toric arrangements (DCP) as the correct, under-named home for
the LRC relation lattice**, with the finite-field method as its engine.

Links: HYP-8830, HYP-8825, HYP-8775, THM-1820, THM-805,
[[the-missing-region-law-is-a-braid-arrangement-and-the-vandermonde-is-its-defining-polynomial-boxeph-S208]].
