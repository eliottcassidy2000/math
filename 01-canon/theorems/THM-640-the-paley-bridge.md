---
id: THM-640
title: The Paley Bridge — the observer-inclusive LRC(p) family, under the QR-mod-p pair-statistic cutoff, IS the Paley tournament T_p; so the density-floor M-minimizer equals the H-maximizer, and step-reversal ↔ tournament-complementation makes the palindromic AP self-complementary
status: PROVED (elementary; the construction is an identity of circulant tournaments; verified exact p=3,7,11,19,23 and SC via x↦gx). Realizes THM-126 (Paley spectral-flatness ↔ density floor) as a literal tournament construction from the owner's binary-cutoff directive.
source: kind-pasteur-2026-07-07-S64 (HYP-4927; owner directive: "pair statistics among integer families can be considered tournaments by deciding a meaningful binary cutoff")
depends_on:
  - THM-126   # Paley spectral-flatness principle: Paley maximizes H (flat |λ|) ⟺ AP at t=1/p minimizes M (same roots of unity)
  - THM-024   # every SC tournament has an involutive anti-automorphism (Moon+Cauchy); Paley is SC
related:
  - HYP-4877  # the step-gauge (reversal = complement involution); palindromic = SC
  - HYP-4897  # heresy that a 2-point tournament invariant sees the floor (T_good/T_witness are transitive — negative)
  - HYP-4937  # monad-explorer perspective-frame runner tournaments (centering/exclusive-attention)
  - THM-639   # mac-mini runner tiling model (steps = path, balanced lattice = tiles)
external: Paley tournaments; Rédei's theorem (H = #Hamiltonian paths); Lonely Runner Conjecture (LRC(14) = first open case).
---

# THM-640 — The Paley Bridge

## Statement

Fix a prime `p ≡ 3 (mod 4)` and the quadratic-residue set `QR_p = {a² mod p}`. Turn the
**pair statistics** of an integer speed family into a tournament by the **QR-mod-p binary
cutoff**:

> `T_p^{QR}(V)`:  for distinct speeds `v_i, v_j`,  **arrow `i → j`  ⟺  `(v_j − v_i) mod p ∈ QR_p`.**

Because `p ≡ 3 (mod 4)`, `−1` is a non-residue, so `QR_p` and `−QR_p` partition the nonzero
residues: the cutoff is a well-defined tournament off the `p | (v_j−v_i)` diagonal.

**(A) The bridge.** For the **observer-inclusive LRC(p) family** `V = {0, 1, 2, …, p−1}`
(the stationary observer `0` plus the `p−1` unit-speed movers — the LRC(p) density-floor
minimizer base),
> **`T_p^{QR}({0,1,…,p−1})` is exactly the Paley tournament `T_p`.**

Consequently `T_p^{QR}(\{1,…,p−1\})` (movers only) is `T_p` with the origin vertex deleted.

**(B) Two extremal principles, one object.** The same family `V = {0,…,p−1}` is:
- the **LRC density-floor / loneliness minimizer**: at `t = 1/p` its orbit is the `p`-th
  roots of unity, the maximally spread configuration that **minimizes `M`** (THM-126,
  Paley spectral flatness); and
- via the QR cutoff, the tournament that **maximizes `H`** (Rédei count of Hamiltonian
  paths): `H(T_7) = 189`, `H(T_11) = 95095`, the canonical maximizers.

So the object that **minimizes loneliness** is, read as a pair-statistic tournament, the
object that **maximizes Hamiltonicity**. THM-126 asserted this correspondence spectrally;
THM-640 realizes it as an explicit construction.

**(C) Reversal ↔ complementation (the SC corollary).** For every family the QR cutoff
**intertwines step-reversal with tournament complementation**: with `V* = (max V + min V) − V`
and the order-reversing relabel `σ`,
> `T_p^{QR}(V*) = σ · (T_p^{QR}(V))^{op}`.
Hence a **palindromic** family (`V* = V` up to translation) maps to a **self-complementary**
tournament. The AP is palindromic, and indeed `T_p` is self-complementary (via `x ↦ gx`,
`g` a non-residue — verified). This fuses the LRC **palindromic-extremizer** phenomenon
(HYP-4877) with the tournament **SC-maximizer** theory (THM-024) into one statement.

## Proof

**(A)** `T_p` is by definition the circulant tournament on `ℤ/p` with connection set `QR_p`
(`a → b ⟺ b − a ∈ QR_p`). The family `{0,1,…,p−1}` has values equal to the residues
`0,1,…,p−1` themselves, and `v_j − v_i = j − i`, so `T_p^{QR}` on it is *literally* the same
adjacency rule as `T_p`. Identity, not isomorphism. ∎

**(C)** Reversal sends `v_i ↦ c − v_{σ(i)}` (`c = max+min`, `σ` the order-reversal), so
`v_{σ(j)} − v_{σ(i)}` under `V*` equals `−(v_j − v_i)`. Since `p ≡ 3 (mod 4)`,
`−QR_p = \overline{QR_p}` (the non-residues), so every arrow flips: complementation, up to
`σ`. The AP is fixed by `σ` (all steps equal), so `T_p^{QR}(AP)` is isomorphic to its own
complement — self-complementary. (Classic: Paley tournaments are SC.) ∎

**(B)** is THM-126 (loneliness/`M` side) combined with `H(T_p) = ` max (the Paley
maximizer, canon) applied to the identity of (A).

## Verification

`lrc_paley_bridge_kps_S64.py` (+out): exact identity `T_p^{QR}({0..p−1}) = T_p` for
`p = 3, 7, 11` (and `T_p^{QR}({1..p−1}) = T_p∖\{0\}` exact for `p = 3,7,11,19,23`);
`H(T_7)=189, c3=14, regular`; `H(T_11)=95095, c3=55`; SC via `x↦gx` confirmed for `p=7,11`;
reversal→complement intertwining confirmed for every tested family (mod-7 and CRT-14 rules).

## The composite-14 shadow (why LRC(14) is hard, in tournament language)

`14 = 2·7` is **not** prime, so there is no residue field and no clean Paley cutoff on all
differences. Under the CRT rule on `ℤ/14 ≅ ℤ/2 × ℤ/7`, the AP `{1,…,13}` maps to a tournament
that is **near-regular but not regular** (scores `5,5,5,6,…,7`, `c3 = 88`, `H = 3 285 381`):
**42 of the 78 pairs have non-unit differences mod 14** (divisible by 2 or 7), and those are
exactly the arcs the CRT cutoff cannot orient by a QR rule. This is klein-S151's
"`k+1` prime" requirement (`14 = 2·7` breaks the SOTA polynomial method) seen from inside the
bridge: **LRC(14) is hard precisely because the QR-cutoff tournament of its AP cannot be a
regular Paley tournament.** The composite obstruction is now a visible tournament defect.

## What it does NOT do (honest negatives, converging with monad-S5 / HYP-4767)

The *geometric* pair-statistic cutoffs are transitive, not informative:
- **Good-set leader** `T_good`: `i→j` iff `i` leads `j` more often on the load-bearing set
  `Good_{1/7}` — this is the **consensus phase-lead total order**, always the transitive
  tournament (`c3 = 0, H = 1`); the Good-conditioning does not break the pair symmetry
  enough at small `k` to create cycles.
- **Witness-semicircle** `T_witness`: at the max-gap witness time `t*` the phases cluster
  (that is *why* the gap is wide), so the semicircle order is transitive there too.

So no *single* orientation invariant read from one time or one average separates tight from
loose families (extends HYP-4767 / monad-S5). The Paley bridge is a different kind of result:
not a separator, but an **identity welding the two projects' opposite extremal principles**.

## Why it matters

The whole project's founding dictionary — "tournament = complete pairwise binary relations;
the Hamiltonian path is the gauge" — meets LRC's founding fact — "the pairwise data is
*forced featureless*" (exact uniformity, S59). The bridge shows these are not merely
analogous: the **arithmetic that makes the pairwise data featureless (roots of unity) is the
arithmetic of the QR field**, and reading that field as a cutoff *is* the Paley construction.
The observer-inclusive family `{0,…,p−1}` — `p` runners for `p` "players", the owner's
"tree on `n` events has `n−1` edges = the Hamiltonian path" — is the shared fixed point.
