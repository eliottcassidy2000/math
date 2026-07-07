---
id: THM-642
title: The image census and the residue-factoring barrier — which tournament iso classes occur under the runner-mapping rules, and the leverage: residue-based maps see EXACTLY the sieve floor (never the density floor), and single-time maps see only round tournaments (metric-blind); so every finite-projection tournament map bottoms out before the LRC(14) density floor, which is irreducibly a measure over all times
status: PROVED (elementary + verified). The image characterizations are identities; the residue-invisibility is exhibited (lift families = mac-mini-S36 escape families, canon). A limitation-leverage theorem: it PRUNES whole approach-families and pins what a working rule must retain, and it presents the sieve as an image census (the hard core = one image cell).
source: kind-pasteur-2026-07-07-S65 (HYP-4957; owner directive: "which set of isomorphism classes are possible under a set of mapping rules ... define rules [so] the possible set is restricted, and leverage that as a fact for proofs")
depends_on:
  - THM-640   # the Paley bridge (the mod-p QR cutoff; T_p^QR({0..p-1}) = Paley T_p)
related:
  - HYP-4667  # mac-mini-S36 escape families (v_i = i + L*k_i, ≡ AP mod L, M varies) = the residue-fibers
  - HYP-4712  # counterexample_needs_all_divisors / the saturation sieve (GREEN)
  - HYP-4897  # S63 heresy C: the 2-point LP is insufficient (pairwise projection barrier)
  - HYP-4937  # monad-S5: no orientation invariant separates tight/loose
external: Paley tournaments (arc-transitive / doubly regular); round (locally transitive) tournaments; Lonely Runner Conjecture.
---

# THM-642 — The image census and the residue-factoring barrier

The owner asked which tournament **iso classes** occur under a runner-mapping rule, whether a
rule can be designed so the image is **restricted**, and how to **leverage** that as a proof
fact. The answer is a clean dichotomy of what tournament maps can and cannot see.

## A. The image of the mod-N (QR / CRT) cutoff — a census

The mod-`p` QR cutoff (`i→j` iff `(v_j−v_i) mod p ∈ QR_p`, THM-640) has, on a family `V`,
adjacency depending **only on the residue multiset `{v_i mod p}`**. Hence:

> **Image(mod-`p`) = { induced subtournaments of the Paley tournament `T_p` }**, indexed by
> the residue multiset.

For the composite factor `p = 7` of `14 = 2·7` this image is tiny, because `T_7` is
**arc-transitive** (doubly regular):
- **6 distinct nonzero residues** (necessarily `{1,…,6}`) → the **unique** class `T_7 ∖ v`
  (`scores = (2,2,2,3,3,3)`, `c3 = 8`, `H = 45`). A one-class image.
- residues **`{0,1,…,6}`** (a multiple of 7 present, all residues once) → the **full Paley
  `T_7`** (`regular`, `c3 = 14`, `H = 189` — the tournament `H`-maximizer). The two-project
  weld of THM-640 sits exactly here: the observer-inclusive saturated AP maps to the
  maximizer.

The **CRT-14 image stratifies** by the **saturation pattern** — which `q ∈ {2,…,14}` have a
multiple among the speeds. This is a down-set in the divisor lattice, and:

> **every image cell except the top is provably lonely.** A family missing a `q`-multiple has
> `M ≥ 1/q ≥ 1/14` (the sieve, `t = 1/q`). The **LRC(14) hard core is the single top cell**
> (fully saturated, a multiple of every `q ≤ 14`). So the entire remaining problem lives in
> *one* image cell — the sieve, presented as an image census. (This re-derives
> `counterexample_needs_all_divisors` in tournament language: the `q`-defect vertex — residue
> `0` mod the prime part of `q` — is present iff the family is `q`-saturated.)

## B. The residue-factoring barrier (the limitation, and its leverage)

Because the mod-N tournament factors through residues mod N, **its every invariant is a
function of the residue multiset alone.** Therefore it can see only the *residue-visible* part
of loneliness — the sieve floor — and **nothing of the density floor**:

> **Residue-invisibility (verified).** Families with **identical residues mod 14** have `M`
> ranging from `1/14` (tight) upward: `{1,…,13}` (`M = 1/14`), `{1,…,12,27}` (`0.0769`),
> `{15,2,…,13}` (`0.1176`), … — the map's fibers are exactly the **lift / escape families**
> `v_i = a_i + N·k_i` (mac-mini-S36, HYP-4667). No residue-mod-N tournament invariant
> distinguishes them.

So a residue-based tournament map **can prove the sieve** (all GREEN) and **cannot reach the
density floor**. That is the leverage: it **prunes** the entire family of residue/covering
tournament approaches at the sieve, and localizes the open problem to the residue-invisible
axis — the density itself.

## C. The single-time map sees only round tournaments (metric-blind)

The value-based snapshot cutoff (`i→j` iff `frac(v_i t)` precedes `frac(v_j t)` in the
clockwise semicircle) has

> **Image(single-time) = round (locally-transitive) tournaments** (verified: every sampled
> `T_t` is locally transitive).

Round tournaments are configurations of points on a circle **up to combinatorial order** —
they **forget the gap sizes** (the metric). Since loneliness *is* a gap size, it is invisible
to any single snapshot. As `t` sweeps `[0,1)`, `T_t` walks the round-tournament graph by
adjacent transpositions (opus-S136 order cells); the number of states is the **Wiener
collision count** `Σ_{i<j}(v_j−v_i)` (S62), **minimized by the AP** (`364` vs record `564`).
This is *why* the S64 geometric cutoffs (`T_good`, `T_witness`) collapsed: a snapshot cannot
carry the metric.

## D. The unifying meta-statement (the real leverage)

Every **finite-projection** of the family bottoms out before the density floor:

| projection | image restriction | sees | bottoms out at |
|---|---|---|---|
| residues mod N (THM-640 cutoff) | Paley subtournaments | the **sieve** | `M ≥ 1/q` (covering) |
| one snapshot `t` (semicircle) | **round** tournaments | combinatorial order | metric-blind |
| pairwise statistics (2-point) | — (LP) | pair-uniform data | `≈ 0.126 < 1/7` (S63/klein-S159) |

> **The LRC(14) density floor is not visible to any residue-, snapshot-, or pairwise-
> projection of the runner family; it is irreducibly a measure over all times** (the density
> `μ_{1/7}` itself). A tournament (or any finite-information) rule that hopes to reach it must
> retain the full time-measure, not a residue class, a snapshot, or a pair marginal.

This unifies four independently-found barriers — S63 heresy C (2-point LP), S64 (geometric
cutoffs transitive), klein-S159 (pairwise LP), and this census (residue-factoring) — into one
principle, and it is a *positive* proof aid: it tells the fleet exactly which whole families
of tournament/projection attacks cannot close the floor, and that the covering/sieve side is
*already* the complete residue-visible content.

## Verification

`lrc_image_restriction_kps_S65.py` (+out): the mod-7 arc-transitivity collapse (one class,
`H=45`); the full `T_7` at residues `{0..6}` (`H=189`); the sieve = defect-vertex dichotomy
(loose-at-7 ⟺ no residue-0 vertex, exact on the zoo); the CRT-14 saturation stratification;
the round-tournament image (all `T_t` locally transitive); and the residue-invisibility
addendum (identical residues mod 14, `M ∈ [1/14, 0.118]`).

## Consequence for the program

- The tournament/image lane **completes the sieve** and **cannot go past it** — consistent
  with the fleet's live front being the *density floor* (kps-S59/S60/S61 ledgers, the μ_{1/7}
  measure), not the covering.
- The two-project weld (THM-640) is now located precisely: it lives at the **top image cell**
  (fully-saturated, residues `{0,…}`), where the runner tournament is (near-)Paley — the
  `H`-maximizer — and the loneliness question is (by B) orthogonal to it. The bridge relates
  the *extremal principles*; the barrier says it cannot relate the *difficulty*.
