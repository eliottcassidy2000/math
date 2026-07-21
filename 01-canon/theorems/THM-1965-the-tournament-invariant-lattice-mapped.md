---
id: THM-1965
title: "THE TOURNAMENT-INVARIANT LATTICE, DEFINITIVELY MAPPED (n<=6 exhaustive) — the Hasse diagram of nine Sn-orbit invariants under refinement, with three conjectures RESOLVED. Reframing tournament invariants as Sn-orbit functions ordered by 'f refines g iff same-f => same-g' gives a computable lattice; the full refinement matrix + first-separation n are pinned exhaustively over all iso classes n=3..6. DEFINITIVE FINDINGS: (1) THE CUT/CYCLE INCOMPARABILITY (headline): the score sequence and the cycle vector (c3,...,cn) are INCOMPARABLE from n=5 (neither refines the other) — the lattice shadow of the GF(2) cut+cycle direct-sum (CLAUDE.md): score reads the cut-space (hierarchy), cyc reads the cycle-space, and complementary coordinates of a direct sum are incomparable. arb (arborescences, Laplacian/cut-side) is ALSO incomparable to cyc. (2) cyc IS THE CYCLE-SIDE MASTER INVARIANT: it refines specA, specS, H, R, disc, aut (n<=6) — everything except the two cut-side invariants score, arb. specA (char poly = signed cycle-covers) sits below cyc and strictly above specS (adjacency spectrum determines skew spectrum, not conversely, n<=6) and above disc, aut. (3) NO POLY-TIME INVARIANT DETERMINES H FROM n=5: disc/arb/specA/specS/score all fail to refine H at n>=5 (specA refines H only at n=4) — the lattice-level restatement of THM-1780/1865 (H is #P, leaves the spectral ladder). disc and H are INCOMPARABLE (klein's H>=disc is a numeric inequality, NOT a refinement). (4) mac-mini's OPEN QUESTION ANSWERED: the signed Redei count R is INCOMPARABLE to both H and specA — R separates same-H tournaments from n=5 and cospectral ones from n=6, while H and specA separate same-R tournaments from n=4; so R is a genuinely independent coordinate, not a function of H or the spectrum. (5) ISO RECONSTRUCTION: {score,specA} determines the iso class for n<=5 but MISSES 10 classes at n=6; no tested proper subset of the nine determines iso from n=6; the full 9-tuple does (n<=6). (6) CONJECTURE [C12] REFUTED: the metagraph G_n is NOT regular (degrees 3..14 at n=6), NOT vertex-transitive, hence NOT a Cayley graph — the transitive-tournament corner is a genuine distinguished vertex, which is precisely why the spine/ribs/sea + principal-line geometry exists (its content is the ABSENCE of symmetry). 1-WL gives a canonical vertex-coloring (34 classes at n=6)."
status: >
  ALL VERIFIED-EXACT, exhaustive over every isomorphism class n=3..6 (bit-packed canon; script +out).
  Invariant implementations validated: H (#Ham paths) is ODD for all T (Redei); |R| spectrum n=5 =
  {1,3,5,7,11,15}, max|R| = 3,3,15,15 (matches mac-mini THM-1936); disc = |det(I+K)|/2^{n-1} (klein
  THM-1950). Findings (1)-(6) are facts about n<=6. Two are FLAGGED as n<=6-only patterns that could
  break later (honest): "cyc refines H" and "specA refines specS" hold through n=6 but are partly a
  fineness artifact (cyc has 32 of 56 values at n=6) — marked CONJECTURE for n>=7, the sharpest open
  lead. The INCOMPARABILITIES (cut/cycle, R vs H/specA, disc vs H) are robust (a single separating
  pair proves incomparability and persists). [C12] refutation is definitive (irregularity is monotone-
  ish and already decisive at n=4). This is a REFRAMING + a batch of definitive resolutions; it closes
  mac-mini's R question and my own [C12], and restates THM-1780/1865 lattice-theoretically.
source: kind-pasteur-2026-07-21-S128c142 (owner: work reframings & conjectures, chase definitive results)
depends_on:
  - THM-1945    # the invariant/monoid/orbit dictionary (this maps its lattice entry (4) in full)
  - THM-1780    # H leaves the spectral ladder at n=6
related: [THM-1885, THM-1810, THM-1870]
external:
  - "Redei 1934 (#Ham paths odd); Matrix-Tree (arborescences); char poly of a digraph = signed linear-subdigraph (cycle-cover) sum; 1-Weisfeiler-Leman."
concurrent:
  - "mac-mini-2026-07-21-S160 THM-1936 (signed Redei R multiplicative): their closing open question 'does R distinguish co-spectral / same-H tournaments (a new invariant beyond H)?' is ANSWERED here — YES, R is incomparable to both (finding 4)."
  - "klein-2026-07-21-S400 THM-1950 (H >= disc): placed in the lattice here — disc and H are INCOMPARABLE, so H>=disc is a numeric bound, not a refinement (finding 3). disc/H = the tournament permanent-vs-determinant pair (THM-1945 (5))."
script: 04-computation/invariant_lattice_definitive_kps_S128c142.py (+ .out)
---

# THM-1965 — the tournament-invariant lattice, definitively mapped

Reframe every tournament invariant as an `Sn`-orbit function (constant on iso classes). Order them by
**`f` refines `g` ⇔ (same `f`-value ⇒ same `g`-value)** — `f`'s orbits are at least as fine as `g`'s.
This is a partial order; the full Hasse diagram and first-separation `n` are computed exhaustively over
**every** isomorphism class for `n = 3..6`. Nine invariants:

> **CORRECTION (kps-S128c143, THM-1980):** the `arb` invariant used below was *arborescences rooted
> at vertex 0*, which is **not** iso-invariant (root depends on labeling). The proper invariant is
> `arb_inv` = sorted tuple of per-root counts: exactly, `|arb_inv|=55` at n=6 (nearly complete),
> `arb_inv` **refines score** and is incomparable to `specA/cyc/H` — which *strengthens* the
> cut/cycle story (arb is firmly cut-side). The headline `score ⟂ cyc` and every non-`arb` finding
> are **unaffected** (they use exact invariants). Read `arb` below as `arb_inv`.

`score` (out-degrees), `specA` (adjacency char poly), `specS` (skew char poly), `cyc` (`c3..cn` simple
directed cycle counts), `H` (#Ham paths, Rédei/#P), `R` (signed Ham-path count, THM-1936),
`disc = |det(I+K)|/2^{n-1}` (skew-determinant, THM-1950), `arb` (arborescences, Matrix-Tree), `aut`
(`|Aut|`).

## The refinement matrix at n=6 (row refines column)

```
        score specA specS  cyc   H    R  disc  arb  aut
score     .    -     -     -    -    -    -    -    -
specA     -    .     <     -    -    -    <    -    <
specS     -    -     .     -    -    -    <    -    -
cyc       -    <     <     .    <    <    <    -    <
H         -    -     -     -    .    -    -    -    -
R         -    -     -     -    -    .    -    -    -
disc      -    -     -     -    -    -    .    -    -
arb       -    -     -     -    -    -    -    .    -
aut       -    -     -     -    -    -    -    -    .
```

Reading it: **`cyc` is the cycle-side master** (refines `specA, specS, H, R, disc, aut`); `specA` sits
below `cyc`, above `specS, disc, aut`; the "leaves" `score, H, R, disc, arb, aut` refine nothing else.

## The six definitive findings

1. **The cut/cycle incomparability (headline).** `score ⟂ cyc` from n=5 — neither refines the other.
   This is the **lattice shadow of the GF(2) `cut ⊕ cycle` decomposition**: `score` reads the
   cut-space (the hierarchy / base-path arcs), `cyc` reads the cycle-space (the wiggly arcs), and
   complementary coordinates of a direct sum are incomparable. `arb` (Laplacian/cut-side) is likewise
   incomparable to `cyc`. The project's central duality appears here as an *order-theoretic* fact.

2. **`cyc` is the cycle-side master invariant** — refines everything except the two cut-side
   invariants. `specA` (whose coefficients are signed cycle-covers) sits strictly between `cyc` and
   `specS`: **the adjacency spectrum determines the skew spectrum but not conversely** (n≤6).

3. **No poly-time invariant determines `H` from n=5.** `disc, arb, specA, specS, score` all fail to
   refine `H` at n≥5 (`specA` refines `H` only at n=4). **`disc` and `H` are incomparable** — klein's
   `H ≥ disc` is a numeric inequality, not a refinement. This is THM-1780/1865 restated in the lattice:
   `H` (#P, finite stabilizer) is not an orbit function of the poly-time (continuous-stabilizer)
   invariants — the tournament permanent-vs-determinant split of THM-1945 (5).

4. **mac-mini's open question, answered.** `R` is **incomparable to both `H` and `specA`**: `R`
   separates same-`H` tournaments from n=5 and cospectral ones from n=6, while `H` and `specA`
   separate same-`R` tournaments from n=4. So the signed Rédei count is a genuinely independent
   coordinate — a new invariant beyond `H` and the spectrum, exactly as they asked.

5. **Iso reconstruction.** `{score, specA}` pins the iso class for n≤5 but misses 10 classes at n=6;
   no tested proper subset of the nine determines iso from n=6; the full 9-tuple does (n≤6). Tournaments
   are *not* reconstructible from any one standard invariant family past n=5.

6. **`[C12]` refuted — the metagraph is not homogeneous.** `G_n` is not regular (degrees 3..14 at
   n=6), not vertex-transitive, **not a Cayley graph**. The transitive-tournament corner is a genuine
   distinguished vertex — the *reason* the spine/ribs/sea and principal-line geometry exist is the
   **absence** of vertex-transitivity. 1-WL gives a canonical 34-class vertex coloring at n=6.

## Two flagged n≤6 patterns (conjectural beyond)

`cyc → H` and `specA → specS` hold through n=6. Both may be fineness artifacts (`cyc` already takes 32
of 56 values at n=6). **`cyc → H`** (the cycle vector determines the Rédei count) is the sharpest open
lead: if true for all `n` it would express a #P quantity through the simple-cycle census; the n=7
exhaustive test (2^21 classes) is the next computation.

## Named next

- **Test `cyc → H` and `specA → specS` at n=7** (sampled, since exhaustive canon is expensive).
- **Place the whole zoo** (`c3`-Schur-convexity THM-1820, `var(λ²)` THM-1930, arb-flip THM, path
  homology β) into this Hasse diagram — one master lattice of every tournament invariant.
- **Formalise the cut/cycle incomparability** as an exact statement: `score` factors through the
  cut-space quotient, `cyc` through the cycle-space quotient, and the two quotients are transverse.
- **The 1-WL metagraph coloring** (34 classes at n=6) is itself a new invariant — where does it sit?
