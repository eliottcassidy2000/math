# Deep archaeology: the OEIS-singleton and uncanonized-HYP layers

*kind-pasteur-2026-07-21-S128c137. Owner: find even more niche forgotten ideas via archaeology;
collaboratively map past ideas. Companion to the theorem-citation atlas (kps-S128c136) and to
klein-S399's tournament-invariant zoo atlas — this mines the layers BELOW theorems, where the
niche ideas hide. Tool: `04-computation/deep_archaeology_kps_S128c137.py`, full dump in
`05-knowledge/results/deep_archaeology_kps_S128c137.out`.*

The theorem-citation graph (c136) finds forgotten *theorems*. But the niche stuff lives in
reflection bodies and session logs, tagged by things a title scan misses: **OEIS numbers, rare
named objects, and confirmed-but-uncanonized hypotheses**. Three passes.

## (A) OEIS singletons — 24 concrete forgotten numerical connections

114 distinct A-numbers are referenced across the corpus; **24 appear in exactly one file** — each
a sequence someone computed, matched, and never returned to. The ones worth reviving:

- **A000182 (tangent numbers)** — "`W(0)` of the transitive tournament = A000182" (era-1 session
  log only). **Revived today as THM-1875:** the transitive skew form is `((x+1)ⁿ+(x−1)ⁿ)/2` with
  cotangent spectrum — the tangent thread is the odd/even shadow of the GIT nullcone vertex.
- **A113077** (`1,3,15,123,1656,36987`) and **A368322** (`1,5,37,389,5413,94085`, EGF
  `exp(2x)/…`) — two *tournament counting sequences* referenced once. Unidentified which tournament
  statistic; a direct lead for the H/moment work.
- **A003049** (`2,3,7,16,54,243,2038`) — connected-Eulerian/even-graph count, the `Eₙ` sibling of
  `A002854` (`V(Eₙ)`); feeds the "`Eₙ` GIT nullcone" gap.
- **A001121 (skew-Hadamard census)** — the doubly-regular-tournament / skew-Hadamard existence
  thread, one mention. **A006125** (`2^{C(n,2)}`) — labelled switching classes.
- **A001175 (Pisano periods)** — the Fibonacci-mod-10-period-60 / "1001 = three sixties" thread
  (a standing owner motif), cited once against THM-481 (Golay).
- **A357242/A357248/A357257/A357266** — a *J. Integer Sequences 27 (2024)* paper's sequences,
  referenced once as "= exactly 2/3/4/5-something"; an external result never chased.
- **A000797** — "241 exceptions, largest 343867, OPEN" — an open OEIS problem the repo brushed.
- **A002324** — `#units(L_t) = 12 + 6·A002324(t)`, the pentagonal/η²⁴-code units thread.

## (B) rare named objects — thin threads

Appearing in `≤ 1` file: **Apéry**, **Macdonald**, **Cayley–Menger**. The last is notable —
Cayley–Menger determinants are the natural tool for the LRC/unit-distance geometry (simplex volumes
from pairwise distances) and were mentioned once and dropped. `≤ 3–4` files: umbral,
quasisymmetric, hypergeometric — the Sheffer/umbral toolkit (which my THM-1620 Pochhammer bridge
leans on) is barely tagged.

## (C) confirmed-but-uncanonized HYPs — 26 proven ideas that never became theorems

Hypotheses marked CONFIRMED/PROVED whose id appears in no theorem file. A proof happened and never
got a THM. The substantive ones:

- **HYP-2058** — "LRC@14 proof-lite: structured counterexamples impossible."
- **HYP-2078** — "anti-automorphism Burnside extends to self-converse oriented graphs."
- **HYP-2189** — "the Cauldron Game as Schur coloring with removal."
- **HYP-2198/2199** — "single-core signature gap: complete for all lengths; density ½, no simple
  closed form."
- **HYP-2210** — "perspective-flip carriers give practical compression for LRC / unit-distance."
- **HYP-2212** — "the π/e quadratic carrier has a discriminant sheet."
- **HYP-2331** — "Erdős 625 (σ, the order-2 involution) on the arc's 2-adic seam."

Each is a proved result sitting only in the HYP log. A canonization pass (promote confirmed HYPs to
LEM/THM with a back-pointer) would recover 26 results the theorem layer cannot see.

## The collaborative map

klein-S399 catalogs the *invariants*; c136 catalogs the *theorems*; this catalogs the *sub-theorem
layer* (OEIS / named-objects / uncanonized proofs). Together they are the repo's index. The three
tools (`corpus_archaeology`, `deep_archaeology`, klein's generators) should be **rerun each era** —
that, not re-reading, is how a 1392-theorem multi-agent corpus keeps from forgetting itself.

Standing revival queue (new, beyond c136's): identify A113077/A368322 (which tournament stat?);
canonize the 26 confirmed HYPs; revive Cayley–Menger for LRC geometry; THM-293 (succession GF `W`)
+ the tangent identity (THM-1875 named next); secant numbers A000364 on the even side.

## Cross-links
`deep_archaeology_kps_S128c137.py` · THM-1875 (tangent revival) · THM-1870 (c136) · klein-S399
(invariant atlas) · THM-293 (forgotten W) · THM-1810/1555 (nullcone / half-dictionary).
