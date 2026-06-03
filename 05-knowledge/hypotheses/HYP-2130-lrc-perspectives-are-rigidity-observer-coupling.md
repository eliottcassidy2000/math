---
id: HYP-2130
status: REINTERPRETATION + EVIDENCE — the human's "perspective" curiosity is about automorphism
  RIGIDITY (vertex-orbits), not chirality; it is the tournament face of LRC's single-corrector
  (observer-coupled) vs multi-sieve (observer-blind) dichotomy. Corrects the HYP-1824/1825
  "chirality residue" reading. Conceptual; the LRC bridge is a conjecture.
source: claudebox-2026-06-03-S590
related:
  - HYP-1824
  - HYP-1825
  - HYP-1981
  - HYP-2063
  - HYP-2075
  - HYP-2120
---

# HYP-2130: the perspective curiosity is rigidity, and it is the LRC single-vs-multi corrector line

The human has stated, in rough wording several times (T075 "perspective conjecture" from the
Tournament Tiling Explorer; T174 / INV-083 "P(n)=A093934"), a curiosity:

> the number of tournament **perspectives** (vertex-orbits under the automorphism group; pointed
> tournaments) on `n` vertices coincides with the number of tournament **structures** on `n+1`.
> n=3: 2 structures, **4** perspectives (3 from the transitive, 1 from the cyclic). n=4: 4
> structures, **12** perspectives (4+4+2+2). 12 = structures on 5.

and noted it as "probably the key to the LRC, where and why it works and doesn't." This reframes it.

## What it is (exact, Burnside — `lrc_perspective_rigidity_observer_s590.py`)

`perspectives(n) = Σ_T (#vertex-orbits of T) = ` A093934 `= 2,4,12,48,296,3040,…`;
`structures(n) = ` A000568 `= 1,2,4,12,56,456,…`.

- The coincidence `perspectives(n) = structures(n+1)` holds **only for n ≤ 4** (gap 0,0,0 then
  **8, 160, 3840, …**). Both sides briefly equal `2(n-1)!`: `perspectives(n) = (n−1)·T(n) =
  2(n−1)!` for n ≤ 5 (each class loses *exactly one* perspective to symmetry); `T(n+1) = 2(n−1)!`
  only for n ≤ 4. The agreement is the overlap of two small-n accidents, not an identity.
- `perspectives(n)/structures(n) → n`: **almost every large tournament is rigid** (trivial Aut, `n`
  distinct perspectives). The deficit `n·T(n) − perspectives(n)` (the "lost" perspectives) is the
  **symmetry mass**, and `deficit/T(n) → 0` — symmetric tournaments are a vanishing fraction that
  nonetheless set the small-n gap.

## The misinterpretation it corrects (HYP-1824 / HYP-1825)

Those notes read the gap `56 − 48 = 8` as a **chirality residue** (`56 = 12 self-converse + 44
chiral`; "8 = projection/stencil layer", matched to the n=14 LRC "8 alpha stencils"). But chirality
is the `T ↔ T^op` involution, a **different symmetry** from the automorphism rigidity that defines
perspectives. The numbers refute the chirality reading: self-converse counts are `2,2,8,12` and
chiral `0,2,4,44`, while the perspective gap is `0,0,8,160` — the gap equals **neither**. The `8`
was a numerical coincidence (`T(6)−Σorbits(T(5))`), matched to a mechanism it does not share.
**Perspectives measure rigidity (vertex-orbits), not chirality.**

## The LRC reading (the actual point): observer-coupled vs observer-blind

A tournament's rigidity is exactly how distinguishable its vertices are *from a marked vertex's
viewpoint* — the **observer** of the LRC dictionary (HYP-1981, THM-381: a runner is lonely iff its
marked observer vertex is a source). Two regimes:

- **Observer-coupled** (rigid, many perspectives, trivial Aut): every runner is distinguishable, so
  the loneliness obstruction **localizes to a single privileged perspective** — the apex /
  co-observer (HYP-2063: the apex `q=n/2` is the unique zero-divisor / the single uncorrectable
  runner). A **single-modulus / single-corrector** argument (Sungkawichai–Trakulthongchai Prop 4.1)
  suffices. This is the LRC-*easy* regime.
- **Observer-blind** (symmetric, few perspectives, nontrivial Aut): the symmetry means **no single
  runner is the unique obstruction** — it is shared across an automorphism orbit. The single
  corrector **fails at the apex**, and one must use the **multi-sieve** (pair-sum moduli — *pairs*
  of perspectives). The repo found this empirically without naming it: **HYP-2075 "multi-sieving
  has no apex"** — the single-perspective (apex) obstruction *dissolves* once multiple perspectives
  (pair-sum sieves) are used. This is the LRC-*hard* regime.

> **Claim.** The break in the perspective↔structure coincidence is the tournament shadow of the
> **single-corrector → multi-sieve transition**: where the loneliness obstruction stops being one
> perspective (the apex) and becomes an orbit of coupled perspectives. "Where and why LRC works and
> doesn't" = whether the observer is coupled to a single critical vertex (works, Prop 4.1) or blind
> across a symmetric orbit (fails, needs the pair-sum sieve). The apex obstruction (HYP-2063) and
> its multi-sieve dissolution (HYP-2075) are the two sides of this line.

This also re-reads the **observer-blind / observer-coupled** dichotomy I have been developing as
circuit-free (Lemma A) vs 3-term (Lemma B) (HYP-2120): a 3-term fusion `χ_a χ_b = χ_c` is an
*additive coupling among runners* that breaks the observer's single-perspective view — the
structured/hard side; circuit-free is the rigid/generic/easy side.

## Open / next

- **Make the bridge precise:** is the symmetry-deficit (the observer-blind mass) of the marked LRC
  tournament a quantitative predictor of where the single corrector fails (n=14 apex)? Compute the
  automorphism/orbit content of the *marked* observer tournaments for the tight families and test
  against the Prop-4.1 failure locus.
- Per-`n` symmetry mass vs the LRC structural→computational threshold (n≤7 structural, n≥8
  computational): does the rigidification `persp/T → n` track the loss of structural proofs?
- The chirality (self-converse) decomposition is real but ORTHOGONAL — both invariants exist; only
  rigidity is the perspective count.

**Artifacts:** `04-computation/lrc_perspective_rigidity_observer_s590.py` (+`.out`),
`07-reflections/lrc-perspectives-are-rigidity-s590.md`. Catalog of prior occurrences: T075, T174
(`00-navigation/TANGENTS.md`), INV-083 (`INVESTIGATION-BACKLOG.md`), HYP-1824, HYP-1825 (the reading
corrected here), the observer/apex thread (HYP-1981, THM-381, HYP-2063), HYP-2075 (multi-sieve has
no apex), HYP-2120 (circuit-free/3-term = observer-blind/coupled).
