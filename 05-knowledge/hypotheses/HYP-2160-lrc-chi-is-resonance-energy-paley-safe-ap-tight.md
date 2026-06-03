---
id: HYP-2160
status: VERIFICATION + ANSWER — among regular (vertex-transitive) tournaments, χ adds beyond
  vertex-transitivity: it is the resonance-energy signature. Flat-χ doubly-regular Paley = SAFE;
  structured-χ rotational AP = TIGHT. No tight config is regular off the AP/AP-lift orbit. Verified
  small-n; the "χ = resonance energy among regulars" claim is conjectural in general.
source: claudebox-2026-06-03-S608
related:
  - HYP-2155
  - HYP-2124
  - HYP-2150
  - HYP-2055
---

# HYP-2160: χ adds beyond vertex-transitivity — it is the resonance energy; Paley is safe, AP is tight

**Question (human):** among the maximally-cyclic (regular) tournaments, does χ add anything beyond
vertex-transitivity — is there a tight config that is regular but not the Paley/AP orbit, and does
its χ differ?

## The answer

**Yes, χ adds — and it is the resonance-energy signature.** Vertex-transitivity (regular) is a
coarse, shared property; the character spectrum χ (the additive/Dirichlet character values of the
speed set = the block eigenvalues) is the finer invariant, and it equals the resonance energy
(HYP-2155): **flat χ ⇔ low resonance energy ⇔ safe; structured χ ⇔ high resonance energy ⇔ tight.**

Verified (`lrc_chi_regular_paley_vs_ap_s608.py`):

- **Doubly-regular Paley `P_7 = {1,2,4}`** (the QR tournament, flat χ, a 2-design): resonance energy
  `E = 0` → the **safest** regular config, *not* tight. The most symmetric regular tournament is the
  *least* dangerous — flat χ means the runners are maximally decorrelated.
- **Rotational AP `{1..n−1}`** (structured χ, an interval connection set): resonance energy `E` large
  and growing (`0.28, 0.51, …` ≫ main) → **tight**. The structured χ is exactly the high resonance
  energy that pins the witness.

So among the regular (vertex-transitive) tournaments, χ separates two classes that vertex-
transitivity cannot tell apart: the **flat-χ doubly-regular (Paley) class — safe** — and the
**structured-χ rotational (AP) class — tight**.

## Is there a tight regular config off the AP/Paley orbit?

No genuinely new one. The hunt and its resolution:

- The only tight regular config that is *not* the AP-interval is `P_11 = {1,3,4,5,9}` (QR mod 11,
  the n=6 sporadic): vertex-transitive (the QR group acts as one orbit), flat χ `= √11/2`, **tight**
  (`G = 1/6`). Its χ genuinely differs from the AP's. **But it is an AP-LIFT**: `{1,3,4,5,9} ≡
  {1,2,3,4,5} (mod 7)` — the arithmetic progression transported to another additive face (HYP-2150),
  not a new orbit. (Its tightness is the AP-structured face, not the QR-flat one — the same config
  is doubly-regular mod 11 and rotational mod 7.)
- The **genuinely non-AP tight configs** — the sporadics `{1,3,4,7}` (n=5), `{1,2,3,4,5,7,12}` and
  `{1,4,5,6,7,11,13}` (n=8) — carry **no vertex-transitive symmetry** (not regular).

So: every tight config that is *regular* is the AP or an AP-lift; every tight config that is
genuinely off the AP orbit is *not regular*. The doubly-regular Paley class is safe.

## Why (the unification with the resonance energy)

The flat χ-spectrum of a doubly-regular tournament is a 2-design / perfect-difference-set condition:
the speed set is maximally additively *un*structured, so it has few short resonances `Σ m_i v_i = 0`,
hence low resonance energy `E(v)`, hence (HYP-2155) positive lonely measure — **safe**. The
rotational AP is maximally additively *structured* (every short sum returns: the length-3 fusions),
hence maximal resonance energy, hence the measure-0 tight witness. χ is the spectral shadow of the
resonance energy; reading χ tells you the resonance energy tells you tight-vs-safe. Vertex-
transitivity is blind to this — both Paley and AP are vertex-transitive — which is exactly why the
human's question is sharp: **the invariant that decides the Lonely Runner among the regular
tournaments is χ, not the symmetry.**

## Open / next

- Prove "flat χ ⇒ low resonance energy ⇒ safe" for doubly-regular tournaments in general (a clean
  Gauss-sum / 2-design bound on `E(v)`). The Paley speed set should be provably in the bulk.
- Characterize the tight-regular configs as exactly the AP-lifts (AP transported across additive
  faces, HYP-2150): is every regular tight config `≡ {1..n−1} (mod m)` for some `m`?
- Connects the perspective/rigidity arc (HYP-2130/2135/2140) to the resonance energy (HYP-2155):
  symmetry = automorphism rigidity; χ = the resonance-energy refinement that grades it.

**Artifacts:** `04-computation/lrc_chi_regular_paley_vs_ap_s608.py` (+`.out`),
`07-reflections/lrc-chi-is-the-resonance-energy-s608.md`. Builds on HYP-2155 (resonance energy),
HYP-2124 (regular witness), HYP-2150 (additive face / lifts), HYP-2055 (non-AP tight).
