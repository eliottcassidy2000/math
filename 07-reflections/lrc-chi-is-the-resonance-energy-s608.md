---
source: claudebox-2026-06-03-S608
status: REFLECTION — χ is the resonance-energy signature; among regular tournaments the flat-χ
  Paley is safe and the structured-χ AP is tight; no tight config is regular off the AP orbit.
tags: [chi, character-spectrum, resonance-energy, Paley, doubly-regular, rotational, AP, tight,
  safe, vertex-transitive, 2-design, LRC]
---

# χ is the resonance energy: Paley is safe, the AP is tight

**Question (human):** among the maximally-cyclic (regular) tournaments, does χ add anything beyond
vertex-transitivity — a tight config that is regular but off the Paley/AP orbit, and does its χ
differ?

The question is sharp because it pits two notions of symmetry against each other. Vertex-transitivity
says all runners look alike — the witness is regular. But there are *two* maximally-symmetric regular
tournaments, the rotational (the AP) and the doubly-regular Paley, and the question is whether the
finer invariant χ tells them apart and whether tightness cares. It does, and it does.

## The two regular tournaments, and which is tight

The doubly-regular Paley tournament `P_7 = {1,2,4}` has a flat character spectrum — `|χ| = √7/2` for
every nontrivial χ, the 2-design condition. Its resonance energy is **zero**: it is the *safest*
config, maximally lonely, nowhere near tight. The rotational AP has a structured χ and a large,
growing resonance energy: it is the tight one. So the most symmetric regular tournament — Paley, the
one with the largest automorphism group and the most uniform spectrum — is the *least* dangerous,
and the merely-rotational AP is the worry. Vertex-transitivity cannot see this; χ is exactly the
difference.

Why: a flat χ-spectrum is a perfect-difference-set / 2-design, which means the speed set is
additively *unstructured* — few short resonances `Σ m_i v_i = 0`, low resonance energy, positive
lonely measure (HYP-2155). The AP is the opposite: every short sum returns, the length-3 fusions are
all present, the resonance energy is maximal, the witness is the measure-0 n-gon vertex. χ is the
spectral shadow of the resonance energy, and the resonance energy is tight-vs-safe. So **χ is the
invariant that decides the Lonely Runner among the regular tournaments — not the symmetry.**

## Is there a tight regular config off the AP orbit?

I chased the one candidate and it dissolved beautifully. `{1,3,4,5,9}` is the n=6 sporadic tight
config; it is QR mod 11, the Paley `P_11`, vertex-transitive with a flat χ — genuinely different
from the AP's χ. For a moment it looks like a tight, regular, non-AP, different-χ config. But it is
an **AP-lift**: `{1,3,4,5,9} ≡ {1,2,3,4,5} (mod 7)`. The same config is doubly-regular mod 11 and the
rotational AP mod 7 — it is the arithmetic progression wearing the Paley face, the AP transported
across additive faces (the S602 two-faces machinery, now showing up as a single config with two
symmetry descriptions). Its *tightness* comes from the AP (structured, mod-7) face, not the Paley
(flat, mod-11) one. And the genuinely non-AP tight configs — the sporadics at n=5, n=8 — are *not*
regular at all. So every tight config that is regular is the AP or an AP-lift; everything genuinely
new is non-regular.

## The closing of the arc

This answers the human's question and ties the whole recent thread into one statement. The
perspective/rigidity sessions said symmetry = automorphism rigidity = the worry-set is regular. This
session says: regularity is necessary but blind; the *grading* within the regular tournaments is χ,
and χ is the resonance energy. Flat χ (Paley) — the most symmetric — is safe; structured χ (AP) is
tight. The Lonely Runner does not live in the symmetry; it lives in the character spectrum, which is
the resonance energy, which you sidestep by construction on the core. One invariant, read three
ways: automorphism (coarse), character (fine), resonance energy (analytic).

**Artifacts:** `04-computation/lrc_chi_regular_paley_vs_ap_s608.py` (+`.out`); new **HYP-2160**.
Builds on HYP-2155 (resonance energy), HYP-2124 (regular witness), HYP-2150 (faces/lifts), HYP-2055.
