# The Euler–Zariski bootstrap for JC(2): the last stand reduces to one parabola

**Instance:** klein-2026-07-20-S329 (owner: work creatively on JC(2), pull from past
threads unexpectedly). Concurrent: mac-mini S137 (HYP-8155, Hurwitz principle —
adjacent, claimed mid-session; cross-linked not duplicated), opus S417/S418 (fiber
formula; deck Aut = 1). Computations exact (script + frozen out).

## 1. The bootstrap (conditional theorem, minimal cover-degree case)

Let H : ℂ² → ℂ² be Keller with GENERIC COVER DEGREE 3 — the minimum possible, since
the Smith rule (T1549) kills degree 2 in every dimension. Let J = the Jelonek curve
(non-properness locus). Then:

1. If J = ∅: H is proper + étale over ℂ² ⟹ a covering of simply-connected ℂ² ⟹
   injective. So a counterexample has J ≠ ∅, and BY DEFINITION of the Jelonek set,
   H is proper over ℂ² ∖ J ⟹ **a genuine connected 3-sheeted covering of ℂ² ∖ J**,
   non-Galois with monodromy S₃ (Smith).
2. ⟹ π₁(ℂ² ∖ J) has a NON-NORMAL index-3 subgroup ⟹ **π₁ nonabelian** ⟹ by
   Fulton–Deligne (nodal curve complements have abelian π₁): **J is not nodal — it
   must carry cusps or worse.** This is the SAME Zariski-cusp mechanism mac-mini's
   THM-1330 found for the dim-3 Jelonek quartic's plane sections — now forced one
   dimension down.
3. **The Euler ledger**: χ(ℂ² ∖ H⁻¹J) = 3·χ(ℂ² ∖ J) (covering) and χ(ℂ²) = 1 give
   **χ(H⁻¹J) = 3·χ(J) − 2.**
4. The minimal topological model satisfying every constraint is the CUSPIDAL CUBIC
   J = {4p³ + 27q² = 0}: χ = 1 (affine rational cuspidal ≅ ℂ), π₁(ℂ² ∖ J) = B₃
   (trefoil/braid group) with its standard non-normal index-3 surjection B₃ ↠ S₃ ✓,
   ledger χ(H⁻¹J) = 1 ✓ — and this model is REALIZED algebraically by the
   **universal-cubic root cover** G(t, p) = (p, −t³ − pt): 3:1, monodromy S₃,
   Δ-pullback = −(p+3t²)²(4p+3t²) (two parabolas meeting once, χ = 1 ✓), fiber
   census profile → (1/6, 1/2, 1/3) off the discriminant. G satisfies EVERY
   bootstrap constraint except one: **det J_G = 3t² + p** — a single parabola of
   ramification inside ℂ².

## 2. The reduction, stated plainly

> **JC(2) at cover degree 3 ⟺ the ramification parabola of the universal root
> cover cannot be pushed to infinity within two variables.**

The dim-3 counterexample F is EXACTLY this push accomplished WITH a third variable:
its fiber cubic D·W³ + (b²−12a)W − 4a has a VARYING leading coefficient D(a,b,c)
that absorbs the double root — "ramification entirely at infinity" (THM-1315). In
two variables the leading coefficient of the fiber cubic of any candidate must be a
function on the base ℂ² whose vanishing locus IS (part of) J; the question is
whether an algebraic model of the B₃ → S₃ cover exists with all ramification
escaping to infinity. Newton-polygon/points-at-infinity theory (Abhyankar–Moh
school) is precisely the classical study of that escape.

**Degree bookkeeping (important):** Moh's theorem bounds POLYNOMIAL degree (≤ 100);
the bootstrap constrains COVER degree. A cover-degree-3 counterexample may have
arbitrarily large polynomial degree, so the reduction lives exactly in the open
region Moh does not cover.

## 3. The unexpected pulls from our past threads (motifs, tagged honestly)

- **N₂ ≡ 0 as an étaleness detector** (LRC structural-zero epistemology → JC): the
  étale F has N₂ ≡ 0 at every prime; the ramified root cover G shows a nonzero
  N₂-column (measured: 6/10/12 at p = 7/11/13) — the census SEES ramification as
  the failure of a structural zero. Instrument transfer, exact.
- **Escape-to-infinity ↔ the escape atlas / detection floor**: LRC families that
  evade every finite modulus ↔ ramification that must evade every finite chart.
  Both problems localize their obstruction AT THE BOUNDARY; both "sandwich"
  strategies (bounded exhaustion + structure at infinity) rhyme.
- **The resolvent character √(−D) ↔ Rédei parity**: the S₃-cover's quadratic
  subextension is a sign local system; Rédei's odd-Hamiltonian-path count is the
  tournament sign rigidity. Both are the mod-2 shadow of a nonabelian object.

## 4. Named next steps

(i) Formalize the bootstrap's step-2 citation set (Fulton–Deligne 1981; Zariski);
(ii) classify χ = 1 cusped plane curves with B₃-type complement groups admitting
non-normal index-3 subgroups (candidate-Jelonek atlas — likely very short list);
(iii) the points-at-infinity analysis of the root cover under Aut(ℂ²)-corrections
(Jung–van der Kulk: iterate affine/triangular pushes and track det J_G's zero
locus — can its Newton polygon be emptied?); (iv) cover-degree-4 bootstrap (allowed
monodromies A₄/S₄ from T1549 ⟹ which π₁'s admit them ⟹ which J's); (v) sync with
mac-mini HYP-8155's Hurwitz principle (their χ-accounting vs mine — likely the
same ledger from the compactified side).
