# Audit D (2026-10-12) of THM-4626 (line fibres / surjective dual fields are locally nilpotent; HYP-9255 ⟺ JC₂)

**Auditor:** independent adversarial audit (Claude subagent, worktree `math-wt-logic4`), requested by
session mac-mini-2026-10-12-logic4. **Scope:** `01-canon/theorems/THM-4626-*.md`. The audit started on the
first version (commit `66c0596e0c`: "surjective ⇒ locally nilpotent" by the leading-coefficient filtration)
and, at the coordinator's request, finished on the strengthened version on disk (commit `1e1c7acc20`: new
condition (b') `μ_j = 0`, direct local-nilpotence argument). Both are covered below. Context read, not
edited: THM-4620, THM-4562, HYP-9255, the session note §4.7, audits A–C. Literature checked from the
sources themselves: Berson–van den Essen–Maubach (Nijmegen report 9918, PDF), Gurjar–Masuda–Miyanishi
(arXiv:1211.0744, PDF), Stein 1989 (Springer abstract), Kambayashi 1975 (citing works). Script:
`04-computation/experiments/logic4_20261012/audit_D_toy_and_padding.py` (+ `.out`, sympy 1.14, ~1 s).
Nothing in the audited files was edited.

## Verdict

**PASS WITH CORRECTIONS.** The mathematics of statement (1) is right in both versions, and the strengthened
form is the correct one: under `μ_j = 0` alone THM-4620(4) gives `K[x] ⊂ K[x]_S = L_j[F_j]` with `V_j = ∂/∂F_j`,
so `V_j` is locally nilpotent on `K[x]` as the restriction of a locally nilpotent derivation, and a cokernel of
generic rank zero is zero. I had reached the same conclusion independently before the strengthening landed
(the first version's "integrating" step and its surjectivity hypothesis were never used after `μ_j = 0`, and
its "Still open" paragraph was closed by its own differentiating step). The remaining corrections are a
wording issue in (b'), a missing citation that (2) needs over a non-algebraically-closed `K`, a wrong journal
in the status block, and stale caveats left in THM-4620 / HYP-9255 / the session note.

## Answers to the coordinator's two questions

1. **Does `K[x]` embed in `L_j[F_j]` with `V_j = ∂/∂F_j` under `μ_j = 0` alone?** Yes. Re-derived: (i) `K[x]`
   is a domain, so `K[x] ⊂ K[x]_S = K[x] ⊗_(K[G_j]) K(G_j) = O(C_η)`; (ii) `L_j ⊂ K[x]_S` because an element of
   `L_j` is integral over the field `K(G_j) ⊂ K[x]_S` and `K[x]_S` is normal (THM-4620(1)); so `C_η` is an
   `L_j`-curve, geometrically integral since `L_j` is algebraically closed in its function field `K(x)` and
   char 0; (iii) `μ_j = 0` means each of the `d_j` geometric components (one per embedding `L_j → L̄_j`) has
   `γ = 0`, `r = 1`, so `C_η ⊗ L̄_j ≅ A^1`; the smooth completion of `C_η` is a genus-0 curve over `L_j` whose
   boundary is one closed point with a single geometric point, hence of degree 1, hence rational; a conic with
   a rational point is `P^1`, so `C_η ≅ A^1_(L_j) = Spec L_j[h]`; (iv) `V_j` extends to `K[x]_S` by the quotient
   rule (`S ⊂ K[G_j] ⊂ R_j`), kills `L_j`, so `V_j = V_j(h)·d/dh`, and `V_j F_j = 1` forces `V_j(h) ∈ L_j^×`,
   `F_j = αh + β`, `L_j[h] = L_j[F_j]`, `V_j = ∂/∂F_j`. No step uses surjectivity.
2. **Is "a cokernel of generic rank zero is zero" right?** Yes: `μ_j = 0 ⇒ V_j` locally nilpotent (restriction of
   `∂/∂F_j` on `L_j[F_j]`) `⇒ K[x] = R_j[F_j]` (slice theorem, `V_j F_j = 1`) `⇒ V_j` surjective (char 0)
   `⇒ coker V_j = 0`. The generic rank is `dim_(K(G_j)) coker V_j ⊗ K(G_j) = μ_j` by THM-4620(3), so the chain is
   exactly (b') ⇒ (d) ⇒ (c) ⇒ (a) ⇒ (b). The parenthetical explicit induction in the file
   (`c_m = V_j^m a/m! ∈ K[x] ∩ L_j = R_j`) is also correct and gives (c) without the slice theorem.

## Corrections

1. **(b'), wording: "disjoint union of lines" is a geometric statement.** Quote: "`μ_j = 0`: the generic fibre of
   `G_j` is a disjoint union of (`d_j`) lines". Over `K(G_j)` the generic fibre `Spec K[x]_S` is irreducible
   (`K[x]_S` is a domain); it is `A^1` over the degree-`d_j` extension `L_j`, and becomes `d_j` disjoint lines
   only after base change to an algebraic closure of `K(G_j)`. Corrected: "every geometric component of the
   generic fibre of `G_j` is a line (equivalently the generic fibre is `A^1_(L_j)`, `[L_j : K(G_j)] = d_j`), i.e.
   the cokernel has generic rank zero" — as THM-4620(4) phrases it. Same phrase in the title and in HYP-9255's
   title ("disjoint union of two or more lines"): add "geometrically".

2. **(2) over a non-algebraically-closed `K`: cancellation needs Kambayashi.** Quote: "surface cancellation
   (Fujita; Miyanishi–Sugie) gives `X_j ≅ A^2`". Fujita and Miyanishi–Sugie are over an algebraically closed
   field; they give `X_j ⊗ K̄ ≅ A^2_(K̄)`. To conclude `X_j ≅ A^2_K` one needs T. Kambayashi, *On the absence of
   nontrivial separable forms of the affine plane*, J. Algebra 35 (1975) 449–456 (every separable `A^2`-form is
   trivial; char 0 ⇒ separable). The alternative "base-change to `K̄` first" does not work as written, because
   the hypothesis `JC_2` is assumed over `K` and `JC_2(K) ⇒ JC_2(K̄)` is not immediate. Corrected: add Kambayashi
   to the status block and the proof of (2), or state (2) for `K` algebraically closed. The converse direction
   (padding) is fine over any `K`. Minor, same paragraph: "a Keller map (its Jacobian is `det JF`)" — in
   coordinates `(u, v)` of `X_j ≅ A^2` the Jacobian of `G_j|_(X_j)` is `det JF / det ∂(u, v, F_j)/∂(x, y, z)`,
   a nonzero constant times `det JF` (the denominator is the Jacobian of the coordinate change), not `det JF`
   itself.

3. **Status block: wrong journal and imprecise precursor statements.** Quote: "Berson-van den Essen-Maubach
   (Nijmegen report 9918 / J. Pure Appl. Algebra 2001)". The paper is Israel J. Math. 124 (2001) 115–124 (as the
   file's own `related` list says); the preprint is Nijmegen Report No. 9918, April 1999 (checked in the PDF).
   Their result is Theorem 3.8: property S(R) — "any surjective `R`-derivation of `R[X,Y]` having divergence
   zero is locally nilpotent" — for every *Noetherian* `Q`-algebra `R`, with the field case attributed to Stein
   (their introduction paraphrases Stein as "any surjective `k`-derivation on `k[X,Y]` is locally nilpotent",
   which is false without divergence zero: `∂_x + y∂_y` is surjective and not locally nilpotent, as the file's
   own (4) illustrates in three variables; the divergence-zero hypothesis is in their Theorem 3.8). Stein's
   actual theorem (J. Analyse Math. 52 (1989) 291–300): for `P ∈ C[x, y]` and `D_P f = P_x f_y − P_y f_x` on the
   ring `E` of entire functions on `C^2`, the image of `D_P` is dense in `E` iff `P = σ(x)` for an automorphism
   `σ` of `C[x, y]` — so "Stein (1989, Jacobian derivations)" is right in substance (surjectivity on polynomials
   gives density, hence `P` is a coordinate, hence `D_P` is locally nilpotent) but it is a density theorem on
   entire functions, not a statement about polynomial surjectivity. Gurjar–Masuda–Miyanishi Cor. 2.4 (over `C`):
   "a non-trivial *reduced* `C`-derivation of `C[x,y]` is locally nilpotent iff it is surjective and
   `Ker D ⊋ C`" (reduced = coprime coefficients); both side conditions hold for a Keller derivation `D_f`
   (`D_f g = 1`, `f ∈ Ker`), so the citation is apt; Remark 2.5 is the `JC_2 ⟺ D_f` surjective statement
   quoted in "Reading". Also note their abstract: Cerveau's classification of surjective derivations of
   `C[x,y]` is "suspended" (error in his Lemma 5.5), so Cerveau should not be cited for the planar result
   beyond what THM-4620 already says.

4. **Housekeeping in the file.** Status: "Not yet independently audited (audit D requested)" — update. The H1
   still reads "surjective dual fields are locally nilpotent" while the title now leads with line fibres
   (cosmetic). The status sentence "The dual-field statement in dimension `n ≥ 3` was not found in the
   literature" is consistent with my search (see below), but it should say where the content lies: the step
   (b') ⇒ (d) is the standard local-slice principle (a derivation that becomes `d/ds` on `L[s] ⊇ A` is locally
   nilpotent on `A`; Freudenburg's book, local slices), and the substance is THM-4620(4) (line fibres ⇒
   `K[x]_S = L_j[F_j]` with coordinate `F_j`).

5. **Stale caveats elsewhere (consistency, not errors in THM-4626).** THM-4620 statement (4), line 103: "a
   cokernel of generic rank zero can be a nonzero torsion module only when `d_j ≥ 2`" — now: never
   (THM-4626(1)); line 109 and status line 13 "dimension three is HYP-9255" — now: HYP-9255 ⟺ `JC_2`; Reading
   bullet line 242 "whether it can be a union of several lines is HYP-9255" — now: a coordinate-pair fibration
   of a Keller map of `A^3` has geometrically `d ≥ 2` line fibres iff `F ≅ G' × id` with `G'` a planar Keller
   map of degree `d` (THM-4626(1) + cancellation), i.e. iff a `JC_2` counterexample exists. HYP-9255 title
   ("equivalent to if the sub-question below has a positive answer") and line 29: the condition is now met, so
   HYP-9255 ⟺ "no coordinate-pair fibration with `d ≥ 2` geometric line components" ⟺ `JC_2`. Session note
   line 235, open list item 5 "Is `μ_j = 0` enough for `coker V_j = 0`?" — settled (yes) by §4.7 of the same
   note.

## Checks that passed (with how they were checked)

* **(1), first version (commit `66c0596e0c`).** The embedding, `K[x] ∩ L_j = R_j`, the `R_j`-module (hence
  `K`-vector space) structure of `M_m`, the differentiating inclusion `m M_m ⊆ M_(m−1)` (zero coefficients
  included, `M_m` being the set of `F^m`-coefficients over `A_(≤m)`), the integrating inclusion (`∂_F a' = a` in
  `L[F]` forces `deg a' ≤ m+1` with `F^(m+1)`-coefficient `c/(m+1)`, zero when `c = 0`), and the induction
  `a − cF^m ∈ A_(≤m−1)` are all correct. But the chain `M_m ⊆ M_(m−1) ⊆ … ⊆ M_0 = R_j` already puts every leading
  coefficient in `R_j`, so the integrating step and surjectivity were redundant — which is what the
  strengthening fixed. No MISTAKE entry seems required (nothing false was asserted; a consequence was missed
  and is now recorded).
* **Abstract version (brief item 4).** True, and stronger than asked: let `A ⊇ Q` be a domain, `D` a
  derivation, `s ∈ A` with `Ds = 1`, `L = Ker(D on Frac A)`. If `A ⊆ L[s]` (equivalently `A_S = L[s]` for some
  `S ⊆ R ∖ 0`, `R = Ker D`), then `D` is locally nilpotent on `A` and `A = R[s]`, with no surjectivity
  hypothesis: `D` acts on `L[s]` as `d/ds`, which kills `L[s]_(≤m)` after `m+1` steps, and the slice theorem
  (or `c_m = D^m a/m! ∈ A ∩ L = R` plus induction on `deg_s`) gives `A = R[s]`. So no counterexample "surjective,
  `A_S = L[s]`, `A ≠ R[s]`" exists; the toy derivations fail the hypothesis `A ⊆ L[s]` (`x∂_x − y∂_y + ∂_z`:
  `L = K(xy)`, `K[x,y,z] ⊄ K(xy)[z]`; `∂_x + y∂_y`: `L = K`). A dual field that is surjective, or merely has
  `μ_j = 0`, but is not locally nilpotent would therefore contradict THM-4620(4), which I re-derived above;
  I found no gap.
* **(4), the non-example** (script (A)): `div D = 0`; explicit preimage of every monomial `x^a y^b z^c` of
  degree `≤ 6` (closed form `p = Σ_k (−1)^k c!/(c−k)! z^(c−k)/(a−b)^(k+1)` for `a ≠ b`, `z^(c+1)/(c+1)` for
  `a = b`), all 84 verified; rank check `P_(≤5) ⊆ D(P_(≤6))`; `D^5 x = x`; `Ker D` on `P_(≤6)` has dimension 4
  with basis `1, xy, (xy)^2, (xy)^3`. Transcendence degree one, so not a dual field. Correct as stated.
* **(3), padding** (script (B)): THM-1300's formula as typed in THM-4562 has `det J = −2` (Keller); the
  padded map's dual field is `V_4 = ∂_w` with `V_j F_i = δ_ij` for all `i, j`; the generic fibre of `F'` over a
  random rational point has 3 points (resultant after eliminating `z`, which is linear in `F'_3`; numeric
  roots at 50 digits; `Res_y = x^12 ·` (irreducible cubic)), so `deg F' = d_4 = 3 > 1`: `F' × id` is a Keller
  map of `A^4` with a surjective dual field and is not an automorphism. (3) is correct, and does not depend on
  anything but THM-1300 being non-injective, which this count re-confirms (matching the `S_3` Galois group of
  THM-4621/audit C).
* **(2), `JC_2 ⇒ HYP-9255` and padding converse.** `G_j ∈ R_j^2` (killed by `V_j`), so `G_j` is constant along
  the `A^1`-factor; `K[x] = R_j[F_j]` gives `F = (G_j|_(X_j)) × id` after the coordinate change
  `(x,y,z) ↦ (u,v,F_j)`; `G_j|_(X_j)` is Keller (Jacobian = constant × `det JF`, see correction 2); the padding
  converse is immediate (`V_3 = ∂_z` for `G' × id`). With Kambayashi (correction 2) the equivalence holds over
  every field of characteristic zero.
* **Mechanism illustration** (script (C)): a non-product automorphism `F = φ∘(G×id)∘ψ` of `A^3` has
  `V_3 = −2z∂_y + ∂_z` with `V_3 F_i = δ_(3i)` and nilpotency indices `1, 3, 2` on `x, y, z` — consistency only.

## Literature (brief item 6)

* **Stein 1989:** correct title and venue, pp. 291–300; theorem as stated in correction 3 (density of
  `D_P(E)` in entire functions on `C^2` iff `P` is a coordinate).
* **Berson–van den Essen–Maubach:** Israel J. Math. 124 (2001) 115–124 = Nijmegen Report 9918 (April 1999);
  Theorem 3.8 (S(R) for Noetherian `Q`-algebras, field case credited to Stein); their closing comment: "if
  the two-dimensional Jacobian Conjecture is true then ... the condition `1 ∈ Im(D)` is equivalent to the
  surjectivity of `D`" (divergence zero), with Questions 1–2 on finite surjectivity certificates.
* **Gurjar–Masuda–Miyanishi 2013:** Theorem 1.2 (a normal affine curve over `C` with a surjective derivation is
  `C[t]`, `D = d/dt`), Theorem 2.3, Cor. 2.4 and Remark 2.5 as quoted above; the paper treats dimension one and
  two only (Example 1.8 is a non-normal surface in three variables); nothing on `n`-variable Keller maps or on
  dual fields.
* **Higher dimension:** no statement of "`V_j` surjective (or line fibres) ⟺ `V_j` locally nilpotent" for dual
  fields of Keller maps in dimension `n ≥ 3` was found (searches on surjective derivations / locally
  nilpotent / Keller / slice / Jacobian, van den Essen, Nowicki, Daigle, Bavula; GMM's and BvdEM's
  reference lists). The ingredients are classical (local-slice principle; `2γ + r − 1`; forms of `A^1`), and
  the equivalence "`F` invertible iff all `V_j` locally nilpotent" is folklore via the slice theorem, but the
  per-`j` fibre criterion appears to be new in this form. I did not consult van den Essen's book directly.

## Scripts

* `04-computation/experiments/logic4_20261012/audit_D_toy_and_padding.py` → `audit_D_toy_and_padding.out`
  (sections (A) toy derivation, (B) padded THM-1300 with the generic-degree count, (C) mechanism illustration;
  all exact except the 50-digit root count in (B)).

## Summary

The strengthened THM-4626 is correct: `μ_j = 0` ⟺ surjective ⟺ `coker = 0` ⟺ `K[x] = R_j[F_j]` ⟺ locally
nilpotent, by THM-4620(4) plus the one-line observation that a derivation acting as `∂/∂F_j` on `L_j[F_j] ⊇ K[x]`
is locally nilpotent; the first version's filtration was correct but over-engineered, and its "Still open"
sub-question was already closed by its own differentiating step. Required edits: say "geometric components"
in (b') and in the titles; cite Kambayashi 1975 for cancellation over non-closed `K` in (2) and fix the
Jacobian phrase; correct the BvdEM venue to Israel J. Math. 124 (2001) 115–124 and state Stein's theorem as a
density theorem on entire functions; update the audit line; and remove the three stale "nonzero torsion
cokernel / dimension three is HYP-9255" caveats in THM-4620, the conditional phrasing in HYP-9255's title and
reduction paragraph, and item 5 of the session note's open list. The dimension-`n` statement was not found in
the literature; its substance is THM-4620(4).
