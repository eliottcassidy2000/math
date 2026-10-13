# Audit C (2026-10-12) of THM-4625 (tangent-sweep Keller maps, fibre by fibre)

**Auditor:** independent adversarial audit (Claude subagent, worktree `math-wt-logic4`), requested by
session mac-mini-2026-10-12-logic4. **Scope:** `01-canon/theorems/THM-4625-*.md` (at commit
`bbc8b997fb`, unchanged on disk during the audit), with THM-4620, THM-4621, HYP-9255 and the three
session scripts in `04-computation/experiments/logic4_20261012/` read as context; Gao, arXiv:2608.00222
(v1, 31 July 2026) read in its HTML version for the citation check. Nothing in the audited files was edited.
All computations below are my own (sympy 1.14, exact unless marked numeric).

## Verdict

**PASS WITH CORRECTIONS.** Every identity and every structural claim that the theorem rests on checks
(for symbolic `p` of degree `d = 2..6`, and `d = 7, 8` with a random integer specialization where the
symbolic discriminant is too slow), the PROVED typing is justified by the general arguments in the file
(the machine checks are consistency only), and the citations to Gao are for the right things. But
statement (3) as written is false at finitely many points (the fibre has points with `bx = 2`, which the
chart excludes: the plane model has a node there), statement (2) overstates the puncture count for
non-generic `a`, the "no coordinate-pair fibre ... union of lines" clause of (4) is only true for generic
fibres, the resultant sentence in (5) misreads the squared factor as coincidence pairs, and — the one
pleasant correction — the hypothesis of (5) is automatically satisfied, so the Galois group is `S_(d+1)`
for every member of the family without any side condition.

## Corrections

1. **(3), the chart misses the `bx = 2` points (false as stated; harmless downstream).** Quote:
   "the fibre over `(a, b)` minus its finitely many points with `F_3 = 0` is isomorphic to the plane curve
   `Ã x^2 + B̃ x + C̃ = 0, x ≠ 0, bx ≠ 2, p(w) ≠ 0`", and in the proof "(1) gives `2c/x = cb - p(w)`, so
   `c = p(w)x/(bx - 2)` (`p(w) ≠ 0`, `bx ≠ 2`)". The equation is `c(bx - 2) = p(w)x`; when `bx = 2` it
   forces `p(w) = 0` but leaves `γ` (hence `c = F_3 = 2γ/b`) free, and `F_1 = a` then reads
   `q(w_0) + γ w_0 = 4aγ^2/b^2`, a quadratic in `γ` for every nonzero root `w_0` of `p`. So the fibre minus
   `{F_3 = 0}` contains up to `2(d - 1)` points with `bx = 2`, generically exactly `2(d - 1)`, all with
   `F_3 ≠ 0` (needs `q(w_0) ≠ 0`). Verified: THM-1300 at `(a, b) = (8, 2)`: two such points, `|F_3| = 0.272`;
   exact instance `(a, b) = (2/3, 2)`: `(1, 0, 2/3)` and `(1, 1, -5/3)` with `F = (2/3, 2, 4/3)` and
   `(2/3, 2, 2/3)`, `w = 4/3`, `p(4/3) = 0` (in THM-4621(1)'s `y`-parametrization these are the roots of
   `2y^2 - by + 6a - b^2 = 0`, i.e. `y ∈ {0, 1}`, confirming the two files are consistent only after the
   points are put back). Gao's `G` at `(2/5, 3/2)`: four points, `w = 3 ± √3` (the nonzero roots of `p`),
   `F_3 ≈ -2.19, 11.06, 1.19 ± 0.57i`. Geometrically: at `w = w_0` the quadratic degenerates to
   `q(w_0)(bx - 2)^2`, so the plane model has a **node** at `(x, w) = (2/b, w_0)` (the discriminant
   `(p/w)^2 h` has a double zero there, `h(w_0) ≠ 0`), and the two fibre points with `bx = 2, w = w_0` are
   its two branches; the smooth fibre is the normalization. Corrected statement: "the fibre minus its
   finitely many points with `F_3 = 0` maps onto the plane curve `Ãx^2 + B̃x + C̃ = 0` minus `{x = 0}` and
   minus its vertical components, bijectively off the nodes `(2/b, w_0)` (`p(w_0) = 0 ≠ q(w_0)`) and
   two-to-one over each node; equivalently the fibre minus `{F_3 = 0} ∪ {bx = 2}` is isomorphic to the chart
   as displayed." Birationality, `h`, the genus and irreducibility are unaffected (finitely many points on a
   smooth curve of pure dimension one), but any puncture count of the `(F_1, F_2)`-fibre read off the chart
   would be wrong by `2(d - 1)`.

2. **(2), `ac ≠ 0` is not enough for `d + 3` punctures.** Quote: "For `ac ≠ 0`, the fibre over `(a, c)`
   is bijectively parametrized by `{w ∈ A^1 : w ≠ 0, q(w) ≠ c^2 a}` ..., so it is `P^1` minus `d + 3`
   points". The bijection holds for all `ac ≠ 0`; the count needs `c^2 a` to be a non-critical value of `q`,
   i.e. `c^2 a ∉ {q(ρ) : p'(ρ) = 0}` (at most `d - 1` values; `a ≠ 0` only excludes the critical value
   `q(0) = 0`). THM-1300: `c^2 a = 4/27` (`ρ = 2/3`) gives `P^1` minus 4 points, not 5 (THM-4621(5) itself
   says "five points generically"); Gao's `G`: `c^2 a = -11/2 ± 4√2` (`ρ = 2 ∓ √2`) gives 5, not 6.
   Corrected: "For `ac ≠ 0` ... bijectively parametrized ...; for generic `a` (`c^2 a` not a critical value
   of `q`) it is `P^1` minus `d + 3` points." Same sentence: "its critical points are exactly the `d + 1`
   finite punctures" — there are `d + 2` finite punctures (`w = 0` is a simple pole, not a critical point);
   corrected: "exactly the `d + 1` finite punctures other than `w = 0`".

3. **(5), the hypothesis is automatic — the theorem is stronger than stated.** Quote: "If
   `w ↦ (p(w), q(w))` is generically injective (a birational parametrization of the sweep curve), then the
   Galois group ... is `S_(d+1)`" (also in the title: "whenever ... is a birational parametrization").
   Generic injectivity holds for every `p` of degree `≥ 1` with `q' = (w/2)p'`: the tangent direction
   `(p', q') = p'(1, w/2)` recovers the parameter. Proof: if infinitely many pairs `r_1 ≠ r_2` had
   `(p, q)(r_1) = (p, q)(r_2)`, the curves `{(p(r_1) - p(r_2))/(r_1 - r_2) = 0}` and
   `{(q(r_1) - q(r_2))/(r_1 - r_2) = 0}` would share an irreducible component `Z`, not the diagonal (the first
   restricts to `p'(r_1) ≢ 0` there) and not a horizontal or vertical line; differentiating along a local
   parametrization of `Z` gives `p'(r_1)ṙ_1 = p'(r_2)ṙ_2` and `(r_1/2)p'(r_1)ṙ_1 = (r_2/2)p'(r_2)ṙ_2`, hence
   `(r_1 - r_2) p'(r_1) ṙ_1 = 0`; generically on `Z`, `r_1 ≠ r_2` and `p'(r_1) ≠ 0` (`r_1` is non-constant on
   `Z`), so `ṙ_1 = 0`, contradiction. Machine check: the resultant in `r_2` of the two difference quotients
   is a nonzero polynomial for SYMBOLIC `p` of degree `2, 3, 4` (degrees `2, 6, 12` in `r_1`), for random
   integer `p` of degree `2..6`, and at a random specialization for `d = 7, 8`; it always carries the
   diagonal factor `p'(r_1)^2`. So (5) holds unconditionally: **the Galois group of `K(x,y,z)/K(F)` is
   `S_(d+1)` for every tangent-sweep Keller map.** Drop the hypothesis in (5) and the "whenever" clause in
   the title.

4. **(5), the resultant sentence misreads the squared factor.** Quote: "the pairs `r_1 ≠ r_2` with
   `p(r_1) = p(r_2)`, `q(r_1) = q(r_2)` are the roots of `9(r_1^2 - 4r_1 - 2)(r_1^2 - 4r_1 + 2)^2/64`,
   finitely many". The resultant is re-derived exactly, but `r_1^2 - 4r_1 + 2 = p'(r_1)/3`: at its roots
   `2 ± √2` (the critical points of `p`, Gao's two cusps) the common root of the two difference quotients is
   `r_2 = r_1`, the diagonal, not a coincidence pair. The genuine pairs are `r_1 = 2 ± √6`, `r_2 = 4 - r_1`
   (`p = -4`, `q = 1/2`; exact), the single node of Gao's sweep curve (Gao, Theorem 3.5: "two cusps and one
   node"). So the only bad value is `cb = -4`: for `b = -4/c` two of the three critical values of `F_1` on the
   `(F_2, F_3)`-fibre coincide (`1/(2c^2)`), and over `(1/(2c^2), -4/c, c)` the fibre of `F` is empty (Gao's
   fibre size `d - 3 = 0` at the node). Corrected: "the pairs are `r_1 = 2 ± √6`, `r_2 = 4 - r_1`; the factor
   `(r_1^2 - 4r_1 + 2)^2 = p'(r_1)^2/9` is the diagonal contribution".

5. **(4), "no coordinate-pair fibre ... union of lines" must read "no generic coordinate-pair fibre".**
   Quote: "no coordinate-pair fibre of any member of the family is a disjoint union of lines (the fibre form
   of HYP-9255)". Special fibres do contain lines: the `(F_2, F_3)`-fibre over `(b, 0)` contains the line
   `{x = 0, y = b}` (THM-1300, `F_2|_{x=0} = y`) resp. `{x = 0, 2y = b}` (Gao's `G`, `F_2|_{x=0} = 2y`), `z`
   free. THM-4620(4) and HYP-9255 concern the generic fibre `C_η`, so the intended statement is "no generic
   coordinate-pair fibre"; the sentence before it already says "generic".

6. **(3) proof, the division by `w^2`.** Quote: "All three coefficients vanish to order two at `w = 0`
   (`p(0) = 0`, `q = O(w^2)`) ... so the equation may be divided by `w^2`." The orders are exactly
   `(2, 3, 2)` when `p'(0) ≠ 0` (`B̃(0) = (p_1 - 4q_2)b = 0` because `q_2 = p_1/4`); when `p'(0) = 0` all
   three vanish to order `≥ 3` and the line `w = 0` survives the division as a component of the plane curve.
   Likewise the lines `w = w_0` for common nonzero roots of `p` and `q` are components (example
   `p = w(w - w_0)(w - w_0/2)`, `d = 3`: `p(w_0) = q(w_0) = P(w_0) = 0`); generically the content of
   `(Ã, B̃, C̃)` is a nonzero constant (checked `d = 2..6`). All such vertical lines lie in `{p(w) = 0}`, so
   the chart condition `p(w) ≠ 0` removes them and the non-vertical part is irreducible iff `h` is not a
   square (Gauss). Wording fix: "the divided equation may still contain vertical lines `w = w_0` with
   `p(w_0) = q(w_0) = 0` (including `w = 0` when `p'(0) = 0`); they are outside the chart, and 'the plane
   curve' means its non-vertical component".

7. **(5) proof, the transfer from one variable to `K(F)` is compressed.** Quote: "this is the Galois
   group of `K(x, y, z)/K(F)` because the fibre of `F` over `(a, b, c)` is the fibre of this polynomial with
   the critical points deleted and `K(x, y, z) = K(F)(w)`." Valid, but two steps are implicit and should be
   written: (i) `K(x,y,z) = K(F)(w)` because `x = 2F_3/(F_2F_3 - p(w))`, `γ = (F_2F_3 - p(w))/2` are
   identities in `K(x,y,z)` (checked for both instances), and `w` is a root of
   `Φ(W) = q(W) + W(F_2F_3 - p(W))/2 - F_1F_3^2 ∈ K(F)[W]` of degree `d + 1 = deg F`, hence `Φ` is the
   minimal polynomial; (ii) the Galois group of a separable specialization `Φ_(b_0, c_0)` over `K̄(a)` is a
   subgroup of the Galois group of `Φ` over `K̄(a, b, c)` (decomposition group with trivial inertia; Lang,
   Algebra VII §2), and the former is `S_(d+1)` by the one-variable monodromy, so the geometric and the
   arithmetic Galois groups are `S_(d+1)`. Independent confirmation by rational specializations (sympy
   `galois_group`): `S_3` (THM-1300), `S_4` (`G`), `S_5` (`p = w^4 - 2w^3 + 3w^2 + w`), `S_6`
   (`p = 2w^5 + w^4 - 3w^3 + w^2 - w`).

8. **Status block / scripts.** `gao_degree4_fibrations.out` prints "polynomial components:
   `[True, True, False]`" for Gao's `G`; this is a false negative of the test `sp.denom(...) == 1` (the
   coefficient denominator 8 of `q`), not a defect of the map: all three components are polynomials
   (degrees `12, 11, 4` in the file's order, `det = -2` in that order, `+2` in Gao's). The status block
   records the script's wrong "genus 2" line but not this one; add a note. Also declare `γ_0, α, β ∈ K`
   constants in the Setting, and state the two side conditions explicitly (`γx | p(w) + 2γ` and
   `(γx)^2 | q(w) + γw` in `K[x,y,z]`; Gao §3.3: finitely many linear conditions on the coefficients of `p`
   and on `(γ_0, α)`). Provenance nit inherited from THM-4621 (flagged by audit A): in (5) "This holds for
   THM-1300 (`S_3`, THM-4621)" should also credit THM-1310 (July 2026), where `deg F = 3`, the `S_3` Galois
   group and the Jelonek set of THM-1300 were first established in-repo.

## What checks (no change needed)

* **(1), (2) parametrizations.** Bijections onto the full fibres (for `c ≠ 0`, every `b`; for `ac ≠ 0`), both
  directions morphisms, so isomorphisms with open subsets of `A^1_w`; excluded points exactly as stated
  (`p(w) = cb`, resp. `w = 0` and `q(w) = c^2 a`, plus `∞`). Derivative identities, `deg_w F_1 = d + 1`,
  leading coefficient `-ℓ/(2(d+1)c^2)`, `F_2`'s pole of order `d` at `∞` with coefficient `ℓ/((d+1)c)` and
  simple pole `2ca/w` at `0`: exact for symbolic `p`, `d = 2..8`.
* **(3) discriminant.** `B^2 - 4AC = p^2(b^2w^2 - 8a(wp - 2q))` and `wp - 2q = P` are identities for all
  `d` (also by hand: `(wp - 4q)^2 - 4q(4q - 2wp) = w^2p^2`); `deg h = d - 1`; `h` squarefree iff `b^2/(8a)`
  is not a critical value of `P/w^2`; genus `⌊(d-2)/2⌋` by Riemann–Hurwitz with the branch point at `∞`
  counted iff `d - 1` is odd (`d` even: `d` branch points, `g = (d-2)/2`; `d` odd: `d - 1`, `g = (d-3)/2`).
  Gao's `G`: `h = -4w^2/5 + 32w/5 - 147/20` at `(2/5, 3/2)`, genus 0 — the session's "genus 2" is indeed
  wrong. The chart-to-fibre map was verified on 10 random chart points per instance.
* **(3) irreducibility of the whole fibre.** The fibre is smooth of pure dimension one (Jacobian rank), meets
  `{F_3 = 0}` in finitely many points (3 for each instance at the tested base points: one with `x = 0`, two
  with `γ = 0`; finiteness also from quasi-finiteness of `F`), so it is irreducible iff the chart is; it cannot
  be reducible "as a whole" since every component is one-dimensional and meets the chart.
* **(4).** `F_j` restricted to a fibre of `G_j` is unramified (`dF` invertible, `T C = ker dG_j`), its degree
  as a map from the generic fibre is the generic fibre size `d + 1` of `F`, and an unramified polynomial on
  `A^1` has constant derivative: valid, if redundant given (1)–(3). THM-4620(1) with `d_j = 1` gives
  `L_j = K(G_j)`; the step `K(G_j) ∩ K[x] = K[G_j]` (THM-4621(7)) is correct and genuinely needs the
  openness of `G_j` (smooth ⇒ flat ⇒ open): for the non-flat `G = (x, xy)` the same intersection contains
  `y`. Cokernel ranks `d`, `d + 2` match `2γ + r - 1` with `(0, d+1)`, `(0, d+3)` and THM-4621's `2, 4` at
  `d = 2`. "HYP-9255 holds on the family" is vacuous-true (no `V_j` surjective, no member invertible).
* **(5) one-variable part.** Simple critical points iff `cb` is not a critical value of `p`; critical values
  `q(r_i)/c^2`; `d` transpositions generating a transitive subgroup of `S_(d+1)`: standard and correct.
* **Citations (Gao, §3, HTML v1).** Lemma 3.1: tangency polynomial `W = q + (w/2)(X - p) - Y`, degree
  `d + 1`, leading coefficient `-lc(p)/(2(d+1))`, `W' = (X - p)/2` — the file's `F_1|_fibre = (W + Y)/c^2`
  at `(X, Y) = (bc, ac^2)`. Proposition 3.2: generic fibre `d + 1` (after the twist `d - 1` at smooth curve
  points, `d - 2` at cusps, `d - 3` at nodes) — cited correctly. §3.3: varying `deg p` gives geometric
  degree `d + 1` for every `d ≥ 2`, which Gao calls Gallagher's family — "by Gao's account" is the right
  hedge (Gao's introduction describes Gallagher's July-20 construction only as a "uniform one-variable
  construction" of geometric degree `d` for `d ≥ 3`). Theorem 3.5: `G` of degrees `4, 11, 12`, `det 2`,
  generic fibre 4, two cusps and one node — consistent with correction 4. The paper has no monodromy,
  Galois, genus or dual-field content (both halves read), so (3)–(5) are not in Gao.

## Scripts (all in `04-computation/experiments/logic4_20261012/`, outputs beside them)

* `audit_C_identities.py` → `.out` (3 min): symbolic `p`, `d = 2..6`: (1)/(2) degrees, leading coefficients,
  derivatives; (3) discriminant identity, `wp - 2q = P`, orders of vanishing at `w = 0` with and without
  `p'(0) = 0`, `deg h`, symbolic discriminant of `h` nonzero, `disc~ = (p/w)^2 h`, Riemann–Hurwitz genus,
  generic content of `(Ã, B̃, C̃)`; the vertical-line example; symbolic resultant for (5), `d = 2, 3, 4`.
* `audit_C_identities_fast.py` → `.out` (1 s): `d = 7, 8`, exact identities plus discriminant and resultant at
  a random integer specialization.
* `audit_C_fibres.py` → `.out` (27 s): THM-1300 and Gao's `G`: polynomiality and `det`, `K(x,y,z) = K(F)(w)`
  identities, the `bx = 2` points (numeric 40 digits + exact instance), the `F_3 = 0` points (resultants,
  exact count over the splitting field), the puncture-count drop at critical values of `q`, chart→fibre
  verification, `h` and genus, the special fibre over `(b, 0)`.
* `audit_C_monodromy.py` → `.out` (1 s): Gao resultant re-derived and its roots sorted into node pairs vs
  diagonal; THM-1300 (no coincidences for any `b`); random `p`, `d = 2..6`; Galois groups of specializations
  `S_3, S_4, S_5, S_6`.

## Summary

The theorem is right in substance and its PROVED typing is earned by general arguments (nothing rests on
the `d ≤ 6` machine checks). Three statements need repair: (3) must exclude, or correctly account for, the
`2(d-1)` fibre points with `bx = 2` (nodes of the plane model, two fibre points each); (2) needs `c^2 a` not
a critical value of `q` for the `d + 3` count (and "`d + 1` finite punctures" is miscounted); (4)'s
union-of-lines clause is for generic fibres only; (5)'s resultant sentence confuses the diagonal factor
`p'(r_1)^2` with coincidence pairs. The main upgrade: the birational-parametrization hypothesis of (5) is a
theorem (the tangent slope `q'/p' = w/2` recovers `w`), so every tangent-sweep Keller map has Galois group
`S_(d+1)`; this was confirmed by specialization for `d = 2..5`. Gao is cited for the right things; nothing
in (3)–(5) is in Gao.

## Git state at hand-off

No commit, push, checkout, reset or rebase was run; no existing file was edited. New untracked files:
the four `audit_C_*.py` scripts with their `.out` files and this report.
