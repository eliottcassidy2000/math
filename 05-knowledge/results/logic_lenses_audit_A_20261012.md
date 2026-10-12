# Audit A (adversarial, independent) — THM-4620, THM-4621, HYP-9255, HYP-9257 (session mac-mini-2026-10-12-logic4)

**Auditor:** independent Claude session, 2026-10-12, worktree `math-wt-logic4`. No existing file was edited.
**Scripts (new):** `04-computation/experiments/logic4_20261012/audit_A_fibrations.py` (+ `.out`),
`04-computation/experiments/logic4_20261012/audit_A_kernel.py` (+ `.out`). sympy 1.14, exact arithmetic; every
identity is asserted (the run aborts on failure; it printed `ALL ASSERTIONS PASSED`).

## Verdicts

| file | verdict |
|---|---|
| THM-4620 | **PASS WITH CORRECTIONS** — the mathematics of (1)–(4), (6) is correct; (5) (the planar criterion) is **KNOWN**, not new (Gurjar–Masuda–Miyanishi 2013, Remark 2.5, verbatim), and the status block's novelty claim must be withdrawn; statement (4) is weaker than its proof and (5)(d)⇒(a) uses the stronger form |
| THM-4621 | **PASS WITH CORRECTIONS** — every exact claim re-derived and confirmed; one displayed identity (the tangency prefactor `9c/2`) is wrong as written and the session's own script output records `False` for it; THM-1310 (July 2026) is under-credited for `c_3`, `deg F = 3`, `S_3` and the Jelonek set; (9) is a correct sketch, not "elementary" |
| HYP-9255 | **PASS WITH CORRECTIONS** — "equivalently" in title/statement is an overstatement (one direction only); Gao's map `G` is transcribed without `u = 1 + xy, w = γu`; plane case is KNOWN |
| HYP-9257 | **PASS WITH CORRECTIONS** (as OPEN) — local structure claim is right; the cusp curve is exactly `Sing{c_3 = 0}` (verified here); the Alexander-polynomial remark is heuristic |
| scripts/outputs | **PASS** — reproduced independently; `thm1300_module_checks.out` line "equal? False" contradicts the text it is cited for |
| session note §4.1–4.2, 5, 6 | **PASS WITH CORRECTIONS** — same items as above; Gao's AI-disclosure quote verified; attribution rule respected |

Owner-ruling check (item 5): none of the four files, nor the note, attributes THM-1300 to Alpöge as discoverer. THM-4621's
provenance note and the note's §1/§6 *report* Gao's external attribution and explicitly do not adopt it. **PASS.**

## 1. THM-4620 — corrections

1. **Novelty claim is false for the plane (status block, title, (5), "Reading").** The sentence "*the planar surjectivity
   criterion [was] not found in the literature*" must go. R. V. Gurjar, K. Masuda, M. Miyanishi, *Surjective derivations in
   small dimensions*, J. Ramanujan Math. Soc. 28A (2013) 221–246 (arXiv:1211.0744) contains, verbatim:
   * Remark 2.5: "*Let {f, g} be a Jacobian pair ... D = f_y ∂/∂x − f_x ∂/∂y. Then the Jacobian conjecture holds if and only if D is surjective.*"
   * Corollary 2.4: a non-trivial reduced derivation of `C[x,y]` is locally nilpotent iff it is surjective and `Ker D ⊋ C`;
     Theorem 1: `D = ∂/∂x_1` in some coordinates iff `D` surjective and `Ker D ⊋ C`.
   * Theorem 1.2: a normal affine curve with a surjective derivation is `A^1`, proved via `dim H^1(X; C) = 2g − 1 + n` —
     i.e. exactly the mechanism of THM-4620(3)–(4) in the plane. Cerveau, J. Algebra 195 (1997) 320–335, Prop. 1.2, is the
     de Rham interpretation of surjectivity; Dimca–Saito, Compositio 85 (1993) 299–309, is behind it. (Stein 1989, cited by
     Berson–van den Essen–Maubach for "surjective derivations of k[X,Y] are locally nilpotent", I could not bibliographically confirm.)
   **Corrected typing:** (5)(a)⇔(b)⇔(e)⇔(f) KNOWN (GMM 2013, Remark 2.5/Cor. 2.4); (b)⇔(c) is the standard identification
   `coker D_f = H^1_dR(A^2/A^1; f)`; (d)⇔(a) also follows from Abhyankar–Moh–Suzuki (generic fibre `A^1` ⇒ general closed fibres
   `A^1` ⇒ `f` a coordinate ⇒ `(f,g)` automorphism). What remains attributable to THM-4620 is the `n`-dimensional dictionary
   ((1)–(4), the `t`-primary lattice (2), the generic-rank formula `μ_j` and the `deg F = d_j` conclusion) and the AMS-free
   re-proof of (d)⇒(a). The title's "(iff D_f is locally nilpotent, KNOWN)" should become "(all KNOWN: GMM 2013)".
2. **Statement (4) is weaker than its proof; (5) uses the stronger form.** After the first clause ("Surjectivity gives
   `coker V_j = 0`, hence `μ_j = 0`") the proof of (4) uses only `μ_j = 0`. (5)'s step "(d)⇒(a): by (4) with d = 1" applies (4)
   under the hypothesis `μ(f) = 0`, which (4) as stated does not allow. **Fix:** restate (4) as "If `μ_j = 0` (in particular if
   `V_j` is surjective) then every geometric component is a line, `C_η ≅ A^1_(L_j)`, `K(x) = L_j(F_j)`, `deg F = d_j`; if moreover
   `d_j = 1`, `F` is an automorphism." **Consequence (strengthening):** in every dimension, `μ_j = 0` and `d_j = 1` force
   `F` invertible, hence `coker V_j = 0`; so the "What is not claimed" caveat ("the cokernel could in principle be a nonzero
   torsion module") and HYP-9255's sub-question concern only `d_j ≥ 2`.
3. **Base change in (2) is justified only for `K[x]`.** "`K̄ ⊗_K K[x] = K̄[x]`" does not cover `K(x)` (`K̄ ⊗_K K(x) ≠ K̄(x)`).
   **Fix:** "`K(x) ⊂ K̄(x)` compatibly with `V_j` (same coefficients), and `F` is Keller over `K̄`, so it suffices to prove
   `Ker (V_j − λ)^k = 0` on `K̄(x)`." The eigenvector lemma and the ODE lemma are correct as written (re-derived: the
   `T^d`-cancellation in `N(T)` and `V_j(a_i) = λ(d − i)a_i` check; the ODE lemma's `q | q'` ⇒ `q' = 0` ⇒ `q = 1` ⇒ `p = 0` is right).
4. **Typo, ODE lemma:** "`c ∈ K^×`" should read "`c ∈ K'^×`" (the constant lies in the coefficient field `K' = K(G_j)`).
5. **Status wording:** "Statement (4)'s converse fails in dimension ≥ 4" — the padding example does not refute a converse of (4);
   it refutes "surjective ⇒ automorphism" with the connectedness hypothesis dropped. Reword.
6. **Normalisation in the Setting/(5):** `V_2 = D_f`, `V_1 = −D_g` hold when `det J(f,g) = 1`; in general `V_2 = c^(−1) D_f`,
   `V_1 = −c^(−1) D_g` with `c = det J`. Harmless but state it.
7. **(6), Bavula:** journal data confirmed — C. R. Math. 362 (2024) 731–738, DOI 10.5802/crmath.556 (the published abstract says
   "cyclic", the arXiv abstract "1-generated"; the explicit length bound is in terms of `deg σ`). Add pages/DOI.
8. **Verified as stated (no change):** (1) `R_j` algebraically closed, `Frac R_j = L_j`, `[L_j : K(G_j)] = d_j` (geometric
   integrality over `L_j` + `L_j ⊗ K̄(G_j) ≅ K̄(G_j)^(d_j)`), the ladder; (3) `Ω^1_rel = K[x] dF_j`, `coker V_j = H^1_dR` relative,
   the `A_(n−1)`-structure via the lifts `V_i` (`L_(V_i)(h dF_j) = V_i(h) dF_j` because `[V_i, V_j] = 0`), holonomicity via
   `G_(j,+)O` (Bernstein), base change to `C_η` by flatness, `2γ + r − 1`; (4) the rational point at infinity (a closed point
   splitting into one geometric point has degree 1), `V_j = p(h) d/dh` with `p F_j' = 1`, ZMT (birational + quasi-finite +
   normal target ⇒ open immersion), codimension-2 complement via units, extension by normality; the padding example
   (`V_n = ∂_n`, generic fibre = `deg F'` lines). (5)'s Nagata–Nowicki step (`R = K[h]`, `f = p(h)`, `p'(h) ∈ K^×`) is right;
   Nagata–Nowicki, J. Math. Kyoto Univ. 28 (1988) 111–118 confirmed (char 0, two variables: ring of constants `= k[f]`).

## 2. THM-4621 — corrections

1. **In-repo prior art under-credited (status, title, (2), (4)).** THM-1310 (kind-pasteur, 2026-07-19) already proved:
   `deg F = 3` ((3)); arithmetic *and* geometric monodromy `S_3` via the resolvent `−L` ((4)); the Jelonek set is exactly
   `{L = 0}` ((6)); and `L(a,b,g) = 27a²g² − 18abg + 16a + b³g − b²` **is** `c_3` (checked: `L(a,b,c) − c_3 = 0`). Gao's
   Theorem 3.4(1) identity `disc_x(h_1) = −c_3 a_2²` is THM-1310's `Δ = −4Q²L` (`a_2 = 2Q`). THM-4621 cites THM-1310 only for
   mod-`p` fibre sizes and presents (2)/(4) as re-derivations of Gao. **Fix:** in (2) and (4) write "re-derivation of
   THM-1310(3),(4),(6) (and of Gao, Theorem 3.4)"; call the polynomial "THM-1310's `L` = Gao's `c_3`". New in THM-4621
   relative to THM-1310: the three fibre parametrizations, the restricted functions with critical points = punctures, the
   empty-fibre curve `C`, rings of constants, nilpotent parts, the residue witnesses, the length bound.
2. **(3) tangency identity is wrong as displayed.** "*the equation `F_3(y) = c` becomes `(9c/2)(w^3 − 2w^2 + bc·w − 2ac^2) = 0`*":
   the exact identity is
   `F_3(b − 3ac/w) − c = −(c/w^3)·(w^3 − 2w^2 + bc·w − 2ac^2)` (audit_A_fibrations.py, asserted). The factor `9c/2` is
   sympy's numerator normalisation (`together`) divided by `−2a`, not a mathematical prefactor. Worse, the cited support,
   `thm1300_module_checks.out`, prints "`equal? False`" for this check because the script compares `(9c/2)·T` with `T`.
   Replace the displayed identity and fix the check (compare after `cancel(.../(-c/w**3))` or test divisibility).
   Gao's text says the tangency polynomial equals the cubic "up to a constant factor"; the proof's parenthetical "(the factor 9c/2)" must go.
3. **(2) stratum `a = 0` is cited to Gao but is now in-repo (audit §[2']).** For `a = 0`: the `u = 0` branch gives exactly one
   point iff `b ≠ 0` (`x = 2/b, y = −b/2`); the `u ≠ 0` branch forces `y = b` and `F_3 = c` becomes
   `b(1 − bc)x² + 2(1 − bc)x − c = 0` (discriminant `4(1 − bc)`). Hence `#F^(−1)(0,b,c) = 3` if `b(1 − bc) ≠ 0`, `= 1` if
   `bc = 1` or `b = 0`; never 0 or 2; `c_3|_(a=0) = b²(bc − 1)`. So "fibre sizes 3/1/0, Jelonek set `{c_3 = 0}`, image
   `A^3 ∖ C`" is PROVED in-repo for all `(a,b,c)`, and the sentence "(Gao, Theorem 3.4, whose stratification is complete)" can
   become a cross-check. Also verified: `{c_3 = 0} ∩ {b² = 12a, a ≠ 0} = C` (`c_3|_(a=b²/12) = b²(3bc − 4)²/48`), and
   `Sing{c_3 = 0} = C` exactly (`∇c_3 = 0 = c_3` has the unique solution family `(4/(27c²), 4/(3c), c)`).
4. **(9) is a sketch, not "elementary" (status: "elementary arguments on top of THM-4620 and standard facts").** The conclusion
   length `≥ 3` is correct, but the argument is analytic/étale-local and uses, unlabelled: restriction to an open set (or étale
   neighbourhood) is exact and carries simple holonomic modules to simple ones or zero, so composition factors of `V_F|_B` lift;
   `j_+ O_(B∖S) = O_B[1/s]`; `j_+ L_(−1) = j_(!*) L_(−1)` (no eigenvalue 1); `O[1/s]/O ≅ H^1_([S])(O) ≅ i_+ O_S` is simple for smooth
   `S`; and the decomposition `F^(−1)(B) = B_1 ⊔ W` with `W → B∖S` finite étale and connected needs `B` small, `p ∈ S ∖ C`
   (`S` is smooth there, see 3.) and the fibre over `S ∩ B` a single point. Retype (9) "PROVED (sketch; standard D-module facts
   CITED)" and list them. "Delta module along the surface" is right at generic points of `S`; globally the constituent is the
   minimal extension of `i_+ O_(S°)`, `S° = S ∖ C`. Local length near a generic point of `S` is 4 (`O`, `O`, `i_+O_S`, `j_+L_(−1)`),
   consistent with the three global constituents `O`, `j_(!*)T` (whose restriction is `O ⊕ j_+L_(−1)`), and one supported on `S`.
5. **(1) proof, cosmetic:** "*`F_3` is a polynomial in `y` because ... `F_3'(y) = 1/c(y)` [has] poles only where `F_3` has poles*"
   is circular; polynomiality is the direct computation (denominator `27a²`). Delete the sentence.
6. **(4), precision:** `disc = −(2/(27a³))²·c_3`, so the condition is "`−c_3` is not a square in `K(a,b,c)`", which holds because
   `c_3` is irreducible (indeed absolutely irreducible: `disc_a(c_3) = 4(4 − 3bc)³`, odd exponent) and appears to the first power.
   Irreducibility of the cubic: it is satisfied by `y`, `K(F)(y) = K(x,y,z)` has degree `deg F = 3` (THM-1310(3)), so it is the
   minimal polynomial — state this chain (the file's "degree 3 = deg F hence irreducible" is correct but terse).
7. **(5)/(6), completeness:** `P^1_X` minus five points for `27ac² ≠ 4` (four when `aX³ − X + c` has a double root);
   `P^1_s` minus three points for `bc ≠ 4/3` (two when `3cs² − 4s + b` has a double root).
8. **Verified as stated (no change):** `det JF = −2`, degrees 7,6,4, `(F_2 − y)u = 3xF_1`; the `(F_1,F_2)` parametrization
   derived by elimination (`(F_2 − b)|_(F_1=a) = (xQ − (b − y))/u`), `u = 3a/Q`, the displayed `z`, `F_3|` cubic with
   `dF_3/dy = κQ`, `κ = −2/(9a²)`, `c_0`; `disc = −4c_3/(729a⁶)`, `c_3` to the first power ⇒ `S_3` and transposition local
   monodromy; `C`: `Q = (y − 2/(3t))²`, `F_3(2/(3t)) = t`; the `(F_1,F_3)` fibre: `ax³ − (1+w)(2+w)x + c(1+w)³ = 0`,
   `u = X/(aX³ − X + c)`, `F_2| = 2aX + 2/X − c/X²`, `dF_2/dX = 2(aX³ − X + c)/X³`, `X = x/(1+xy)` recovered; the `(F_2,F_3)`
   fibre: conic `bx² − (6+4w)x + 3c(1+w)² = 0`, `x = 2/(3cs² − 4s + b)`, `F_1| = s(cs² − 2s + b)/2`, `dF_1/ds = (3cs² − 4s + b)/2`,
   `u| = 2s/(3cs² − 4s + b)`; residues: `x²dF_3` at `r_1, r_2` = `−2/(r_1²(r_1 − r_2))`, `2/(r_2²(r_1 − r_2))` (also `x²y`, `x³`);
   `x dF_2` at `X = 0` = 2 (also `y`: `6a`, `z`: `18a²c`; `x²` via `res_ρ = 2ρ/(3aρ² − 1)`); `x²dF_1` at `σ` = `1/(3cσ − 2)` (also
   `x³`, `x²y`). Controls: `q = V_j(h)|_C` (for `h = x, x²z, x²y`) has all residues 0, as it must. The residue logic is sound: only
   the direction "`q ∈ V_jK[x]` ⇒ `q dF_j` exact on `C` ⇒ residues vanish" is used; `V_j|_C = d/dF_j` because `V_j` is tangent to
   `C`, `Der K(C)` is one-dimensional over `K(C)` and `V_jF_j = 1`. (7): `R_j = K[G_j]` argument correct (openness of the smooth
   `G_j`, `q | P` by density; base-change to `K̄` implicit); `z|` not quadratic in `F_3|` (no solution over `K(a,b)`), `x|` has
   poles, `y|` has degree 1. Kernel control reproduced independently mod `p = 1000003` (`rank_p ≤ rank_Q` pins
   `dim Ker_Q = 4` at `D = 12`; `= 3` at `D = 8, 10`), iterate degrees `1 + 7m, 1 + 8m, 1 + 10m` (`y`: two less) for `m ≤ 4`.
   `μ = 2, 4, 2` ✓. Gao's Theorem 3.4 confirmed from the PDF: `c_3`, `a_2 = ∂c_3/∂v_1`, `σ(v) = (v_2v_3, v_1v_3²)`,
   `v_3²c_3 = E∘σ`, sizes 3/1/0, "no fiber has exactly 2 points", image `C³ ∖ C`, Jelonek set `{c_3 = 0} = σ^(−1)(K)`;
   twisted form `F = ((q+γw)/(γx)², (p+2γ)/(γx), γx)` with `f_3 = γx`, `f_2f_3 = p + 2γ`, `f_1f_3² = q + γw`, i.e. `(X,Y) = (F_2F_3, F_1F_3²)` ✓.

## 3. HYP-9255 — corrections

1. **"equivalently" is one-directional.** "*no Keller map of A^3 has a coordinate-pair fibration whose generic fibre is a disjoint
   union of two or more lines*" **implies** HYP-9255 (surjective ⇒ `μ_j = 0` ⇒ lines, `deg F = d_j`; `d_j ≥ 2` excluded ⇒
   `d_j = 1` ⇒ automorphism), but HYP-9255 implies the fibration statement only if `μ_j = 0` forces `coker V_j = 0` for
   `d_j ≥ 2` (the file's own sub-question). Replace "equivalently" by "implied by (equivalent iff the sub-question below has a
   positive answer)" in the title, the Statement and the "Reduction" paragraph ("So the hypothesis is equivalent to" → "is implied by").
2. **Sub-question, sharpened:** by the proof of THM-4620(4), `μ_j = 0` with `d_j = 1` already gives an automorphism, hence
   `coker V_j = 0`; the torsion question is open only for `d_j ≥ 2`.
3. **Gao's map `G` is incompletely transcribed:** `w` is undefined. Gao §3.5: `u = 1 + xy`, `w = γu`, `γ = 2 − 4xy − x²z`,
   `p = w³ − 6w² + 6w`, `q = (3/8)w⁴ − 2w³ + (3/2)w²`, `G = (γx, (p(w)+2γ)/(γx), (q(w)+γw)/(γx)²)`; Theorem 3.5: polynomial,
   degrees 4, 11, 12, `det = 2`, generic fibre 4; fibre sizes `{4,3,2,1,0}` (size 3 over the hyperplane `v_1 = 0`, Theorem A.1).
4. **Status:** "PROVED in the plane (THM-4620(5))" → "KNOWN in the plane (Gurjar–Masuda–Miyanishi 2013, Remark 2.5; THM-4620(5) re-proves it)".

## 4. HYP-9257 — corrections

1. Keep OPEN. "PROVED: length ≥ 3" stands with the sketch caveat of §2.4. The local model `O ⊕ O[1/s] ⊕ j_+L_(−1)` near a
   generic point of `S` is correct (analytic/étale-local).
2. Add the verified fact `Sing{c_3 = 0} = C` (the cusp curve is exactly the singular locus of the Jelonek surface), which is
   what makes "generic point of `S`" = "point of `S ∖ C`" and supports the "transversal slice is cuspidal" picture (the latter
   is not verified here).
3. The Alexander-polynomial sentence is heuristic (correct arithmetic: `Δ_trefoil(−1) = 3 ≠ 0`), label it so.
4. Bavula's length bound: "explicit upper bound in terms of `deg σ`" confirmed from the abstract.

## 5. Session note (§4.1–4.2, 5, 6) and scripts

* §4.1: the bold "**in the plane, (f, g) is an automorphism iff D_f is surjective iff H^1 = 0**" must be labelled KNOWN
  (GMM 2013 Remark 2.5; Cerveau 1997). §4.2: "independent re-derivation of his Theorem 3.4" → add "and of THM-1310(3),(4),(6)".
  §5(ii): "reproduces Gao's `c_3` exactly" → "reproduces THM-1310's `L` = Gao's `c_3`". §6 Gao line: the quote
  "Claude Fable 5 assisted" is verified (title-page AI disclosure: "*The main idea and framework are due to the author, and
  Claude Fable 5 assisted in the proofs and in the writing up of the paper*"); attribution to Alpöge (July 19), Gallagher
  (July 20), Speyer (July 23) is reported, not adopted — consistent with MISTAKE-205. Bavula: add pages 731–738 / DOI.
* `thm1300_module_checks.py/.out`: the tangency check is mis-normalised and prints `False` (§2.2). Everything else in both
  scripts reproduces. Chini arXiv:2609.20184 and Zhang arXiv:2609.10180 exist with the described content (abstract pages).

## 6. Commands run

    cd 04-computation/experiments/logic4_20261012
    python3 audit_A_fibrations.py | tee audit_A_fibrations.out     # all assertions passed
    python3 audit_A_kernel.py     | tee audit_A_kernel.out         # dim Ker V_3 = 3,3,4 at D = 8,10,12; iterate degrees
    pdftotext arXiv:2608.00222 (Gao) and arXiv:1211.0744 (Gurjar–Masuda–Miyanishi), grep for Theorem 3.4 / 3.5 / Remark 2.5
    web: arXiv abstracts 2112.03177, 2609.20184, 2609.10180; numdam record CRMATH_2024__362_G7_731_0; Project Euclid Nagata–Nowicki 1988

## 7. Summary

All exact claims of THM-4621 are correct and were reproduced from scratch (parametrizations derived by elimination, restricted
functions, derivatives, discriminant `−4c_3/(729a⁶)`, curve `C`, residues, kernel dimension, iterate degrees), and the `a = 0`
stratum is now proved in-repo, closing the only citation gap in the stratification. Two things must change: the tangency
identity's prefactor (`−c/w³`, not `9c/2`; the session's own output says `False`), and the credit for `c_3`/`S_3`/Jelonek set,
which THM-1310 established in July. THM-4620's general-dimension dictionary is correct (one statement/proof mismatch in (4),
one base-change sentence, one typo), but its planar corollary (5) — "JC_2 iff every Keller derivation `D_f` is surjective" —
is a verbatim known result (Gurjar–Masuda–Miyanishi 2013, Remark 2.5, via Cerveau 1997), proved there by the same
`2g − 1 + n` mechanism; the novelty claim must be withdrawn and (5) retyped KNOWN. HYP-9255's "equivalently" is a
one-way implication; HYP-9257 stands as OPEN with its length-`≥ 3` input correct but sketched. No attribution-rule violation.
