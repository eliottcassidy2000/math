# Second independent audit: the parts written after the first audit (kind-pasteur-2026-11-11-S24; THM-5138, HYP-9668)

Audited: `README.md` at commit `c440c4a4b7` (the working tree differs only in the header lines 10-13): Proposition 2.4;
Theorem 4.1 as rewritten and Corollary 4.2; Lemma 5.2, Lemma 5.3, Theorem 5.4; Corollaries 5.5, 5.6, 5.7 and the table of
5.3; the `a = 0` statements of 6.1; "What is proved on the slice" in 6.2; section 6.3 (Proposition 6.3, added to the scope
by the coordinator during the audit); section 8, first bullet; the credit paragraphs; the application of the first
audit's corrections. Not audited: `THM-5138-...md`, `HYP-9668-...md` (being rewritten), and the numerical work already
done by the first audit. Date: 2026-11-11. Every program in `audit/audit2/` is mine (`eng2.py` is a third engine,
written from scratch: charge basis `x^c q^n`, `y^c q^n` with `q = 1 + xy`, filtration `|c| + 2|n|`, prime `2^61 - 1`;
exact series arithmetic for the places at infinity). The author's two new scripts were read for bugs and re-run only
to check that the stored outputs reproduce (they do, byte for byte). Nothing outside `audit/AUDIT-2.md` and
`audit/audit2/` was written.

## 0. Verdict

**PASS WITH CORRECTIONS.** I tried to break every new statement, by re-deriving each proof and by independent
computation, and found no mathematical error in Proposition 2.4, Theorem 4.1, Corollary 4.2, Lemma 5.2, Lemma 5.3,
Theorem 5.4, Corollaries 5.5-5.7 (apart from a hypothesis missing in the last sentence of 5.6, m1), Proposition 6.3,
or the identities of section 8. One MAJOR finding concerns the
status of the Image Conjecture (the note treats its failure in characteristic zero as possibly unrecorded and quotes a
survey of 2010; canon has refuted it since July 2026, in the theorem the note itself cites). Sixteen MINOR findings
(one false sentence by a missing hypothesis, labels, credit, stale counts). Five strengthenings came out of the attempts
to break the statements; they are proved in section 3.

| statement | verdict |
|---|---|
| Proposition 2.4 (`M_z = (E_a - a) B + (E_b - b) B`, the Weyl relations, the heat operator for the standard pair) | PASS |
| Theorem 4.1 (trace argument; `S_z = 0` or isomorphism) | PASS |
| Corollary 4.2 (i), (ii) | PASS |
| Lemma 5.2 (valuation at infinity) | PASS (m12: one expository step) |
| Lemma 5.3 (exact shift `w(P(E_a, E_b) phi) = w(phi) - n`; torsion-freeness) | PASS |
| Theorem 5.4 (injectivity of the Stein transform) | PASS (S4: the exceptional set is in fact a proper algebraic subset, CITED) |
| Corollary 5.5 | PASS (m8: label; S3: a sharper reading) |
| Corollary 5.6 | PASS WITH CORRECTION (m1: last sentence false without "`h` in `B_0`") |
| Corollary 5.7, items 1-3 | PASS; "does not use the comparison theorem": true |
| Corollary 5.7, item 4 | PASS (m7: one more CITED fact is used) |
| table of 5.3 | PASS (m6, m15) |
| 6.1, `S_z = 0` at the centres `(0, b)` | PASS |
| 6.2, "What is proved on the slice" (i)-(iii) | PASS; (ii) holds for all three generators and both directions (S1) |
| Proposition 6.3 (dimension 2, `[1]_z != 0`, SN and CV hold, finite chaos fails) | PASS; every line of the proof checked, machine output reproduced with my own code |
| 6.3, "for a Keller pair on `C^2` this cannot happen" | PASS for the geometric configuration only (m6) |
| section 8, first bullet: the identity, the functional, the claim for `x^p y^q` | PASS |
| section 8, first bullet: what it says about Zhao's Conjecture 1.3 | the example does contradict the conjecture as stated; the account of what is known is WRONG BY OMISSION (M1) |
| credit to Session269 | accurate and fair (m9: two smaller credits missing) |
| first audit's M1-M6, m1-m16 | applied in the README, with the residues m3, m5, m13, m14 |

## 1. MAJOR finding

### M1. The Image Conjecture is already refuted in characteristic zero, in canon and in public sources; the README says the opposite by omission

*Location.* README section 8, first bullet, lines 576-581 ("so this is also an instance where the statement of
Zhao's Conjecture 1.3, read for an arbitrary potential, fails in characteristic zero. I do not claim this is new: the
one search made (first audit) found no record of it, and van den Essen's survey (arXiv:1006.5801) reports that no
counterexample in characteristic zero was known in 2010; it needs its own sourced note before anyone relies on it");
Corollary 5.7, line 400 ("the Image Conjecture itself fails for some other potentials (section 8)"); section 1, line 63
("Conjecture 1.3 (Image Conjecture: ...)", no status); section 0 item 4 and section 5.1 ("Zhao's Image Conjecture").
The same gap is in the first audit (`AUDIT.md`, M1, "Side remark").

*Finding.* (a) I read the primary source (arXiv:0902.0210v2, text extracted with pdftotext). Conjecture 1.3 is stated
for any commutative algebra `A` over a field of characteristic zero and any commuting set `C` of operators
`h(z) + sum c_i d_i` with `h` in `A[z]`, `c_i` in `A`, with `Im C` the sum of the images: `f^m` in `Im C` for all
`m >= 1` implies `f^m g` in `Im C` for `m` large. With `A = C`, `n = 2`, `C = {d/dx - 2xy, d/dy - x^2}`, `f = 1`,
`g = x` the hypothesis holds and the conclusion fails, so the example **does contradict Conjecture 1.3 as stated**
(it is also of the reduced form of Zhao's Theorem 2.5, and in the finite-codimension case of his Corollary 2.11).
The mechanism (`1` in the image, image proper) is the one Zhao uses himself in characteristic `p` (Example 2.7), and
van den Essen calls it the 1-property. So far the README is right.

(b) What the README then says about the state of knowledge is stale and misses canon. The special case of
Conjecture 1.3 with `A = C[xi_1, ..., xi_n]` and the operators `d/dz_i - xi_i` (the Special Image Conjecture `SIC(n)`;
an instance of the form with a regular sequence to which van den Essen-Wright-Zhao restricted the name "Image
Conjecture" in 2011, Conjecture 2.1 of arXiv:1008.3962) is **false for every `n >= 2`: THM-2801** (PROVED, VERIFIED-EXACT,
independently audited; in canon since 2026-07-28; `ACTIVE-GUARDRAILS.md` item 58; `CURRENT-FRONTIER.md` line 92;
`CORE-PAPERS.md`, entry Derksen-van den Essen-Zhao). The README cites THM-2801 in section 1, for another purpose.
Public sources of 2026 say the same (found by web search today; I read summaries only and did not check the
mathematics): the MathWorld entry "Image Conjecture" lists it among refuted conjectures and cites Z. Zhang, "Direct
Consequences of the Three-Dimensional Counterexample to the Jacobian Conjecture" (20 July 2026); A. Dvorsky,
arXiv:2608.07338 (7 August 2026), gives an explicit counterexample in five variables to the Generalized Vanishing
Conjecture, which the Image Conjecture implies. Van den Essen's sentence of 2010 (arXiv:1006.5801, section 4,
"General Image Conjecture", Comment 1: "no counterexamples to this conjecture are known") is quoted correctly but is
sixteen years old.

(c) So a reader of section 8 is led to think that a characteristic-zero counterexample to the Image Conjecture may be
unrecorded, and a reader of section 0 item 4 that "K-GMC is an instance of the Image Conjecture" lends it support.
Neither is the case. What the example `x^2 y` does add: it is a failure in the smallest setting (a field of
coefficients, as many operators as variables, quotient of finite dimension), and of the trivial kind. An observation,
not a repaired conjecture: `W_x = 2xy` and `W_y = x^2` vanish together on a line, while for `W = (f - a)(g - b)` of a
Keller pair the critical points are the finitely many preimages of the centre. The counterexamples of THM-2801 are of
another kind (`n` operators in `2n` variables, infinite codimension, the regular sequence `xi_1, ..., xi_n`). I found
no record of the instance `x^2 y` (or of monomial potentials) in the five papers read (section 6); that search does
not cover the 2021 book of van den Essen, Kuroda and Crachiola.

*Correction.* Section 8, replace the two sentences by: "The two operators commute and have order one with constant
leading coefficients, so the statement of Zhao's Conjecture 1.3 fails for this potential (`A = C`, `f = 1`, `g = x`).
That the conjecture fails in characteristic zero is KNOWN: its special case `SIC(n)` is false for every `n >= 2`
(THM-2801; guardrail 58), and public sources of 2026 record the failure (MathWorld, "Image Conjecture", citing Zhang
2026; Dvorsky, arXiv:2608.07338, refutes the Generalized Vanishing Conjecture, which the Image Conjecture implies;
neither re-checked here); van den Essen's survey of 2010 still had no counterexample. What
`x^2 y` adds is a failure over a field with two operators in two variables, by the 1-property (Zhao's Example 2.7 in
characteristic `p`); there the critical set of `W` is a line, while for `W = (f - a)(g - b)` it is finite. No record
of this instance was found in Zhao 2010 (both papers), van den Essen 2010, van den Essen-Wright-Zhao 2011 (both
papers); the search is not exhaustive." Line 400: "Statement 2 is the instance `q = fg` of Zhao's Image Conjecture, which is false
in general in characteristic zero (THM-2801: `SIC(n)`, `n >= 2`; section 8: `x^2 y`), so statement 2 inherits no
support from it." Line 63: add "(REFUTED in general: THM-2801)". Log the missed canon theorem with MISTAKE-764.

## 2. MINOR findings

- **m1. Corollary 5.6, last sentence, is false as stated.** Line 384: "More precisely, if `d_u^m h != 0` for all `m`
  then `[d_u^m h]_z != 0` for all `m`." The proof says "`h` in `B_0`"; the statement does not. Without it the sentence
  fails already at `m = 0`: on the double cover, `h = x^(-2) - 1/a` is nonzero, `d_u^m h != 0` for all `m`, and
  `lambda_+(h) = lambda_-(h) = 0`, so `[h]_z = 0` at every centre with `a != 0` (`b3`, part 4). *Correction:* "More
  precisely, for `h` in `B_0` (defined over `k_0`): if `d_u^m h != 0` for all `m`, then `[d_u^m h]_z != 0` for all `m`."
- **m2. K-GMC is an instance of Zhao's conjecture only for `X = C^2`.** Line 272-273 defines K-GMC "for an etale `F`"
  and says "This is the instance `q = (u - a)(v - b)` of Zhao's Image Conjecture". Zhao's conjecture is about
  polynomial algebras `A[z]`; for the surfaces of section 6 the ring is not one, and Zhao notes himself that the
  straightforward generalisation to Laurent polynomials fails (arXiv:0902.0210, end of 1.1; arXiv:0902.0212,
  Proposition 2.6). *Correction:* "For a Keller pair (`B = C[x, y]`) this is the instance ...; on other surfaces it is a
  generalisation, which fails (Proposition 6.1)."
- **m3. Section 0 item 6 and the title of section 6 blend labels and are out of date.** Line 49: "Examples
  (FINITE-EXACT, sections 6-7). On the two planar non-proper etale maps available ... finite chaos fails, and SN fails
  at every centre in general position." Section 7 is HEURISTIC and contains no example; there are now three maps
  (6.3); "finite chaos fails" and "SN fails at every centre in general position" are PROVED (Theorem 4.1, Corollary
  5.6), Propositions 6.1 and 6.3 are PROVED by hand, only 6.2 is FINITE-EXACT; and the item does not say that at
  special centres SN can hold (6.3). Line 428 labels all of section 6 FINITE-EXACT. *Correction:* "6. **Examples
  (section 6).** On three planar non-proper etale maps (a double cover on `C^* x C`; a slice of the `C^3` Keller map;
  the inclusion of `C^2` minus a hyperbola) the count is confirmed (dimensions 2 and 3: PROVED by hand and by machine;
  11 and the linear potentials 2, 4, 6: FINITE-EXACT), finite chaos fails wherever `S_z != 0` (PROVED), and SN fails at
  every centre in general position (PROVED); at special centres SN and CV can hold (6.3)." Section 6 title: drop the
  label, or "(6.1, 6.3: PROVED by hand; 6.2: FINITE-EXACT)".
- **m4. Section 0 item 3 hides a citation.** Line 32-33: "For a Keller pair `S_z != 0`, so finite chaos at any single
  centre forces an automorphism." `S_z != 0` at every centre is the CITED count (Corollary 4.2 (i) says so).
  *Correction:* "For a Keller pair `S_z != 0` at every centre (CITED count; without citation at centres in general
  position, Theorem 5.4), so ...".
- **m5. The audit numbers quoted after Proposition 2.2 are wrong.** Line 149-150: "on the slice of section 6 (15 and
  21 cases, 0 failures; the audit: 21 and 27)". The first audit's 21 for (b) is on the double cover (`i + j <= 5`); on
  the slice it checked 15 of 15 for (b) and 27 for (d) (`AUDIT.md`, section 3). *Correction:* "(15 and 21 cases, 0
  failures; the audit: 15 and 27, and 21 of 21 for (b) on the double cover)".
- **m6. 6.3, last two sentences: scope and completeness.** Line 540-543: "The centre (0, 0) is special because the
  missed curve `uv = -1` is itself a member of the pencil ... For a Keller pair on `C^2` this cannot happen (...)".
  (i) THM-4122 (4) does justify that for a Keller pair no component of the non-properness curve is a member of the
  pencil or is missed by it, so the count of Theorem 3.1 has no special centre on `C^2`. It does not justify "SN and
  CV cannot hold at a special centre of a Keller pair that is not an automorphism", which is how "this" can be read and
  which is OPEN (table of 5.3). (ii) The proof of Proposition 6.3 uses the grading by the torus
  `(x, y) -> (lambda x, y / lambda)`, which preserves both `W = xy` and the missed curve; say so. (iii) (0, 0) is not
  the only special centre of this map: on the two lines `ab = 0` the dimension is 2 and one-sided Stein nilpotence
  holds (S2). *Correction:* "The centre (0, 0) is special: `W = xy` and the missed curve `uv = -1` are invariant under
  the torus `(x, y) -> (lambda x, y / lambda)`, which grades `M_z`, and the missed curve is a member of the pencil
  `uv = t`, so the count drops. All centres on the lines `ab = 0` are special for this map (dimension 2, `SN_u` on
  `b = 0`, `SN_v` on `a = 0`: second audit). For a Keller pair on `C^2` the count has no special centre (Theorem 3.1;
  THM-4122 (4)); whether SN or CV can hold at some centre of a Keller pair that is not an automorphism is OPEN."
- **m7. Corollary 5.7, item 4, uses one more cited fact.** Line 397-398 lists "the CITED count and the CITED fact that
  a holonomic module is a vector bundle over a dense open set". It also needs that `B` is holonomic over the Weyl
  algebra generated by `E_a, E_b, d_u, d_v` (Bavula's theorem for the twisted module, CITED in THM-4620 (6), and the
  linear automorphism of section 2); with it, the generic rank is `dim S_z` for generic `z` although `B` is not
  finitely generated over `C[E_a, E_b]`: the fibre of `B[1/delta]` at a point of `{delta != 0}` is `B / m_z B = S_z`.
  *Correction:* add "the CITED holonomicity of the twisted module (Bavula; THM-4620 (6))".
- **m8. Labels of Corollaries 5.5, 5.6, 5.7.** "PROVED" without source, while `1 => 2` is THM-2022 and `6 => 1` is the
  KNOWN criterion; Proposition 2.3 is labelled "PROVED from THM-2022". The proofs say it; the labels should.
  *Correction:* "PROVED from THM-2022 and the KNOWN criterion THM-4620 (5)".
- **m9. Two credits.** (i) The reduction of one-direction detection to the torsion-freeness of `B` over
  `C[E_a, E_b]` is the first audit's (`AUDIT.md`, section 7, (iii): "if `K' = 0`, `d_u` is locally nilpotent and the
  pair is an automorphism", verified there on the double cover and for open immersions, left without proof in
  general). MISTAKE-764 says so; the README credits the audit only for the module `G = B` and for Theorem 4.1. Add
  to 5.3: "the route (first hits in one direction lie in the torsion of `B`) was sketched by the first audit;
  Lemmas 5.2-5.3 supply the missing torsion-freeness". (ii) Testing at one centre algebraically independent over the
  field of coefficients is Session269's Theorem 2a (two directions, every monomial, conditional); the README credits
  the rank step inside the proof of Theorem 5.4 but not the formulation in Corollary 5.5. One clause there suffices.
- **m10. Section 8: "`M` is proper because `H^2` has dimension `gcd(p, q)`"** rests on the CITED comparison theorem
  inside an otherwise elementary remark. It is elementary (S5): a linear form `m(i, j)` on monomials kills `M` iff it
  is supported on `q (i + 1) = p (j + 1)` and satisfies `m_(n+g) = (n / g) m_n` along that line (`g = gcd(p, q)`), so
  the codimension is `g` and `1` lies in `M` iff `p != q`. Machine: `1 <= p, q <= 4` (`b7`).
- **m11. Section 8: "(`x + x^2 y`: the functional takes the value -1 on 1)"** depends on a normalisation that is not
  stated (value 1 on `x`). Confirmed with my engine: `[1] != 0` and `[x] = -[1]` (`b8`). *Correction:* "the class of 1
  is not zero (`[x] = -[1]`)".
- **m12. Lemma 5.2, last step.** "on `K` the extensions of `d/du, d/dv` are `d_u = ...`": the derivation
  `-tau D~_1 - s tau D~_2` lives on the completion, `d_u` on `K`. The step is: two derivations of `K` with values in
  the field `K_w` that agree on `k` are equal, `K / k` being finite separable. One clause. No hidden assumption was
  found: the residue extension need not be trivial (on the slice `f = 2`), `K / k` need not be Galois, `w` need not be
  unique (the slice has two places), and `c` need not be constant (`b1`, `b1c`).
- **m13. Residue of the first audit's m12.** Line 7-8: "predictions S2-S6 all CONFIRMED". S4 is confirmed relative to
  the truncation (the pre-registration file now says so in its audit notes). *Correction:* "(S4 relative to the
  truncation)".
- **m14. Residue of the retracted framing in the title.** Line 1: "Gaussian integration by parts between GMC(2) and
  JC(2)", while section 9 says that nothing here lies strictly between them. *Correction:* "... in the coordinates of a
  Keller pair: GMC(2) and JC(2) in one ring".
- **m15. Table of 5.3, third row.** Correct, but a reader may take 6.3 as an example where K-GMC holds. K-GMC at the
  centre (0, 0) of 6.3 is not decided (section 4, J). Add "(K-GMC there: not decided; with `S_z = 0` it holds
  trivially)". The row can also cite the lines `ab = 0` of S2 for the one-direction properties.
- **m16. Section 5.2, "Checks".** "every monomial tested ... has a nonzero class": on the slice a nonzero normal form
  is non-membership relative to the truncation (the README's own convention, section 6). Add "(relative to the
  truncation on the slice)".

## 3. Strengthenings found while trying to break the statements

- **S1 (PROVED; answers "which generator" in 6.2 (ii)).** With `r = x w` the slice is `C^2_(r,u)` minus the curve
  `Omega = {om = 0}`, `om = 1 - r u + 3 r^2 = w`, and `v = 2 r^3 - u r^2 + 2 r`: the slice map is the restriction of the
  finite map `(r, u) -> (u, v)` to the complement of its critical curve, whose image is the quartic of 6.2 (its
  discriminant in `r` is `-4 (27 v^2 + u^3 v - 18 u v - u^2 + 16)`). Here `x = r / om`, `y = u - 3r`,
  `d_v = (1 / (2 om)) d/dr`, `d_u = d/du + (r^2 / (2 om)) d/dr`. *Lemma.* If `phi = psi / om^n` with `n >= 1` and `om`
  not dividing `psi`, then `om^(n+2) d_u phi = -(n/2) psi r (r u - 2)` and `om^(n+2) d_v phi = -(n/2) psi (6r - u)`
  modulo `om`; on `Omega`, `r (r u - 2) = r (3 r^2 - 1)` and `6r - u = (3 r^2 - 1) / r` are not zero and `om` is
  irreducible; so `d_u` and `d_v` raise the order of pole along `Omega` by exactly 2. Since `x` has a simple pole and
  `d_u y = 1 - 3 r^2 / (2 om)`, `d_v y = -3 / (2 om)`, `d_u w = r (r u - 2) / (2 om)`, `d_v w = (6r - u) / (2 om)` have
  simple poles: **`d_u^m h` and `d_v^m h` have a pole of order `2m + 1` (`h = x`) or `2m - 1` (`h = y, w`) along
  `Omega` and are never zero.** By Corollary 5.6, at every centre algebraically independent over `Q`, all classes
  `[d_u^m h]_z`, `[d_v^m h]_z` (`h = x, y, w`; all `m`) are nonzero. (`b5`: lemma symbolic; pole orders for `m <= 6`.)
  Replace "for at least one of the generators" by "for each of the generators".
- **S2 (PROVED by hand, dimension FINITE-EXACT; a new family of special centres for 6.3).** For the inclusion of
  `C^2` minus `{1 + xy = 0}` and a centre `(a, 0)`, `a != 0`: `T_u (x psi) = (d/dt - 1)(t psi)` (charge 0) and
  `T_v (y phi) = (d/dt - 1)(t phi) + a y phi`, so `y R` lies in `M_z`; then
  `T_v (y^(j+1) phi) = y^j [(j + 1) phi + t phi' - t phi] + a y^(j+1) phi` gives `y^j R` in `M_z` for all `j >= 1` by
  induction. So `M_z` contains every element of negative charge, `d_u^m h` lies in `M_z` as soon as `m` exceeds the
  charges of `h`: **`SN_u` holds at every centre of the line `b = 0`** (and `[h U^k V^(m+k)]_z = 0` for `m` large),
  although `d_u = d/dx` is not locally nilpotent; symmetrically `SN_v` on `a = 0`. Machine (`b6b`, `b9`): dimension 2
  at (0, 3), (0, -7), (2, 0), (5, 0) (Remark 3.3: `2 - 1 + 1`), 3 at (1, 1), (2, 3) and at three pseudo-random
  centres; the classes vanish exactly as the charges predict, and in the other direction nothing vanishes beyond
  what the charge forces (`m <= 7`). So on surfaces the exceptional set of Corollary 5.6 can be a whole curve, with
  `[1]_z != 0` (machine).
- **S3 (PROVED from THM-4620 (1); a sharper reading of Corollary 5.5).** For a Keller pair the union of the kernels of
  the powers of `d_u` is `C[g][f] = C[f, g]` (THM-4620 (1)). With Theorem 5.4: at a centre algebraically independent
  over `k_0`, a polynomial `h` in `k_0[x, y]` has a first hit in the direction `d_u` iff `h` lies in `k_0[f, g]`. So
  item 5 of Corollary 5.5 says "`x` and `y` are polynomials in `f, g`", and the implication to item 1 is immediate.
  This is the plainest form of "the ladder is not intermediate".
- **S4 (PROVED from CITED holonomicity; sharpens the last sentence of Theorem 5.4).** `B` is holonomic over the Weyl
  algebra of the plane of centres, so there is `delta != 0` in `C[a, b]` with `B[1/delta]` a vector bundle on
  `{delta != 0}`; a nonzero `h` is a nonzero section (Lemma 5.3), and it vanishes in the fibres only on a proper
  closed subset. So for each `h != 0` the set of centres with `[h]_z = 0` lies in a **proper algebraic subset** of
  `C^2` (the curve `delta = 0`, which depends only on `F`, and a proper closed subset depending on `h`), not only in a
  countable union. The countable union is still needed in Corollary 5.5 (all `r`).
- **S5 (PROVED).** The elementary codimension count for `x^p y^q` of m10.

## 4. Proof checks

**A. Lemma 5.2.** Every step holds. `nu` is the `tau`-adic valuation of `C(s)(tau)` (`p = tau^(-n) (p_n(1, s) + ...)`,
`p_n(1, s) != 0`); the formulas `d/dv = tau D_2`, `d/du = -tau D_1 - s tau D_2` are right (check on `tau`, `s`); the
completion at `w` is a finite extension of `C(s)((tau))` of degree `e f`, its residue field lifts uniquely (Hensel,
residue characteristic zero), the extension over `k'((tau))` is totally and tamely ramified, and a uniformizer with
`pi^e = c tau`, `c` in `k'`, exists (write `Pi^e = tau c_0 (1 + m)` and divide `Pi` by an `e`-th root of the 1-unit);
`D~_1 = (1/e) pi d/dpi` and `D~_2` (the derivation with `D~_2 | k' = '` and `D~_2 pi = (1/e)(c'/c) pi`) are continuous
derivations that preserve the valuation, with `D~_1 tau = tau`, `D~_2 tau = 0`; they extend `D_1, D_2`; uniqueness of
extension gives the restriction to `K` (m12). Conceptually: `d/du`, `d/dv` are `tau` times logarithmic vector fields
along the line at infinity, and logarithmic fields lift to logarithmic fields. Tests (`b1`, `b1b`, `b1c`): the double
cover (`e = 2`, monomial valuation `w(x) = -1/2`, `w(y) = -3/2`, which restricts to `-deg` on `C[u, v]`); the double
cover `(u, eta^2)` (`c = 1/s`); an automorphism; the slice, where the places over the line at infinity are
`w_1`: `r = rho + (rho^2 + 1) tau + ...`, `rho^2 = -s` (`e = 1`, `f = 2`) and `w_2`: `r = 1/(2 tau) + ...`
(`e = f = 1`), with `w_1(x, y, w) = (1, -1, -1)`, `w_2(x, y, w) = (1, -1, -2)`; and a place with `e = 2`, `f = 2`,
`c = 1/s` on `zeta^2 = u eta^2 + 1` over `(u, eta^2)`, where the valuation is not monomial. No violation in 1202
inequalities; equality is frequent, so the bound 1 is sharp; the formulas of the proof hold on the series (662
comparisons, 0 mismatches); the growth exponents along rays `t (1, s_0)`, `t = 10^8, 10^12`, computed in the README's
coordinates on the three sheets, drop by at least 1.0000.

**B. Lemma 5.3.** Correct: a word with `k >= 1` derivative letters shifts `w` by at least `-(i + j) + 2k`, the word
`u^i v^j` by exactly `-(i + j)`, and `w(P_n(u, v)) = -n` because `P_n(1, s)` is a nonzero element of the residue
field; non-commutativity is harmless because each word is estimated separately. Brute force (`b2`): the rank of
`h -> P(E_a, E_b) h` on boxes of 961 (double cover), 729 (`C^* x C`), 841 (`C^* x C^*`), 703 (plane minus hyperbola,
inclusion) and 507 (slice) basis monomials is full for 18 polynomials each (`E_a - a_0` with `a_0 = 0, 1, 2, 12345`;
`E_b`, `E_b - 3`, `E_b + 7/2`; `E_a E_b - 1`; `E_a E_b`; `E_a^2 + E_b`; `E_a - E_b`; three more; four random
quadratics); a full rank modulo a prime is a full rank over `Q`. No torsion. On the double cover the exact shift
`w(P(E) phi) = w(phi) - deg P` holds in 720 random cases.

**C. Theorem 5.4.** Correct. The linear system has entries affine in `(a, b)` over `k_0`; ranks descend from `C` to
`k_0(a, b)`, which is a rational function field because `a, b` are algebraically independent; clearing denominators is
legitimate; `Phi` is well defined on `C[a, b] (x) B` (a basis of monomials `a^i b^j`), kills both kinds of generators
because `E_a E_b = E_b E_a` and nothing else, and sends `P h` to `P(E_a, E_b) h`; the countable union is the union of
the curves `P = 0`, `P` in `k_0[a, b]`. Independent tests (`b3`): (1) on the double cover
`lambda_+-(x^k y^l) = (+-1)^(k+l) (1/2) sum_j C(l, j) ff((k+l)/2, j) b^(l-j) alpha^(k+l-2j)` (`alpha^2 = a`; against the
residue for `|k| <= 4`, `l <= 3`; kills the generators for `|k| <= 6`, `l <= 5`, symbolic centre); the leading terms
`b^l alpha^(k+l) / 2` are pairwise distinct, so `lambda_+` alone is injective on `C[x, 1/x, y]`: this proves the
theorem on the double cover without Lemmas 5.2-5.3. Zero loci: `[y]_z = 0` iff `2ab + 1 = 0`; `[xy - 1]_z = 0` iff
`b = 0`; `[x^2 - y^2]_z = 0` iff `a - a b^2 - 2b = 0` (machine membership at seven centres agrees): curves over `Q`, as
the theorem requires. (2) Standard pair: `h = (e^(d_x d_y) h)(E_a, E_b) 1` (25 random polynomials), so the zero locus
of `[h]_z` is the curve `e^(d_x d_y) h = 0`. (3) Open immersion `C^* x C -> C^2`: `L(h) = (e^(d_x d_y) h)(a, b)` kills
`M_z`, `L(1/x) = 1/a`, and `S_z = 0` exactly on the line `a = 0`. I could not construct a counterexample; by the proof
there is none.

**D. Corollaries 5.5-5.7.** 5.5: `1 => 2 => 3 => 4 => 5` is Proposition 5.1 (2.2 (b) holds for every etale `F`);
`5 => 6`: `d_u = {., g} / {f, g}` has coefficients in `k_0`, Theorem 5.4 applies, and the elements killed by a power of
a derivation form a subalgebra; `6 => 1`: THM-4620 (5) gives `Ker D_g = C[g]` (KNOWN: Cheng-McKay-Wang; Nagata-Nowicki)
and the slice theorem gives `C[x, y] = C[g][f]`; this is the last step of the proof of Theorem 4.3 in
arXiv:1004.0521. The version for a set of centres is right (the set contains a centre algebraically independent over a
countable field). 5.6: Theorem 5.4 on `B_0`, linearity for `B = C (x) B_0`; the slice theorem twice
(`d_v` preserves `Ker d_u`, `v` is a slice there) and `Ker d_u cap Ker d_v = C`; "K-GMC at `z` iff isomorphism" uses
`[1]_z != 0` implicitly (Theorem 5.4) and THM-2022 for the converse. 5.7: `2 <=> 1` by translation of the pair and a
centre in general position; item 3 the same; item 4 with m7. "Does not use the comparison theorem": true for
`1 <=> 2-3`.

**E. Theorem 4.1, Corollary 4.2.** The module `P` (six commutators), surjectivity of `Theta` (induction on the top
degree; the constant term gives `S' = S_z`), the extension of `tau` to a `k`-linear form on `K`, nondegeneracy of the
trace form, `Tr` commuting with the derivations (unique extension to the Galois closure), `kappa` constant
(`Ker d_u cap Ker d_v = C`), `kappa = l([1]) / d`, integrality by Newton's identities, and "finite etale connected
cover of `C^2`" are all right. 4.2 (ii) needs only Theorem 5.4.

**F. Sections 6.1, 6.2.** `S_z = 0` at `(0, b)`: `x^(i+1) y^j = j x^i y^(j-1)` modulo `M'` for all `i` in `Z` kills the
powers of `x`, then everything. 6.2 (ii): "a locally nilpotent derivation of a domain kills the units" is right, and
"at least one of `x, y, w`" follows (the elements killed by a power of `d_u` form a subalgebra; if it contained
`x, y, w` it would be `B`). A single unit can be killed by a power of a derivation that is not locally nilpotent
(`t` for `d/dt` on `C[t, 1/t]`), so the README's weak form is the right inference from that argument; S1 gives the
strong form.

**G. Section 8.** `2 (d/dy - W_y)(y) - (d/dx - W_x)(x) = 1`; the functional `x^k y^l -> l! [k = 2l + 1]` kills both
families of generators (`k < 12`, `l < 8`) and is 1 on `x`; the Euler identity for `x^p y^q`; `1` in `M` iff `p != q`
and codimension `gcd(p, q)` (`p, q <= 4`). See M1 for the conjecture.

**H. Proposition 6.3 (coordinator's questions).** (1) The two formulas for `k >= 1`, the case `k <= -1`, the charge-0
identity `M_z cap R = (d/dt - 1)(t R)`, the codimension (`R / (d/dt - 1) R` is spanned by `w`; `1 = (d/dt - 1)(-1)`;
with `q = 1 + t` the generators are `g_n = -q^(n+1) + (n + 2) q^n - n q^(n-1)`, and `<1, 1/q>` is a complement), and
`1` not in `M_z` are right (`b6`, A and B). The two functionals are `lam_G(phi) = int_0^inf phi e^(-t) dt` (the
Gaussian expectation, which exists on `B` because `1 + |z|^2 > 0`) and `lam_res(phi) = Res_(t = -1) phi e^(-t) dt`.
(2) Own code: `c(N) = 2` for `N <= 20` at `N' = 28, 36`; classes of `d_u^m h`, `d_v^m h` (seven elements, `m <= 7`)
zero exactly when the charge is not zero; `(d_u d_v)^m w` nonzero for `m <= 4`; the cone classes vanish for `m` above
the charge. (3) Centre (2, 3): 3, nothing vanishes. (4) Consistent with Theorem 5.4 and Corollary 5.6: (0, 0) is not
in general position, the exceptional sets found (`a = 0`, `b = 0`) are curves over `Q`, `d_u = d/dx` is not locally
nilpotent (`d_x^m w = (-1)^m m! y^m w^(m+1)`), and at three pseudo-random centres no basis monomial has a zero
class. (6) See m6.

**I. First audit's corrections.** Checked one by one against the README: M1-M6 are applied; of the minors, m1-m3,
m5-m11, m14 and m16 are applied in the README, m4, m13 and m15 concern the canon files (not audited; the README has
the condition `a != 0` of m4), and m12 is handled by audit notes appended to the pre-registration. Residues and one
new slip made while applying them: m3, m5, m13, m14 of this audit. No sentence of the README contradicts another,
apart from the counts in m3 and the title (m14). Status labels: m3, m4, m8.

**J. K-GMC at the centre (0, 0) of 6.3: not decided.** `M_z` is the sum of all `B_k`, `k != 0`, and of
`M_0 = Ker lam_G cap Ker lam_res` in charge 0. Elements with all charges of one strict sign are null and harmless
(`h P^m` has one-signed charges for `m` large). Proved here: (i) an element `phi` of `R` with a pole at `t = -1`
cannot have `lam_res(phi^m) = 0` for all large `m` (the generating function is `-z d/dz` of a sum of exponentials
`e^(-s_i(z))` over the small roots of `s^n = z N(s)`; exponentials of non-constant algebraic functions with
non-constant differences are linearly independent over algebraic functions), and a polynomial `phi` with
`lam_G(phi^m) = 0` for all `m` is zero (the one-variable factorial conjecture: van den Essen-Wright-Zhao,
arXiv:1008.3962, Theorem 4.9); so the only null element of charge 0 is 0, and a null element with charges `>= 0` has
no charge-0 part; (ii) a two-term element `x^l a(t) + y^k c(t)` is null only if `a c = 0`. What remains is whether a
null element with charges of both signs (necessarily with at least three charge components) exists and, if it does,
whether `h P^m` returns to `M_z`; this is the analogue of NC2 (THM-2022) for a ring with a pole. I did not prove or
refute it, so K-GMC at this centre is open in both directions.

## 5. What was recomputed

| quantity | README / author | second audit |
|---|---|---|
| `w(d_u phi), w(d_v phi) >= w(phi) + 1`, double cover | proof only | 230 inequalities, 0 violations (214 equalities) |
| the same, `(u, eta^2)` (`c = 1/s`) | not considered | 230, 0 violations; formula `d_v = tau D~_2` exact |
| the same, automorphism `(x + y^2, y + (x + y^2)^3)` | proof only | 80, 0 violations |
| places of the slice over the line at infinity | not computed | `w_1` (`e = 1`, `f = 2`), `w_2` (`e = f = 1`); valuations of `x, y, w`: (1, -1, -1), (1, -1, -2) |
| the same inequalities on the slice, exact series | proof only | 424 on 108 elements, 0 violations; formulas of the proof: 0 mismatches of 424 |
| growth exponents along rays, three sheets | not computed | 9 elements, smallest drop 1.0000 |
| a place with `e = 2`, `f = 2`, `c = 1/s` | not considered | 238, 0 violations; formulas: 0 mismatches of 238 |
| torsion of `B` over `C[E_a, E_b]` | Lemma 5.3 | rank full for 18 polynomials on boxes of 961, 729, 841, 703, 507 monomials (five maps) |
| `[E_a, E_b] = 0` | 7 elements (slice) | samples on all five boxes, 0 failures |
| exact shift `w(P(E) phi) = w(phi) - deg P`, double cover | Lemma 5.3 | 720 cases, 0 failures |
| `lambda_+-` in closed form; injectivity of `lambda_+` | not stated | 72 residues, 312 generators, 119 leading terms |
| zero loci of `[y]`, `[xy - 1]`, `[x^2 - y^2]` on the double cover | not computed | `2ab + 1`, `b`, `a - a b^2 - 2b`; engine agrees at 7 centres |
| standard pair, `h = (e^(d_x d_y) h)(E_a, E_b) 1` | stated | 25 of 25 |
| open immersion `C^* x C` | 1 at (2, 3), 0 at (0, 3) (first audit) | functional `L` on 110 generators; `L(1/x) = 1/a`; `S_z = 0` iff `a = 0` by hand |
| iterates of `d_u, d_v` on the slice | "N" for `m <= 4` | pole orders `2m + 1` (`x`), `2m - 1` (`y, w`) along `Omega`, `m <= 6`, and the lemma for all `m` |
| Proposition 6.3, formulas of the proof | by hand | 6 identities symbolic; `g_n`, `lam_G`, `lam_res` for `|n| <= 6`; codimension 2 in the window `|n| <= 40` |
| 6.3, centre (0, 0): dimension | 2 (`N' = 34`) | 2 (`N' = 28, 36`) |
| 6.3, centre (0, 0): classes | 5 elements, `m <= 6` | 7 elements, `m <= 7`; cone classes for 3 elements, `m <= 4`, `k <= 3` |
| 6.3, centre (2, 3) | 3; no vanishing | 3; no vanishing |
| 6.3, centres on `ab = 0` | not considered | 2 at (0, 3), (0, -7), (2, 0), (5, 0); one-sided vanishing (proved by hand) |
| `C^* x C^* -> C^2` | not considered | 2 at generic centres, 1 on `ab = 0`; no vanishing at pseudo-random centres |
| `runs/special_centre.txt`, `runs/inject.txt` | stored | reproduce byte for byte |
| `x^2 y`: identity, functional | stated | exact; functional on 96 monomials |
| `x^p y^q`, `1 <= p, q <= 4` | "`gcd(p, q)`; 1 in `M` iff `p != q`" | same in all 16 cases; elementary proof |
| `x + x^2 y` | "value -1 on 1" | `[1] != 0`, `[x] = -[1]`, `[x^2 y] = [1]` |

No number of the new parts was found wrong.

## 6. Literature and credit

*Sources read (text extracted from the arXiv PDFs; arXiv:0902.0212 and arXiv:1008.3962 in the relevant parts only).*
W. Zhao, arXiv:0902.0210v2 (J. Algebra 324 (2010)): Conjecture
1.3; conditions (C1), (C2) and Theorem 2.5; Example 2.6 (non-constant leading coefficient), Example 2.7
(characteristic `p`, by `f = 1`); Lemma 2.9, Theorem 2.10, Corollary 2.11; Theorems 3.5-3.7 (JC for all `n` as special
cases for `d/dz_i - xi_i`). A. van den Essen, arXiv:1006.5801v1: `IC(n, A)` with a regular sequence; the 1-property;
"General Image Conjecture" with Comments 1 and 2. A. van den Essen, D. Wright, W. Zhao, arXiv:1008.3962v2 (J. Algebra
340 (2011)): Conjecture 2.1 and its footnote (which records that the formulation of 2010 allows more general sets of operators);
Remark 2.3; Theorem 4.9. The same authors, arXiv:1004.0521: Lemma 2.2, Example 2.4 (`d_x - y^2 d_y`: `1` in the image,
image proper), Question 4.2, Theorem 4.3; the README's description of Theorem 4.3 is accurate. W. Zhao,
arXiv:0902.0212: Lemma 4.5 (1-property), Proposition 2.6 (the conjecture fails for Laurent polynomials). Summary only:
van den Essen, arXiv:1907.06107; MathWorld, "Image Conjecture"; Dvorsky, arXiv:2608.07338 (abstract). Canon: THM-2801,
guardrail 58.

*Is the content of 5.2-5.3 in the literature?* Not found: the family of Stein kernels over the plane of centres, the
torsion-freeness of `B` over `C[f - d_v, g - d_u]`, or "JC(2) iff `Im(d/dx - (fg)_x) + Im(d/dy - (fg)_y)` is a Mathieu
subspace for every Keller pair". Zhao's Theorems 3.6-3.7 are statements in all dimensions for the operators
`d/dz_i - xi_i`; Theorem 4.3 of arXiv:1004.0521 is planar and untwisted (`Im D`, `D = {., g}`, by the 1-property).
Torsion-freeness is equivalent to "an algebraic function times `e^(-uv)` is not annihilated by a nonzero
constant-coefficient operator in `(u, v)`", which is elementary by asymptotics and may well be folklore; for the
standard pair it is the invertibility of the heat operator on polynomials. One targeted web search found nothing; that
is not a certificate. The README claims no priority, and I found no sentence that overclaims, apart from M1.

*Session269.* The README's description is accurate: Theorem 1 (first hits in both directions for every monomial on an
open set give mixed finite chaos outside a countable union of curves; constructible membership, Baire, the connection
`d/da + V`, `d/db + U`), Theorem 2a (added at 13:08: one centre algebraically independent over the coefficient field),
Conditional Corollaries 3-4 (with the count and the finite-chaos criterion of THM-5138). They were on `origin/main` at
13:03 and 13:08, built on the pin `1ca9c5f83b`, which contains the checkpoint `86ab0707d6`; the statements of 5.2-5.3
were pushed at 14:58. The first-pusher credit is stated fairly; see m9 for two smaller credits. Session269 has no
torsion-freeness and no one-direction statement (only the existential remark that some pure ray is not killed at a
very general centre).

*Later commits* (`git fetch origin` at 15:30; `git log origin/main --since 13:00`). After `5384be8b8f` (Session269,
13:08) the only commits that touch Stein modules, Keller pairs or detection are the author's own (`c440c4a4b7`, 14:58;
`182fd02979`, 15:01). Session272 (`26fa4d8d52`, 15:03) mentions the audit as a guardrail and uses none of it.
Sessions 270, 271, 273 and opus-S92/S93 do not touch the subject.

## 7. Not verified

- The proofs of the cited theorems (Dimca-Saito, Sabbah, Nguyen Van Chau, Bavula) and the canon theorems used as
  given: THM-2022, THM-4620 (1), (5), (6), THM-4122 (4), THM-2801 (statements read, proofs not re-checked).
- The public sources of 2026 quoted in M1 (summaries only); the 2021 book of van den Essen, Kuroda, Crachiola and
  Zhao's "Mathieu subspaces of associative algebras" were not consulted. "No record found" is bounded by section 6.
- K-GMC at the centre (0, 0) of 6.3 (section 4, J); whether `SN_v` fails at the centres `(a, 0)` of S2 (no vanishing
  for `m <= 7`, not proved); `[1]_z != 0` at those centres by machine only.
- Machine statements of non-membership are relative to the truncation; memberships and dimensions are computed
  modulo one prime (`2^61 - 1`); only the full ranks of section 4, B are certificates over `Q`.
- The novelty of torsion-freeness (one search).
- THM-5138 and HYP-9668 were not read.

## 8. Files (`audit/audit2/`)

`eng2.py` (engine; self-test against sympy); `b1_lemma52_valuation.py`, `b1b_lemma52_slice_places.py` (about 100 s),
`b1c_lemma52_ramified_place.py` (Lemma 5.2); `b2_lemma53_torsion.py` (Lemma 5.3; 35 s); `b3_thm54_tests.py`
(Theorem 5.4, m1); `b5_slice_which_generator.py` (S1; 45 s); `b6_prop63_special_centre.py`, `b6b_hyp_axis_centres.py`
(Proposition 6.3, S2); `b7_section8_image_conjecture.py`, `b8_broughton_class_of_one.py` (section 8);
`b9_generic_centres_other_surfaces.py`; the `.out` files next to them; `repro_special_centre.txt`,
`repro_inject.txt` (the author's scripts re-run). Reproduce with `python <script>` in that folder.
