# Independent adversarial audit: supersingular counts, divisor shapes, six q-series families (audit_ssj)

*The main session (opus-2026-10-22-ssj) saved this file. The auditing subagent's own write of REPORT.md was refused by its session, so it returned the text instead. Content is as the auditor returned it, lightly formatted. Its scripts and outputs are in this directory. The downloaded sources are kept local (`downloads/` is git-ignored); their sha256 prefixes are listed at the end.*

**Audited:** `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md` (commit 33f9f99615, "pre-audit"), with the scripts and outputs in `04-computation/experiments/supersingular_fus_qseries_20261022/`.
**Auditor:** independent subagent, temporary identity `audit_ssj`, 2026-10-22.
**Rules kept:** the note and the original scripts were not edited, and nothing was committed. PARI stacks were at most 512 MB, one heavy job ran at a time, and `al_genus_check.gp` was rerun only at reduced size. Downloads are stored as data under `audit_ssj/downloads/` and were never executed.
**Method:** each claim was checked by an independent route:
* point counting instead of Hasse invariants;
* curves and Tate normal forms instead of Hauptmodul fibres;
* Hecke traces instead of Atkin–Lehner operators;
* numerical evaluation instead of q-expansion linear algebra;
* my own genus formula instead of `mfdim`.

The committed scripts were also rerun in a scratch copy and their outputs compared.

## Summary verdict: PASS WITH CORRECTIONS

All of these pass every independent test:
* F, U, S, their level sets and the recursion;
* the cusp and Atkin–Lehner dictionary;
* the supersingular census and Ogg's list;
* the 17 DS15 identities;
* the Hauptmoduln;
* the mass-formula counts;
* the Frobenius = ⟨p⟩ proposition;
* the §6.3 class-number theorem and all its local numbers;
* the level-N Ogg lists;
* the Deligne–Rapoport genus identities.

Four statements are factually wrong (F1–F4), several overclaim or rest on an unstated convention (F5–F9), and the rest are wording or reproducibility points. No error reaches a proved result.

| # | Severity | Location | Verdict |
|---|---|---|---|
| F1 | high (false; the note's own output contradicts it) | §6.3, last paragraph | ERROR |
| F2 | medium (bibliographic) | §1 | ERROR |
| F3 | medium (formula wrong as written) | §6.3 table and last sentence | ERROR |
| F4 | low (NUMEROLOGY row) | §7, C8 | ERROR |
| F5 | medium | §0 item 1, §2 | OVERCLAIM |
| F6 | medium (convention) | §5 "Diamond operators", §6.2 T4 | WORDING |
| F7 | low–medium | §6.2 "Reading", §0 item 4 | OVERCLAIM |
| F8 | low | §0 item 2 versus §9 | WORDING |
| F9 | low | §6.3 (c) | WORDING |
| F10 | low | §6.2 proof, step 1 | WORDING (gap) |
| F11 | low | §6.3 Ogg lists, C11 | WORDING |
| F12 | low | assorted | WORDING |
| F13 | low | README step 4 | reproducibility |

## Findings and repairs

**F1. "At the identity levels the split is F_{p^2}-heavy for every p." (§6.3) — ERROR.**
* *Checked:* curve-level counts for p < 104 (`audit_fp_points.gp`, Part 1), and a scan of the confirmed formulas for p < 2000 (`audit_misc.gp`, M2).
* *Finding:* X_0(8) has more F_p-rational points than not at two primes:
  * p = 7: 4 of 6 points are F_7-rational;
  * p = 23: 12 of 22 points are F_23-rational.
* The note's own `al_genus_check.out` shows both, in rows `[8,7,5,1,6,4]` and `[8,23,21,5,22,12]`.
* X_0(7) and X_0(9) at p = 5 split exactly evenly (2 of 4).
* There are no other exceptions below 2000. Above 2000 none can occur, because the F_p count is at most 6n_p = O(√p log p), while the total grows linearly in p.
* *Repair:* "…F_{p^2}-heavy for every p except p = 5, 7, 23: X_0(8) is F_p-heavy at p = 7 and 23, and X_0(7), X_0(9) split evenly at p = 5."

**F2. "The most recent arXiv papers on supersingular curves are 2609.11037, 2609.03249, 2608.13371" (§1) — ERROR.**
* *Checked:* an arXiv API search, `all:supersingular` sorted by date. The three titles themselves are described correctly.
* *Finding:* later papers exist, for example:
  * 2610.21488 (15 Oct): minimal isogeny degree from E to E^(p);
  * 2610.10139 (7 Oct): maximal orders and isogeny graphs;
  * 2609.22018 (18 Sep): the supersingular isogeny problem in p^{1/3+o(1)};
  * 2609.17819 (15 Sep): J_0(N) is supersingular mod p for only finitely many N.
* *Repair:* "three recent papers (not exhaustive)". 2610.21488 and 2609.17819 fit the brief's F_p/F_{p^2} and level themes better, so they could be offered to the user as candidate sources.

**F3. "For ℓ = 7 and 9 … (1 + (−p/ℓ)) n_p", and the table row "X_0(ℓ), X_0(ℓ^2), ℓ odd (here 7, 9)" — ERROR (notation).**
* *Problem:* (−p/9) = 1 for every p ≠ 3, so the formula as written gives 2n_p. The scripts correctly use (−p/3).
* *Checked:* X_0(9) counted independently as Frobenius-stable ⟨2⟩-orbits on the Tate-normal-form X_1(9), for all p < 400. The count is (1+(−p/3)) n_p; at p = 7 it is 0, not 2.
* *Repair:* "N = 7 (ℓ = 7) and N = 9 (ℓ = 3, e = 2)".

**F4. C8: "X_0(27) (27a: CM by discriminant −27)" — ERROR.**
* *Checked:* `ellisomat` and `ellweilcurve` on the class of y² + y = x³ − 7.
* *Finding:* the strong Weil curve is X_0(27) itself (genus 1): [0,0,1,0,−7], with j = 0 and CM by discriminant −3. Discriminant −27 belongs to the class members [0,0,1,−30,63] and [0,0,1,−270,−1708], which have j = −12288000.
* *Repair:* correct the parenthetical; the NUMEROLOGY verdict is unchanged.

**F5. "for all primes 5 ≤ p < 3000 by three independent routes" / "Three routes compute ss_p" — OVERCLAIM.**
* *Finding:* the routes cover different ranges:
  * route A covers all 428 primes p < 3000;
  * route B covers only the 123 primes p < 700;
  * route C covers only the 301 primes p < 2000, and only counts F_p-roots — it does not compute ss_p.
* *Repair:* state these ranges in §0 and §2.

**F6. "⟨d⟩ permutes the Klein forms f_a ↦ f_{da}; on M7, ⟨2⟩ is x ↦ (x−1)/x" — WORDING (unstated convention).**
* *Checked:* numerical values at τ = 0.1 + 0.3i.
  * For γ = [[4,1],[7,2]], whose lower-right entry is 2, M7(γτ) = 1/(1 − M7(τ)) = μ²(M7), to 1e−91.
  * For γ = [[2,1],[7,4]], whose upper-left entry is 2, M7(γτ) = μ(M7).
  * M9 behaves the same way.
* *Consequence:* μ = ⟨2⟩ holds when ⟨d⟩ is taken from the upper-left entry, i.e. the Q-rational Z/N-model, where Frobenius = ⟨p⟩. In the Diamond–Shurman lower-right convention, μ = ⟨2⟩^{−1} — this is the X_1(N) versus X_μ(N) ambiguity.
* T4 does not depend on this choice and is correct, and the Proposition was verified intrinsically (item 7 below).
* *Repair:* state the convention in one sentence.

**F7. "The identity index set is the Galois group of the F_{p^2}-structure… which is why every supersingular j lies in F_{p^2}"; "Frobenius… permutes each family's product sides" — OVERCLAIM.**
* The Galois group is ⟨p⟩, of order k(p), and it is trivial when p ≡ ±1 mod N.
* F_{p^2}-rationality is an input to the proof, not a consequence of it.
* Frobenius acts on points, while the diamond operators act on the Klein forms (functions).
* *Repair:* "Frobenius acts through the image of p in (Z/N)^×/±1, the group that permutes the product sides; the Galois group over F_{p^2} is ⟨p⟩."

**F8. "A new prime-adjunction recursion" versus §9 "No novelty" — WORDING.**
* The recursion is a two-line consequence of DB1; it is proved below.
* *Repair:* drop "new".

**F9. 6.3(c): "the two F_p-twists of j correspond to the classes of gπ" — WORDING.**
* This holds only at j = 1728, where G = Z/2 gives two classes for the two twists.
* For generic j, G is trivial; at j = 0, G = Z/3 has one twisted class.
* In both of those cases the two twists differ by −1, which fixes every level structure, so the factor 1/2 is still correct. Curve-level fibre counts confirm this.
* *Repair:* reword along these lines.

**F10. 6.2 proof, step 1, presupposes a curve E/F_{p^2} with π = −p — minor gap.**
* *Repair:* add "e.g. the base change of any supersingular E_0/F_p (a_p = 0, so π_p² = −p)".
* Optional: for j ∉ {0, 1728}, any F_{p^2}-model already has π = ±p.

**F11. Level-N Ogg lists — WORDING.**
* The lists are proved only for p < 2000, so C11 should say "(p < 2000)", or the note should add an explicit bound on h(−p) and h(−4p) that closes them for all p.
* "For N ≥ 2 the list is finite": it is finite for N = 1 too.
* The genus-relation sentence applies only to the X_0(N) lists.

**F12. Minor wording:**
* §4 table, "non-reduced": at p = 2 the middle multiplicity is φ(2) = 1, so the fibre is reduced.
* "affine in p with slope ψ/12": write "grows linearly".
* "bounded functions of the exponent profile": write "functions of the profile only".
* "The only genuine meeting point": soften "only".
* Ogg's list "reproduced for p < 3000": the primes 2 and 3 are hard-coded, so write "5 ≤ p < 3000".
* §0, "RS at q², Bailey at q³, GG at −q": add "(numerators)".
* §9, "a theorem once the local numbers are computed": add that s(floor) depends only on p mod N and s(surface) only on p mod 4N, so the finite enumeration covers every case.
* §1: list the 16 neighbouring arXiv IDs that were checked.

**F13. Reproducibility:** `level_hauptmoduln.gp` writes `level_relations.gp` with PARI `write`, which appends. Rerunning README step 4 in place doubles the committed file (3188 → 6376 bytes, same values). *Repair:* delete the file before writing.

## Confirmations, by brief item

1. **§1 source:**
   * The arXiv API and HTML confirm Bizeul–Klartag–Lehec, 2026-10-04, math.FA; the HTML sha256 is identical to the note's (`ebf2bf23e635b236…`).
   * The text has no supersingular, j-invariant or isogeny content; both "elliptic" hits are elliptic regularity.
   * The summary of Theorem 1.1 and the ingredients is accurate: (11) with the Hilbert–Schmidt norm, Lemma 4.2, Theorem 4.1, Eldan's process with inverse covariance, Proposition 7.1, Theorem 3.1 = [34, Thm 5.1] giving C_P ≤ CR², and [34] = arXiv:2610.01447.
   * The acknowledgement sentence matches verbatim, and the list of prior bounds is faithful.
   * None of the 15 neighbouring IDs I checked is about elliptic curves.
2. **§3:** confirmed and proved.
   * Values: rechecked with sympy divisor lists.
   * Level sets: δ = τ − 2^r − r − 1 + [sqfree] + [prime], and τ ≥ 3·2^{r−1} for non-squarefree N. This is DB2, re-derived, and matches an exhaustive run up to 3·10^5.
   * Recursion: F(Nr) = 2F + 2, S(Nr) = 2S + 1 + [N squarefree], U(Nr) = U + 1 + [N prime]; subtracting gives the recursion. It also holds symbolically for r ≤ 6, e ≤ 6.
   * Klein quadrant: matches THM-4733's conventions.
   * Cusps: counted by explicit Γ_0(N)-equivalence for 20 levels, equal to Σφ(gcd(d, N/d)).
   * Atkin–Lehner, with explicit W_Q matrices: each W_Q toggles the Q-part of d, the Fricke involution sends d to N/d, and ⟨C,Q⟩ meets the Atkin–Lehner group only in ⟨C⟩. Correctly typed as CITED/classical.
3. **§2 census:**
   * Full counts for 39 primes, covering every requested prime:
     * by numpy point counting over F_{p^2} for p ≤ 101;
     * by `ellcard` over F_{p^2} for p up to 383.
   * n_p by own point counting for all p < 400.
   * The Delfs–Galbraith formula holds in both its forms.
   * Degree law = g(X_0(p)) + 1 holds.
   * #pairs = #{newforms with a_p = −1}, from trace forms; Ogg's list to 3000 is confirmed via these traces.
4. **§5 identities:**
   * All 17 transcriptions match the DS15 PDF (2nd edition, Dec 2022, 61 pages).
   * All 17 hold to O(q^1200) in PARI, using product forms; the hostile controls fail, as they should.
   * The theta-family relations and sum-side corollaries hold.
   * The S55 note does state both relations the note attributes to it, and its §7 does type Ogg's primes as NUMEROLOGY.
5. **Hauptmoduln (beyond F6):**
   * j = A/B holds at the three requested τ for N = 5, 7, 8, 9 and for X_0(8), with errors ≤ 4e−95.
   * The norm maps t7, t9 and the classical level-7 and level-9 j-maps hold.
   * Γ_1(N)-invariance holds.
   * The factorizations of A and B are as stated, the field discriminants are 49, 81, 8, and the cubic fields are Q(ζ_7)^+ and Q(ζ_9)^+.
6. **§6.1:**
   * The ν values, the X_0(7) formula and d_N are correct.
   * #ss(X_1(N)) was recounted in Tate-normal-form coordinates.
   * The Eichler counts were recounted as ⟨2⟩-orbits.
7. **§6.2:**
   * Frobenius = ⟨p⟩ holds intrinsically in Tate normal form, in 5280 tests (N = 5, 7, 8, 9, p < 100).
   * The exceptional fibres are exactly the claimed ones; the j = 0 points on X_1(7) are ⟨2⟩-fixed, i.e. they lie over the elliptic points of X_0(7).
8. **§6.3:**
   * Delfs–Galbraith (a) checked from curves for 166 primes.
   * Brute-force F_p counts from E[N] ⊂ E(F_{p⁴}) on X_0(2), X_0(4), X_0(8) and X_1(8), for 25 primes, including per-fibre counts.
   * Tate-normal-form counts on X_1(4), X_1(5), X_1(7), X_1(8) and X_1(9) for p < 400, including the X_1(8) rule mod 16.
   * X_0(5), X_0(7) and X_0(9) via orbit counts.
   * Local numbers over complete residue systems.
   * Ogg lists for X_0(2) and X_0(4).
   * Genus relation via newform traces on 129 pairs.
   * Reduced Atkin–Lehner run: 66 pairs, 0 mismatches.
9. **§4:**
   * Deligne–Rapoport identities with my own genus formula: on 2070 pairs for pq, 541 pairs for p²q, and 4597 pairs for the general ℓM form.
   * The ℓ = 2, 3 Kronecker conventions are correct.
   * The example rows are correct, and the Katz–Mazur multiplicities equal the cusp weights.
10. **Ledger and numbers:**
    * The typing is correct apart from C8 (F4) and C11 (F11).
    * "No HYP … proved or refuted" is justified in substance.
    * Every number quoted in the note matches the `.out` files.
    * The committed outputs reproduce exactly, apart from the PARI stack-warning line.

**Not checked:**
* The Katz–Mazur / Deligne–Rapoport structure, Tate's theorem and Ogg's theorem beyond 3000 are cited; I checked only their genus consequences.
* The Hauptmodul relations remain FINITE-EXACT, as the note says; no formal proof was attempted.
* The KLS mathematics was not audited.
* The literature form of the §6.3 formula was not located.

## Checks run

All run from `audit_ssj/`:

| Check | Command | Result |
|---|---|---|
| F, U, S; cusps; Atkin–Lehner | `python3 -B audit_fus_cusps.py` | passed |
| Census by point counting | `python3 -B audit_census_pointcount.py` | passed |
| Census over F_{p^2} | `gp -q -f audit_census_ffp2.gp` | passed |
| Ogg's list to 3000 | `gp -q -f audit_ogg_traces.gp` | passed |
| Identities to O(q^1200) | `gp -q -f audit_qseries.gp` | passed |
| Hauptmoduln, numerical | `gp -q -f audit_hauptmoduln_numeric.gp` | passed; see F6 |
| Curve-level counts, Frobenius, genus relation | `gp -q -f audit_fp_points.gp` | passed |
| X_0(5), X_0(7), X_0(9) | `gp -q -f audit_x0_7_9.gp` | passed |
| Eichler counts as orbits | `gp -q -f audit_x0_eichler.gp` | passed |
| Local numbers | `python3 -B audit_local_numbers.py` | passed |
| Delfs–Galbraith classes | `gp -q -f audit_dg_classes.gp` | passed |
| Deligne–Rapoport genus | `python3 -B audit_dr_genus.py` | passed |
| Spot checks | `gp -q -f audit_misc.gp` | evidence for F1, F4 |
| arXiv | `curl` into `downloads/`, then `html_to_text.py` and `atom_titles.py` | evidence for F2 |
| Reproducibility | scratch rerun and diff | identical; F13 |

Each script writes a matching `.out` file in `audit_ssj/`. The downloads are in `downloads/{arxiv_api, arxiv_api_neighbours, arxiv_html, ds15}/`.

## Downloaded sources (kept local; sha256 prefixes)

| file | bytes | sha256 |
|---|---|---|
| `downloads/arxiv_html/2610.05474v1.html` | 765655 | ebf2bf23e635b236 |
| `downloads/arxiv_html/arxiv_2610.05474v1.txt` | 82968 | 6bd2edbd53019c34 |
| `downloads/arxiv_api_neighbours/neighbours.xml` | 38923 | 2f4ac09509cb32c0 |
| `downloads/arxiv_api_neighbours/three_ss.xml` | 6929 | de7c8cd534312d74 |
| `downloads/arxiv_api_neighbours/search_supersingular.xml` | 24414 | 5f92dbffad48ccaf |
| `downloads/arxiv_api/query_2610.05474.xml` | 2059 | 11a26c120191ab10 |
| `downloads/ds15/DS15.pdf` | 1093167 | aad948531839498a |
| `downloads/ds15/ds15_text.txt` | 231213 | 21a78cff8af7c32d |

## Corrections applied by the main session

All of F1–F13 are applied in the results note and scripts; section 10 of the note has the list, and the lineage is MISTAKE-666.

Three items go beyond the report:
* `level_fp_counts.gp` now prints the F_p-heavy and even cases.
* It enumerates the local numbers over complete residue systems mod 4N (92 classes, 0 deviations).
* `level_hauptmoduln.gp` now deletes `level_relations.gp` before writing it.
