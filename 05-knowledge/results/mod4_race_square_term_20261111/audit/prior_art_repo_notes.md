# Audit part D - prior art inside the repository (working notes)

Scope: claims typed NEW or "not found in the repository" in `README.md`, THM-5131, THM-5132 of the packet
`mod4_race_square_term_20261111`. Read-only sweep of `01-canon/theorems`, `05-knowledge/results`,
`05-knowledge/hypotheses` (plus `01-canon/MISTAKES.md`, `04-computation`, `agents/*/inbox` where a grep led there).
No computation was done here. Judgments: IDENTICAL / IMPLIES IT IN ONE LINE / PARTIAL / RELATED ONLY.

## 1. Overlaps, by seriousness

### O1 (MAJOR). Theorem S (2), the mod-8 refinement `H = 1 + 4s (mod 8)` - IDENTICAL, stated at least twice before

- `05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md`
  (committed 2026-10-22, 25a4f6ba21), section 3 "The reflection-overlap congruence":
  - l.33-36: "Let A be the number of two-element cycle orbits. Partition these A orbits into: D admissible pairs: the
    two cycle supports are disjoint; I overlapping pairs: their supports meet."
  - l.60-62: "For every self-converse tournament, H(T) ≡ 1+2F+4I (mod 8)."
  - l.64-66: "Indeed OCF gives H≡1+2 alpha_1+4 alpha_2 modulo eight. Here alpha_1=F+2A. ... so alpha_2≡D modulo two.
    Hence H≡1+2F+4(A+D)≡1+2F+4I (mod 8)."
  - l.68-73: "At even order this specializes to H(T) ≡ 1 (mod 4), H(T) ≡ 1+4I (mod 8). Thus the next binary digit
    measures the parity of overlapping mirror-pairs of odd cycles".
  Same statement (`s = I`), same proof (`alpha_1 = 2t`, `alpha_2 = t_0`), plus independence of the chosen
  anti-involution (section 1, conjugacy), which the note does not have.
- Independent audit of that report: `.../session29-congruence-audit/AUDIT.md` l.15 (ACCEPTED): "alpha_1=F+2A and
  alpha_2≡D mod2. Substitution in H≡1+2alpha_1+4alpha_2 mod8 gives H≡1+2F+4I."
- `05-knowledge/results/session130_blue_black_theorems_20261031/hamiltonian/PROOFS.md` (committed 2026-10-31,
  41a9d96041), Theorem 1, l.61-79: "For every finite tournament with an involutory anti-automorphism,
  `H(T) = 1 + 2q + 4i (mod 8)`. ... For even order, `q=0`, hence `H(T)=1+4i (mod 8)`. This refines, but does not
  replace, inherited THM-791." (also README section 4 of that packet, l.51-53).
- The pointer was in the author's own lane three days ago:
  `05-knowledge/results/eternal_colours_audit_K1_20261108.md` l.257-259 (audit of mac-mini's THM-4930): "Session29's
  converse-congruences report (§§1, 3) proves the same parity, with a mod-8 refinement `H ≡ 1 + 4I (mod 8)` at even
  order. The THM-4922 audit observed it on all 8,784 self-converse classes of order 10."
  Also `first_creators_genealogy_20261107/audit/AUDIT.md` l.326: "Session29 gives H = 1 + 2F + 4I (mod 8) without
  the balance." THM-4930 (cited by the note) cites "Session29 converse-congruences report, section 3" four times;
  MISTAKE-738 (01-canon/MISTAKES.md l.23747) is the process lesson about exactly these Session29 reports.
- What the note says: THM-5131 status "NEW here: (T); the mod-8 refinement (S) (2); ..."; README outcome 3
  "(PROVED + FINITE-EXACT)"; Theorem S (2) carries no credit.
- New in the note on this point: only the class split `1/1, 6/6, 82/94, 4448/4336` (H = 1 vs 5 mod 8 over
  self-converse classes, n = 4, 6, 8, 10). (Session130 has the split over blue tilings, a different count:
  `1/3, 28/36, 1852/2244` at n = 4, 6, 8.)

### O2 (MAJOR as a credit/typing matter). Theorem T (3), `E(-1)^(c_3) = 2^(-floor((n-1)/2))` - PARTIAL (same closed form displayed in canon since June; other ensemble, verified not proved)

- `01-canon/theorems/THM-555-tiling-cycle-moment-hierarchy.md` (kind-pasteur-2026-06-20-S20, f727806f8b), l.49:
  "E[(-1)^c3] = 1 / 2^{floor((n-1)/2)}   (0 at n=3) => Pr[c3 odd] = 1/2 - 1/2^{floor((n-1)/2)+1}  [VERIFIED n<=10]";
  l.52-53: "always biased toward EVEN `c3`, the bias halving every two steps"; l.78-79: "THM-552 (c3-parity
  dichotomy): the ensemble even-bias `E[(-1)^c3]=1/2^{floor((n-1)/2)}>0` is the statistical SHADOW of the pointwise
  dichotomy." Ensemble there: the tiling model (tournaments containing a fixed base path = labelled tournaments
  weighted by H), hence the exception at n = 3; status VERIFIED n <= 10, not proved.
- `05-knowledge/hypotheses/HYP-2707-tournament-clifford-magic-hierarchy.md` (kind-pasteur-2026-06-20-S22),
  l.44-52: "The c3-parity `q(b) = c3(b) mod 2` is a GF(2) quadratic ... function ... Its associated bilinear form ...
  has GF(2) rank `r = 2*floor((n-1)/2)` ...: r = 2,4,4,6,6 for n=4..8. The standard Boolean Gauss sum then gives
  `E[(-1)^c3] = ... = 2^{-r/2} = 2^{-floor((n-1)/2)}` — which is EXACTLY THM-555's c3-parity formula."
  (mechanism = Gauss sum of the quadratic form; rank checked n = 4..8 only). Also HYP-2705 l.261, HYP-2711 l.83,
  HYP-2702 l.319, `07-reflections/the-eight-external-ideas-...md` l.24, and the broadcasts
  `agents/mac-mini/inbox/MSG-166-...`, `MSG-171-...` (both delivered to the author's machine).
- THM-559 (mac-mini-2026-06-20-S7, the author's machine): `c3` is an exact 2-body form on the arc spins
  (`c3 = n(n^2-1)/24 - E_score/2`), the same quadratic form.
- What the note says: README section 3, Typing: "Elementary. Not found in the repository (THM-478 has the quadratic
  form, THM-552 the self-converse case)."; THM-5131 "NEW here: (T)".
- What is actually new: the proof for every n; the uniform labelled ensemble (THM-555's is H-weighted; the two agree
  only because conditioning on the base path does not change the Gauss sum, which nobody proved); the pointwise
  statement (T) (2) (every vertex cube has sum `2^ceil(m/2)`) - not found anywhere.
- Suggested correction: credit THM-555 / HYP-2707 (formula and Gauss-sum mechanism, tiling ensemble, verified
  n <= 10) and THM-559; retype (T) (3) as "PROVED here for all n and the uniform ensemble; closed form first recorded
  in THM-555"; delete "Not found in the repository".

### O3. Theorem S (3), the count bound - IMPLIES IT IN ONE LINE; exact values only partly new

- `.../session29_framed_profiles_20261022/session29-halfpath/REPORT.md` l.72: "Fix sigma and choose uniformly among
  all labelled tournaments admitting it. There are 2^(m^2) choices at n=2m and 2^(m(m+1)) at n=2m+1: these are the
  independent unordered-edge orbits under sigma."; l.82: "A self-converse class T contributes exactly
  2^m m!/|C_G(sigma)| distinct labelled tournaments in the fixed-sigma ensemble. Consequently ...
  sum_T 1/|C_G(sigma)| = 2^(o-m)/m!". With `|C_G(sigma)| <= |G|` this is
  `SC_n = sum_T n!/|G| <= n! 2^(o-m)/m! = (n-1)!! 2^(n^2/4)` (even n), `n!! 2^((n^2-1)/4)` (odd n).
- THM-4743 l.78-79: "Every self-converse tournament has an anti-involution with at most one fixed point, so
  enumerating the tournaments anti-invariant under one fixed involution covers all of them. At `n = 11` that is
  `2^30` labelled cases". THM-4930 l.175-177: "all such involutions are conjugate. So the `2^25` labelled
  tournaments anti-invariant under one fixed involution contain every self-converse class of order 10".
- `.../session29-selfconverse-asymptotics/REPORT.md` l.14-15, 30-37: exact partition formula and
  `b_(2m) = 2^[m(m-1)]/m! [1 + ...]`, `b_(2m+1) = 2^(m^2)/m! [1 + ...]` for the number of self-converse classes
  (`n! b_n ~ (n-1)!! 2^(n^2/4)`: the note's bound is the leading term); cites Eplett 1979, Harary-Palmer 1966,
  A002785, and THM-283 (`01-canon/theorems/THM-283-sc-burnside-formula.md`, class counts
  `1, 2, 2, 8, 12, 88, 176, 2752, 8784, ...`).
- Numbers: `05-knowledge/results/selfk_type_ledger_fixed_opus_S312.out` (opus S312, July): "(1, 2, 2): 960",
  "(2, 2, 2): 7680", "(1, 2, 2, 2): 430080" = number of pairs (labelled tournament, anti-involution) at n = 5, 6, 7
  = the note's bounds 960, 7680, 430080. The superseded `selfk_type_ledger_opus_S312.out` l.8, 16, 26 prints
  "704", "7040", "380416" (= SC_5, SC_6, SC_7) as an unlabelled by-product of a mis-typed cycle type.
- Not found: SC_n for n = 3, 4, 8, 9, 10 (48, 6737920, 958578688, 31596429312), the shares of D_n, "majority
  exactly for n <= 5" (the 8 : 4 split itself is in `cpq_defect_family_mirror_folds_20261104.md` l.192-195, typed
  NUMEROLOGY there).
- THM-5131 types "the count (S) (3)" NEW; it should be "one-line corollary of the fixed-involution ensemble
  (Session29 half-path section 5, THM-4743 (6), THM-4930); exact values new".
- Theorem S (1): the note credits THM-791 / THM-552. The first sentence ("if T has an anti-automorphism it has one
  of order 2", at most one fixed vertex) is not in THM-791; it is THM-4930 statement 2 and proof 2, and Session29
  converse-congruences section 1 (which also proves conjugacy).

### O4. Theorem P (4), first sentence, and the citation "KNOWN: THM-1470"

- THM-1470 (klein-2026-07-20-S338, not "opus S409" as README section 1 says; THM-1460-babai-cameron is opus S409)
  proves the bijection only for odd n inside a switching class (l.73-76: "For `n` odd the `2^{n−1}` members of a
  switching class realise `2^{n−1}` distinct score-parity vectors"; l.86: n even "the parity vector is pinned only
  up to global complement. No analogue.").
- The all-n statement is `.../session72_even_graph_unifying_carriers_20261026/session72-even/SCORE_PARITY_QUOTIENT.md`
  l.5-17: "Reversing an edge uv toggles exactly p_u and p_v ... The kernel of ∂ is Z, and for n≥1 its image is all
  even-weight vectors ... Tour/Z ≅ P_n as affine S_n-sets", l.31: "for n≥1 the probability normalization is
  binom(n,r)/2^(n−1) on the allowed parity sector." This is the fact the proof of (4) uses
  (`2^(C(r,2) - (r-1))` tournaments per parity vector). Judgment: the KNOWN typing is right, the citation is wrong
  for even n.
- The joint law of `(s, t)` itself: not found.

### O5. Objects of Theorem P that the note types "not found"

- The matrix `A` over `F_2`, its rank and kernel: THM-4857 (mac-mini-2026-11-01-lenses, the author's machine),
  "The GF(2) automaton of a tournament ... `Phi_T(x) = Ax` over GF(2)", with "de Caen (1991) for rank >= (n-1)/2"
  (status block) and the Paley/QR-code cases (THM-481). No statement there about weights mod 4 of kernel elements.
- The vector `A^2 1`: over the integers it is the Kendall-Wei second-order score ("score-of-scores"), THM-896
  (boxeph S30) title: "self-composition R∘R counts 2-paths (the score-of-scores)"; HYP-7029 (INDEX-HISTORICAL
  l.15247): "per-vertex joint (s^(1),...,s^(6)) = rowsums of A..A^6"; THM-555 l.60-63: "`c5` and `alpha_2` are NOT
  score-determined ... Adding `A^2` row-sums does not rescue the distribution — fails already at `n=6`";
  THM-4911 l.225-226 "classical second-order score quantities". README 6.2 Typing says "Not found in the
  repository: the vector `t = As`": true for its reduction mod 2 and its weight mod 4, not for the vector.
  Judgment: RELATED ONLY; credit lines missing.
- (2), first equality `x^T (A + A^3)(1 + x)`, and (T) (1): the layers "paths on 2 and 4 vertices" of THM-5127
  (4.4), `q_R(X) = sum_{i in X, j notin X} P_ij (mod 2)` (tilted README l.232). Credited in substance.

### O6. Values already on record

- `D_4 .. D_7`: THM-466 l.149 "P(v_2(H−1) = 1) = P(α_1 odd) climbs 0.250 → 0.297 → 0.375 → 0.448", i.e.
  `D = 1/2, 13/32, 1/4, 853/8192` (also THM-5127, credited).
- `E_5[{5}] = 7/32`: THM-466 l.164 "n = 5: 400/1024 tournaments have c5 odd".
- Class race at n = 6: THM-1980 l.66 "n=6: 32 classes `≡1`, 24 `≡3`".
- Partial law for `H mod 4` at n = 5, 6: HYP-6935 (INDEX-HISTORICAL l.919, opus S316) [ii] "for ODD c3,
  H = 3 mod 4 iff c3 = 1 mod 4 holds exactly at n = 5, 6". The corollary of Theorem P is the complete formula.
- Mod-8 masses at n = 10, 11: first and second audit of THM-5127 (credited; values agree, see section 2).

### O7. Deck / subset remarks - RELATED (k = 4 case)

- THM-1455 (kind-pasteur) title (III): "k₄ is EVEN whenever n is ODD ... double counting — each 4-subset lies in
  exactly n−4 five-subsets, so (n−4)·k₄ is even". Since `(n-3) c_3 = k_4 + 2 (#strong 4-sets)`, this is the case
  `k = 4` of the note's `Y_k = sum_l C(n-l, k-l) c_l`. Nothing found for `chi(T) = prod_v chi(T - v)` or for
  `p_k mod 4`.

### O8. Proposition 5.1 - IMPLIES IT IN ONE LINE (classical identity with a multiplicative weight)

- Unsigned identity: THM-4816 l.75, THM-4817 l.36 "`T_lab = 1/(1 - S_lab)`".
- Same-day codex note `spectrum_multiplicity_signed_tables_20261111/WEIGHTED_ATLAS.md` l.81-94, eq. (5):
  "N_n(h)-N_(n-1)(h) = sum_(q|h,q>1) sum_(r>=3) a_r(q) N_(n-r)(h/q) ... Replacing a_r(q) by the sum of reciprocal
  automorphism orders gives the same identity for P." Summing (5) against the multiplicative sign gives
  `T_chi = 1/(1 - S_chi)`. The note credits the codex sign coordinates. The signed strong counts
  (-2, 24, 176, 4080, ...) were not found anywhere.

## 2. Citation accuracy

| item | result |
|---|---|
| THM-466: `H = 1 + 2 alpha_1 + 4 alpha_2 (mod 8)`, alpha_2 = pairs of vertex-disjoint odd cycles | OK (l.22-28: "H ≡ Σ_{k<m} α_k 2^k (mod 2^m)", α_k = "collections of k pairwise vertex-disjoint directed odd cycles") |
| THM-478: `c_3 mod 2` quadratic in the arcs | OK (l.67-69) |
| THM-791: fixed-point-free involutive anti-automorphism => c_odd even, H = 1 mod 4 | OK (l.49-50). It does not contain the reduction from an arbitrary anti-automorphism (see O3) |
| THM-552: phi-self-converse, n even => c_3 even | OK (l.25) |
| THM-4922 (4): odd-order even split, 1/1, 4/4, 44/44, 1376/1376 | OK (l.85-87) |
| THM-1470: score-parity vector uniform on its coset | PARTIAL MISMATCH: odd n, inside a switching class only (O4); author is klein S338, not opus S409 |
| THM-1980: H mod 4 spectral for n <= 6, not n = 7; d(6) = 2 | OK (l.52-59) |
| THM-500: `alpha_1 = tr A^3/3 + tr A^5/5` for n <= 6; `tr A^7 = 7(c_7 + TQ)` | OK (l.33-35, l.73) |
| THM-118: `tr A^k = k c_k`, k = 3, 4, 5 | OK (l.13) |
| THM-5127: D_n fractions for n <= 11 | OK (l.189-196; tilted README l.297-303, 391) |
| THM-5127: `abs(D_n) <= 1/2`, Corollary 4.6, gate (C) (7) | OK (tilted README l.424-427 "equality holds at n = 3, 4"; THM-5127 l.168-170) |
| mod-8 masses at n = 10, 11 | OK up to rounding: audit l.398 "0.25082, 0.24901, 0.25039, 0.24978" (x4 = 1.00328, 0.99604, 1.00156, 0.99912); audit2 l.372-373 "0.25013, 0.24987, 0.25011, 0.24989" (x4 = 1.00052, 0.99948, 1.00044, 0.99956). The fifth digits at n = 10 (1.00330, 0.99603, 1.00155) are not derivable from the five-digit fractions on file |
| HYP-9355: "absent at n = 11"; rows n = 7-9; 1.0030 +- 0.0009; 0.9996 pooled | OK ("resolved through n = 10, absent at n = 11"; table to 3 decimals: 1.179 0.985 1.029 0.807 / 1.064 0.979 1.019 0.938 / 1.015 0.992 1.000 0.993) |
| THM-1370: Busch floor, max H | OK (f = 3, 5, 9, 15, 25, 45; 189 at n = 7, 661 at n = 8). Its crevasses are (f(4), f(5)) = (5, 9) and (f(6), f(7)) = (15, 25) |
| THM-4889 (3): windows of 1 to 4 orders for h < 225 | OK ("length 1–4 (14, 27, 58, 10 values)") |
| "its length grows like 1.86 ln h" (README section 8, given as part of an "exact form") | MISMATCH: not in THM-4889. Source is HYP-9520 l.65, section "Heuristic (HEURISTIC)": "`B(e) ≈ 1.86 ln e` (Busch)" for the upper end of the window; that windows reach B(h) or B(h)-1 in general is HYP-9481 (OPEN) |
| THM-4026 / THM-4028 | OK ("N=896315812331399"; "V X^(25/24)+O(X^(11/12))") |
| HYP-9659 (B''): D_n -> 0 open | OK (l.89) |

Outside part D, noticed in passing: README section 8 row 1 says 7 and 21 are "carried by" bands
`[f(n), max H(n)]` that do not overlap (`f(n+1) > max H(n)` for n <= 4). That covers 7 only; 21 lies inside the
band `[15, 45]` of order 6 and is excluded by "21 ∉ strong(6)" plus `21 = 3·7` (THM-1370 (4)-(5)).

## 3. Claims for which no prior statement was found (patterns used)

- (T) (2), the pointwise cube sum `2^ceil(m/2)`: `2^ceil`, `c_?3\((R|T) ?\+ ?v`, `triangles through`,
  `e_?R?\(X`, `\(-1\)\^c3`, `E\[\(-1\)\^`, `c3-parity`, `triangle parity`.
- `c_3 = C(q_1,2) + floor(n/4) (mod 2)`: `floor\(n/4\)`, `⌊n/4⌋`, `odd[- ]scores?`, `number of odd`,
  `score[- ]parit`, `parity vector`, `evenly many odd`, `mod 4` near odd scores. Nearest: KBS formula (THM-118
  corollary, THM-466 (v), THM-559 (3)), "c3 mod 2 (classical): const" on score classes (THM-466 l.155), the residue
  argument of THM-867 l.32 ("forces evenly many odd e_v"), `x = n(n^2-1)/3 - 8 C3` (THM-790 l.126). One line from
  these, never written.
- (P) (1) for `c_3 + c_5` and `c_5`, (2) second equality, (3), (4) joint law, (5) closed forms and the values
  67/256, 263/2048, 1007/16384, 11/256, 23/2048, 49/16384: `c_?5 (mod|parity)`, `parity of c_?5`, `c5 odd`,
  `tr ?\(?A\^5`, `A ?\+ ?A\^3`, `A\^2 ?(1|row)`, `second-order score`, `ker ?A\^T`, `weight mod 4`,
  `doubly[- ]even`, `Arf`, `de Caen`, `2-rank`, `rank over F_2/GF(2)`, the fractions as strings. Only THM-466
  mentions `c5 mod 2` (non-constant on score classes; 400/1024 odd at n = 5).
- Three strata, class race 36, 264, 1268, 22600, signed strong counts, table `E_n[K]` (811/8192, 1553/8192,
  461/8192, 533081): number strings; `(rigid|non-trivial automorphism).{0,160}mod 4`.
- "Square term" reading: `square term`, `chebyshev` near `self-converse` / `tournament`.
- Persistent lead `D_n > 0`, the `4^(-n)` plateau: `D_n ?> ?0`, `positive (lead|discrepancy)`, `4\^\(-n\)` in
  HYP-9659 and the tilted packet.
- Same-day packets (2026-11-10/11, 50 folders): no statement on triangle parity, self-converse mod 8 or the mod-4
  lead.
- MISTAKE-753 / 755 lane (THM-4097..4184): no mod-4 or odd-cycle statement in the tournament files of that lane
  (counted per file); nothing there overlaps this packet.
