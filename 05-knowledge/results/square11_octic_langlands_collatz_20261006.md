# The 11-square packing octic has Galois group S8 and its only unconditional Langlands piece is a Dirichlet character of conductor 16128215; the formalization's root isolation, reused on log2 3, gives a Lean-checked certificate of the S9 expense records; the 3SUM paper read against S18 (Artin, 20/19) and S19 (Mazur, Lean)

2026-10-06, session opus-2026-10-06-S16 (worktree `codex/session-square-packing-langlands-20261006`).

Owner's seed (verbatim): "consider the below formalization of the 11 square packing in lean checked, no need to double verify it yourself, but leverage similar ideas and expand on them toward our open problems like collatz: https://github.com/Queuingtheorydotcom/11SquaresFormalized Think of how the minimal polynomial of degree 8 here fits in with prior work on langland's program-adjacent ideas. s8−20s7+178s6−842s5+1923s4−496s3−6754s2+12420s−6865=0 [pasted: T = (6u+4)/(1+2u-u^2), u the unique root in (9/25, 37/100) of 5u^8-10u^7-2u^6+14u^5+12u^4-6u^3+2u^2+2u-1] also think about how many of the themes in this paper connect with our work on collatz, especially regarding 18 and 19 https://arxiv.org/html/2610.06783v1"

**Status.**

- **PROVED:**
  - `Gal = S8`, by Dedekind and Jordan with explicit primes;
  - `Q(T) = Q(u) = Q(tan(theta/2))`;
  - the elementary non-compactness proposition of section 5;
  - the glue lemma of the certificate.
- **FINITE-EXACT:** the field discriminant and ramification (Dedekind criterion), the class number of the quadratic resolvent, the Chebotarev and trace-moment statistics, the exact S8 entanglement table, and the S9 record list for `l <= 320`.
- **LEAN-CHECKED:** `04-computation/lean/standalone/collatz_expense_records_certificate_20261006.lean`. It uses core Lean 4.30 with `native_decide`, has no `sorry`, and rests on one native-decide axiom per theorem.
- **AUDITED:** an independent adversarial audit re-derived every computation with its own code and confirmed all computational and literature claims; its corrections to the Collatz statements and wording are applied in place (section 7) and logged as MISTAKE-566.
- **CITED:**
  - the packing formalization, whose optimality claim is not re-verified here (Friedman's 2005 survey still lists `s(11)` as unproven);
  - Stromquist, and Goebel's `s(5)`;
  - Trump's 11-square packing;
  - Eliahou 1993;
  - Hicks-Mullen-Yucas-Zavislak 2008 and Alon-Behajaina-Paran 2024;
  - Aramata-Brauer, Uchida, van der Waall;
  - Chebotarev, Jordan, Dedekind, class field theory;
  - Alman-Vassilevska Williams, arXiv:2610.06783.
- **ANALOGY / NUMEROLOGY:** typed in sections 6 and 7.

Nothing here proves a Collatz statement beyond the certified finite table.

Scripts:

- [part 1](../../04-computation/experiments/square11_octic_field_20261006.py): link, discriminants, Chebotarev, blocks;
- [part 2](../../04-computation/experiments/square11_octic_field_part2_20261006.py): S8 certificate, Dedekind criterion, resolvent, entanglement;
- [part 3](../../04-computation/experiments/square11_octic_field_part3_20261006.py): trace moments, tilt angle;
- [expense records](../../04-computation/experiments/collatz_expense_records_exact_20261006.py), with outputs in [square11_octic_langlands_collatz_20261006.out](square11_octic_langlands_collatz_20261006.out);
- the Lean file above.

## 0. Answers

1. **The octic is as generic as a degree-8 polynomial can be.**
   - The side length `T = 3.8770835900228...` and `u` generate the same field `K`. In fact `u = tan(theta/2)`, where `theta = 40.1819373 degrees` is the tilt of the rotated squares in Trump's packing, and `T = (2 + 2 cos theta + 3 sin theta)/(sin theta + cos theta)`.
   - `Gal(K^gal/Q) = S8`.
   - `d_K = -2^8 * 5 * 31 * 104053`.
   - Signature `(2, 3)`.
   - The earlier proven optima `s(5) = 2 + 1/sqrt 2` (Goebel) and `s(10) = 3 + 1/sqrt 2` (Stromquist, EJC 2003) live in `Q(sqrt 2)`, with Galois group `C2`. S8 is the generic Galois group of a degree-8 polynomial. The contrast is that those optima use only 0 and 45 degree tilts, while `n = 11` needs Trump's angle, which is neither a rational multiple of `pi` nor constructible (Corollary in section 2).
2. **Langlands, honestly.**
   - `zeta_K = zeta * L(s, rho_std)` with `rho_std` the 7-dimensional standard representation of S8, of Artin conductor `|d_K| = 4128823040`. Its traces `a_p = #{roots of P mod p} - 1` have the predicted moments `0, 1, 1, 4`; observed to `10^6`: `-0.003, 0.998, 0.975, 3.89`.
   - The only unconditionally automorphic piece is GL(1). The sign character of S8 is the Kronecker character of `Q(sqrt(-16128215))`, class number 5444, verified as the Frobenius parity on 2258 primes.
   - Automorphy of `rho_std` on GL(7) is open, and Dedekind's conjecture for `K` is not covered by the known theorems (Galois case or solvable closure).
   - This sits on the repo's existing Langlands ladder. Every exact algebraic number found in the Collatz notes generates an abelian field (GL(1)); this rests on a search, not a proof. The level-22 oldform bridge used GL(2) (`X_0(11)`, modularity known). The octic is a GL(7) Artin case (unknown).
   - The oldform bridge's lesson holds again: Frobenius data see the field and its conjugacy classes, never which conjugate (side length versus `-1.853` versus `u`) or the geometry.
3. **S18 revisited.**
   - Artin's `20/19` correction for base 5 is entanglement through `Q(sqrt 5)`, which lies inside `Q(zeta_5)`.
   - The octic has the same phenomenon through the quadratic subfield of its Galois closure (`K` itself has no proper subfields, since S7 is maximal in S8). For example, P has exactly one root mod p with density `31/84` when `(D/p) = +1` and `11/30` when `(D/p) = -1`, against `103/280` overall.
   - S18's Collatz verdict stands: no 2-adic correction exists, and all entanglement is 3-adic.
4. **What the formalization's method buys Collatz.**
   - The packing proof has four parts: (a) a finite closed-cell cover, (b) certified slack on every cell, (c) one exactly isolated algebraic endpoint (`u` in `(9/25, 37/100)`), (d) `native_decide` certificates.
   - (a) has no Collatz analogue for descent. Stopping times are unbounded, so no finite set of descent cylinders covers all large positive integers (Proposition 5.1). The coefficient-survival set `E_inf` (dimension `h* = 0.950`) is uncountable, not one rigid point. Even it does not capture actual descent, because coefficient exit ignores the carry (section 5).
   - (c) and (d) transfer. Isolating `log2 3` once, between the consecutive convergents `16785921/10590737` and `301994/190537`, decides all 320 S9 expense values and their eleven records. The decidable part is checked in Lean in seconds; the monotonicity glue is a two-line paper proof. This independently confirms S9's table; S9's 60-digit evaluations were already safe, since the closest value to an integer is `0.0019` at `l = 76`.
5. **The 3SUM paper.** Of its seven themes:
   - **five are typed analogies:**
     - computing only wanted entries (mac-mini's Proposition 5.1 is a proved Collatz statement in that spirit, but the link to the paper is an analogy: there is no product or identity);
     - identity-driven recursion with pruning (the S14 two-sheet fusion relation);
     - amortized encoding across products (banks);
     - lopsidedness (2-adic sources against 3-adic landings, S15 and S18);
     - refuting standing hypotheses;
   - **two are meta-themes:** AI discovery plus verification, shared with S19's Mazur paper and the 11-square formalization; and the algorithmic or finite reduction those successes rely on.

   On "18 and 19":
   - S18 is the Artin directive with "19/20", where 19 = `5^2 - 5 - 1`.
   - S19 is the Mazur Lean-checked positive-density paper.
   - The paper's 18, in `N >= D^18` and the `D^(1/18)` saving, is unrelated to S18 as a number (NUMEROLOGY).
   - Thematically, both S18's directive and the paper are about a small gain unlocked by restructuring a recursion. The paper's gain is a polynomial saving, which is a different kind of object from an Artin density factor. S18's evidence (2-adic statistics exactly Haar under population measures) suggests that a Collatz analogue would have to be 3-adic; no theorem says so.

## 1. Inheritance and portfolio

- **Prior Langlands-adjacent work** (searched as `langlands` across canon, results and navigation):
  - [S8 crossings / coincidence atlas](collatz_coincidence_atlas_20260927.md): `139 = 3^7 - 2^11`, `X_0(11)`; verdict NUMEROLOGY.
  - [level-22 oldform bridge](collatz_level22_oldform_bridge_20261005.md), section 4: Tate module of `X_0(11)`; Hecke data are commutative and lose route order (`F_2 F_1 (n) = (9n+5)/8` against `F_1 F_2 (n) = (9n+7)/8`).
  - [Pascal boundary note](collatz_pascal_boundary_leaf_section_20261005.md): only the GL(1) layer, HEURISTIC.
  - [S18 Artin corrections](collatz_artin_corrections_20260927.md): Chebotarev, entanglement, `A * 20/19`.
  - [fs_e8_finiteness](fs_e8_finiteness_20260917.md): local Langlands, not Collatz.
- **S18 and S19 decoded.**
  - S18 = [collatz_artin_corrections_20260927.md](collatz_artin_corrections_20260927.md): owner directive "think of Artin's conjecture and 19/20".
  - S19 = [mazur_positive_density_20260928.md](mazur_positive_density_20260928.md): Lech Mazur's AI-assisted, Lean-checked positive-density log-time Collatz convergence, with harmonic mass and negative cycles as 3-adic resonances.
- **Concurrent reading of arXiv:2610.06783.** The mac-mini [five-papers note](collatz_five_papers_synthesis_20261006.md), section 5, already transferred the paper's "wanted entries" theme: stopping-time records empty the last-dip residual windows below `2^33`.
- **Certificate culture already in the repo:** [LrcExtremalCert.lean](../../04-computation/lean/LrcExtremalCert.lean) (the WOWII-103 `native_decide` template) and the standalone Collatz descent certificates. The new piece here is root isolation of an irrational endpoint.
- **Function-field side** (Weil's Rosetta stone, the Langlands-adjacent habit): the `F_2[x]` Collatz map is solved; see the [sibling dimension ladder](collatz_procgen_20260922_sibling_dimension_ladder.md) and the HMYZ and ABP references there.

Portfolio.

- **Anchor:** the octic's arithmetic and its exact place in the Langlands ladder.
- **Niche:** reuse the formalization's root isolation on `log2 3`, giving a Lean certificate.
- **Wildcard:** the 3SUM themes against S18 and S19.

## 2. The octic field (PROVED + FINITE-EXACT)

**Link.** `Res_u(U(u), (1 + 2u - u^2) s - (6u + 4)) = -1792 P(s)`. So P is the minimal polynomial of `T`, and since both octics are irreducible, `Q(T) = Q(u)`. The real roots of U are `-0.54609...` and `u = 0.365769307604677293388545...`. With `u = tan(theta/2)`:

- `cos theta = (1-u^2)/(1+u^2)`;
- `sin theta = 2u/(1+u^2)`;
- `T = (6u+4)/(1+2u-u^2) = (2 + 2 cos theta + 3 sin theta)/(sin theta + cos theta)`;
- `theta = 40.18193729032971646...` degrees, the tilt of Trump's packing. So `cos theta`, `sin theta` and `T` all lie in `K`.

**Theorem (Gal = S8; PROVED).** The Galois group of P over Q is the full symmetric group S8.

*Proof.* P is irreducible, so G is transitive. P is unramified at 7 and 8017 and factors there as degrees `(7, 1)` and `(2, 1, 1, 1, 1, 1, 1)` respectively. By Dedekind, G contains a 7-cycle and a transposition. A transitive group containing an `(n-1)`-cycle is 2-transitive, hence primitive, and a primitive group containing a transposition is `S_n` (Jordan). QED.

Other first appearances: an 8-cycle at `p = 29`, type `(5, 3)` at `p = 23`. A cross-check found no block system among all 35 + 105 candidates (80-digit roots).

**Discriminants and ramification (FINITE-EXACT).**

- `disc P = -2^8 * 5 * 31 * 367^2 * 104053 * 614321^2`.
- `disc U = -2^32 * 5 * 31 * 104053`.
- Dedekind's criterion says `Z[T]` is p-maximal at `2, 5, 31, 104053` and not at `367, 614321`.

Hence `d_K = -2^8 * 5 * 31 * 104053 = -4128823040` and `[O_K : Z[T]] = 367 * 614321`.

| p | `P mod p` | ramification |
|---|---|---|
| 2 | `(x^2+x+1)^4` | `e = 4, f = 2` (wild); `T mod 2` is a primitive cube root of unity in `F_4` |
| 5 | `(lin)(lin)^2(quintic)` | one prime with `e = 2` |
| 31 | `(lin)(lin)^2(quintic)` | one prime with `e = 2` |
| 104053 | `(lin)^2(sextic)` | one prime with `e = 2` |

**Signature** `(2, 3)`: real roots `3.8770835900...` (the packing) and `-1.8530324790...`.

**Corollary (PROVED). Trump's tilt angle is neither a rational multiple of `pi` nor ruler-and-compass constructible, and neither is the optimal side `T`.**

- If `theta = r pi` with `r` rational, then `tan(theta/2)` would lie in a cyclotomic field. Its Galois closure would then be abelian, but `U` has Galois group S8.
- Constructible numbers have Galois closures of 2-power degree, and `|S8| = 40320` is not a power of 2.
- `T` generates the same field `K`.

**Quadratic resolvent.**

- The squarefree part of `disc P` is `D = -16128215 = -5 * 31 * 104053`, which is `1 mod 4`, so `Q(sqrt D)` has fundamental discriminant `D`.
- `h(D) = 5444`, counted from reduced forms. The analytic class number formula agrees: `sqrt|D|/pi * L(1, chi_D)` with the Euler product over `p < 2*10^6` gives `5444.9`, `L(1, chi_D) = 4.259`. The value is large because `2` and `3` both split (`D = 1 mod 8`, `D = 1 mod 3`).
- The Frobenius at `p` is even iff `(D/p) = +1`: 0 mismatches on 2258 primes.

**Chebotarev (FINITE-EXACT, `p < 2*10^5`, 17979 primes).** Selected factorisation types, observed against the S8 class proportion:

| type | observed | S8 |
|---|---|---|
| `(7,1)` | `0.1414` | `1/7 = 0.1429` |
| `(8)` | `0.1236` | `1/8` |
| `(5,2,1)` | `0.0963` | `1/10` |
| `(6,2)` | `0.0843` | `1/12` |
| `(4,3,1)` | `0.0842` | `1/12` |
| has a root | `0.6305` | `1 - 14833/40320 = 0.6321` |

**Trace moments.** `a_p = #{roots mod p} - 1 = chi_std(Frob_p)`. For `p < 10^6` (78492 primes) the moments are `-0.0028, 0.9977, 0.9751, 3.8919`, against `0, 1, 1, 4`.

## 3. The Langlands ladder (CITED + typed)

The permutation character of S8 on the roots is `1 + chi_std`, so

`zeta_K(s) = zeta(s) * L(s, rho_std)`,

and by the conductor-discriminant formula `cond(rho_std) = |d_K| = 4128823040`. The Galois closure has 22 irreducible representations (dimensions 1 to 90). What is known:

| layer | object | status |
|---|---|---|
| GL(1) | sign character = `chi_D` | **automorphic** (class field theory): `L(s, sgn) = L(s, chi_{-16128215})` |
| GL(7) | `rho_std` (conductor 4128823040) | Langlands reciprocity predicts a cuspidal `pi` on GL(7) with `a_p(pi) = #roots mod p - 1`. **Unknown.** Brauer gives meromorphy. Entireness of `zeta_K/zeta` (Dedekind's conjecture) is a theorem for Galois `K` (Aramata-Brauer) and for solvable closure (Uchida, van der Waall); S8 is neither. |
| other S8 irreducibles | dimensions 7..90 | unknown |

The ladder across the repo:

| object | field / representation | Langlands layer | status |
|---|---|---|---|
| Collatz exact arithmetic (cycle points in Q; Pell unit `3 + 2 sqrt 2`; base-phi Collatz in `Q(sqrt 5)`; transfer spectra mod `3^k` in cyclotomic fields; the cyclic quintics of the Brocard note) | abelian | GL(1) | known (class field theory); the [Pascal note](collatz_pascal_boundary_leaf_section_20261005.md) already says "only GL(1)" |
| `X_0(11)` ([oldform bridge](collatz_level22_oldform_bridge_20261005.md)) | 2-dimensional l-adic | GL(2) | known (modularity) |
| the 11-square octic | `S8` Artin, 7-dimensional | GL(7) | open |

*Check of the first row.* Searching the Collatz notes for Galois groups, splitting fields and minimal polynomials found only abelian ones. The census is VERIFIED by search, not exhaustive.

**Loss ledger.**

- **Map:** (packing constant) to (its Galois representation / Frobenius statistics).
- **Preserved:** the field `K`, conjugacy classes, L-functions, conductors.
- **Lost:** which of the eight conjugates (`3.877`, `-1.853`, three complex pairs), the generator (`T` and `u` give identical L-functions), and all geometry.
- **Sidecar that restores it:** the real embedding plus the root-isolation interval, i.e. exactly the formalization's `(9/25, 37/100)`.
- **Hostile control:** `s(5)` and `s(10)` share the single field `Q(sqrt 2)` while being different optima. The Galois data cannot even tell `n = 5` from `n = 10`.

This is the oldform bridge's lesson in a new setting. Commutative, conjugacy-invariant data forget order and geometry.

## 4. S18 revisited: entanglement, exactly (FINITE-EXACT)

Artin's base-5 correction `A * 20/19` arises because "`p = 1 mod 5`" implies that `p` splits in `Q(sqrt 5)`, and "`5` is a square mod `p`" says exactly that; and `Q(sqrt 5)` sits inside `Q(zeta_5)`. Its `19` is `q^2 - q - 1` at `q = 5`.

The octic has the same structure through `Q(sqrt D)`, the quadratic subfield of its Galois closure. It is not a subfield of `K`, which has no proper subfields because S7 is maximal in S8. Densities of `#roots of P mod p`, exact from S8:

| roots | all p | `(D/p) = +1` (even) | `(D/p) = -1` (odd) |
|---|---|---|---|
| 0 | 2119/5760 | 353/960 | 53/144 |
| 1 | 103/280 | 31/84 | 11/30 |
| 2 | 53/288 | 13/72 | 3/16 |
| 3 | 11/180 | 1/15 | 1/18 |
| 4 | 1/64 | 1/96 | 1/48 |
| 5 | 1/360 | 1/180 | 0 |
| 6 | 1/1440 | 0 | 1/720 |
| 8 | 1/40320 | 1/20160 | 0 |

The analogue of "`20/19`" for the event "P has a root mod p" given `(D/p) = +1` is `(12747/20160)/(25487/40320) = 3642/3641 = 1 + 1/3641`. (Even derangements number 7413 and odd ones 7420; the common factor 7 cancels, so it carries no meaning here.)

Principle: Chebotarev conditions entangle exactly through the intersection with the maximal abelian subextension, which for an S8 closure is the quadratic resolvent alone. `K^gal` meets `Q(zeta_m)` in `Q(sqrt D)` when `16128215 | m`, and in `Q` otherwise. So for `m` prime to the conductor, conditions on `p mod m` are independent of the Frobenius class, the root count included.

*Interpretation (not FINITE-EXACT).* S18's Collatz answer is the opposite extreme. Under population measures the 2-adic side shows no correction, and all observed entanglement is 3-adic, attached to one orbit's past. That suggests a Collatz "Artin correction" would have to be a 3-adic, orbit-coupled factor (THM-4506) rather than a population constant; nothing proves it.

## 5. The formalization's method, mapped onto Collatz

| packing proof | role | Collatz counterpart | transfers? |
|---|---|---|---|
| finite closed-cell cover of configuration space | compactness | 2-adic cylinders mod `2^k` | **no** for descent (Proposition 5.1) |
| certified positive slack on every non-optimal cell | quantitative gap | descent margin `A - l log2 3 > 0` on a first-descent segment | yes, segment by segment (S8-S10) |
| one rigid zero-slack configuration, an algebraic endpoint of degree 8 | the optimum | coefficient-survival set `E_inf` (`h* = 0.950`); its periodic points are the negative cycle points `x_w = c_w/(2^A - 3^l)` with `3^l > 2^A` (S6, S19 Theorem B), and it also contains preperiodic rationals such as `-11/3 -> -5` | **no**: uncountable, not rigid, and it ignores the carry |
| root isolation `u in (9/25, 37/100)` | exactness | isolation of `log2 3` between convergents | **yes** (section 5.1) |
| `native_decide` certificates, trust = kernel + native compiler | checking | repo template `LrcExtremalCert.lean`, standalone descent certificates | yes, already practised |

**Proposition 5.1 (PROVED, elementary).** No finite set of 2-adic cylinders, each certified to descend within its own bounded number of steps, covers all sufficiently large positive integers.

*Proof.* A finite family has a largest certified descent time `K`. Take `n = 2^(K+1) - 1`. Under `T(n) = (3n+1)/2` its first `K+1` iterates are `3^j 2^(K+1-j) - 1` (`j = 1..K+1`), all above `n`. So `n` does not descend within `K` steps, and no member of the family covers it. QED.

The packing proof's compactness has no descent analogue. The valid reformulations are those of the current [S7](collatz_measurement_independence_20261005.md) text: Collatz is equivalent to a computable deadline for every source, and to "every odd `n > 1` has an actual strict first descent".

- "No positive integer in `E_inf`" is necessary (positive integers in `E_inf` cannot descend), but not known to be sufficient. Coefficient exit ignores the carry: a nontrivial cycle's minimum has whole-cycle coefficient below 1, so it lies outside `E_inf`. (An earlier equivalence with `E_inf` was withdrawn in S7.)
- The packing proof's step "the zero-slack configuration is unique and algebraic" has no finite Collatz counterpart.

### 5.1 The transferred piece: a Lean certificate of the S9 expense records (LEAN-CHECKED)

S9's segment expense is `q(l, A) = ceil(l/(A - l log2 3))`; its worst type at length `l` is `A = A_min(l)`. S9 computed the record lengths by exact power comparisons. For `q` beyond `l = 94` it used 50-digit floats, exact "unless within `10^-50` of an integer".

The certificate decides all of these values from one isolation:

- **(B)** `2^16785921 < 3^10590737` and `3^190537 < 2^301994`, i.e. `16785921/10590737 < log2 3 < 301994/190537`. These are consecutive convergents; the width is `4.96e-13`.
- **(A)** For each `l`: `2^(A-1) < 3^l < 2^A`.
- **(Q)** `f(t) = l/(A - l t)` is increasing for `t < A/l`. With `dlo = A*QLO - l*PLO` and `dhi = A*QHI - l*PHI` (both positive, checked), the inequalities `(q-1) dlo <= l QLO` and `l QHI <= q dhi` give `q - 1 < l/(A - l log2 3) < q`.

Glue: (B) and (Q) imply `q = q_max(l)`, by monotonicity. This is PROVED in the comment and not formalized. Everything decidable is in Lean: `bracket_lo`, `bracket_hi`, `all_lengths_certified` (all `l = 1..320`), `records_320`, and a theta-free cross-check `direct_small_records`, which proves `q = least m with 3^(ml) < 2^(mA-l)` for the first six records. The check takes 12 s. The axiom audit shows only `X._native.native_decide.ax_1_1` per theorem.

Certified records `(l, A, q_max)`: `(1,2,3), (3,5,13), (5,8,67), (17,27,306), (29,46,804), (41,65,2480), (94,149,6951), (147,233,13984), (200,317,26668), (253,401,56382), (306,485,207489)`. This is S9's table exactly.

**Relation to classical cycle bounds (CITED, corrected after audit).** Eliahou (1993): a nontrivial cycle with minimum above `2^40` has length `301994a + 17087915b + 85137581c` (`b >= 1`, `ac = 0`). His generators are the numerators of the convergents `c13 = 301994/190537`, `c15 = 17087915/10781274` (the mediant of this bracket) and `c16 = 85137581/53715833` (below `log2 3`). `16785921 = c14` is not among them.

This bracket cannot reproduce his theorem. At `(l, A) = (190537, 301994)` it only gives `0 < A - l log2 3 < 9.44e-8`, while excluding that shape for minima above `2^40` needs a certified lower bound above `8.33e-8` (the true value is `9.31e-8`). Certifying Eliahou's Diophantine step needs a sharper bracket built from `c15` and `c16`.

## 6. arXiv:2610.06783 (Alman-Vassilevska Williams) against the Collatz work

The paper: truly subquadratic 3SUM (`O(n^1.9992)`) and subcubic APSP (`O(n^2.9995)`).

- It computes sparse wanted entries of thin products `X (N x D) * Y (D x N)` by recursing on Schonhage's 10-term identity (`3x1 * 1x3` outer product plus `1x4 * 4x1` inner product).
- It prunes recursion leaves that feed no wanted entry.
- It shares the encoding across many products (`N >= D^18`), for a `D^(1/18)` saving.
- Its core algorithm was discovered by Claude.

| theme | Collatz counterpart | type |
|---|---|---|
| compute only wanted entries | stopping-time records empty the residual windows ([five papers](collatz_five_papers_synthesis_20261006.md), Prop. 5.1, a proved Collatz statement) | ANALOGY (no product or identity in Collatz) |
| an algebraic identity with one contributing term per output enables pruning | the two-sheet fusion relation (one relation for the 3x-1-rooted multipliers) and the carry cocycle | ANALOGY (no saving shown) |
| amortized encoding (`N >= D^18`) | banks of certified excursions buy coverage, not deadline (S8/S9) | ANALOGY |
| lopsided thin products (`D <= N^0.12`) | sources are 2-adic, landings are 3-adic (S15, S18): two short lists against one long hovering list (five-papers section 5) | ANALOGY |
| refuting standing hypotheses | the repo's refutation culture (S7 measurement independence, MISTAKE ledger) | meta |
| AI discovery, then human simplification and verification | S19 (Mazur: AI-assisted, Lean-checked); the 11-square formalization; this note's Lean certificate | meta |
| a finite or algorithmic reduction behind each success | packing: compactness; 3SUM: an algorithm; Mazur: positive density with explicit constants; Collatz: no finite reduction (Proposition 5.1) | the dividing line |

**On 18 and 19.**

- **S18** was the owner's "Artin's conjecture and 19/20" directive. `20/19 = 1 + 1/(5^2 - 5 - 1)` is the entanglement correction for base 5 (section 4).
- **S19** was the Mazur digest.
- **The paper's 18** is the amortization exponent (`N >= D^18`, saving `D^(1/18)`; `eps < 0.1204`). As a number it has no relation to S18 or S19 (NUMEROLOGY).
- **As a theme it resembles the S18 directive** ("correction factors that may unlock recursion"). The paper's gain is a polynomial saving bought by pruning a recursion to the leaves that matter. That is a different kind of object from an Artin density factor, so this is a loose ANALOGY.
- **Where a Collatz version could live (heuristic).** S18 found the 2-adic side exactly Haar under population measures, so a pruning that lowers THM-4499's exponent would presumably have to be orbit-coupled and 3-adic (THM-4506). No theorem forces this. Cheapest decisive test (open): does the landing tree (THM-4506) admit a Schonhage-type sparsity, i.e. one ancestor class contributing to each wanted landing class mod `3^k`?

## 7. Verification, hostile controls, open threads

**Independent paths.**

- `Gal = S8` by an explicit-prime proof and by Chebotarev frequencies on 17979 primes;
- `d_K` by Dedekind's criterion; `disc U` independently confirms only its odd part (`disc U / d_K = 2^24`); the audit's Eisenstein check confirms `v_2(d_K) = 8`;
- the sign character by Kronecker symbols against Frobenius parity;
- the records by Python bracket arithmetic, by S9's table, by Lean, and by the audit's 200-digit recomputation;
- the trace moments against exact S8 values.

**Hostile controls.**

- `s(5)` and `s(10)` share `Q(sqrt 2)`: the Galois data cannot separate different optima.
- `T` and `u` give the same field, so L-functions cannot see the generator.
- The bracket's failure condition (a ceiling range wider than 1) never fires for `l <= 320`.
- The theta-free direct check agrees on six records.

**Numerology guards.**

- The paper's 18 versus S18 is NUMEROLOGY.
- "11" squares versus the conductor 11 of `X_0(11)` (S8 crossings) is NUMEROLOGY.
- The 19 of `20/19` versus 19 in S19 is NUMEROLOGY.

**Audit (2026-10-06, independent code).**

- *Confirmed:* every computation (irreducibility, resultant, `Gal = S8` primes, discriminants, Dedekind criterion, `d_K`, class number by two Dirichlet-sum formulas, Frobenius parity on 17979 primes, the S8 table, Chebotarev and moment figures, the Lean file including rejection of false variants, a 200-digit recomputation of the 11 records), Proposition 5.1, the Corollary, and every citation.
- *Corrected (MISTAKE-566):*
  - the withdrawn `E_inf` equivalence restated as fact;
  - `E_inf`'s rational points misdescribed;
  - the Eliahou link overclaimed;
  - a cosmetic "7" in the entanglement factor;
  - `Q(sqrt D)` called a subfield of `K`;
  - an overstated float caveat in S9;
  - several over-typed connections in sections 0, 4 and 6;
  - Proposition 5.1 needed "sufficiently large".

**Open threads.**

1. Extend the Lean certificate to the Diophantine core of Eliahou's theorem (all `(l, A)` with `0 < A - l log2 3 < delta`, `l <= 10^7`), using a sharper bracket from the convergents `c15`, `c16`.
2. A Collatz-derived polynomial with non-abelian Galois group. Candidates: generalized Collatz maps on rings of integers of non-Galois fields. None is known here, and the Langlands ladder predicts nothing beyond GL(1) without one.
3. The Schonhage-sparsity test of section 6 on THM-4506's landing tree.
4. A polredabs-type small model of `K`, and its class group (PARI was not available in this session).

## 8. Reproduction

```bash
python3 04-computation/experiments/square11_octic_field_20261006.py 200000
python3 04-computation/experiments/square11_octic_field_part2_20261006.py
python3 04-computation/experiments/square11_octic_field_part3_20261006.py 1000000
python3 04-computation/experiments/collatz_expense_records_exact_20261006.py 320
lean 04-computation/lean/standalone/collatz_expense_records_certificate_20261006.lean
```

Requires `sympy`, `mpmath`, `python-flint` (0.9), Lean 4.30 core. Each Python script runs in under 5 s; Lean takes 12 s.
