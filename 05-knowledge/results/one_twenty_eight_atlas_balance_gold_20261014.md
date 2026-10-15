# 128 across the corpus, three levels down: the 2-adic Mersenne law of a level-32 newform, perfect QR codes, a tournament 1/3–2/3 theorem, the gold partition, Mahler's middle third and Borsuk in dimension four

2026-10-14, session opus-2026-10-14-S35 (worktree `/tmp/math-wt-fractal-recursion-23`). Companion to [paley_least_resistance_runs_one_modular_form_20261014.md](paley_least_resistance_runs_one_modular_form_20261014.md).

**The owner's message (verbatim).** "look back for other work in different session regarding 128 and deeply explore up to a few levels of recursion any related themes those sessions touched."

**Pasted material** that came with it:
- a request to "merge in a genuine attempt at also proving the 1/3 -2/3 conjecture and peczarski's gold partition conjecture";
- the conjectures.io Borsuk-four page;
- the instruction to treat the hexagonmath.org works as solid and "create bold ideas instead of reverifying computations";
- an older golden-lane snippet: period-5 rational cycles 0, −5 and 1/13; BS(1,3) mod 11; the temperaments 5, 12, 41, 53; the apotome and the −17 cycle; periods {1, 2, 3, 11}; `3⁵ − 1 = 2·11²` and the ternary Golay code.

**Status.**
- **PROVED:**
  - the tournament component reduction and the strong-majority lemma (HYP-9277);
  - the equality-class count `A000930(n) − 1`;
  - planar Fourier duality and the Jacobi-sum wheels (companion note).
- **FINITE-EXACT:**
  - the square law of the `K_6 − P_6` terminal at every prime `q ≡ 3 (mod 4)` below 200,000;
  - the 2-adic depth table, the Mersenne values and the valuation law;
  - the tournament balance sweep to `n = 7` [N8];
  - the closed-family counterexamples;
  - the Collatz sieve count.
- **CITED:**
  - McCarthy–Papanikolas (2015);
  - Akiyama–Frougny–Sakarovitch (2008);
  - Flatto–Lagarias–Pollington (1995);
  - Dubickas;
  - Brightwell–Felsner–Trotter (1995);
  - Kahn–Saks (1984);
  - Peczarski (2006, 2008);
  - Pless and MacWilliams–Sloane (QR codes);
  - Kahn–Kalai (1993), Bondarenko (2014), Jenrich–Brouwer (2014).
- **UNREFEREED, treated as solid as instructed:**
  - arXiv:2607.23926 (balance constants through 14 elements);
  - Grinsztajn (2026, 63-dimensional Borsuk counterexample);
  - Tolmachev–Voronov (2026, `b(4) <= 8`);
  - hexagonmath.org 2610.00249 (Collatz positive lower density) and 2610.00225 (BMV lattice method).
- **NUMEROLOGY:** typed in place.
- **Not yet audited** (section 9).

## 0. Answers in brief

1. **Every 128 in the corpus is `2^7` read in one of three ways:**
   - as a 2-adic level, the order of 3 modulo `2^9`;
   - as the Boolean cube `F_2^7` on the seven Fano points;
   - as `127 + 1`, a Mersenne prime plus one.

   This session adds two more:
   - the `K_6 − P_6` terminal is `≡ 1 (mod 128)` at Paley primes;
   - exactly 128 classes of 14-element posets attain the balance constant `1/3`, i.e. `A000930(14) − 1`, Narayana's cows.
2. **Level 2: the Mersenne law.**
   - The `mod 128` signature comes from a congruence `|b(q)| ≡ ±(q+1) (mod 2^m)` of the weight-3 level-32 newform `b` (mac-mini's HYP-9276).
   - Followed to Mersenne primes, where `q + 1 = 2^p` and `q ≡ −1` 2-adically, the newform's coefficient is 2-adically tiny: `b(127) = 0`, `b(8191) = 0`, `|b(131071)| = 2^14`, `|b(524287)| = 63·2^13`.
   - Over all primes `q ≡ 3 (mod 4)` below 200,000, `v_2(b(q)) >= min(v_2(q+1), [K])`.
   - This is the 2-adic continuity of Frobenius traces near `−1`, in the same topology where `⟨3⟩` is dense in `{1, 3 mod 8}` (the Collatz theme).
3. **Level 2: the codes.**
   - The three perfect nontrivial codes are the quadratic-residue codes of the Paley primes `7, 11, 23`: Hamming `[7,4,3]`, ternary Golay `[11,6,5]` and binary Golay `[23,12,7]`.
   - `128 = 2^7 = 16·8` is the Hamming code's perfect packing. E8, built from the extended Hamming code, has `240 = 112 + 128` roots, and `128` is the half-spin representation of `Spin(16)` in `e8 = so(16) ⊕ 128`. That is the honest Bott-periodicity content of 128.
   - The sphere-packing identities contain the Paley affine orders `C(q,2) = T_{q−1}`: `3^5 = 1 + 2·11 + 4·55` with `55 = |Aff(P_11)| = T_10 = F_10`.
4. **Level 3: the 1/3–2/3 conjecture, attempted through tournaments** (HYP-9277).
   - Every non-transitive tournament on `n <= 7` vertices [N8] has a pair ordered with probability in `[1/3, 2/3]` in a uniform Hamiltonian path.
   - Equality holds exactly at transitive compositions of singletons and cyclic triangles, the exact analogue of the 128 poset classes.
   - The proof reduces, by proved lemmas, to strong tournaments.
   - Peczarski's gold partition property **fails** once precedence answers are recorded. So the gold partition conjecture, if true, is a convexity phenomenon of order polytopes, not a property of "sorting" in general.
5. **Mahler's middle third.**
   - Akiyama–Frougny–Sakarovitch: countably many `z > 0` have `{z(3/2)^n}` in `[0, 1/3) ∪ [2/3, 1)` for all `n`, avoiding the middle third forever.
   - The 1/3–2/3 conjecture says linear-extension measures always hit the middle third.
   - Flatto–Lagarias–Pollington's spread bound for `3/2` is `1/p = 1/3`.
   - The common source is the three-way branch: three digits for `×3/2`, three linear extensions of `2+1`, three Hamiltonian paths of `C_3` (ANALOGY).
6. **Borsuk in dimension 4.**
   - Open. Counterexamples exist from dimension 63 or 64 up, all built from strongly regular graphs.
   - The hexagonmath BMV method turns SRG existence into a discriminant-form condition, a mod-8 Milgram constraint. So "how low can SRG counterexamples go" is a mod-8 lattice question.
   - Nothing here touches dimension 4 itself.

## 1. Level 0: the inventory

A read-only search of the corpus (1,258 raw hits for `\b128\b`, most of them run sizes) left six substantive themes.

| theme | what 128 is | representative sources |
|---|---|---|
| A. Collatz exponent phases | `ord_{2^9}(3) = 2^7`: `3^k mod 2^b` has period `2^(b−2)` | `paid_guard_cover_tree_20261005.md`, `mersenne_joint_guards_20261014c.md`, `mersenne_switch_parity_f21_compression_20261006.md` (S18: "`128 = 2^(K−2)` at `K = 9`"), `collatz_five_mirrors_20260929.md` (`3^128 ≡ 1 mod 2^9`), HYP-9175 |
| B. Collatz residues mod `2^7` | the 7-bit parity-vector sieve | `collatz_boundary_compiler_20261004.md` (`7 mod 128` certified), `collatz_guards_20260921_valves.md`, THM-4601, THM-4257 (`123 = −5 mod 128`) |
| C. Paley and Mersenne 127 | `2^7 ≡ 1 (mod 127)`, so `±⟨2⟩` has 14 elements | THM-4529 (all-odd 14-vertex tournament from `QR_127`), S33 (`H(P_31)/465 ≡ 41 mod 256`) |
| D. LRC14 | `F_2^7` over the 7 Fano/`χ_7` labels | THM-2504, THM-1247, THM-2591, HYP-3581, HYP-3527 (Kraft masses `/128`) |
| E. codes and lattices | Fano sign choices `2^7 = 16 + 112`; E8's 128 half-integer roots | `arithmetic_braids2_20260917_fano_code.md`, `collatz_mod6_20260922_paley_fano_octonion_design.md` |
| F. AMM 12592 | `C(127, j)` all odd (Lucas, `127 = 1111111₂`) | `amm12592-r128-lp-parity-rounding-boxeph.md`, THM-3028 |

**This session** adds:
- `Y_q ≡ 1 (mod 128)` for the `K_6 − P_6` terminal;
- the 128 extremal 14-element posets.

## 2. Level 1: one mechanism per theme

- **A, B.** `Z_2^× = {±1} × (1 + 4Z_2)`, and `1 + 4Z_2` is pro-cyclic, generated by 5. The closure of `⟨3⟩` is the index-2 subgroup `{x ≡ 1, 3 (mod 8)}` (HYP-9175), so `3` has order `2^(b−2)` modulo `2^b`. The mod-8 quotient `Z_2^×/(Z_2^×)² ≅ (Z/8)^×` is the 2-adic shadow of everything "mod 8" in the corpus. A residue class modulo `2^k` is the same datum as a parity vector of length `k` (Terras–Everett), so theme B is theme A's sieve.
- **C, F.** `127 = 2^7 − 1` is `−1` in `Z/128`, and multiplicatively `2^7 ≡ 1 (mod 127)`. The first reading is additive and gives Lucas' theorem, with all binomials odd. The second is multiplicative and gives `±⟨2⟩` of order 14.
- **D, E.** `F_2^7` is the space of functions on the Fano plane. The Hamming code (16 words, 8-point balls) tiles it perfectly: `16·8 = 128`. Sign gauges on the Fano lines form a Hamming coset, and the 16 octonion tables are the code.

## 3. Level 2: recursion into each theme

**(a) The `K_6 − P_6` signature, down to the newform** (companion note §4; mac-mini's HYP-9276 for the identification).
- At Paley primes, `Y_q − 1 = (|b(q)| − (q+1))(|b(q)| + (q+1))`, with `b` the level-32 weight-3 newform.
- The data give `|b(q)| ≡ ±(q+1) (mod 2^m)` with `m = 4, 5, 6, 6` for `q ≡ 11, 3, 7, 15 (mod 16)`. Hence `v_2(Y_q − 1) >= 7, 8, 10, 11`, with the minima attained.
- The square law and these minima hold at all 7,373 primes `q ≡ 3 (mod 4)` in `[30,000, 200,000]` as well (`paley_k6p6_two_adic_mersenne_20261014.out`).

**(a′) Mersenne primes (level 3, FINITE-EXACT).** At `q = 2^p − 1`, `q + 1 = 2^p` and `q` is 2-adically next to `−1`:

| `p` | `q` | `|b(q)|` | `M_q = Y_q − 1` |
|---|---|---|---|
| 3 | 7 | `8 = 2^3` | 0 |
| 5 | 31 | `32 = 2^5` | 0 |
| 7 | 127 | **0** | `−2^14` |
| 13 | 8191 | **0** | `−2^26` |
| 17 | 131071 | `2^14` | `−63·2^28` |
| 19 | 524287 | `63·2^13` | `−127·2^26` |

- Two of the six Mersenne coefficients vanish. Below 700 the newform has no other zero.
- The table suggests a **valuation law**: `v_2(b(q)) >= min(v_2(q+1), K)`. Its test over all primes below 200,000 (`paley_k6p6_two_adic_valuation_law_20261014.out`): [VALLAW]
- The reading is that `q ↦ b(q)` behaves like a 2-adically continuous function of `q ∈ Z_2^×`, carrying an Eisenstein factor `(1 + q)` up to a bounded depth. Mersenne primes are where that factor is 2-adically smallest. That is the same 2-adic topology in which theme A's `⟨3⟩` is dense in `{1, 3 mod 8}`. Themes C and F, the prime `127`, enter the arithmetic of a modular form exactly through `127 ≡ −1 (mod 2^7)`.
- The mechanism is proposed, not proved. McCarthy–Papanikolas' trace-formula expression for `b(p)` runs over Legendre curves with `Z/2 × Z/8` torsion, i.e. `16 | #E(F_p)`.

**(b) Perfect codes are the QR codes of Paley primes (CITED; the Paley tie-in is exact).**
- The Hamming `[7,4,3]`, ternary Golay `[11,6,5]` and binary Golay `[23,12,7]` codes are the quadratic-residue codes of lengths 7, 11, 23.
- Their defining sets are the out-neighbourhoods of the Paley tournaments `P_7`, `P_11` and `P_23`, the squares modulo `q`.
- Their sphere-packing identities are:
  - `2^3 = 1 + 7`;
  - `3^5 = 1 + 2·11 + 4·C(11,2)`;
  - `2^11 = 1 + 23 + C(23,2) + C(23,3)`.

  Here `C(q,2) = T_{q−1} = |Aff(P_q)|`: 55 for `q = 11` (the owner's triangular-Fibonacci `55 = T_10 = F_10`) and 253 for `q = 23`.
- `P_3`, `P_7` and `P_11` are the orders where `T_{q−1}` is Fibonacci and Paley maximises `H` (companion note §5).
- E8 is Construction A of the extended Hamming code. Its 240 roots split as `112 + 128`. The Leech lattice comes likewise from the extended binary Golay code.
- The golden-lane snippet's `3^5 − 1 = 2·11²` is the ternary Golay identity rewritten: `1 + 22 + 220 = 243`. Its binary partner is Ramanujan–Nagell, `2^7 = 11² + 7`. That is the sphere-packing equation `(2n+1)² + 7 = 2^(r+3)` of 2-error-correcting binary perfect codes at `n = 5`, the repetition code. Both 11's are NUMEROLOGY, but each equation is exact.

**(c) E8's 128 and tournament scores (THM-868, REFEREED there).**
- THM-868 places the score half-vectors `s_i − 7/2` of 8-vertex tournaments in E8's half-integer coset. Its near-regular shell is exactly the `C(8,4) = 70` roots with coordinate sum 0.
- Those 70 lie among the 128 half-integer roots `(±1/2)^8` (even number of minus signs): they are the ones with exactly four minus signs.
- THM-868 has a further 128 that the search missed. Landau's condition leaves every shell `x = Σd² <= 128` untouched and first bites at `x = 136`.
- Bott: half-spin representations of `Spin(8m)` have dimension `2^(4m−1)`, which is 8 (triality) and then 128 (`e8 = so(16) ⊕ 128`). The real Clifford algebras repeat with period 8, `Cl_(n+8) ≅ Cl_n ⊗ M_16(R)`. This is where 128 and "mod 8" genuinely meet, and the owner's "Bott periodicity" points there (ANALOGY for the tournament side).

**(d) The Collatz sieve passes through 128 (FINITE-EXACT; the apotome link is NUMEROLOGY).**
- The number of residue classes mod `2^k` whose `T`-trajectory has not dropped below its start within `k` steps, i.e. parity vectors with `3^(o_j) > 2^j` for all `j <= k`, is `1, 1, 2, 3, 4, 8, 13, 19, 38, 64, 128, 226, 367, 734, …` for `k = 1, 2, …`.
- It is 128 at `k = 11`. At that step `3^7 = 2187` first exceeds `2^11 = 2048`: the apotome, the ratio of the snippet's `−17` cycle (7 odd steps, 11 halvings, denominator `2^11 − 3^7 = −139`).
- The limma `2^8/3^5 = 256/243` gives the `1/13` cycle (`2^8 − 3^5 = 13`), and `2^3/3^2` gives the `−5` cycle.
- The temperaments 5, 12, 41, 53 are convergent and semiconvergent denominators of `log_2 3`. So is the sieve's doubling pattern (doublings at `k = 6, 9, 11, 14, 17, 19, …`): least-resistance lattice paths above the line of slope `log_3 2`.

**(e) Narayana's 128 and the golden ratios.**
- Compositions into parts 1 and 3 grow like `ψ^n`, with `ψ = 1.46557…` the supergolden ratio (`x³ = x² + 1`, a Pisot number).
- The gold partition conjecture's `t_0 >= t_1 + t_2` is a Fibonacci recursion, giving sorting in `log_φ e(P)` comparisons.
- Brightwell–Felsner–Trotter's bound `δ(P) >= (5 − √5)/10 = 1/(1 + φ²) = 0.2764` is best possible for infinite posets.
- So the golden ratio governs the infinite and sorting sides, and the supergolden ratio counts the finite extremal posets. Mahler's `3/2` is the archetypal non-Pisot base, where fractional parts should not converge.

## 4. Level 3: a genuine attempt at 1/3–2/3 and the gold partition (HYP-9277)

**What a counterexample must look like (PROVED, for any measure on linear orders).** If no pair is balanced, `x ≻ y ⟺ p(x,y) > 2/3` is a linear order `σ*`, since `p(x,z) >= p(x,y) + p(y,z) − 1`.
- For posets, `σ*` is a linear extension.
- A majority 3-cycle always contains a balanced pair. That ties the 1/3–2/3 problem to arXiv:2607.23926's "linear-extension majority cycles": their longest is 8 at order 14, and a counterexample can have no majority 3-cycle.

**Both problems reduce to indecomposables (PROVED).**
- Posets reduce to ordinally indecomposable posets, since `δ(P ⊕ Q) = max`.
- Tournaments reduce to strong tournaments, since `δ(T) = max` over components.
- The extremal blocks are `2 + 1` and `C_3`, each with exactly 3 orders.

**Tournaments (FINITE-EXACT).**
- `δ >= 1/3` for every non-transitive labelled tournament with `n <= 7` [N8].
- Equality holds exactly at the predicted compositions, with labelled counts 2, 16, 120, 1040, 10080.
- The least strong balance is 2/5, then 4/11 (`n = 5, 6, 7`). The poset analogues through 14 elements are 37/106 and 254/725.

**The gold partition does not survive generalisation (FINITE-EXACT counterexamples).**
- Record precedence answers, so that the surviving set is `HP(T) ∩ Ext(Q)`. Then the balance can drop to `1/4` (tournaments at `n = 6`; semicomplete digraphs at `n = 4`), and `t_0 >= t_1 + t_2` fails for 176 family members at `n = 6`.
- So the poset conjectures are not instances of a general "sorting" principle. They need the convexity of linear-extension sets, the order polytope with its Stanley/Kahn–Saks log-concavity.
- **Bold reading.** The gold partition conjecture should be attacked as a statement about mixed volumes of order polytopes, a "golden Grünbaum theorem" for braid-hyperplane cuts. Kahn–Linial already got `1/(2e)` from Grünbaum's centroid theorem. Pure combinatorics of permutation sets cannot prove it.

**Where the attempt stops.** Proving `δ > 1/3` for strong tournaments on `>= 4` vertices would make HYP-9277 a theorem. The strong-majority lemma turns it into the claim that no strong tournament has a linear order `σ*` with all earlier-before-later probabilities `> 2/3`. That claim is open here.

## 5. Mahler's middle third (CITED + ANALOGY)

- **Akiyama–Frougny–Sakarovitch** (Israel J. Math. 168, 2008, Corollary 4): the set of `z > 0` with `{z(3/2)^n} ∈ [0, 1/3) ∪ [2/3, 1)` for all `n >= 0` is countably infinite. These are the reals with two `3/2`-expansions.
- **Dubickas**: `||x(3/2)^n|| < 0.238117…` for all `n` forces `x = 0`.
- **Flatto–Lagarias–Pollington**: the spread of `{ξ(p/q)^n}` is at least `1/p`, which is `1/3` for `3/2`.
- **Mahler's own question** (Z-numbers, `[0, 1/2)`) sits between the Dubickas bound and the FLP bound.

The 1/3–2/3 conjecture is the opposite assertion for a different dynamics: linear-extension measures cannot avoid the closed middle third. In all three settings the `1/3` is the reciprocal of the smallest nontrivial branching: three digits, three extensions of `2+1`, three paths of `C_3` (ANALOGY; no transfer of proofs is claimed).

## 6. Borsuk's conjecture in dimension 4 (status only)

- The conjectures.io page lists dimension 4 as the smallest open case. Its references are:
  - Kahn–Kalai 1993;
  - Bondarenko 2014 (dimension 65, two-distance sets from the `G_2(4)` strongly regular graph);
  - Jenrich–Brouwer 2014 (dimension 64);
  - Grinsztajn 2026 (dimension 63, UNREFEREED);
  - Tolmachev–Voronov 2026 (`b(4) <= 8`, UNREFEREED, improving Lassak's 9).
- All known counterexamples are spherical two-distance sets from strongly regular graphs. Their dimension is an eigenvalue multiplicity.
- hexagonmath's BMV paper (2610.00225, treated as solid) gives a uniform discriminant-form condition for SRGs with integral eigenvalues. It is a Milgram/Gauss-sum constraint, the lattice-theoretic "mod 8", and it rules out 70 open parameter sets.
- **Bold idea:** the smallest SRG-based Borsuk counterexample dimension is a question about which eigenvalue multiplicities survive the mod-8 discriminant condition. That is the same mod-8 that ties 128 to E8 and `Spin(16)`.
- No statement about dimension 4 is made. Two-distance sets in `R^4` are far too small to be counterexamples.

## 7. The hexagonmath works as inputs (UNREFEREED; treated as solid per the owner)

- **2610.00249 (Collatz).** Claims a positive lower density of `n < X` reaching 1 within `⌊10.46 log n⌋` steps, with a tower-type constant. If solid, this is the density statement that this repository's renewal and deadline notes (S8–S10: density zero for every excess-rate family) could not reach. Its step constant sits between the typical `6.95 log n` and the S8 deadline `1.051 q log₂ n` families.
- **2610.00225 (BMV).** See section 6. It is also the natural home for Paley graphs at square `q`, where eigenvalues are integral.
- **2610.00248 (kissing numbers in 19 and 21 dimensions).** Context for section 3(b): E8's kissing configuration is `112 + 128`.
- **2610.00236 (Haglund's monotone zeros fail at `N = 27`).** No contact with this note.

## 8. Hostile controls

1. The square law and the 2-adic tables are FINITE-EXACT through 200,000. The identification is mac-mini's HYP-9276. The congruence mechanism is a proposal.
2. The Mersenne table is exact, but six primes are a small sample, and the valuation law is tested only below 200,000.
3. The QR-code facts are classical. The appearances of 11, 55 and 253 inside sphere-packing identities are exact identities whose link to Paley Hamiltonian-path counts is NUMEROLOGY.
4. HYP-9277 is FINITE-EXACT to `n = 7` [N8]. The analogy table with arXiv:2607.23926 uses an unrefereed source.
5. The gold-partition counterexamples refute only the generalised (closed-family) statement, not Peczarski's conjecture for posets.
6. Nothing here proves 1/3–2/3, the gold partition conjecture, Borsuk-4 or Mahler's conjecture.

## 9. Audit

Pending.

## 10. Reproduction and references

**Scripts** (`04-computation/experiments/`, outputs in `05-knowledge/results/`):
- `paley_k6p6_two_adic_mersenne_20261014.py`;
- `paley_k6p6_two_adic_valuation_law_20261014.py` (writes `paley_k6p6_cq_table_20261014.json`);
- `tournament_one_third_balance_20261014.py`;
- `tournament_balance_labelled_sweep_20261014.py`;
- `tournament_balance_spectrum_gold_20261014.py`.

The Collatz sieve count is a ten-line enumeration of parity vectors, reproduced in section 3(d)'s description.

**Literature:**
- D. McCarthy and M. A. Papanikolas, arXiv:1205.1006;
- S. Akiyama, Ch. Frougny and J. Sakarovitch, *Israel J. Math.* 168 (2008);
- L. Flatto, J. C. Lagarias and A. D. Pollington, *Acta Arith.* 70 (1995);
- A. Dubickas (on `||x(3/2)^n||`);
- G. Brightwell, S. Felsner and W. T. Trotter, *Order* 12 (1995);
- J. Kahn and M. Saks, *Order* 1 (1984);
- J. Kahn and N. Linial, *Combinatorica* 11 (1991);
- M. Peczarski, *Order* 23 (2006) and 25 (2008);
- arXiv:2410.12494;
- arXiv:2607.23926;
- F. J. MacWilliams and N. J. A. Sloane, *The Theory of Error-Correcting Codes*, ch. 16 (QR codes);
- J. H. Conway and N. J. A. Sloane, *SPLAG* (E8, Leech);
- J. Kahn and G. Kalai, *Bull. AMS* 29 (1993);
- A. Bondarenko, *DCG* 51 (2014);
- T. Jenrich and A. E. Brouwer, *EJC* 21 (2014);
- conjectures.io (Borsuk-four page);
- hexagonmath.org (2610.00249, 2610.00225, 2610.00248, 2610.00236);
- OEIS A000930 and A076227.

**This repository:** HYP-9175, HYP-9273, HYP-9276 and THM-4653 (mac-mini), HYP-9277, HYP-9278, THM-868, THM-4529, THM-3028, and the theme sources in section 1.
