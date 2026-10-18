# Binary information, class groups and largest anomalies: idoneal 1848, the class-number ranking of Paley tournaments, the exponent ladder, and the cutoffs at sixteen

2026-10-17, opus-2026-10-17-S40. Owner prompt: imaginary quadratic class groups and imaginary primes; the 65 idoneal
numbers as a model for proving that other "large number anomalies" do not exist; 1848 as the largest, its sixteen-fold
self-similarity, cutoffs at powers of two and the Cayley-Dickson tower; merged versus solo classes of the tournament
metagraph; past work on 463; binary relational information.

**Typing.**
- **PROVED (elementary):**
  - THM-4666: class number = ranking excess; arithmetic orders; local optimality of the natural Paley order = Kalmynin's
    L+; a real zero forces an improving move; randomised moves = Fekete positivity.
  - The involution dictionary (section 3).
  - The binary-sufficiency reading (section 2, from classical CITED theorems).
  - Euler's criterion for 463 (section 7).
- **FINITE-EXACT:**
  - the identity checks of THM-4666 and exact MAS for `p <= 31`;
  - the exponent ladder to `4·10^6` (section 5), which reproduces Euler's 65 idoneal numbers and the
    Elsenhans-Klüners-Nicolae lists;
  - the idoneal class-number histogram.
- **REFUTED in its threshold form:** HYP-9287. Two pre-registered tests passed as executed, but certified rotation orders
  beat the natural order at high-L primes; a revised conjecture is stated.
- **CITED:** Weinberger, Kani, Heegner-Baker-Stark, Rabinowitsch, Watson, Lorch-Kirschmer, Elsenhans-Klüners-Nicolae,
  Siksek, Ono-Soundararajan, Fekete-Polya, Kalmynin (via Baker-Montgomery), Joshi, Chindelevitch-Harutyunyan.
- **ANALOGY / NUMEROLOGY:** the sixteen atlas and the 463 atlas (verdicts per row).

**Audit.** Independently audited 2026-10-17 (`04-computation/experiments/s40_audit_20261017/REPORT.md`, MISTAKE-626).
All computations reproduced. Text fixes applied: Weinberger/Kani, Gauss's cubic and quartic criteria, the repository's
adopted idoneal status, EKN coverage, credits.

**Repository prior art used here.**
- clock_duality_oct15's [escape_renewal_idoneal_20261017d.md](escape_renewal_idoneal_20261017d.md) section 7 already has:
  - the orbit law `(|S| + |Fix|)/2` for converse and for class inversion;
  - the sixteen rotational classes at `n = 17` set against `Cl(-7392)`;
  - the Cayley-Dickson zero-divisor contrast.
- [prime463_conductor_clock_20261017e.md](prime463_conductor_clock_20261017e.md) has the 1848 tower and both
  class-number identities at 463.
- The merged count `(A000568 + A002785)/2` is in [trees5_tournaments4_converse_20261005.md](trees5_tournaments4_converse_20261005.md).
- THM-871's self-converse sentence is clock_duality's scope repair (90913a62ff).
- The exponent-2 row (101 discriminants, 65 fundamental, 5460, 7392) is in THM-4566's oai3 reader (`exp2_analysis.out`).
- The Paley half-row imbalance (statement 2 of THM-4666) appears in ogg_triangular_chowla_20260927.md section 2 and
  special_primes_curling_genus_zero_langlands_20261015.md section 4.

## 0. Answers in brief

1. **"Only 65 idoneal numbers": two layers.**
   - *Classically,* Weinberger (1973) proved Euler's list complete under GRH. Unconditionally at most two further idoneal
     numbers can exist: Kani (Ann. Sci. Math. Québec 35 (2011), Cor. 23 and Rem. 24) shows that the often-repeated
     "at most one" is unjustified, and OEIS A000926 also says two. The obstruction is the possible Siegel zero: the class
     number of a binary form is an `L(1, χ)`-value, and lower bounds for it are ineffective. In that literature the
     idoneal problem is the one rank of its family where "no large anomaly" is conditional.
     - Class number one (largest `|D| = 163`) is unconditional (Heegner, Baker, Stark).
     - One-class genera of rank `>= 3` are effectively and completely classified (Watson: none in rank `> 10`;
       Lorch-Kirschmer 2013), because the mass formula is effective there.
   - *In the repository's adopted status* (THM-4566, owner directive 2026-10-07; golden_tail_idoneal_20261017c.md
     section 6), completeness is PROVED modulo openai/math #003. That is the quasi-Riemann zero-free half-plane
     `Re s > 7/8`, which excludes exactly the Siegel zeros. So "truly only 65" holds there, and the #003 input removes
     precisely the obstruction that keeps the classical statement conditional.
2. **What the knowledge can do elsewhere.** A finiteness proof of this kind needs three ingredients:
   - an involution (or local-global structure) whose fixed points are the anomalies;
   - a growth gap between the whole set and the fixed set;
   - an **effective** form of that gap.

   Section 1 sorts the classical "largest anomaly" problems by which ingredient is missing. None of the open ones
   follows by analogy alone. The transferable part is the diagnosis: when the gap is a value `L(1, χ)`, the obstruction
   is a Siegel zero.

   THM-4666 turns that obstruction into a tournament statement. A real zero of `L(s, χ_(-p))` forces the natural order
   of the Paley tournament `P_p` to be improvable by moving one vertex, and Fekete positivity is the same statement for
   randomised moves. So combinatorial certificates (local optimality of a ranking) exclude the Siegel zero prime by
   prime. The deterministic certificates exist only on a set of primes of density zero (Kalmynin Thm 1, via
   Baker-Montgomery). The randomised ones cover about 80% of the primes below 3000 (numerically).
3. **1848 and sixteen.**
   - The class numbers of the 65 idoneal discriminants `-4n` are exactly `1, 2, 4, 8, 16`, with multiplicities
     `5, 15, 24, 17, 4` (FINITE-EXACT). Sixteen is the maximal 2-rank 4.
   - In the conductor tower `Cl(-1848·4^k) = C_(2^k) × C_2^3`, sixteen solo (self-inverse) classes persist at every
     level `k >= 1` (eight at `D = -1848`) (clock_duality_oct15, prime463_conductor_clock_20261017e.md).
   - The repository's sedenion rung has the same count in a different predicate: sixteen rotational classes at
     `n = 17`, all self-converse (THM-871). The two rank-16 even unimodular lattices are THM-869 (S). clock_duality's
     escape_renewal section 7 made the first comparison.
   - The shared numeral is NUMEROLOGY. The shared shape, a group of exponent 2 on which the involution acts trivially,
     is ANALOGY (section 6).
4. **Merged and solo (PROVED dictionary).** Each piece has an exact counterpart:
   - flipping all edges ↔ inverting a form class ↔ complex conjugation of ideals ↔ improper equivalence;
   - merged pairs `{T, T^op}` ↔ pairs `{g, g^-1}`;
   - solo (self-converse) classes ↔ ambiguous classes (the 2-torsion);
   - the merged metagraph ↔ Lagrange's `GL_2(Z)` classes, counted as `(h + |Cl[2]|)/2`.

   "Idoneal" means that every class is solo. The tournament analogues are all tournaments on `n <= 3` vertices, and
   every rotational family (section 3).
5. **Binary relational information.** By genus theory, `n` is idoneal iff the primes `p ∤ 2n` with `p = x^2 + ny^2`
   are decided by finitely many quadratic residue symbols, i.e. by Paley-tournament arcs and Paley-graph edges.
   - So 1848 is the largest `n` for which binary information suffices: under GRH, or on #003 in the adopted status.
   - Beyond it, cubic or higher residuosity is always needed. For example (Gauss), `p = x^2 + 27y^2` iff `p ≡ 1 (mod 3)`
     and `2` is a cubic residue mod `p`, and `p = x^2 + 64y^2` iff `p ≡ 1 (mod 4)` and `2` is a quartic residue mod `p`.
   - This is the arithmetic form of clock_duality's observation that eliminating hidden vertices creates higher-arity
     relations (marked_contacts_root_current_20261017e.md, sections 8-9).
6. **Imaginary class groups as tournament invariants (THM-4666, HYP-9287).**
   - The natural order of `P_p` has excess `p h(-p)`, the area under the Legendre walk.
   - The inverse order has excess `3h(-p)`. This is proved by a marked-coefficient argument, the marked-cubic mechanism
     in a tournament.
   - Pair-block orders have excess `(1 - 2(2/p))h(-p)`. The residue sub-tournament has `(p - 3)h(-p)/4`.
   - The arcs at a vertex are the split primes of `Q(sqrt(-p))`.
   - Class number one, the least possible excess, holds exactly for `P_7, ..., P_163`, which correspond to Euler's lucky
     numbers.
   - The natural order is not the best ranking in general. Rotation orders beat it at high-L primes (167, 263, 383, 503,
     599, 647, ...), which refutes the threshold form of HYP-9287. It is never optimal for `h(-p) = 1`, `p > 7`.
     Whether it is optimal at only finitely many primes (a tournament largest anomaly) is open.
   - Paley tournaments exist exactly on the rational primes that stay prime in `Z[i]`.
7. **463 (section 7).** Three occurrences are structural:
   - `463 = 1 + 462` is prime **because** 462 is idoneal (Euler's criterion);
   - the norm-463 generator of the 1848 conductor clock (clock_duality);
   - the ranking excess `463·7 = 3241` of `P_463` (THM-4666 (1); both identities first computed in
     prime463_conductor_clock P4).

   The rest are coincidences of the numeral, typed NUMEROLOGY.

## 1. The family of "largest anomaly" problems, by mechanism

| problem | largest | status | what proves it / what blocks it |
|---|---|---|---|
| class number one, imaginary quadratic | 163 | PROVED unconditionally (Heegner; Baker 1966; Stark 1967) | effective: modular functions, linear forms in logarithms |
| Euler's lucky numbers (`n^2+n+q` prime, `0 <= n <= q-2`) | 41 | PROVED (Rabinowitsch 1913 ⟺ class number one) | as above |
| class number two | 427 | PROVED (Baker; Stark) | linear forms in logarithms |
| idoneal = one class per genus, `D = -4n` | 1848 | GRH (Weinberger 1973); at most two unknown unconditionally (Kani 2011); repository: complete modulo openai/math #003 (THM-4566) | Siegel ineffectivity of `L(1, χ)` |
| class group of exponent 3 | 4027 | ERH (Elsenhans-Klüners-Nicolae 2020); unconditional for `|D| < 3.1·10^20`; Boyd-Kisilevsky 1972 | same |
| class group of exponent 5 | 37363 | ERH (EKN; Heath-Brown 2008 for finiteness) | same |
| one-class genera of positive forms, rank `>= 3` | finite list, none in rank `> 10` | PROVED (Watson; Lorch-Kirschmer, LMS JCM 2013) | mass formula, effective |
| even unimodular lattices with one class | rank 8 (`E_8`); rank 16 has two | PROVED (mass formula; Witt) | effective |
| sums of at most seven cubes | 454 | PROVED (Siksek 2016) | explicit Linnik + computation |
| sums of at most six cubes | 8042 | OPEN (Jacobi's 1851 observation from Dase's table) | no effective method |
| odd numbers not of the form `x^2+y^2+10z^2` | 2719 | GRH (Ono-Soundararajan 1997); ineffective unconditionally (Duke-Schulze-Pillot) | Siegel (ternary forms) |
| decimal `2^n` without the digit 0 | `2^86` | OPEN (heuristic) | no growth law to make effective |
| ternary `2^n` without the digit 2 (Erdős) | `n = 8` | OPEN | 3-adic dynamics (Lagarias 2009) |
| Fermat primes `2^(2^k)+1` | `k = 4` (65537) | OPEN (heuristic) | none |
| Collatz cycles with `m` local minima | none beyond the trivial cycle for `m <= 91` | PROVED for `m <= 91` (Simons-de Weger: `m <= 68`; Hercher 2023: `m <= 91`); effective for each fixed `m` | linear forms in logarithms in `log2 3` (THM-4632: the boundary growth exponent of the trivial-character transfer law is `2(μ(log2 3) - 1)`) |

- The proved rows share an **effective** growth gap.
- The GRH rows have a gap whose constant is ineffective.
- The open rows have no gap at all, only heuristics.

The Collatz cycle question is on the effective side for every fixed number of excursions, but not uniformly. Its
uniform version is a "no anomaly at all" statement and not a finiteness statement.

## 2. Binary sufficiency: idoneal numbers decide primes by quadratic residues

Genus theory (Gauss) assigns to each form class its vector of genus characters, which are quadratic residue symbols of
the represented primes. The following are equivalent (Cox, *Primes of the form x² + ny²*, §3.C, Thm 3.22, and §6):

- `n` is idoneal;
- every genus of discriminant `-4n` contains one class;
- the ring class field of `Z[sqrt(-n)]` equals the genus field, a multiquadratic field;
- `p = x^2 + ny^2` (`p ∤ 2n`) iff `p` satisfies congruence conditions modulo `4n`, i.e. conditions on finitely many
  Legendre symbols.

The genus characters evaluated at `p` are `(-1/p)`, `(2/p)` and the symbols `(p/q)` for the odd primes `q | n`. Each
`(p/q)` is binary relational data: the arc between `0` and `p mod q` in the Paley tournament `P_q` when
`q ≡ 3 (mod 4)`, and the edge or non-edge in the Paley graph when `q ≡ 1 (mod 4)`. So **for idoneal `n` the
representation question is decided by binary relational data**.

For `n = 1848` the ring class field is `Q(i, sqrt2, sqrt3, sqrt7, sqrt11)`, of degree `32 = 2·16`. The audit checked
that `p = x^2 + 1848y^2` iff the five Legendre symbols equal 1, for all 20,639 such primes below `10^7`.

When `Cl(D)` is not elementary abelian, the ring class field is a non-abelian generalised dihedral extension. No finite
set of quadratic conditions then decides representation. The obstruction is the principal genus `Cl(D)^2`, the
"non-binary part", of size `h/2^(μ-1)`. Its size tends to infinity (Siegel: `h >> |D|^(1/2-ε)`, while
`2^(μ-1) = |D|^(o(1))`), which gives finiteness, but ineffectively.

## 3. The involution dictionary (PROVED, elementary)

| tournaments (repository metagraph) | binary quadratic forms / ideal classes |
|---|---|
| converse `T ↦ T^op` (flip every arc) | inversion `(a, b, c) ↦ (a, -b, c)`; complex conjugation of ideals; improper equivalence |
| merged node `{T, T^op}`, `T ≇ T^op` | inverse pair `{g, g^-1}`, `g^2 ≠ 1` |
| solo node (self-converse) | ambiguous class (`g^2 = 1`), counted by genus theory `2^(μ-1)` |
| merged metagraph (classes up to converse) | Lagrange's classes under `GL_2(Z)` |
| number of metagraph nodes `(T(n) + SC(n))/2` | `(h + |Cl[2]|)/2` |
| all classes solo | idoneal (`Cl` of exponent dividing 2) |

Both counts are Burnside's lemma for an involution; clock_duality's escape_renewal section 7 states the orbit law. For
tournaments `T(n) = 1, 1, 2, 4, 12, 56, 456, 6880` and `SC(n) = 1, 1, 2, 2, 8, 12, 88, 176` (OEIS A000568, A002785;
CITED; recomputed by the audit with Burnside). Every class is solo exactly for `n <= 3`, and the solo fraction falls to
`0.026` at `n = 8` and tends to 0.

For class groups the solo part is thinner still: `|Cl[2]| = 2^(μ-1) = |D|^(o(1))`, while `h >> |D|^(1/2-ε)` (Siegel).
In both settings the fixed points of the involution become a vanishing fraction, but only the class-group version has
an ineffective constant. Cox's Theorem 3.22 lists "two forms of discriminant `-4n` are equivalent iff properly
equivalent" among the equivalent forms of idoneality. In the dictionary this says Lagrange's count equals Gauss's:
there are no merged nodes.

**Families with every class solo.**
- Rotational (circulant) tournaments are always self-converse, since vertex negation sends `S` to `-S`
  (clock_duality's THM-871 scope repair). In the dictionary they are the ambiguous forms of the tournament world.
- At the Fermat rung `n = 17` they form exactly sixteen classes (THM-871).
- For `Cl(-7392)` (`n = 1848`) all sixteen classes are ambiguous.

The two sixteens arise from different mechanisms: a free `Z/16` multiplier action on `2^8` connection sets, against five
genus characters with one product relation. That makes them NUMEROLOGY as numbers and ANALOGY as shapes.

## 4. Class numbers as rankings (THM-4666) and the refuted threshold conjecture HYP-9287

See [THM-4666](../../01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md)
for statements and proofs. The points for this note are the following.

* **Excess = class number.** The natural order of `P_p` has (forward − backward) `= Σ_t S(t) = p h(-p)`, where `S` is
  the Legendre walk, and every affine order has excess `±p h(-p)`. The identity is checked at 9006 primes below
  `2·10^5` and at 6079 fundamental discriminants (the oriented graphs `Γ_D`).
* **Arithmetic orders.**
  - The discrete-log order has excess 1 (Jacobsthal).
  - The inverse order has excess `3h(-p)`. The proof is a marked-coefficient argument: substitute `b = <at>`, make the
    order condition linear in least residues, and get `N(t) = h(1 + χ(t - 1) - χ(t))`.
  - Every pair-block order (each pair `{x, -x}` adjacent) has excess `(1 - 2(2/p))h(-p)`.
  - The residue sub-tournament has `(p - 3)h(-p)/4`.
  - Non-affine Moebius orders give only small odd excesses.
* **Imaginary primes.** `0 → q` in `P_p` iff `q` splits in `Q(sqrt(-p))`. The class number is the harmonically weighted
  out-flow of a vertex: large `h(-p)` means small primes split, so short steps point forward. The tournament exists
  exactly when `p` is inert in `Z[i]`.
* **The largest poor ranking.** Excess `p` (the minimum) holds iff `h(-p) = 1`, iff `p ∈ {7, 11, 19, 43, 67, 163}`,
  iff `p = 4q - 1` with `q` one of Euler's lucky numbers `2, 3, 5, 11, 17, 41` (HYP-2226, HYP-3730 record the lucky
  numbers). This is unconditional: the proved member of the anomaly family, restated as a tournament theorem.
* **Siegel zeros as improvable rankings.** Moving one vertex over `t` neighbours changes the excess by `-2 S(t)`. So the
  natural order is insertion-locally optimal iff `p ∈ L+` (all partial sums nonnegative). Then Fekete's polynomial is
  positive on `(0, 1)` and `L(s, χ_(-p))` has no real zero, so a Siegel zero makes the natural order improvable by a
  single move.
  - Geometrically randomised moves change the excess by `-(2/x) f_p(x)` in expectation. Fekete positivity is local
    optimality in mixed strategies.
  - `L+` contains 2937 of the 9005 primes `p ≡ 3 (mod 4)` in `[7, 2·10^5]`.
  - Fekete positivity holds at 173 of the 217 primes below 3000 (numerical), against 97 in `L+`.
  - The `L+` share rises with `L(1, χ)`: 0.000 for `L < 1`, 0.451 on `[1.5, 2)`, 0.996 for `L >= 3`. Its share among
    primes decays slowly (0.54 to 0.32 across dyadic ranges), consistent with density zero. Kalmynin's upper bound is
    `<< π(x)/(log log x)^0.0368`, and it is not known whether `L+` is infinite.
* **Optimal rankings (HYP-9287: threshold form REFUTED).**
  - Exact MAS for `p <= 31`: the natural order is optimal at 3, 7 and 23 only. At 43 a witness of excess 183 beats the
    natural 43 (543 is the certified maximum: Chindelevitch-Harutyunyan, unrefereed preprint).
  - Pre-registered local search on `131 < p <= 500` left 9 of the 10 primes with `L >= 2.2` unbeaten. But certified
    rotation orders `x -> {αx}` beat 167, 263 and 383 (4 of 10 beaten with 439).
  - On `500 < p <= 1000`, rotation orders beat 503, 599, 647 and 983 (audit) and 887, 911 (S40 scan). Both "unbeaten"
    predictions fail; both "low L is beaten" predictions stand.
  - Unbeaten by any known order: 47, 71, 191, 239, 311, 359, 431, 479 (audit, full rotation family). 503 (`L = 2.942`) is
    beaten while 191 (`L = 2.955`) is not, so a pure `L`-threshold fails, and the boundary rises quickly with `p`.
  - Revised conjecture: the natural order is optimal only on a thin set inside `L+`, where `L(1, χ_(-p))` dominates the
    small-conductor twists that govern the rotation orders. Whether that set is finite is open.
* **Spectral ceiling.** Every order has excess `<= (p^(3/2)/π)(log p + γ + log(4/π) + o(1))` (trace duality). Exact
  optima at `p <= 43` sit at 0.40-0.50 of this ceiling.

## 5. The exponent ladder (FINITE-EXACT to `X = 4·10^6`)

The [script](../../04-computation/experiments/class_group_exponent_ladder_20261017.py):
- takes a reduced-form census of all discriminants `-X <= D < 0`;
- computes genus counts;
- computes exponents by Dirichlet composition, with a randomised pre-test and an exact cross-check against full orders
  for `|D| <= 20000`.

Validation: the `e = 2` list restricted to `D = -4n` is exactly Euler's 65 idoneal numbers. The lists agree with
Elsenhans-Klüners-Nicolae (Acta Arith. 193 (2020); quoted by Klüners-Komatsu, Math. Comp. 90 (2021)) wherever EKN give
them. The audit reproduced rows `e = 1, 2, 3, 4, 5, 7, 9` with independent composition code.

| e | discriminants | fundamental | largest `|D|` | largest fundamental `|D|` | `h` there | new members in `(10^6, 4·10^6]` | completeness / literature |
|---:|---:|---:|---:|---:|---:|---|---|
| 1 | 13 | 9 | 163 | 163 | 1 | none | PROVED (Heegner-Baker-Stark) |
| 2 | 101 | 65 | 7392 | 5460 | 16 | none | GRH (Weinberger); repository: modulo #003 (THM-4566); idoneal `n = 1848` |
| 3 | 42 | 26 | 4027 | 4027 | 9 (`C_3^2`) | none | ERH (EKN); 26 = 9 + 17 matches EKN |
| 4 | 586 | 268 | 887040 | 435435 | 256; 128 (`C_2^3 × C_4^2`) | none | ERH-complete (EKN Thm 2); 268 = 9 + 56 + 203 matches EKN |
| 5 | 44 | 36 | 37363 | 37363 | 25 (`C_5^2`) | none | ERH (EKN); 36 = 9 + 27 matches EKN |
| 7 | 53 | 42 | 118843 | 118843 | 49 (`C_7^2`) | none | unconditional for `|d| <= 3.1·10^20` (EKN Thm 1); 9 + 33 = 42 matches |
| 9 | 112 | 67 | 430411 | 430411 | 81 (`C_9^2`) | none | not in the literature checked |
| 6 | 1338 | 512 | 3892980 | 3892980 | 576 | yes (still growing) | EKN Thm 1: 514 fields below `3.1·10^20` (512 below `4·10^6`) |
| 8 | 2727 | 1028 | 3953280 | 3608760 | 512 | yes (still growing) | ERH (EKN: 1046 fields), beyond this range |
| 10, 12, 16 | 1513, 8784, 8062 | 907, 2330, 2975 | near `4·10^6` | | | yes | open in range |

`e` means "exponent dividing `e`". The `10^6` run is kept as `class_group_exponent_ladder_1e6_20261017.out`.

The odd ladders stop at `|D| ≈ (50-65) e^4`. The largest member is a rank-two homocyclic group `C_e × C_e` with a small
`L(1, χ)` (`h = e^2 ≈ sqrt(|D|) L(1, χ)/π`). The ladders `e = 6, 8, 10, 12, 16` keep growing within the range (`e = 2, 4`
stop), because the 2-rank `μ - 1` grows with the number of prime factors.

## 6. Sixteen: an atlas with verdicts

| object | the sixteen | mechanism | verdict |
|---|---|---|---|
| idoneal discriminants | class numbers `1, 2, 4, 8, 16` (5, 15, 24, 17, 4 discriminants); the 2-rank stops at 4 | genus count `2^(μ-1)` against `h ≈ sqrt(|D|) L(1, χ)/π`: four odd primes and two 2-adic characters is the last balance | FINITE-EXACT (+ GRH, or #003, for completeness) |
| 1848 conductor tower | `|Cl[2]| = 16` for `k >= 1` (8 at `k = 0`), `Cl = C_(2^k) × C_2^3` | marked ideal of norm 463, `v_2(B_n) = v_2(n)` (clock_duality P1-P3) | PROVED (clock_duality) |
| Cayley-Dickson | sedenions (dimension 16) lose division | Hurwitz; Adams (Hopf invariant one) | CITED |
| even unimodular lattices | rank 16: two classes, `E_8 ⊕ E_8` and `D_16^+`, isospectral | mass formula; THM-869 (S) already ties this to the sedenion rung | PROVED/CITED; `E_8 ⊕ E_8` is the Cayley-Dickson double of the octavian lattice, `D_16^+` is not (elementary remark) |
| rotational tournaments at `n = 17` | sixteen classes, all self-converse | free `Z/16` action on `2^8` connection sets (THM-871); compared with `Cl(-7392)` in escape_renewal section 7 | PROVED (repository) |
| Fermat primes | `65537 = 2^16 + 1` is the last known | none known | OPEN |
| ring class field of `Z[sqrt(-1848)]` | `Q(i, sqrt2, sqrt3, sqrt7, sqrt11)`, degree `32 = 2·16` | genus theory | PROVED (classical; audit-checked below `10^7`) |

The one common structure is an elementary abelian 2-group with a quadratic refinement: genus characters on `Cl(D)`,
the sign cocycle of a Cayley-Dickson basis on `(Z/2)^k`, glue codes of unimodular lattices. Its size `2^4` is the cutoff
in each row, for unrelated reasons:
- a growth balance for class groups;
- topology for division algebras;
- the mass formula for lattices;
- multiplier freeness for tournaments.

Treating the numeral as a law is NUMEROLOGY. "Exponent-2 structure is where the binary (quadratic) description is
complete" is the ANALOGY that survives.

Recursions in the repository that use values around sixteen are the same kind of cutoff:
- the Cayley-Dickson rungs `n = 2^k + 1` (THM-1935 lists 2, 3, 5, 9, with the quaternion wall at `n = 5`; 17 is the next
  rung);
- the rank-16 isospectral pair (THM-869 (S));
- the two-octave energy transport `E_4X = 4E_X + 8C_X + 4D_X + R_X` of the Collatz census (THM-4658);
- the sixteen solo classes in every level `k >= 1` of the 1848 tower.

## 7. 463, typed

| occurrence | predicate | verdict |
|---|---|---|
| `463 = 1^2 + 462·1^2` is prime | Euler's idoneal criterion: for idoneal `n`, an odd `m` prime to `n` with a unique proper representation by `x^2 + ny^2` is a prime or a prime power; `m = n + 1` has the unique representation `1 + n` | structural (PROVED, classical; all 65 idoneal `n` checked: `n + 1` is a prime, a prime power or twice one) |
| `N(1 + sqrt(-462)) = 463` generates the conductor kernel of `Cl(-1848·4^k)` at every level `k >= 1` | marked ideal; `462 = 2·3·7·11` idoneal | structural (clock_duality, PROVED) |
| `h(-463) = 7`; `P_463` has natural excess `463·7 = 3241` and half-window imbalance `7` (`463 ≡ 7 mod 8`) | THM-4666 (1), (2); first computed in prime463_conductor_clock_20261017e.md P4 | instance of a theorem |
| `463 ∉ L+`, `L(1, χ_(-463)) = 1.022`; the natural order of `P_463` is beaten (excess 9553 against 3241) | THM-4666 (5), HYP-9287 | instance |
| `QR_463[μ14]` all-odd | procgen_petersen_20261001_petersen_heawood_paley.md section 6.4 | NUMEROLOGY |
| EI17 embeds in `G_463` (THM-4665 twisted controls) | unit-distance graph over `F_(463^2)` | NUMEROLOGY |
| THM-144: `r_9 = 463` for the Interval tournament `T_11` | `Ω_9 - r_10 = 611 - 148` under an exactness hypothesis | NUMEROLOGY |
| HYP-2880: `H_count = 463` at `w = 3` | number of distinct values | NUMEROLOGY |
| THM-463 | a theorem identifier (two-cube divisor criterion) | not about 463 |

Euler's criterion has a further consequence ([check](../../04-computation/experiments/idoneal_euler_criterion_463_20261017.py)).
In all fifteen cases `n = 2m` with `m` odd squarefree (clock_duality's P1 family), `n + 1` is prime: 3, 7, 11, 23, 31,
43, 59, 71, 79, 103, 131, 191, 211, 331, 463. So the marked ideal `(1 + sqrt(-n))` of the conductor clock always has
prime norm. Also `1849 = 43^2` for `n = 1848`. Euler's own certificate `18518809 = 197^2 + 1848·100^2` (a unique
representation, hence prime) is replayed. **463 is prime because 462 is idoneal.**

## 8. Open

- An exact MAS for `P_47` or `P_71`, the first unbeaten primes beyond the DP range. The method would be the exact DP
  modulo `Aut(P_p)` of Chindelevitch-Harutyunyan. A plain MIP solves `p = 23` in 72 s but stops at `p = 31` with a gap
  of 22 after 1200 s.
- The largest prime where the natural order is optimal, if there is one: full rotation scans with pre-registration.
- An explicit formula for rotation-order excesses in terms of twisted `L`-values (character analogues of Dedekind sums),
  and a pre-registered test of the revised HYP-9287 on `p <= 2000`.
- A combinatorial route to the absence of Siegel zeros beyond Fekete positivity. In the tournament reading these would
  be block moves, i.e. weighted positivity of the Legendre walk.
- The `e = 9` list (112 discriminants, 67 fundamental, largest 430411): is it complete under ERH? (`e = 4` is
  ERH-complete by EKN Thm 2, and `e = 7` is unconditional to `3.1·10^20` by EKN Thm 1.)

## Reproduction

```
python 04-computation/experiments/paley_class_number_ranking_20261017.py 200000 20000 23
python 04-computation/experiments/paley_mas_exact_20261017.py 11 --brute ; ... 19 ; ... 23 ; ... 31
python 04-computation/experiments/paley_legendre_walk_positivity_20261017.py 200000
python 04-computation/experiments/paley_fekete_geometric_moves_20261017.py 3000
python 04-computation/experiments/paley_arithmetic_orders_20261017.py 3000
python 04-computation/experiments/paley_mobius_orders_20261017.py 71 300
python 04-computation/experiments/paley_ranking_local_search_20261017.py 131 30 150
python 04-computation/experiments/paley_ranking_local_search_20261017.py 132:500 12 80    # pre-registered
python 04-computation/experiments/paley_ranking_local_search_20261017.py 501:1000 6 40   # pre-registered
python 04-computation/experiments/paley_rotation_orders_20261017.py 12 167 191 239 263 311 359 383 431 479
python 04-computation/experiments/class_group_exponent_ladder_20261017.py 4000000
python 04-computation/experiments/idoneal_euler_criterion_463_20261017.py
```
