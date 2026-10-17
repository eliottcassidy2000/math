# Binary information, class groups and largest anomalies: idoneal 1848, the class-number ranking of Paley tournaments, the exponent ladder, and the cutoffs at sixteen

2026-10-17, opus-2026-10-17-S40. Owner prompt: imaginary quadratic class groups and imaginary primes; the 65 idoneal
numbers as a model for proving that other "large number anomalies" do not exist; 1848 as the largest, its sixteen-fold
self-similarity, cutoffs at powers of two and the Cayley-Dickson tower; merged versus solo classes of the tournament
metagraph; past work on 463; binary relational information.

**Typing.** PROVED (elementary): THM-4666 (class number = ranking excess; local optimality of the natural Paley order
= Kalmynin's L+; a real zero forces an improving move), the involution dictionary (section 3), the binary-sufficiency
reading (section 2, from classical CITED theorems). FINITE-EXACT: the identity checks of THM-4666, exact MAS for
`p <= 31`, the exponent ladder to `4·10^6` (section 5; it reproduces Euler's 65 idoneal numbers and the published
exponent-3 and exponent-5 lists), the idoneal class-number histogram. CONJECTURE: HYP-9287 (pre-registered test
PASSED). CITED: Weinberger, Heegner-Baker-Stark, Rabinowitsch, Watson, Lorch-Kirschmer, Elsenhans-Klüners-Nicolae,
Siksek, Ono-Soundararajan, Fekete-Polya, Baker-Montgomery, Kalmynin, Chindelevitch-Harutyunyan. ANALOGY / NUMEROLOGY:
the sixteen atlas and the 463 atlas (verdicts per row).

## 0. Answers in brief

1. **"Only 65 idoneal numbers" is a theorem under GRH, not unconditionally.** Weinberger (1973) proved that Euler's
   list is complete under GRH, and that unconditionally at most one or two further idoneal numbers can exist (the
   sources differ: OEIS A000926 says one; Kani and Wikipedia say two). The obstruction is the possible Siegel zero: the class number of a binary form is an
   `L(1, χ)`-value, and lower bounds for it are ineffective. Its neighbours fare differently. Class number one (largest
   `|D| = 163`) is unconditional (Heegner, Baker, Stark). One-class genera of rank `>= 3` are effectively and completely
   classified (Watson: none in rank `> 10`; Lorch-Kirschmer 2013), because the mass formula is effective there.
   **The idoneal problem is the one rank in its family where "no large anomaly" is still conditional.**
2. **What the knowledge can do elsewhere.** A finiteness proof of this kind needs three ingredients:
   - an involution (or local-global structure) whose fixed points are the anomalies;
   - a growth gap between the whole set and the fixed set;
   - an **effective** form of that gap.

   Section 1 sorts the classical "largest anomaly" problems by which ingredient is missing. None of the open ones
   follows by analogy alone. The transferable part is the diagnosis: when the gap is a value `L(1, χ)` the obstruction
   is a Siegel zero.

   THM-4666 turns that obstruction into a tournament statement. A real zero of `L(s, χ_(-p))` forces the natural order
   of the Paley tournament `P_p` to be improvable by moving one vertex. So a combinatorial certificate (local
   optimality of a ranking) excludes the Siegel zero prime by prime. Such certificates exist only on a set of primes
   of density zero (Baker-Montgomery).
3. **1848 and sixteen.** The class numbers of the 65 idoneal discriminants `-4n` are exactly `1, 2, 4, 8, 16`, with
   multiplicities `5, 15, 24, 17, 4` (FINITE-EXACT). Sixteen is the maximal 2-rank 4. In the conductor tower
   `Cl(-1848·4^k) = C_(2^k) × C_2^3`, sixteen solo (self-inverse) classes persist at every level (clock_duality_oct15,
   prime463_conductor_clock_20261017e.md). The repository's sedenion rung has the same count in a different
   predicate: sixteen rotational classes at `n = 17`, all self-converse (THM-871), and the two rank-16 even unimodular
   lattices (THM-869 (S)). The shared numeral is NUMEROLOGY. The shared shape, a group of exponent 2 on which the
   involution acts trivially, is ANALOGY (section 6).
4. **Merged and solo (PROVED dictionary).** Each piece has an exact counterpart:
   - flipping all edges ↔ inverting a form class ↔ complex conjugation of ideals ↔ improper equivalence;
   - merged pairs `{T, T^op}` ↔ pairs `{g, g^-1}`;
   - solo (self-converse) classes ↔ ambiguous classes (the 2-torsion);
   - the merged metagraph ↔ Lagrange's `GL_2(Z)` classes, counted as `(h + |Cl[2]|)/2`.

   "Idoneal" means that every class is solo. The tournament analogues are all tournaments on `n <= 3` vertices, and
   every rotational family (section 3).
5. **Binary relational information.** By genus theory, `n` is idoneal iff the primes `p ∤ 2n` with
   `p = x^2 + ny^2` are decided by finitely many quadratic residue symbols, that is, by Paley-tournament arcs. So 1848
   is the largest `n` for which binary information suffices (under GRH). Beyond it, cubic or higher residuosity is
   always needed: for example `p = x^2 + 27y^2` iff `2` is a cubic residue, and `p = x^2 + 64y^2` iff `2` is a quartic
   residue. This is the arithmetic form of clock_duality's observation that eliminating hidden vertices creates
   higher-arity relations (marked_contacts_root_current_20261017e.md, sections 8-9).
6. **Imaginary class groups as tournament invariants (THM-4666, HYP-9287).**
   - The natural order of `P_p` has excess `p h(-p)`, the area under the Legendre walk.
   - The arcs at a vertex are the split primes of `Q(sqrt(-p))`.
   - Class number one, the least possible excess, holds exactly for `P_7, ..., P_163`, which correspond to Euler's lucky
     numbers.
   - The natural order is globally optimal only when `L(1, χ_(-p))` is large (pre-registered test PASSED).
   - Paley tournaments exist exactly on the rational primes that stay prime in `Z[i]`.
7. **463 (section 7).** Two occurrences are structural:
   - the norm-463 generator of the 1848 conductor clock (clock_duality);
   - the ranking excess `463·7 = 3241` of `P_463` (THM-4666 (1)).

   The rest are coincidences of the numeral, typed NUMEROLOGY.

## 1. The family of "largest anomaly" problems, by mechanism

| problem | largest | status | what proves it / what blocks it |
|---|---|---|---|
| class number one, imaginary quadratic | 163 | PROVED unconditionally (Heegner; Baker 1966; Stark 1967) | effective: modular functions, linear forms in logarithms |
| Euler's lucky numbers (`n^2+n+q` prime, `0 <= n <= q-2`) | 41 | PROVED (Rabinowitsch 1913 ⟺ class number one) | as above |
| class number two | 427 | PROVED (Baker; Stark) | linear forms in logarithms |
| idoneal = one class per genus, `D = -4n` | 1848 | GRH (Weinberger 1973); at most one (two) unknown unconditionally | Siegel ineffectivity of `L(1, χ)` |
| class group of exponent 3 | 4027 | ERH (Elsenhans-Klüners-Nicolae 2020); unconditional for `|D| < 3.1·10^20`; Boyd-Kisilevsky 1972 | same |
| class group of exponent 5 | 37363 | ERH (EKN; Heath-Brown 2008 for finiteness) | same |
| one-class genera of positive forms, rank `>= 3` | finite list, none in rank `> 10` | PROVED (Watson; Lorch-Kirschmer, LMS JCM 2013) | mass formula, effective |
| even unimodular lattices with one class | rank 8 (`E_8`); rank 16 has two | PROVED (mass formula; Witt) | effective |
| sums of at most seven cubes | 454 | PROVED (Siksek 2016) | explicit Linnik + computation |
| sums of at most six cubes | 8042 | OPEN (Jacobi's 1851 observation) | no effective method |
| odd numbers not of the form `x^2+y^2+10z^2` | 2719 | GRH (Ono-Soundararajan 1997); ineffective unconditionally (Duke-Schulze-Pillot) | Siegel (ternary forms) |
| decimal `2^n` without the digit 0 | `2^86` | OPEN (heuristic) | no growth law to make effective |
| ternary `2^n` without the digit 2 (Erdős) | `n = 8` | OPEN | 3-adic dynamics (Lagarias 2009) |
| Fermat primes `2^(2^k)+1` | `k = 4` (65537) | OPEN (heuristic) | none |
| Collatz cycles with `m` local minima | none beyond the trivial cycle for `m <= 91` | PROVED for each bounded `m` (Simons-de Weger; Hercher) | linear forms in logarithms in `log2 3` (THM-4632: Collatz growth exponents are irrationality measures of `log2 3`) |

The proved rows share an **effective** growth gap. The GRH rows have a gap whose constant is ineffective. The open
rows have no gap at all, only heuristics.

The Collatz cycle question sits on the effective side for every fixed number of excursions, but not uniformly. Its
uniform version is a "no anomaly at all" statement and not a finiteness statement.

## 2. Binary sufficiency: idoneal numbers decide primes by quadratic residues

Genus theory (Gauss) assigns to each form class its vector of genus characters, which are quadratic residue symbols of
the represented primes. The following are equivalent (Cox, *Primes of the form x² + ny²*, §3.C, Thm 3.22, and §6):

- `n` is idoneal;
- every genus of discriminant `-4n` contains one class;
- the ring class field of `Z[sqrt(-n)]` equals the genus field, a multiquadratic field;
- `p = x^2 + ny^2` (`p ∤ 2n`) iff `p` satisfies congruence conditions modulo `4n`, i.e. conditions on finitely many
  Legendre symbols.

The genus characters evaluated at `p` are `(-1/p)`, `(2/p)` and the symbols `(p/q)` for the odd primes `q | n`.
Each `(p/q)` is binary relational data: the arc between `0` and `p mod q` in the Paley tournament `P_q` when
`q ≡ 3 (mod 4)`, the edge or non-edge in the Paley graph when `q ≡ 1 (mod 4)`. So **for idoneal `n` the representation
question is decided by binary relational data**. For `n = 1848`
the ring class field is `Q(i, sqrt2, sqrt3, sqrt7, sqrt11)`, of degree `32 = 2·16`.

When `Cl(D)` is not elementary abelian, the ring class field is a non-abelian generalised dihedral extension. No
finite set of quadratic conditions then decides representation: the obstruction is the principal genus `Cl(D)^2`, the
"non-binary part", of size `h/2^(μ-1)`. Its size tends to infinity (Siegel: `h >> |D|^(1/2-ε)`, while
`2^(μ-1) = |D|^(o(1))`). That gives finiteness, but ineffectively.

## 3. The involution dictionary (PROVED, elementary)

| tournaments (repository metagraph) | binary quadratic forms / ideal classes |
|---|---|
| converse `T ↦ T^op` (flip every arc) | inversion `(a, b, c) ↦ (a, -b, c)`; complex conjugation of ideals; improper equivalence |
| merged node `{T, T^op}`, `T ≇ T^op` | inverse pair `{g, g^-1}`, `g^2 ≠ 1` |
| solo node (self-converse) | ambiguous class (`g^2 = 1`), counted by genus theory `2^(μ-1)` |
| merged metagraph (classes up to converse) | Lagrange's classes under `GL_2(Z)` |
| number of metagraph nodes `(T(n) + SC(n))/2` | `(h + |Cl[2]|)/2` |
| all classes solo | idoneal (`Cl` of exponent dividing 2) |

Both counts are Burnside's lemma for an involution. For tournaments `T(n) = 1, 1, 2, 4, 12, 56, 456, 6880` and
`SC(n) = 1, 1, 2, 2, 8, 12, 88, 176` (OEIS A000568, A002785; CITED). Every class is solo exactly for `n <= 3`. The solo
fraction falls to `0.026` at `n = 8` and tends to 0. For class groups the solo part is thinner still:
`|Cl[2]| = 2^(μ-1) = |D|^(o(1))`, while `h >> |D|^(1/2-ε)` (Siegel). In both settings the fixed points of the involution
become a vanishing fraction, but only the class-group version has an ineffective constant. Cox's Theorem 3.22 lists
"two forms of discriminant `-4n` are equivalent iff properly equivalent" among the equivalent forms of idoneality. In
the dictionary this says Lagrange's count equals Gauss's: there are no merged nodes.

**Families with every class solo.**
- Rotational (circulant) tournaments are always self-converse, since vertex negation sends `S` to `-S` (THM-871 scope
  repair). In the dictionary they are the ambiguous forms of the tournament world.
- At the Fermat rung `n = 17` they form exactly sixteen classes (THM-871).
- For `Cl(-7392)` (`n = 1848`) all sixteen classes are ambiguous.

The two sixteens arise from different mechanisms: a free `Z/16` multiplier action on `2^8` connection sets, versus five
genus characters with one product relation. That makes them NUMEROLOGY as numbers and ANALOGY as shapes.

## 4. Class numbers as rankings (THM-4666) and the conjecture HYP-9287

See [THM-4666](../../01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md)
for statements and proofs. The points for this note are the following.

* **Excess = class number.** The natural order of `P_p` has (forward − backward) `= Σ_t S(t) = p h(-p)`, where `S` is
  the Legendre walk. Every affine order has excess `±p h(-p)`. The identity is checked at 9006 primes below `2·10^5`
  and at 6079 fundamental discriminants (the oriented graphs `Γ_D`).
* **Imaginary primes.** `0 → q` in `P_p` iff `q` splits in `Q(sqrt(-p))`. The class number is the harmonically weighted
  out-flow of a vertex: large `h(-p)` means small primes split, so short steps point forward. The tournament exists
  exactly when `p` is inert in `Z[i]`.
* **The largest poor ranking.** Excess `p` (the minimum) holds iff `h(-p) = 1`, iff `p ∈ {7, 11, 19, 43, 67, 163}`,
  iff `p = 4q - 1` with `q` one of Euler's lucky numbers `2, 3, 5, 11, 17, 41`. This is unconditional: the proved
  member of the anomaly family, restated as a tournament theorem.
* **Siegel zeros as improvable rankings.** Moving one vertex over `t` neighbours changes the excess by `-2 S(t)`, so
  the natural order is insertion-locally optimal iff `p ∈ L+` (all partial sums nonnegative). Then Fekete's polynomial
  is positive on `(0, 1)` and `L(s, χ_(-p))` has no real zero. A Siegel zero therefore makes the natural order
  improvable by a single move.
  - `L+` contains 2937 of the 9005 primes `p ≡ 3 (mod 4)` in `[7, 2·10^5]`.
  - Its share rises with `L(1, χ)`: 0.000 for `L < 1`, 0.451 on `[1.5, 2)`, 0.996 for `L >= 3`.
  - Its share among primes decays slowly (0.54 to 0.32 across dyadic ranges), consistent with density zero.
* **Optimal rankings follow the class number (HYP-9287).**
  - Exact MAS for `p <= 31`: the natural order is optimal at 3, 7 and 23 only. At 43 the certified MAS 543
    (Chindelevitch-Harutyunyan) beats the natural 473.
  - Pre-registered local search on `131 < p <= 500`: the natural order is unbeaten at 9 of the 10 primes with
    `L >= 2.2` (beaten only at 439, `L = 2.249`). All 20 primes with `L <= 1.7` are beaten by certified orders.
  - The threshold drifts from about 2 at `p ≈ 50` to between 2.25 and 2.52 at `p ≈ 450`.
  - The natural order wins exactly when its excess `L p^(3/2)/π` exceeds the generic optimum, about
    `0.75-0.96 p^(3/2)` there.
* **Spectral ceiling.** Every order has excess `<= (p^(3/2)/π)(log p + γ + log(4/π) + o(1))` (trace duality). Exact
  optima at `p <= 43` sit at 0.40-0.50 of this ceiling.

## 5. The exponent ladder (FINITE-EXACT to `X = 4·10^6`)

[Script](../../04-computation/experiments/class_group_exponent_ladder_20261017.py) (reduced-form census of all
discriminants `-X <= D < 0`, genus counts, exponents by Dirichlet composition with a randomised pre-test and an exact
cross-check against full orders for `|D| <= 20000`). Validation: the `e = 2` list restricted to `D = -4n` is exactly
Euler's 65 idoneal numbers. The fundamental lists for `e = 3` (26 fields, largest 4027) and `e = 5` (36 fields, largest
37363) are exactly those of Elsenhans-Klüners-Nicolae (Acta Arith. 193 (2020); quoted by Klüners-Komatsu, Math. Comp.
90 (2021), Cor. 22, Thm 32).

| e | discriminants | fundamental | largest `|D|` | largest fundamental `|D|` | `h` there | new members in `(10^6, 4·10^6]` | completeness |
|---:|---:|---:|---:|---:|---:|---|---|
| 1 | 13 | 9 | 163 | 163 | 1 | none | PROVED (Heegner-Baker-Stark) |
| 2 | 101 | 65 | 7392 | 5460 | 16 | none | GRH (Weinberger); idoneal `n = 1848` |
| 3 | 42 | 26 | 4027 | 4027 | 9 (`C_3^2`) | none | ERH (EKN) |
| 4 | 586 | 268 | 887040 | 435435 | 256; 128 (`C_2^3 × C_4^2`) | none | ERH (EKN list `E = 4`; numbers not compared) |
| 5 | 44 | 36 | 37363 | 37363 | 25 (`C_5^2`) | none | ERH (EKN) |
| 7 | 53 | 42 | 118843 | 118843 | 49 (`C_7^2`) | none | not in the literature checked |
| 9 | 112 | 67 | 430411 | 430411 | 81 | none | not in the literature checked |
| 6 | 1338 | 512 | 3892980 | 3892980 | 576 | yes (still growing) | open in range |
| 8 | 2727 | 1028 | 3953280 | 3608760 | 512 | yes (still growing) | ERH (EKN), beyond this range |
| 10, 12, 16 | 1513, 8784, 8062 | 907, 2330, 2975 | near `4·10^6` | | | yes | open in range |

(`e` means: exponent dividing `e`. The `10^6` run is kept as `class_group_exponent_ladder_1e6_20261017.out`.)

The odd ladders stop at `|D| ≈ (50-65) e^4`: the largest member is a rank-two elementary group `C_e × C_e` with a
small `L(1, χ)` (`h = e^2 ≈ sqrt(|D|) L(1, χ)/π`). The even ladders keep growing within the range, because the 2-rank
`μ - 1` grows with the number of prime factors.

## 6. Sixteen: an atlas with verdicts

| object | the sixteen | mechanism | verdict |
|---|---|---|---|
| idoneal discriminants | class numbers `1, 2, 4, 8, 16` (5, 15, 24, 17, 4 discriminants); the 2-rank stops at 4 | genus count `2^(μ-1)` against `h ≈ sqrt(|D|) L(1, χ)/π`: four odd primes and two 2-adic characters is the last balance | FINITE-EXACT (+ GRH for completeness) |
| 1848 conductor tower | `|Cl[2]| = 16` at every level, `Cl = C_(2^k) × C_2^3` | marked ideal of norm 463, `v_2(B_n) = v_2(n)` (clock_duality P1-P3) | PROVED (clock_duality) |
| Cayley-Dickson | sedenions (dimension 16) lose division | Hurwitz; Adams (Hopf invariant one) | CITED |
| even unimodular lattices | rank 16: two classes, `E_8 ⊕ E_8` and `D_16^+`, isospectral | mass formula; THM-869 (S) already ties this to the sedenion rung | PROVED/CITED; `E_8 ⊕ E_8` is the Cayley-Dickson double of the octavian lattice, `D_16^+` is not (elementary remark) |
| rotational tournaments at `n = 17` | sixteen classes, all self-converse | free `Z/16` action on `2^8` connection sets (THM-871) | PROVED (repository) |
| Fermat primes | `65537 = 2^16 + 1` is the last known | none known | OPEN |
| ring class field of `Z[sqrt(-1848)]` | `Q(i, sqrt2, sqrt3, sqrt7, sqrt11)`, degree `32 = 2·16` | genus theory | PROVED (classical) |

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
- the Cayley-Dickson rung convention `n = 2, 3, 5, 9, 17` (THM-1935, the quaternion wall at `n = 5`);
- the rank-16 isospectral pair (THM-869 (S));
- the two-octave memory of the Collatz census, a factor 4 per two octaves (THM-4658);
- the sixteen solo classes in every level of the 1848 tower.

## 7. 463, typed

| occurrence | predicate | verdict |
|---|---|---|
| `N(1 + sqrt(-462)) = 463` generates the conductor kernel of `Cl(-1848·4^k)` at every level | marked ideal; `462 = 2·3·7·11` idoneal, `463 = 1 + 462` | structural (clock_duality, PROVED) |
| `h(-463) = 7`; `P_463` has natural excess `463·7 = 3241` and half-window imbalance `7` (`463 ≡ 7 mod 8`) | THM-4666 (1), (2) | instance of a theorem |
| `463 ∉ L+`, `L(1, χ_(-463)) = 1.022`; the natural order of `P_463` is beaten (excess 9553 found against 3241) | THM-4666 (5), HYP-9287 | instance |
| EI17 embeds in `G_463` (THM-4665 twisted controls) | unit-distance graph over `F_(463^2)` | NUMEROLOGY |
| THM-144: `r_9 = 463` for the Interval tournament `T_11` | `Ω_9 - r_10 = 611 - 148` under an exactness hypothesis | NUMEROLOGY |
| HYP-2880: `H_count = 463` at `w = 3` | number of distinct values | NUMEROLOGY |
| THM-463 | a theorem identifier (two-cube divisor criterion) | not about 463 |

## 8. Open

- An exact MAS for `P_47` or `P_71`, the first high-`L` primes beyond the DP range: a proof or refutation of HYP-9287
  where it is most informative.
- Whether `L(1, χ_(-p))` is unbounded along `p ≡ 3 (mod 4)` (Joshi's limsup over prime moduli). If it is, THM-4666 (6)
  makes the best Paley rankings exceed the random-tournament scale `Θ(n^(3/2))` by an unbounded factor.
- A combinatorial route to the absence of Siegel zeros beyond `L+`, which has density zero. The tournament reading
  suggests looking at 2-vertex and block moves, i.e. weighted positivity of the Legendre walk.
- The rungs `e = 4, 7, 9` show no new member between `10^6` and `4·10^6`; `e = 7, 9` are not in the literature
  checked. Are these lists complete (under ERH)? The rungs `e = 6, 8` still grow at `4·10^6`.

## Reproduction

```
python 04-computation/experiments/paley_class_number_ranking_20261017.py 200000 20000 23
python 04-computation/experiments/paley_mas_exact_20261017.py 11 --brute ; ... 19 ; ... 23 ; ... 31
python 04-computation/experiments/paley_legendre_walk_positivity_20261017.py 200000
python 04-computation/experiments/paley_ranking_local_search_20261017.py 131 30 150
python 04-computation/experiments/paley_ranking_local_search_20261017.py 132:500 12 80   # pre-registered
python 04-computation/experiments/class_group_exponent_ladder_20261017.py 4000000
```
