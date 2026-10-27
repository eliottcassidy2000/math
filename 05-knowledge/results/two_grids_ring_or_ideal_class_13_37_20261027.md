# The owner's second grid is a ring or an ideal class: Z[√−d] plus its (½,½)-centred copy is the Eisenstein ring at d = 3 and the order-two class over 2 for d ≡ 1 (mod 4), and the two grids are the whole class group exactly for d = 5, 13, 37; which grid carries a norm shell; the cubic character of 2; the centred units complete the quartic designs; four Eisenstein fourth powers miss n ≡ 5 (mod 8)

**opus-2026-10-27-S66.**

- Scripts:
  - `04-computation/experiments/two_grids_class_group_20261027.py` (+ `.out`, ALL CHECKS PASSED, about 10 s);
  - `04-computation/experiments/eisenstein_four_fourth_powers_20261027.c` (+ `.out`);
  - `04-computation/experiments/eisenstein_fourth_power_obstructions_20261027.py` (+ `.out`).
- Canon: THM-4825, THM-4826 and HYP-9429, all filed this session.
- Concurrent lanes on the same directive. This note builds on them and does not repeat them.
  - **mac-mini-2026-10-27-slab**: [THM-4824 / HYP-9428](eisenstein_two_grids_13_37_fourth_powers_20261027.md). Merge
    debt by lattice shape (the Eisenstein grid as the first grid glued with its (½,½) shift), parity-coset shells, the
    zero-sum quartic identity, 13 and 37 as primes ≡ 1 (mod 12) with typed coincidences, and session69's square-later
    (100, 37, 13).
  - **session78**: [FOURTH_POWERS.md](session78_fourth_power_aliquot_20261027/FOURTH_POWERS.md) and
    [ALIQUOT.md](session78_fourth_power_aliquot_20261027/ALIQUOT.md). The four-fourth-powers window below 81, the gaps
    11, 12, 13 and the lone 50, the aliquot fibres of 11, 12, 13 and 37, and s(2^a·25) = 37·2^a − 31.

**Status.**
- **The owner's construction, in general (THM-4825, PROVED).** Take the integer grid Z[√−d] (d odd, squarefree) and lay
  over it its copy centred at (½, ½) in the coordinates of u + v√−d. Either the second grid completes the first to a
  ring, or it is an ideal class:
  - **d ≡ 3 (mod 4): a ring.** The union is the ring of integers. At d = 3 it is the Eisenstein ring, the owner's
    example.
  - **d ≡ 1 (mod 4): an ideal class.** The union is not a ring. It is half the prime ideal over 2, (2, 1 + √−d), whose
    class has order exactly 2 (d > 1).
- **Where 13 and 37 come from (THM-4825 (3)–(5)).** The two grids are the whole class group exactly when h(−4d) = 2,
  which happens exactly for d = 5, 13, 37. This rests on the Baker–Stark class-number-2 list, and was checked exactly
  for all 40,520 squarefree d ≡ 1 (mod 4) below 2·10⁵.
  - **Three grids agree.** For those three d, a split prime is a first-grid norm iff it is ≡ 1 (mod 4), iff it is a
    norm from the square grid Z[i]. Otherwise it is a norm from the centred grid. This is the owner-directive THM-4566's
    H(K) = K(i) read on primes, and it singles out {1, 5, 13, 37} among squarefree d ≡ 1 (mod 4).
  - **The owner's numbers.** In Z[√−5], 5, 25, 100 lie on the first grid and 10, 50 on the second: each factor 2 flips
    the grid.
- **Which grid carries a norm shell (THM-4826 (1)–(2), PROVED).**
  - In Z[ω] and in the Hurwitz quaternions, exactly one third of every odd-norm shell lies on the first grid (Z[√−3];
    the Lipschitz quaternions Z⁴), and all of every even shell does. The mechanism is the unit index 3 = |F₄*|.
  - Among the rings Z[(1 + √−d)/2], only d = 3 spreads every odd shell over both grids.
- **The grid of a prime is the cubic character of 2 (THM-4826 (3), classical in grid language).** The primary
  Eisenstein prime over p ≡ 1 (mod 3) lies on the first grid ⟺ 2 is a cubic residue mod p ⟺ p = x² + 27y².
  - Checked for all 39,231 such p below 10⁶, with no exception.
  - The primes over Mersenne primes 2^q − 1 (q > 3) always lie on the first grid.
  - 13 and 37 lie on the centred grid.
- **Fourth powers (THM-4826 (4)–(6)).**
  - The units of each grid alone do not form a spherical 4-design. Their unions, the hexagon and the 24-cell, are
    5-designs. So the quartic identities need both grids: Candido's for Z[ω], and Liouville's (via (1 + i)·Hurwitz =
    D₄), which gave the first Waring bound g(4) ≤ 53.
  - Eisenstein fourth powers are ≡ 0, u or 9u (mod 16), where u ∈ {1, ω, ω²} are the merge charges of THM-4820.
  - So a rational n ≡ 5 (mod 8) is never a sum of four Eisenstein fourth powers. The owner's 13 and 37 stay excluded on
    the Eisenstein grid.
- **HYP-9429 (OPEN, three pre-registered tests).**
  - (1) CONFIRMED: 6 has no such representation with bases ≤ 150, though it is locally representable everywhere.
  - (2) REFUTED: eight other values in 51..200 lack one with bases ≤ 60, against the predicted three.
  - (3) REFUTED: only two of those eight (102 and 186) gain one with bases ≤ 120.

  So 6, 54, 60, 75, 78, 124 and 140 are open. A post-hoc check finds them mostly in low local-density classes.
- **NUMEROLOGY.** {5, 13, 37} are all ≡ 5 (mod 8), the class excluded in (6). The mechanisms differ: Rédei's theory
  forces p ≡ 5 (mod 8) for h(−4p) = 2, while (6) is a 2-adic count of fourth powers.
- **Independent audit: OWED** (pre-audit checkpoint).

## 0. The directive

"consider the idea that the eisenstein integers are just another coordinate grid on top of the first grid, but centered
at (1/2,1/2), or something along those lines for inspiration. consider prior work on 13 and 37 and consider that of the
nonnegative integers up through 50, 13 natural numbers can be written as the sum of 4 4th powers. 0 to4,16 to 19, 32 to
34,48 to 49 =5+4+3+2, while 37 cannot. the lone 50 + ranges of 11,12,13. consider previous work connecting
5,10,25,50,100 and 11,12,13 in the aliquot related work."

**What the concurrent lanes settled, and what this note adds.**
- Session78 settled the fourth-power window (13 values, 37 gaps, lengths 11, 12, 13, the lone 50) and the aliquot side.
- Mac-mini made the two grids exact for the merge debt and typed the 13/37 coincidences as numerology. It also found
  session69's square-later census (100, 37, 13), which is the likeliest referent of "prior work on 13 and 37".
- This note asks what the owner's construction *is* arithmetically, for every d. The answer gives 13 and 37, with 5, a
  structural role:
  - they are the d for which "the grid and the grid with its centre" are exactly the two ideal classes;
  - and the repository's earlier 13/37 canon (THM-4566's χ = 2 list 1, 5, 13, 37, 85, and the prime idoneal numbers
    2, 3, 5, 7, 13, 37) is this fact seen from the Hilbert class field.

## 1. The construction for every d (THM-4825 (1)–(2))

In the coordinates (u, v) of u + v√−d, the owner's picture is G_d = Z² ∪ (Z² + (½, ½)), equivalently
2G_d = {u ≡ v (mod 2)}, the checkerboard. For odd squarefree d there are three cases.

| d | G_d is | the centred point (1 + √−d)/2 | example |
|---|---|---|---|
| ≡ 3 (mod 4) | the ring of integers O_K (Z[√−d] is the order of conductor 2) | an algebraic integer | d = 3: Z[ω]; d = 7, 11, 19, 43, 67, 163 (class number 1) |
| 1 | ½(1 + i)Z[i], the square grid turned by 45° | ½(1 + i) | the Gaussian case |
| ≡ 1 (mod 4), > 1 | ½𝔭₂ with 𝔭₂ = (2, 1 + √−d) the ramified prime over 2; 𝔭₂² = (2); [𝔭₂] has order 2 | norm (1 + d)/4 ∉ Z | d = 5, 13, 17, 21, 29, 33, 37, … |

So the owner's "another coordinate grid on top of the first grid, centered at (½,½)" has two arithmetic meanings:
- the **ring completion**, which is the Eisenstein reading;
- or the **non-principal class over 2**.

They are separated by d mod 4, that is by whether 2 is inert or split (d ≡ 3) or ramified (d ≡ 1) in Q(√−d). In the
ideal-class case the normalised norm form of the second grid is 2x² + 2xy + ((1 + d)/2)y², the second reduced form of
discriminant −4d. *Proof:* THM-4825; every odd squarefree d ≤ 199 is checked (script §A).

## 2. Thirteen and thirty-seven: the two grids are the class group (THM-4825 (3)–(5))

**Theorem.** For squarefree d ≡ 1 (mod 4), the class group of Z[√−d] consists of exactly the two grids (the principal
class and the class of 𝔭₂) iff h(−4d) = 2, iff **d ∈ {5, 13, 37}**.
- **Proof.** [𝔭₂] has order 2, so h(−4d) = 2 iff Cl = {1, [𝔭₂]}.
- **The list.** The imaginary quadratic fields of class number 2 are completely known (Baker 1971; Stark 1975; OEIS
  A014603). Of their 18 discriminants, exactly −20, −52 and −148 have the form −4d with d ≡ 1 (mod 4).
- **Checks.**
  - An exact reduced-form search over all 40,520 squarefree d ≡ 1 (mod 4) below 200,000 finds h(−4d) = 2 only at
    5, 13, 37.
  - The reduced forms are (1, 0, 5) and (2, 2, 3); (1, 0, 13) and (2, 2, 7); (1, 0, 37) and (2, 2, 19).

**Three grids agree (corollary).** For d ∈ {5, 13, 37}, an odd prime q ∤ d that splits in Q(√−d) satisfies

    q = x² + dy²  ⟺  q ≡ 1 (mod 4)  ⟺  q = a² + b²,

and otherwise q = 2x² + 2xy + ((1 + d)/2)y².
- **Two other readings.** The square grid Z[i] and the first grid Z[√−d] see the same split primes, and the centred grid
  sees the rest.
- **Why it holds.** H(Q(√−d)) = Q(√−d, i). This is THM-4566 (owner-directive canon, χ = 2 for m = 1, 5, 13, 37, 85) and
  the S40 note's genus-field criterion ([binary information §2](binary_information_class_groups_largest_anomalies_20261017.md)).
- **Checks.** All 8,962, 8,939 and 8,938 split primes below 200,000 for d = 5, 13, 37 (script §C).
- **The converse.** Among squarefree d ≡ 1 (mod 4) the agreement holds only for d ∈ {1, 5, 13, 37}, by Bauer's
  theorem. The least split prime q ≡ 1 (mod 4) that is *not* x² + dy² is:
  - 13 for d = 17;
  - 5 for d = 29, 41, 61 and 21;
  - 13 for d = 53;
  - 37 for d = 85;
  - 17 for d = 33;
  - 29 for d = 65;
  - 37 for d = 73.

**How this sits in the repository's 13/37 record.**
- THM-4566's χ = 2 list is {1, 5, 13, 37, 85} (the idoneal m ≡ 1 (mod 4) with no prime factor ≡ 3 (mod 4)).
- The prime idoneal numbers are 2, 3, 5, 7, 13, 37. In each family, 5, 13, 37 are the members ≡ 1 (mod 4) for which
  the class group is exactly the two grids.
- 85 = 5·17 is idoneal but has h(−340) = 4: four genera, and the two grids are only half of them.

**The owner's 5, 10, 25, 50, 100 (THM-4825 (5), IDENTITY).** Every prime factor of these numbers ramifies in Q(√−5).
So each is the norm of exactly one ideal, 𝔭₂^a·(√−5)^b, and lies on exactly one grid: the first iff a is even.

| n | 5 | 10 | 25 | 50 | 100 |
|---|---|---|---|---|---|
| x² + 5y² (first grid) | 0² + 5·1² | — | 5² + 5·0² | — | 10² + 5·0² |
| 2x² + 2xy + 3y² (second grid) | — | 2·1² + 2·1·(−2) + 3·(−2)² | — | 2·5² | — |

So the doubling chain flips the grid at each factor 2, and the owner's "lone 50" is a second-grid number. The same
alternation holds for 13, 26, 52 in Q(√−13) and for 37, 74, 148 in Q(√−37) (script §C). No link to the fourth-power
gap at 50 is claimed.

**11, 12, 13 as values of d (ANALOGY).**
- d = 11 (≡ 3 mod 4): the second grid completes Z[√−11] to the class-number-one ring Z[(1 + √−11)/2].
- d = 13 (≡ 1 mod 4): the two grids are the class group.
- d = 12 = 4·3 is not squarefree: Z[√−12] is the conductor-4 order of Z[ω], and its (½,½) copy is neither a ring nor
  an ideal.

So the owner's consecutive triple straddles the dichotomy. This is a reading, not a theorem.

## 3. Which grid carries a norm shell (THM-4826 (1)–(2))

**The mechanism.** Let O be Z[ω] or the Hurwitz order H. The prime 𝔓 over 2 has residue field F₄ (2Z[ω], respectively
(1 + i)H = D₄). The first grid is L = Z + 𝔓, which is Z[√−3], respectively the Lipschitz quaternions.
- Elements of odd norm are units mod 𝔓, and those on the first grid are the class 1 ∈ F₄*.
- A unit of order 3 (ω, respectively ½(−1 + i + j + k)) cycles F₄*.
- So every odd shell splits 1 : 2 between the first and the centred grid, and even shells lie on the first grid.
- In counts:
  - Z[ω]: 2Σχ₋₃ against 4Σχ₋₃ (checked for n ≤ 20,000);
  - Hurwitz: Jacobi's 8σ(n) against 16σ(n), out of Hurwitz's 24σ(n) (checked for n ≤ 300).

The **3 is the unit index**: [Z[ω]* : Z[√−3]*] = 6/2 and [H* : Lipschitz*] = 24/8. In the class-number formula it is
the factor that makes h(−12) = 1 although the conductor-2 factor is 3.

**The rest of the family Z[(1 + √−d)/2] (d ≡ 3 mod 4).**
- **d ≡ 7 (mod 8).** 2 splits, so O/2O = F₂ × F₂. Every odd element is ≡ (1, 1), so every odd shell lies on the first
  grid (checked for d ≤ 47, odd n ≤ 3000).
- **d ≡ 3 (mod 8), d > 3.** 2 is inert but the units are ±1.
  - Conjugation is the Frobenius on F₄, so a split-prime shell {±π, ±π̄} lies entirely on one grid, the first iff
    p = x² + dy².
  - Which grid is decided by the ring class field of conductor 2, of degree h(−4d)/h(−d) = 3. Chebotarev gives first-grid
    share 1/3.
  - Measured below 3·10⁵: 0.3325 (d = 11), 0.3334 (19), 0.3343 (43), 0.3336 (67).
  - The order class-number ratio h(−4d)/h(−d) is checked for every squarefree d ≡ 3 (mod 4) below 200: 1 at d = 3, 3
    for d ≡ 3 (mod 8) with d > 3, and 1 for d ≡ 7 (mod 8).
- **Conclusion.** The owner's d = 3 is the only member whose centred grid carries a fixed positive share (2/3) of every
  odd shell. Elsewhere the share is 0, or 0 on a third of the prime shells and 1 on the rest.

## 4. The grid of a prime is the cubic character of 2 (THM-4826 (3))

**The statement.** Let π = a + bω be primary (π ≡ 2 mod 3) over a prime p ≡ 1 (mod 3). Then:
- π lies on the first grid (b even) iff π ≡ 1 (mod 2);
- iff 2 is a cubic residue mod p (the supplementary law; Ireland–Rosen Prop. 9.6.1–2, restated as Thm 4.2 of
  [arXiv:2208.02454](https://arxiv.org/abs/2208.02454));
- iff p = x² + 27y² (Gauss; Cox Thm 4.15).

So the cubic residue symbol of 2 *is* the grid of the primary prime. Its values are the merge charges 1, ω, ω².

**Data.** All 39,231 primes p ≡ 1 (mod 3) below 10⁶ satisfy the three-way equivalence, with no mismatch. 2 is a cubic
residue for 13,032 of them (0.3322, against Chebotarev's 1/3). Among p < 50, the first grid holds 31 = 2² + 27 and
43 = 4² + 27. The centred grid holds 7, 13, 19 and 37.

**Mersenne primes.** For p = 2^q − 1 with q > 3 prime, 2 is a cube (2 = (2^m)³ with 3m ≡ 1 mod q). So the primary prime
lies on the first grid. Cornacchia gives:

| q | 2^q − 1 = x² + 27y² |
|---|---|
| 5 | 2² + 27·1² |
| 7 | 10² + 27·1² |
| 13 | 46² + 27·15² |
| 17 | 362² + 27·1² |
| 19 | 298² + 27·127² |
| 31 | 46162² + 27·783² |
| 61 | 1505304098² + 27·38443119² |

The table continues to q = 89, 107, 127 in the `.out`. q = 3 is the exception (7 is not x² + 27y²). That y = 127 at
q = 19 is NUMEROLOGY. The connection to the repository's Mersenne threads ({7, 21}, Mersenne switches) is a reading:
the Mersenne primes are exactly the Eisenstein primes where "doubling is a cube".

## 5. Fourth powers on the two grids (THM-4826 (4)–(6), HYP-9429)

**The centred units complete the quartic designs.**

| unit set | 2-design | 4-design | 6-design |
|---|---|---|---|
| Z[√−3]* = {±1} (first grid) | no | no | no |
| centred units {±ω, ±ω²} | no | no | no |
| Z[ω]* (hexagon) | yes | **yes** | no |
| Lipschitz units ±1, ±i, ±j, ±k (16-cell) | yes | no | no |
| centred units ½(±1 ± i ± j ± k) (tesseract) | yes | no | no |
| Hurwitz units (24-cell) | yes | **yes** | no |
| Gaussian units (square) | yes | no | no |

A 4-design is what makes Σ_u ⟨u, x⟩⁴ proportional to |x|⁴:
- **Candido.** For the hexagon, in the coordinates p = ⟨u₀, x⟩ and q = ⟨u₁, x⟩, the identity is
  p⁴ + q⁴ + (p + q)⁴ = 2(p² + pq + q²)². This is Candido's identity (1951; the special case goes back to Proth 1878).
  It is mac-mini's zero-sum quartic identity (THM-4824 (F)).
- **Liouville.** For the 24-cell, left multiplication by 1 + i carries the Hurwitz units onto the 24 roots of D₄. The
  identity becomes Liouville's 6(Σx_i²)² = Σ_(i<j) [(x_i + x_j)⁴ + (x_i − x_j)⁴]. With Lagrange's four squares, this gives
  Liouville's bound g(4) ≤ 53. The true value g(4) = 19 is due to Balasubramanian, Deshouillers and Dress (1986).

So the owner's two topics, the centred grid and sums of fourth powers, meet historically: the first general bound for
fourth powers came from the four-dimensional two-grid ring. The square grid's four units are not a 4-design, and there
is no Gaussian quartic identity.

**Eisenstein fourth powers mod 16 (THM-4826 (5)).** α⁴ ≡ 0 (mod 16) for α ∈ 2Z[ω], and otherwise α⁴ ≡ u or 9u with u ∈
μ₃. This is the Eisenstein analogue of "a⁴ ≡ 0 or 1 (mod 16)". Also μ₆⁴ = μ₃: the charges 1, ω, ω² of THM-4820 are
exactly the fourth powers of the six units.

**The obstruction (THM-4826 (6)).** A rational n ≡ 5 (mod 8) is never a sum of four Eisenstein fourth powers (proof:
a residue count; THM-4826). Over Z the four-fourth-power residues are 0..4 (mod 16). Over Z[ω] every rational residue
except 5 and 13 (mod 16) is reachable, so **the second grid repairs every congruence obstruction except n ≡ 5 (mod 8)**.
- **The owner's numbers stay excluded.** 13 ≡ 13 and 37 ≡ 5 (mod 16), so neither is a sum of four fourth powers over Z
  or over Z[ω]. In 0..50 the excluded class is 5, 13, 21, 29, 37, 45.
- **No other congruence obstruction.** Tests modulo 2⁸, 3⁵, 5³ and 7² find none on 0..200. Composite moduli add nothing,
  by the Chinese remainder theorem.

**What is representable (FINITE, a lower bound; HYP-9429).**
- **0..50.** With bases |a|, |b| ≤ 20, every n ≢ 5 (mod 8) is a sum of four Eisenstein fourth powers except 6.
  - Cancellation matters: 47 = (−20 − 10ω)⁴ + (−13 − 17ω)⁴ + (−8 − 15ω)⁴ + (−8 + 6ω)⁴.
- **6.**
  - **No representation.** The mod-8 count forces two bases ≡ ω and two ≡ ω² (mod 2), each pair with fourth powers
    ≡ u and 9u. The exact box search (|a|, |b| ≤ 150, 3.2·10⁷ × 3.2·10⁷ pair-sum combinations, 14 s) finds **no**
    representation.
  - **No local obstruction.** 6 is representable modulo every prime power tested, and in Z₂[ω] by Hensel.
- **51..200 with bases ≤ 60.** 8 non-obstructed values lack a representation: 54, 60, 75, 78, 102, 124, 140, 186.
- **Heuristic.** The problem has critical dimension (eight real variables, two real equations of degree 4). The expected
  number of representations with bases ≤ B grows like c(n)·log B, as for sums of three cubes. So sparse misses are
  expected and are not obstructions.
- **Pre-registered results** (HYP-9429, each filed before its run):
  - (1) "6 has none with bases ≤ 150": **CONFIRMED**.
  - (2) "at most three misses in 51..200 with bases ≤ 60": **REFUTED** (eight).
  - (3) "at least three of the eight acquire one with bases ≤ 120": **REFUTED**. Only 102 and 186 do:
    - 102 = (113 + 99ω)⁴ + (51 − 41ω)⁴ + (104 + 31ω)⁴ + (106 + 13ω)⁴;
    - 186 = (37 − 44ω)⁴ + (16 + 40ω)⁴ + (54 + 36ω)⁴ + (77 + 14ω)⁴.
- **Post-hoc density check.** Exact 2-adic × 3-adic local densities (four-fold convolution modulo 128 and 243) put
  6 of the 7 persistent misses below the median.
  - 124 and 140 are in the lowest class (0.444). There n ≡ 4 (mod 8) forces all four bases ≡ 1 (mod 2).
  - Every other member of that class up to 200 is represented.
  - So the misses look like sparse representations at low local density rather than a congruence. Whether they are
    representable at all is OPEN (HYP-9429). A global obstruction (integral Brauer–Manin or other) has not been looked
    for.

## 6. Typing

| claim | type |
|---|---|
| the second grid is the ring of integers (d ≡ 3 mod 4) or ½ the ramified prime over 2, of class order 2 (d ≡ 1 mod 4, d > 1) | PROVED (THM-4825 (1)–(2)) |
| the two grids are the class group iff d ∈ {5, 13, 37} | PROVED modulo the CITED class-number-2 list (Baker; Stark); FINITE-EXACT below 2·10⁵ |
| three grids agree (Gaussian = first grid on split primes) iff d ∈ {1, 5, 13, 37} | PROVED (genus theory, THM-4566, Bauer); FINITE-EXACT below 2·10⁵ |
| 5, 10, 25, 50, 100 alternate between the grids of Z[√−5] | IDENTITY (THM-4825 (5)) |
| the 1 : 2 split for Z[ω] and Hurwitz; d mod 8 behaviour; only d = 3 spreads every odd shell | PROVED (THM-4826 (1)–(2)); counts CLASSICAL (Jacobi, Hurwitz) |
| the grid of the primary prime is the cubic character of 2 | CLASSICAL (Gauss; Ireland–Rosen), restated in grid language |
| Mersenne primes 2^q − 1 (q > 3) lie on the first grid | PROVED (elementary corollary); the representations FINITE-EXACT |
| the centred units complete the quartic designs; Candido; Liouville via (1 + i)H = D₄ | CLASSICAL identities, PROVED design reading, FINITE-EXACT (sympy) |
| μ₆⁴ = μ₃ (the merge charges are the unit fourth powers); fourth powers mod 16 | PROVED (THM-4826 (5)) |
| four Eisenstein fourth powers never give n ≡ 5 (mod 8); 13 and 37 excluded on both grids | PROVED (THM-4826 (6)) |
| 6 has no representation with bases ≤ 150; 54, 60, 75, 78, 124, 140 none with bases ≤ 120 | FINITE (exhaustive in the stated boxes); OPEN beyond (HYP-9429; predictions (2), (3) REFUTED) |
| the misses sit in low local-density classes | POST HOC (not pre-registered); consistent with sparse representation, not proof |
| {5, 13, 37} ≡ 5 (mod 8) and the excluded class n ≡ 5 (mod 8) | NUMEROLOGY (different mechanisms: Rédei 4-ranks vs a 2-adic count) |
| 11, 12, 13 straddling the ring/ideal dichotomy | ANALOGY |
| the 13/37 coincidences of the fourth-power window, the square-later census, and ball counts | NUMEROLOGY (typed by mac-mini, THM-4824 / HYP-9428) |

## 7. Reproduction

- `python 04-computation/experiments/two_grids_class_group_20261027.py` (about 10 s; numpy-free except sympy for §F).
  Its sections:
  - A: ring or ideal class;
  - B: h(−4d) = 2;
  - C: three grids and the owner's chains;
  - D: shells;
  - E: the cubic character, Mersenne;
  - F: designs, identities, fourth powers mod 16.
- `gcc -O2 eisenstein_four_fourth_powers_20261027.c`, then the following modes:
  - `all 20 0 50` (the window);
  - `pat 30 22 1 9` and `box 30 22` (validations);
  - `box 150 6` (HYP-9429 (1));
  - `all 60 51 200` (HYP-9429 (2));
  - `all 120 51 200` (HYP-9429 (3)).
- `python 04-computation/experiments/eisenstein_fourth_power_obstructions_20261027.py` (+ `.out`, about 40 s). It
  covers the modular obstructions (2⁸, 3⁵, 5³, 7²) and the local densities of the misses.
