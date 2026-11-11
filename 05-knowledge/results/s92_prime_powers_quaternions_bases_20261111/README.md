# S92 (2026-11-11): F1 at every odd modulus, the four-square series, and the algebras, bases and Collatz tables around it

Session `opus-2026-11-11-S92`. Canon: THM-5134, THM-5135. Hypotheses: HYP-9665, HYP-9666. Scripts and outputs:
`04-computation/experiments/s92_prime_powers_20261111/`. The Chan-Pak conjecture (every odd number except 7 and 21
is a Hamiltonian path count) is NOT proved here and remains OPEN: every local density below is positive, at 7 and
21 as elsewhere, and decides nothing about it.

The owner's directive had seven parts. Short answers, with the section that holds the detail:

| asked | answer | section |
|---|---|---|
| the odd prime-power case of THM-5128, the series `sigma(n)/n` | PROVED at every odd modulus (THM-5134); pre-registered test confirmed | 1, 2 |
| quaternions versus bicomplex numbers, dihedral groups | prior work listed; the eight extensions of `F_2^2` put the four algebras in one table (THM-5135 (5)) | 3 |
| `(8/pi^2) sigma(n)/n` and the octonions / Fano plane | the series is the quaternion (rank 2) member of a ladder whose rank-4 member is E8, `240 sigma_3(n)`; the repository's octonion work lives at the prime 2, where F1 has no factor; no tournament family is known to carry the rank-4 law | 4 |
| Tao's toy example, all bases | it is the money-plant puzzle; Session266 recorded and generalised it the same morning; the F1 joins are such quadratics, and their local law holds in every base (THM-5134 (H)) | 5 |
| multiplication tables modulo odd values and arcs of `3n + 1`, `3n - 1` | the dictionary exists (credited); added: the F1 density is the additive square of the table, and three exact relabellings | 6 |
| three negative cycles against one positive, and the 4-vertex tournaments | exact pieces and one coincidence, typed; the asymmetry is archimedean and no local density can see it | 7 |
| the pasted note on strong components as primes | profiles and generating function confirmed; all six strong classes of order 5 are self-converse | 8 |

A sweep of the repository by eight independent readers (quaternions / bicomplex, dihedral, octonions / Fano,
tables and Collatz, negative cycles and 4-vertex tournaments, bases and digits, divisor sums, the packets of the
day) preceded the write-up; what it found is credited in section 10. Section 9 types every bridge.

## 1. The law at every odd modulus (THM-5134)

F1: `N = x1 x2 + 2^s y1 y2 - [s = 0]` over pairs of one-apex states `(x, y) = (2 psi(h) + 1, 2h + 1)`; sheet 0 is
`s = 0`, sheet 1 is `s >= 1`. THM-5128 proved the class law by value at odd primes and squarefree moduli and left
odd prime powers numerical. THM-5134 closes that step. For an odd modulus `m`:

- the type of a triple is `d = gcd(x1, y1, m)` for its smaller state; a triple of type `d` has `d | N - c_sigma`
  exactly (`c_1 = 0`, `c_0 = -1`), and its value is uniform on those classes with a power saving;
- the types have limiting shares `c_d = d^-2 prod_{p | m/d} (1 - p^-2)`;
- so the class factor tends to `L_m(e) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)`, `e = N - c_sigma`, which
  is multiplicative, with `L_(p^k) = (1 - p^-2) sigma(p^v)/p^v` at valuation `v < k`;
- the product over the odd primes is `S(N) = (8/pi^2) sigma(n)/n`, `n` the odd part of `N` (sheet 1) or `N + 1`
  (sheet 0).

The factor for uniform states and the series were written down by the audit of THM-5128; the fixed-source law
without a rate is Session206's. New: the power saving, the sum over the smaller state, hence the by-value law.

Data (exact censuses to `2^28`; 854,415,478 and 1,083,261,684 triples on the two sheets). The remainder of the
first-order law times `m`, at `2^28`: `0.02` at 9, `0.03` at 27 and 81, `0.004` at 25, `0.005` at 125, `0.10` at 49,
`0.04` at 45, `0.03` at 75; it decreases slowly (7 is the slow prime, as in THM-5128). The factors themselves are
far from their limits at this size (class 0 modulo 9: `1.45` against `35/27 = 1.296`), as the `1/ln X` term
predicts. Pre-registered before the run (commit 73c8437aac): window `[2^27, 2^28)`, eight moduli; P1 (a-priori
factors from the smaller states alone, within tolerance), P2 (flatness inside each group), P3 (closer to the
first-order value than to the limit, 14 of 14) CONFIRMED; P4 (remainder smaller at `2^28` than at `2^26` in at
least 8 of 12 pairs) CONFIRMED with exactly 8.

What it does not give: an asymptotic for the multiplicity of one `N`. On `[2^24, 2^25)` the multiplicity correlates
with `sigma(n)/n` (`0.40`, `0.43`) and doubles across its quintiles (exploratory; HYP-9665 (B)).

## 2. Three classical faces of `sigma(n)/n`

1. *The multiplication table.* The number of `(a, b)` modulo `m` with `a b = t` is
   `N_m(t) = sum_{d | gcd(t, m)} d phi(m/d)` (in the modular-centres packet, P4). The F1 seam is a sum of two table
   entries (`2^s` is a unit), and the additive square of the table is `(N_m * N_m)(c) = m^3 L_m(c)` (checked for
   eleven odd `m <= 105`). So `L_m` is the density of a sum of two entries of the table modulo `m`.
2. *Ingham (1927, CITED).* `#{a1 a2 + b1 b2 = n} = sum_{m < n} d(m) d(n - m) ~ (6/pi^2) sigma(n) ln^2 n`. Its series
   `sigma(n)/(zeta(2) n)` without the factor at 2 is `(8/pi^2) sigma(n)/n`. F1 is a thin analogue; the theorem is not
   used. (Illustration: on `[20000, 20400)` the correlation of `D(n)/n` with `sigma(n)/n` is `0.979`.)
3. *Jacobi (CITED).* `r_4(n) = 8 sigma(n) = pi^2 n * S` for odd `n`: `S` is the product of the densities of the
   four-square form at the odd primes. Proof of the link: `a1 a2 + u b1 b2` and the sum of four squares are
   unimodular of the same discriminant over `Z_p`, `p` odd (class counts equal modulo 9, 27, 25, 49, 121).

Divisor by divisor: `sigma(n)/n = sum_{d | n} 1/d`, and `d` is the type of the smaller tournament (its count and
its address share the factor `d`), present with density `d^-2` and concentrated on a share `d` of the classes.

## 3. Quaternions, split quaternions, bicomplex numbers, dihedral groups

Prior work found (all credited in section 10):

- THM-2887: 16 normalized cocycles on `F_2^2`, 8 with determinant commutator, three Arf-zero `D_8` classes and one
  Arf-one `Q_8` class. The four commutative classes are counted there and not named.
- THM-5133 (D1) (kind-pasteur S23, same day, PRE-AUDIT): an antisymmetric sign table on `Z_2^n` is a tournament;
  on a line, cyclic = Hamilton's quaternions, transitive = split quaternions with the unit of square `-1` in the
  middle.
- THM-024, HYP-8205, THM-4818, THM-4743: the group `Aut+-(T)` of automorphisms and anti-automorphisms is a split
  extension of an odd group by `Z_2` ("no quaternionic option"); anti-automorphisms have order `2 (mod 4)`;
  `D_5` for the regular 5-tournament, `D_3 = S_3` for the 3-cycle. THM-127, THM-131, THM-586: the dihedral group
  `D_(2p)` on the Hamiltonian paths of the Paley tournament.
- THM-4826, Session116: Jacobi's `8 sigma(n)` and Hurwitz's `24 sigma(n)`, for grids of the Hurwitz order and for
  sedenion zero divisors; never for path counts. MISTAKE-506: Hamilton's algebra ramifies only at 2 and infinity.
- The bicomplex numbers (tessarines, `C` tensor `C`) appear nowhere: in tracked files "bicomplex" always means a
  double complex.

Added (THM-5135 (5), all by enumeration; classical mathematics):

| units of square `-1` | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| algebra | `R^4` | split quaternions `M_2(R)` | bicomplex `C + C` | Hamilton |
| group of signed units | `C_2^3` | `D_4` | `C_4 x C_2` | `Q_8` |
| classes | 1 | 3 | 3 | 1 |
| sign table | signed triangle | transitive tournament | signed triangle | cyclic tournament |
| norm | - | determinant = F1 seam | product of two sums of two squares | sum of four squares |

- The class is the set of units of square `-1`; the algebra is commutative iff that set is even. So quaternions
  and bicomplex numbers share two imaginary units and differ in the third; the dihedral class sits one step
  below the quaternion class, with one imaginary unit.
- The non-trivial coboundary reverses the three arcs: a class is a tournament up to converse, or a signed
  triangle up to a global sign. Commutative algebras have no tournament.
- `r_4 = r_2` additively convolved with itself (quaternion norms, `8 sigma(n)` at odd `n`); the bicomplex norm
  count is the multiplicative convolution of `r_2` with itself and vanishes for `n = 3 (mod 4)`. In the join
  calculus the ordered join multiplies counts (rank 1, divisor-function law) and the one-arc join adds two
  products (rank 2, the law `sigma(n)/n`).
- The split-quaternion norm (F1, transitive triangle, `H = 1`) and the Hamilton norm (cyclic triangle, `H = 3`)
  have the same law at every odd prime (THM-5135 (2)); they differ at 2 and at infinity. Orientation, i.e. what
  THM-5133 calls flip against spin, is invisible to the odd local densities; so is the signature of a two-port
  kernel (THM-5135 (4)).
- Two different "D versus Q" statements should not be confused: as symmetry groups of a tournament both `D_4`
  and `Q_8` are impossible (`|Aut|` is odd); as sign tables `D_4` is `TT_3` and `Q_8` is `C_3`.

## 4. The octonions and the Fano plane

Prior work: the octonion sign table is the Paley tournament `P_7` on the Fano plane (2026-09-17 / 09-22 notes: 16
alternative orientations with `H = 189`, 112 others with `H = 171`; THM-448 and THM-450: the doubling tower;
THM-5133 (E)-(G) and Session265 the same day: Fano and `PG(3,2)` parities, the line cube). All of it is about
signs, i.e. the prime 2. THM-5128 (E) shows that F1 has no factor at 2. So the octonion material and the series
`(8/pi^2) sigma(n)/n` do not meet on the tournament side. They meet in the norm forms:

- rank 1 (`a b`): the divisor function, no series;
- rank 2 (`a1 a2 + b1 b2`, split quaternions): `sigma(n)/(zeta(2) n)`; definite model Hamilton / Hurwitz, ramified
  at 2, whence `8 sigma(n)` and `24 sigma(n)` at odd `n` only;
- rank 4 (split octonions): `sigma_3(n)/(zeta(4) n^3)`; its definite model is E8, the integral octonions, unramified
  at every prime: `(8 pi^4/3) n^3 * sigma_3(n)/(zeta(4) n^3) = 240 sigma_3(n)` for every `n` (THM-5135 (3); the E8
  form and the rank-4 seam have equal class counts modulo 2, 4, 8; shells `240, 2160, 6720` recounted). `Z^8` is
  not that model (MISTAKE-484).

Ranks 1, 2, 4 are the split composition algebras (Hurwitz); the laws hold for every rank and the join calculus
also has rank 3, so the ladder is a reading, not a constraint. Rank-4 forms do occur (the K-form of THM-4916, the
fusion corner of Session264, both unimodular), but the one with a proved law (THM-4916) has factor
`1 + (3p - 2)(p - 1)/(2p^3)` at the class 0, not `1 + p^-3 - p^-4`: its states are confined. Whether a tournament
family carries the E8 law is HYP-9666 (A).
`240 = 7!/21` labelled Paley tournaments against 240 roots of E8: both are `16 x 15`; the repository already
records that no equivariant bijection is known (fano_code, section 4). COINCIDENCE.

## 5. Tao's toy example, and every base

The example was not attached to the message. It is the money-plant puzzle (The Guardian, 26 October 2026), as
recorded by Session266 (`session266_fiber_radix_20261111/radix/PROOF.md`, committed while this session was
running; the newspaper page could not be fetched from here): `C_t = a t^2 + b t + c`, and `100 a + 10 b + c = 427`
with digits gives `(4, 2, 7)`. Session266 proves the all-base version (injective iff the day is at least the
base), the nilpotent algebra `Q[eps]/(eps^3)` behind the recursion and "carries are the kernel `(x - B)`".

What this session adds:

- *The F1 joins are such quadratics.* The odd-cycle packing polynomial of `J(h1, h2, s)` is
  `1 + alpha_1 z + alpha_2 z^2`, `alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2`, `alpha_2 = psi1 psi2` (THM-5070 (2) in
  these coordinates; 314 joins recounted). The path count is its value at 2, the base-`b` packing count
  `N_b = (1 + b psi1)(1 + b psi2) + b (2^s y1 y2 - [s = 0])/2` its value at `b`.
- *The digit bound.* For `b > max(alpha_1, alpha_2)` the base-`b` numeral of `N_b` is `(alpha_2)(alpha_1)(1)`: the
  count displays the packing numbers (the block lemma of THM-4940 is the case `alpha_2 = 0`). In base 2 this
  almost never holds and carries collide: 17 is the count of joins with `(alpha_1, alpha_2) = (6, 1)` and `(8, 0)`,
  45 of three vectors. Of 10,400 small joins the share with digits at most `b - 1` is `0.0002, 0.0004, 0.0014, 0.022`
  for `b = 2, 3, 10, 100`.
- *The local law in every base* (THM-5134 (H)). At every odd modulus prime to `b` the base-`b` counts of F1 have
  the same class law, with the distinguished classes `0` (sheet `s >= 1`) and `-b/2` (sheet 0); `N_b = 1 (mod b)`.
  The rates are those of the binary digit chain and do not depend on `b`: 7 is the slow prime in every base.
  Apart from `h = 0`, balanced states `X = y` exist only for `b = 2` (the fans) and `b = 3` (`h = 3, 45, 117, 309, ...`;
  `X = y = 7` at `h = 3`, which is why 7 is boosted early in base 3).
- Two radices are in play and are independent: the base `b` of the numeral and the binary radix of the address
  `h` (also noted by Session266 and visible in THM-5085).

## 6. Multiplication tables modulo odd `m` and the arcs of `3n + 1`, `3n - 1`

Prior work (credited): THM-4792 (the centre labels `+-1/2` of the table are the fixed points of `3x + 1`, `3x - 1`;
alternation iff `m | 2^t + 1`, Moree); Session255 (every arc of either map is the centre block of one table,
`m = 4x + sigma`); the central-modular packet (`m -> 3m + 2 sigma`, `m -> (m + sigma)/2`); modular centres (cycle labels
`{5, 17, 9}`, `{3, 7}`, `{19, 55, 27, 79, 39}`); codex WEIGHTED_ATLAS and PUNCTURED_ATLAS (the product law, the triangle
`J(x, sigma) = (3x + sigma, -sigma)`, 7 and 21 as the arcs `2 -> 1` of `3n - 1` and `5 -> 16` of `3n + 1`).

One coordinate puts the two sheets together (a restatement; bridges_exact.out S7). Every odd `m >= 1` is `|4n + 1|`
for exactly one integer `n = (chi_4(m) m - 1)/4`; `m` labels the arc of `3n + 1` on `Z` leaving `n`. Then: an odd
step is `m* -> 3 m* + 2` and an even step `m* -> (m* + 1)/2` for `m* = 4n + 1`; multiplying labels is
`n * k = 4nk + n + k`; tripling a count is `n -> -3n - 1`. The sheet is the sign of `n`, that is `chi_4(m)`.

Added:

- *The density is the additive square of the table* (section 2, item 1).
- *The favoured classes in arc coordinates.* For a count `N = |4n + 1|` and an odd prime `p`: `p | N` iff
  `n = -1/4 (mod p)`; `N = -1 (mod p)` iff `n = -1/2` (for `n >= 0`) or `n = 0` (for `n < 0`). So the tail sheet of F1
  favours the arcs leaving the class `-1/4`, and the head sheet the arcs leaving the class of the fixed point
  `-1/2` of `3n + 1` (positive side) or the class 0 (negative side). The argument of the head series is `2n + 1`
  for `n >= 0` and the odd part of `n` for `n < 0`. A relabelling; no dynamics.
- *Fan primes are the alternation primes.* A prime has a balanced `p`-divisible state iff `p | 2^t + 1` for some
  `t` iff `ord_p(2)` is even: the condition of THM-4792 (c). The least `t` is `ord_p(2)/2` and the first such state
  is the fan `h = 2^(t-1)`; it is the first `p`-divisible state for `p = 3, 5, 11, 13, 17, 43, 113` (`p < 200`).
- *Lines.* `c_-(N)` = odd part of `N - 1` (line `1 + c 2^s`, THM-5126) and `c_+(N)` = odd part of `N + 1`
  (anti-line, THM-5130): `c_-(3n) = T_-(n)` and `c_+(3n) = T_+(n)`, the Syracuse maps of `3n - 1` and `3n + 1`. The
  join with a triangle moves a count along a Syracuse arc of its line label (a variant of codex's `J`).

## 7. Negative against positive cycles, and the four tournaments on four vertices

Prior work: THM-4472 (the two converse classes with `H = 3` are the forward map and the inverse tree of one pair:
converse is time reversal composed with negation, and the reading "the two classes are the two sheets" is
REFUTED there); the sheet conjugacy `x -> -x` between `3n + 1` on the negatives and `3n - 1` on the positives, under
which every residue-level datum is invariant (THM-4474, THM-4521, the functional-uniqueness note); THM-4849
(circulant tournaments are rational cycles; reversal fixes the cycles of `1`, `-1`, `{-5, -7}` and moves that of
`-17`); THM-5131 (the 48 labelled self-converse 4-tournaments have `H = 1` or `5`, the 16 others `H = 3`).

Exact statements (bridges_exact.out S7, small_tournaments_collatz_lines.out):

- The five known cycles of `3n + 1` on `Z` (`{0}`, `{1, 4, 2}`, `{-1, -2}`, the cycle of `-5`, the cycle of `-17`) have 29
  arcs with 29 distinct labels; exactly one label, 7 (the arc `-2 -> -1`), is not a path count. The positive
  trivial cycle has labels `5, 17, 9`, all attained; the negative trivial cycle has labels `3, 7` and is cut at 7.
  The other missing count, 21, is the arc `5 -> 16`; the trunk arcs `(4^k - 1)/3 -> 4^k` have labels
  `(4^(k+1) - 1)/3 = 5, 21, 85, ...`. (The cut at 7 and the arc of 21 are in PUNCTURED_ATLAS; the count over all five
  cycles, including the 18 labels of the cycle of `-17`, is new.)
- The path counts `1, 3, 5` of the 4-vertex tournaments are the labels of the arcs leaving `0, -1, 1`: `0 -> 0`,
  `-1 -> -2`, `1 -> 4`, one in each of the three cycles through `0, -1, 1`. At order 5 the integers covered are
  `[-4, 3]` without `-2`; at order 6, `[-11, 11]` without `-10, -9, -2, 5`.
- At order 4: self-converse iff `H = 1 (mod 4)` iff `n >= 0` iff the number of 3-cycles is even. The converse pair
  (`H = 3`) is the only class on the negative side. For even order in general only one direction is true:
  self-converse implies `H = 1 (mod 4)` (KNOWN, THM-4930, THM-5131).

What is not there: a map from cycles to tournaments. "Three cycles against one" and "four tournaments on four
vertices" agree only as counts: COINCIDENCE. And no congruence statement can explain three against one: the two
sheets are conjugate under negation at every finite place, so the asymmetry lives at the real place. The local
densities of this session are of that blind kind. The unit equations behind the cycles (`4 - 3`, `3 - 2`, `9 - 8`)
and the sign-flip criterion `2^delta - 3^(delta-1) = 1` of THM-5130 (B) are the same equations; that link is
kind-pasteur's.

## 8. The pasted note: strong components as primes

- Confirmed by enumeration: at order 4 the arc-usage profiles (number of Hamiltonian paths through each arc) are
  `{1,1,1}`, `{2,2,2,1,1,1}` twice (a converse pair), `{4,3,3,3,1,1}`; the classes by strong components have generating
  function `1/(1 - (x + x^3 + x^4 + 6x^5))`, giving `1, 1, 2, 4, 12` at orders 1 to 5.
- The closing question: the six strong classes of order 5 have `H = 9, 9, 11, 13, 15, 15` and automorphism groups
  of order `1, 1, 1, 1, 3, 5`; all six are self-converse. No converse pair is strong at order 5; of the 12
  classes 8 are self-converse and the other 4 are two converse pairs (non-strong).
- Prior work: the arc-usage function is Session265's `W(u, v)` (it determines the strong components and their
  order up to reversal); Session265 also corrects two statements of the note (a general vertex inflation does not
  multiply `H`; the total weight of a boundary, not each boundary arc, is `H`) and shows the limit of the analogy:
  63 is irreducible among path counts. Multiplicativity over strong components is THM-1370; the product law on
  table codes is codex's.
- In the language of section 4, "strong components as primes" is the rank-1 seam: for uniform factors its class
  law is the divisor-function law `(1 - 1/p)(v + 1)`. Strong counts are not uniform (7 and 21 never occur), so this
  is a comparison, not a theorem.

## 9. The bridges, typed

EXACT = a proved or enumerated identity; MECHANISM = a common cause with a stated loss; ANALOGY = same shape, no
map; COINCIDENCE = equal numbers only.

| bridge | type | evidence |
|---|---|---|
| F1 seam = four-square form at every odd prime | EXACT | THM-5134 (F); class counts modulo 9, 27, 25, 49, 121 |
| `S(N) = r_4(n)/(pi^2 n)` | EXACT (Jacobi CITED) | `r_4 = 8 sigma` for odd `n <= 400` |
| `S` = Ingham's series without the prime 2 | EXACT (CITED); F1 as a thin Ingham problem: ANALOGY | section 2 |
| `L_m` = additive square of the multiplication table | EXACT | eleven moduli |
| rank-4 law over all primes = E8 = integral octonions | EXACT for uniform states (CITED + enumerated at 2) | THM-5135 (3) |
| a tournament family with the rank-4 law | OPEN | HYP-9666 (A); THM-4916 is a negative instance |
| octonion sign tournaments against norm forms | MECHANISM (orientation at 2 and infinity, rank at odd primes) | THM-5135 Remark 1 |
| 240 Paley tournaments, 240 roots | COINCIDENCE | no equivariant bijection known |
| eight extensions = subsets of imaginary units; quaternion / bicomplex / dihedral | EXACT (classical) | THM-5135 (5) |
| quaternion norm additive, bicomplex norm multiplicative in `r_2` | EXACT | `n <= 400` |
| ordered join = rank 1, one-arc join = rank 2 | MECHANISM (laws for uniform states) | THM-5135 (1) |
| `Aut+-(T)` never contains `D_4` or `Q_8` | KNOWN | THM-024, THM-4818 |
| `Out(Q_8) = S_3 = Aut+-(C_3)` | EXACT (classical; from the sweep, not recomputed) | - |
| the sheet of a table label is `chi_4` of the count | EXACT (restatement) | S7 |
| Tao's quadratic = packing polynomial of a join; digit bound = no carry | EXACT | S11; THM-4940 |
| base-`b` law | EXACT (THM-5134 (H)) | censuses, six bases |
| tripling = `n -> -3n - 1`, products = `4nk + n + k` | EXACT, KNOWN in `(x, sigma)` form | Session255, codex |
| one missing label among the 29 arcs of the known cycles | EXACT (uses THM-1370 and THM-4875) | S7 |
| `1, 3, 5` = arcs at `0, -1, 1`; self-converse iff `n >= 0` at order 4 | EXACT (small) | S7 |
| three negative cycles and one positive against four 4-tournaments | COINCIDENCE | THM-4472; sheet conjugacy |
| fan primes = alternation primes of THM-4792 (c) | EXACT (arithmetic) | S8 |
| `sigma(n)/n` against the tournament abundancy `Pi_2/H` of kind-pasteur S22 | ANALOGY | `a(C_3) = 4`, `sigma(3)/3 = 4/3` |
| Session265's prime-power lifts, exceptions 5 and 761 | unrelated (a name clash) | THM-5134 Remark 6 |
| 7 at the minimum of the head series | COINCIDENCE (true of every `2^k - 1`) | THM-5134 (E) |

## 10. Credits and prior work

- THM-5128 with its audit (the factor for uniform states, the series, the enumeration at 9, 25, 27, 49), MISTAKE-761.
- Session206 (weighted parameter capacity), Session193 modular_sync, Session 25: the fixed-source law and the
  digit chain. Session223, THM-4916: valuation-resolved factors for root changes; the K-form.
- kind-pasteur: THM-5133 and HYP-9664 (S23, PRE-AUDIT; cited at that status), THM-5130 and HYP-9662 (S22),
  THM-5070, THM-5085 (base `x`), THM-5126 (fans, lines), the modular-centres packet.
- Codex lanes: Session255, Session260 (the norm form `x x' + 2 y y'`), Session261 (the kernel with determinant 329),
  Session263, Session264 (fusion law), Session265 (`W(u, v)`, the Fano cube, 63), Session266 (Tao's puzzle, radix),
  the central-modular, odd-tables and signed-tables packets (WEIGHTED_ATLAS, PUNCTURED_ATLAS, PHASE_CLOCKS).
- mac-mini: THM-5131, THM-5132, HYP-9663 (the mod-4 race, self-converse tournaments); THM-448, THM-450, THM-480.
- Older canon: THM-2887, THM-024, THM-4818, THM-4743, THM-127, THM-131, THM-586, THM-4472, THM-4792, THM-4826,
  THM-4849, THM-4930, THM-4940, THM-4756, THM-1370, THM-4875, Session112 (norm signs), Session116, the 2026-09-17
  and 2026-09-22 Paley / Fano / octonion notes, MISTAKE-484, MISTAKE-506.
- External (CITED): Jacobi (four and eight squares); Hurwitz (integral quaternions; `1, 2, 4, 8`); A. E. Ingham,
  J. London Math. Soc. 2 (1927), on `sum d(m) d(n - m)`; Coxeter's integral octonions and the theta series of E8
  (OEIS A004009); Gronwall; the classification of unimodular forms over `Z_p`.

An observation for kind-pasteur's audit, not tested here: `H(CD_4) = 198,335,025` with 21 automorphisms (THM-5133
(G)) equals `H(T_15)` with `|Aut| = 21` of THM-448, and the split-octonion value 45 and the split-sedenion value
24,453,597 equal values of THM-450 (7); the Cayley-Dickson sign tournaments are probably the June doubling tower.

## 11. Audit

Pending at the time of this draft; the outcome is recorded below when the audit is applied.

## 12. Reproduce; next obligations

    cd 04-computation/experiments/s92_prime_powers_20261111
    python3 odd_modulus_census.py cum 24            # exact remainder of the first-order law (2^28: about 25 minutes)
    python3 odd_modulus_census.py eval 28 27 prereg_window_28_27.json
    python3 seam_rank_exact.py; python3 bridges_exact.py; python3 every_base.py
    python3 small_tournaments_collatz_lines.py; python3 series_vs_counts.py 25 24

Next: (1) even moduli (HYP-9665 (A)); (2) uniformity in the modulus (HYP-9660 (B)), which is what a statement
about single values needs; (3) a rooted family with equidistributed rank-4 states, or a proof that none of the
known ones is (HYP-9666 (A)); (4) the second digit in bases 3 and 6 (HYP-9665 (C)); (5) the cross-reference
between THM-5133 (G) and THM-448 once THM-5133 is audited.
