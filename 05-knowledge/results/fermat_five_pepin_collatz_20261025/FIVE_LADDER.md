# The three fives as one ladder: {0,1,2,3,4} → {1,2,4,8,16} → {3,5,17,257,65537}

2026-10-25, klein. **EXACT** classical structure (Gauss–Wantzel, Lucas' theorem, Paley graphs),
verified in `fermat-normal.json`; the repository relations in §4 are typed STRUCTURE / COINCIDENCE
as in the earlier ledgers. No novelty claim for the classical parts.

## 1. The ladder

The three sets are three floors of one tower: an index `i in {0,1,2,3,4}`, its power `2^i in
{1,2,4,8,16}`, and the Fermat prime `F_i = 2^{2^i} + 1 in {3,5,17,257,65537}`. Each floor carries a
different structure:

- `{0,1,2,3,4}` is `Z/5` as a set of exponents; it is the discrete logarithm of the second floor.
- `{1,2,4,8,16}` is the subgroup of order five in `F_31^x` (`2` has order five modulo `31 = 2^5 - 1`),
  so the second floor is literally a cyclic group of order five, isomorphic to the first by `log_2`.
- `{3,5,17,257,65537}` has product `2^{32} - 1 = 4294967295` (`fermat-normal.json`), and its `2^5 - 1
  = 31` nonempty subset products are the odd constructible polygon orders below `2^{32}`.

The second floor also counts **square roots**: `[Q(zeta_{F_i}) : Q] = F_i - 1 = 2^{2^i}`, so the
regular `F_i`-gon is built from a tower of exactly `2^i` nested quadratic extensions — one square root
for the triangle, two for the pentagon, four for Gauss's 17-gon, eight for the 257-gon, sixteen for the
65537-gon. In the language of this session's root-surface packets, each quadratic step is a sign
torsor (a Möbius band over a loop around its branch locus), so the heptadecagon is a four-storey Möbius
tower. The owner's list `3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,…` is exactly the constructible
polygon orders `n = 2^a · (distinct Fermat primes)` (OEIS A003401), recomputed in the receipt.

## 2. Products of Fermat primes are Sierpiński rows (EXACT)

For `S ⊆ {0,…,4}` with bit-sum `s = sum_{i in S} 2^i in {0,…,31}`,

    prod_{i in S} (2^{2^i} + 1) = sum_{T ⊆ S} 2^{sum_{i in T} 2^i} = sum_{k : C(s,k) odd} 2^k,

by expanding the product and Lucas' theorem (`C(s,k)` is odd iff the bits of `k` are a subset of the
bits of `s`). So the 32 odd constructible orders below `2^{32}` are rows `0..31` of Pascal's triangle
modulo 2 read as binary numbers, and the five Fermat primes are the rows `s = 1, 2, 4, 8, 16`
(binary `11, 101, 10001, 1 0000 0001, 1 0000 0000 0000 0001`). The receipt checks the identity on all
32 subsets. Completeness of the five is the statement that no 33rd row type exists among odd
constructible orders. (The repository's Sierpiński lane, THM-4754, is 3-adic; this is the 2-adic
Sierpiński, a different carrier with the same name — recorded as COINCIDENCE of names, STRUCTURE
within each.)

## 3. The Fermat Paley graphs (EXACT)

Fermat primes are `≡ 1 (mod 4)`, so the Legendre relation on `Z/F_i` is symmetric: a Paley **graph**,
not a tournament (contrast the primes `≡ 3 (mod 4)` of HYP-9375). `P(F_i)` is self-complementary,
strongly regular of degree `(F_i - 1)/2 = 2^{2^i - 1}`, with spectrum `{(F_i-1)/2, (-1 ± sqrt F_i)/2}`
(Gauss sums). Computed: `P(5) = C_5` with eigenvalues `(-phi, 1/phi, 2)` — the golden ratio of the
Antihydra note returns as the pentagon's spectrum — and `P(17)` with eigenvalues `(-2.5616, 1.5616, 8)`
and clique number `3`: the classical `R(4,4) = 18` witness (self-complementary, so also no independent
4-set). `P(257)` and `P(65537)` continue the family. Completeness of the five says this
Fermat–Paley family has exactly five members.

## 4. Relations to the repository's other fives

Filled from the bounded survey (see the appended section).

### 4a. The ladder has a tournament floor (STRUCTURE, inherited)

THM-4748 (KNOWN, Moon 1972) gives `H(NT_{k+2}) = 2^k + 1` for the near-transitive tournaments, and
THM-4749 statement 4 records `F_i = H(NT_{2^i + 2})`: the five Fermat primes are the Hamiltonian-path
counts of `NT_3, NT_4, NT_6, NT_{10}, NT_{18}`, with `k = 2^i in {1,2,4,8,16}` the second floor. THM-4749
also has `2^{2^j} - 1 = F_0 ⋯ F_{j-1}`, the product identity of §1. THM-871 and THM-4644: Fermat primes are
exactly the primes at which every circulant tournament is rigid (the multiplier group is a 2-group), and,
CITED Akin 2018 there, the multipliers act freely on circulants of order `2n+1` iff `2n+1` is a squarefree
product of Fermat primes — i.e. iff `2n+1` is an odd constructible polygon order. So the owner's list is,
on the tournament side, the list of odd orders with free multiplier action. (That reading is not stated in
the repository; it is an inference from THM-4644's citation.)

### 4b. Other exact appearances (inherited, typed)

| object | file | type |
|---|---|---|
| von Staudt–Clausen: `denom(B_{2^r}) = 2·prod_{2^i <= r} F_i`; completeness of the five ⟺ the denominators stabilise at `2(2^{32}-1)` for `r >= 16` | `forest_20260926_bernoulli.md` | STRUCTURE (a third reformulation of completeness) |
| for `p = 2^k + 1` the Collatz parity graph modulo `2^k` is the `y -> ±sqrt y` graph on `F_p^x` | `collatz_sqdbl_20260925_squares_doubles_foundry.md` Lemma 1(b) | STRUCTURE: the square-root tower of §1 seen inside the parity graph |
| Gauss period at 17: `eta^2 + eta - 4 = 0` | `prime_shells_20260921_shells.md` | STRUCTURE (the first storey of the 17-gon tower) |
| rows `0..31` of Pascal mod 2 in binary are the 32 Fermat products; row 32 is the composite `F_5` | `gilbreath_fermat_platonic_20260926.md` §1–2 (audited) | STRUCTURE — §2 above recovers this; new here only the link to THM-4756's mod-2 forward-step row `(1+x)^{n-1}`, whose value at `x = 2` is the Fermat product with support the bits of `n - 1` |
| doubling-cycle counts `1,1,2,16,2048` and rotational class counts `1,1,16,2^{120},2^{32752}` at the Fermat primes | `fermat_wieferich_order_20261014e.md`, `number_tournament_bridges_20261014.md` | STRUCTURE |
| `P_5, P_17` are the `R(3,3), R(4,4)` extremal graphs; "that they are Fermat primes" | `number_tournament_bridges_20261014.md` §2 | typed NUMEROLOGY there; §3 above adds only the spectrum |

### 4c. The repository's other fives, and whether they relate

- The five integer Terras cycles `{0}, {-1}, {1,2}, {-5,-7,-10}, {-17,…}` (THM-4484, THM-4591): four are
  "free" by Gersonides' `|2^p - 3^a| = 1`, one sporadic. No map to the Fermat five; the only shared
  mechanism is "a power of two plus or minus one", which is COINCIDENCE at the level of the sets.
- The five groups of order 27 (this session, THM-4764): a classification by `p^3`; unrelated.
- The five Niemeier types with prime Coxeter number (`A1^24, A2^12, A4^6, A6^4, A12^2`, THM-4758): the
  primes `2,3,5,7,13` with `h = p`; unrelated to `F_i` (COINCIDENCE of cardinality).
- Platonic solids: the repository already ruled "not the same five" as the Fermat primes
  (`collatz_doubling_tower_brocard_20260927.md` §4 types six fives as having different sources); the
  only common object is `{3,4,5}`.
- `QR_11 = {1,3,4,5,9}` (THM-4592), Ellison's five shapes (THM-4591), the Lehmer five, the `-5` cycles:
  each explained by its own structure; none by the ladder.

Verdict: the three owner fives are one object (a tower with a tournament floor); the repository's other
fives are COINCIDENCES of cardinality, each with its own explanation, and the one earlier attempt to
identify two of them (Platonic/Fermat) was already rejected.
