# Definability lenses on the Hamiltonian-path spectrum: the "Z is Diophantine in Q" architecture, the Carmichael and Korselt themes, and what they say about the {7, 21} programme

**Session:** opus-2026-11-05-S81. Worktree `C:/tmp/wt-s81`. Canon: THM-4884. Hypothesis: HYP-9477. Pre-registration:
`04-computation/experiments/definability_lenses_prereg_20261105.md` (commit before the runs). Independent audit:
section 9.

**Owner prompt (verbatim).** "[attached: diophantine_zq.pdf] merge in prior themes around carmicheal's conjecture and
asking is every ring of integers of an algebraic number field Diophantine in its field of fractions? prioritize not
verification of outside results i mention, but instead reasoning about how their ideas can meaningfully extend our
work"

**Scope.** The attached paper ("Z is Diophantine in Q", dated 2026-11-05, building on an OpenAI preprint on Hilbert's
tenth problem over Q) is used only for its architecture; nothing in it is verified or relied on, as the owner asked.
The repository's Carmichael themes are two: the Korselt-type congruences of strong-tournament counts (THM-4746,
"Camion pseudoprimes") with the supplied Carmichael-counting paper behind `carmichael_affine_transfer_20261019`, and,
by analogy with Carmichael's totient conjecture, the multiplicity of a count. The concurrent Session 200 (active,
same directive) built a Diophantine compiler for physical source words (quadratic systems of size O(log n) whose
characteristic-zero solutions are integral, so C1/C2 membership is in NP) and height-bounded source lifting; this
note works the other lenses and cites that lane where it meets them.

## 0. Summary

- **Architecture map (reasoning, section 1).** The paper defines Z inside Q with four ingredients: integer indices
  modelled as multiples nP_0 of a point on a rank-one elliptic group (G), a guarded multiplication template (M), a
  height estimate (H), and local conditions at primes (pole predicate PC, Heegner/Kummer descent), assembled into one
  polynomial. The spectrum programme already contains the same four ingredients in a different arrangement: the
  one-apex orbit is a Z-indexed family (word value v, count g(v+1)); ordered joins are an exact multiplication
  template (H(T => T') = H(T) H(T')); Moon's floor and the sink-run bound are the height estimates (Northcott: finitely
  many tournaments below a given count); residues modulo every modulus are the local conditions (Session 193's
  modular synchronization, Session 198's modular seeds). What the spectrum lacks is the group law: no operation adds
  counts, and THM-4882 (4) shows no repeated block lowers a count. Session 198's signed theorem is exactly a
  "Diophantine definition of Z by tournament statistics" (every integer is aU - bD on a one-apex source), obtained
  with group structure (binary words), height control (the clock) and local conditions (residues); the positive
  spectrum problem is the same statement without a group, which is why it is a numeration problem (S79, S80).
- **Korselt lens (THM-4884 (1), PROVED, elementary).** For the spectrum polynomial P_n(x) = sum over labelled
  tournaments of x^H and every prime p | n, P_n = F_{n,p} (mod p), where F_{n,p} sums over tournaments fixed by the
  rotation by n/p; for prime p the fixed tournaments are the circulants, whose counts are multiples of p (the
  rotation acts freely on Hamiltonian paths), so **the number of labelled tournaments of prime order p with a count
  not divisible by p is itself divisible by p**. Checked exactly for n <= 7 (2,097,152 labelled 7-tournaments: the
  only coefficients not divisible by 7 are at h = 175 and 189, the circulant counts).
- **Carmichael-multiplicity lens (THM-4884 (2), FINITE-EXACT from THM-4744's inventories).** Let m(h) be the number
  of strong classes, over all orders, with exactly h Hamiltonian paths. Below 625, m(h) = 0 exactly for h = 7, 21,
  m(h) = 1 exactly for h = 3, 5, 11, 13 (each realised by a single self-converse class, as the converse pairing
  P_d = SC_d mod 2 of THM-4744 forces), m(h) = 2 exactly for h = 9, 17, 19, 35, 39, 47, 63, m(h) = 3 only for h = 25,
  and m(h) >= 6 for every attained h in (63, 625) (the minimum 6 at h = 67; m(h) >= 73 on [375, 623]). The
  multiplicity function omits the values 9, 13, 21, 22, 34, 37 below 625 while taking every other value up to 38:
  the analogue of Ford's theorem (every multiplicity >= 2 occurs for Euler's function) fails for the strong
  spectrum, and the analogue of Carmichael's conjecture (no multiplicity 1) fails exactly at the four smallest primes
  other than 7. Pre-registered Q1: contained in [3, 100], at most 12 elements, contains 3, 5, 11, 13: CONFIRMED;
  "contains 41": REFUTED (m(41) = 6); "every multiplicity 1..10 occurs": REFUTED (9 is omitted).
- **Diophantine complexity lens (THM-4884 (3), FINITE-EXACT).** The paper measures a definition by its number of
  variables; the natural measure here is the Hankel rank of the N-rational series whose range is to be cofinite in
  the odd numbers. Sub-languages of the cutwidth-2 grammar G(0,2) (THM-4878, rank 17): forbidding simultaneous double
  events gives rank 14, additionally forcing first-in-first-out closing (crossing chains) keeps rank 14, forcing
  last-in-first-out closing (nested chains) gives rank 10. None of the three is complete to 2047 (pre-registered Q2
  REFUTED for the first two): without double events the values 43, 47, 63, 79, 103, 113, 235, 709 are lost, so the
  shared-endpoint configurations (a vertex beaten by two later vertices, or beating two earlier ones) are
  indispensable pieces of any proof through cutwidth 2. The one-apex series (rank 4) attains half of the odd numbers;
  the minimal rank of a complete candidate lies between 5 and 17, and the two natural rank-14 sub-languages are
  excluded.
- **What transfers to a proof (section 7).** The paper's assembly ("six statements, then a lemma") is a template for
  a proof of Chan-Pak through a grammar: (i) a Z-indexed family with explicit counts (the one-apex orbit, KNOWN),
  (ii) a multiplication template (ordered joins, KNOWN), (iii) a height estimate (counting lemma and floors, KNOWN),
  (iv) local coverage (modular synchronization, KNOWN), and the missing (v): an additive or numeration step on
  realizations that is total above a bound. The repository's candidates for (v) are the C2 digit steps 3h + 2, 3h + 4
  (klein, Session 194: not total with the parent retained as a module), the one-apex Z-step h -> h - E + U (total but
  two-coordinate), and the cutwidth-2 arcs as digits (S79). The Korselt and Carmichael lenses do not supply (v); they
  supply congruence and multiplicity invariants that any proposed (v) must respect.

## 1. The architecture of the attached paper, read as a template

(Section 1 of the paper, taken as stated, not verified.) The polynomial P in 292,890,231,905 variables is assembled
from six statements: **G** canonical point codes form the infinite cyclic group 2E(Q(sqrt 2)) = Z P_0, with a
membership polynomial and raw addition describing the group law; **M** a guarded multiplication template whose
rational witnesses are exactly the index triples (a, b, ab); **H** a six-contact height theorem; **PC** an explicit
pole polynomial, soluble for every input; **C** a five-input natural polynomial with a power template; **B** the
literal template specification. Lemma 1.2 (Assembly) combines them. The consequences claimed are the refutation of
Mazur's conjecture and the many-one equivalence of HTP(Q) with the halting set.

The same skeleton describes what the repository has built for the spectrum, with one slot empty:

| paper | spectrum programme | status |
|---|---|---|
| G: integer indices as multiples of P_0, group law = addition | the one-apex orbit: chain words of binary value v, profile (g(v+1), 1 - last bit, v + 1, 1); letters O/Z = append a bit (THM-4878 G(1,0); Session 197/198) | KNOWN; a free monoid, not a group |
| value relation extracting n from nP_0 | the count g(v + 1) and the signed functional aU - bD (Session 198, THM-4882) | the signed functional is onto Z; the count is onto half of the odd numbers |
| M: guarded multiplication (a, b, ab) | ordered join: H(T => T') = H(T) H(T'); the spectrum is a multiplicative monoid generated by its strong values (THM-4094, session193 prime_induction) | exact, KNOWN |
| H: height estimate | Moon's floor f(n) (THM-4744), the sink-run and counting lemmas (THM-4878), Session 200's height-bounded lifting | KNOWN |
| PC: pole predicate (divisibility of denominators by primes along indices) | divisibility of counts: p | H(circulant of order p); THM-4747's Fermat-Wilson congruences; the Korselt congruence below | elementary |
| local conditions (Heegner descent, Kummer at 2) | residues: every admissible residue modulo every m is attained by one-apex words (Session 193), modular seeds (Session 198), the modulo-4 confinements of THM-4882 | KNOWN |
| assembly into one polynomial | Session 200's quadratic systems of size O(log n) for C1/C2 membership | KNOWN (concurrent) |
| the group law itself | **missing**: no additive operation on counts; THM-4882 (4) forbids lowering by repetition | the obstruction |

Reading: the paper proves a universality statement (every integer is a point code) by a group structure plus size
and local control. The spectrum's universality statement (every odd number except 7, 21 is a count) has size and
local control but no group; its surrogate is a numeration system, and the three candidate digit systems in the
repository are listed in section 7.

## 2. The ring-of-integers question as the physical/integral distinction

The owner's second question - is every ring of integers O_K Diophantine in its fraction field K - is the
number-field form of "what distinguishes integral solutions from fractional ones". Session 200 identified the
spectrum analogue precisely: a profile (D, E, U, d) can be integral and satisfy every congruence test without being
physical (their hostile (1, 0, 6, 1), which would give the C2 count 67 formally), and their compiler forces
physicality by carrying the construction word and its height. In the grammar language of THM-4878 the same
distinction is the difference between the reachable cone of cut vectors (physical) and the lattice they span
(integral): the Hankel rank measures the lattice, the reachable set is the "ring of integers". The seam ranks 3, 13,
62, 310 of THM-4878 are the dimensions of the lattices across one, two, three, four reversed arcs; the physical
reachable sets inside them are what a Diophantine definition of the spectrum would have to carve out. No theorem is
claimed here; the observation fixes the vocabulary: *integral* = in the span, *physical* = reachable, and the
"definability of O_K in K" question becomes "is the reachable set of a grammar existentially definable inside its
lattice by equations of bounded size" - which Session 200 answered positively for C1/C2 words (O(log n) variables).

## 3. The Korselt lens: orbit-counting congruences of the spectrum polynomial (PROVED)

Let P_n(x) = sum over the 2^C(n,2) labelled tournaments T on Z_n of x^H(T). The cyclic group Z_n acts on labelled
tournaments by i -> i + 1; H is invariant. For a prime p | n, the subgroup generated by the rotation by n/p has order
p, so its orbits have size 1 or p; the orbits of size 1 are the tournaments invariant under i -> i + n/p. Hence

    P_n(x) = F_{n,p}(x) (mod p) coefficientwise,   F_{n,p}(x) = sum over the invariant tournaments of x^H.

For n = p prime the invariant tournaments are the circulants (arc i -> j iff j - i in S, S a half-system of Z_p^*),
and the rotation acts freely on the Hamiltonian paths of a circulant (a path's first vertex moves), so every
circulant count is a multiple of p. Corollary: **for prime p and every h with p not dividing h, the number N(p, h) of
labelled p-tournaments with exactly h Hamiltonian paths is divisible by p.** This is the spectrum's version of the
Korselt mechanism of THM-4746 (there: the strong count s(p) = (2/p) mod p by the same action). It gives attainment
information only for circulant counts and is recorded as an invariant, not as a tool for Chan-Pak.

Exact check (`spectrum_polynomial_congruence_20261105.cpp`, all 2^C(n,2) labelled tournaments, n <= 7): the congruence
holds for every (n, p) with p | n <= 7; at n = 5 the only coefficient not divisible by 5 is h = 15 (64 = 4 mod 5: the
four circulants), at n = 7 the only coefficients not divisible by 7 are h = 175 (720 = 6 mod 7) and h = 189 (240 = 2
mod 7), the eight circulants on Z_7 having counts 175 (six of them) and 189 (two: the Paley tournament and its
converse); at n = 6 and p = 3 the thirty-two rotation-by-2-invariant tournaments give the residue 2 at h = 9 only.
The full coefficient tables for n <= 7 are in the program's output.

## 4. The Carmichael-multiplicity lens (FINITE-EXACT)

Carmichael's conjecture says Euler's function takes no value exactly once; Ford proved it takes every multiplicity
>= 2. For the spectrum define m(h) = the number of strong tournament classes (all orders) with H = h; it is finite
because H >= f(n). THM-4744's complete inventories (`strong_ears_375.out`, `strong_ears_625.out`: P_d(1) for odd
d < 625) give, read off here:

| m(h) | odd h < 625 |
|---|---|
| 0 | 7, 21 |
| 1 | 3, 5, 11, 13 |
| 2 | 9, 17, 19, 35, 39, 47, 63 |
| 3 | 25 |
| 4, 5 | 15, 23, 31, 43, 55 (4); 27, 61 (5) |
| >= 6 | every other attained h < 625; the minimum above 63 is m(67) = 6; m(h) >= 73 for all odd h in [375, 623] |

The full table for odd h <= 101 (h:m): 3:1 5:1 7:0 9:2 11:1 13:1 15:4 17:2 19:2 21:0 23:4 25:3 27:5 29:6 31:4 33:8
35:2 37:7 39:2 41:6 43:4 45:12 47:2 49:8 51:12 53:10 55:4 57:7 59:12 61:5 63:2 65:14 67:6 69:18 71:11 73:19 75:12 77:16
79:11 81:26 83:11 85:10 87:11 89:23 91:29 93:15 95:17 97:23 99:24 101:18.

Observations. (a) The unique-class values are the four smallest primes other than 7; each is realised by a
self-converse class (C_3, the strong 4-tournament, and the order-5 realisations of 11 and 13), as THM-4744's converse
pairing P_d = SC_d (mod 2) forces for any odd multiplicity. (b) The multiplicities 9, 13, 21, 22, 34, 37 do not occur
below 625 although every value up to 38 other than these does; since m(h) >= 73 on [375, 623] and the counts grow
with h, these six values are plausibly never multiplicities of strong counts (HYP-9477 (A)). (c) The largest h with
m(h) <= 3 is 63, the largest with m(h) <= 5 is 63 as well (m(63) = 2, m(61) = 5), and m(h) >= 6 beyond; the
"Carmichael exceptional set" {3, 5, 11, 13} is complete below 625 and, by (b), presumably in total. Pre-registered
Q1: CONFIRMED in its containment, size and membership clauses; REFUTED in the clause "contains 41" (m(41) = 6) and
in the clause "every multiplicity 1..10 occurs".

The Korselt theme of THM-4746 and the HHKK Carmichael-counting method (subset products in (Z/L)^* with character
control) are counting methods for congruence-defined composites; their direct spectrum analogue is the enumeration
of ordered joins with a prescribed product, which is not where the conjecture's difficulty lies (the prime atoms).

## 5. The Diophantine-complexity lens: minimal rank of a complete sub-grammar (FINITE-EXACT)

The paper's size statement (variables, degree, coefficients) suggests measuring a candidate definition of the
spectrum by the dimension of its linear representation. Sub-languages of G(0,2) (THM-4878; the engine of S79 with a
letter filter, `spine_grammar_sublanguages_20261105.cpp`; ranks with the empty tournament counting one path):

| sub-language | letters | cut states | Hankel rank | seam c = 2 | attained odd <= 2047 (products added) |
|---|---|---|---|---|---|
| G(0,2) | all | 83 | 17 | 13 | 1022 of 1024 (misses 7, 21; S79) |
| mode 1: no double events | b <= 1 and at most one slot terminates per vertex | 64 | 14 | 10 | 1014 of 1024: misses 7, 21, **43, 47, 63, 79, 103, 113, 235, 709** |
| mode 2: crossing chains | mode 1 and only the older slot may terminate | 64 | 14 | 10 | 998 of 1024: 24 misses besides 7, 21 (37, 43, 47, 63, 67, 79, 103, 113, 127, 139, 163, 173, 235, 431, 491, 519, 559, 563, 709, 751, 797, 1373, 1427, 1499) |
| mode 3: nested chains | mode 1 and only the newer slot may terminate | 58 | 10 | 6 | 440 of 1024: 582 misses besides 7, 21 (first 11, 13, 19, 23, 35) |
| G(1,0) one-apex | | 15 | 4 | | 763 of 1024 (THM-4878) |

Pre-registered Q2: mode 1 and mode 2 complete, mode 3 misses at least 10; ranks mode 2 < mode 1 <= 17 and mode 3 <
mode 1. Outcome: "mode 1 complete" REFUTED, "mode 2 complete" REFUTED, "mode 3 misses >= 10" CONFIRMED, "mode 2 <
mode 1" REFUTED (both ranks are 14: first-in-first-out closing removes no dimension but loses sixteen more values),
"mode 3 < mode 1" CONFIRMED. The structural reading is the useful part: **the double events are indispensable.** A
vertex that opens two slots (two later vertices beat it: a shared head) or terminates two (it beats two earlier
vertices: a shared tail) is needed already for 43, 47, 63, 79, 103, 113 below 256. S79's exact minimal-arc
witnesses (to 511) confirm this for the seven lost values they cover, each having a shared endpoint: 43 = chain 0..6
with 2 -> 0, 4 -> 2, 6 -> 2, 6 -> 4 (head 2 and tail 6 shared); 47 = 0..6 with 6 -> 0, 6 -> 1 (shared tail; THM-4748's
shared-endpoint case); 63 = 0..7 with 3 -> 0, 6 -> 3, 7 -> 3 (shared head 3); 79 = 0..7 with 7 -> 0, 7 -> 2 (shared
tail); 103 = 0..8 with 3 -> 0, 6 -> 3, 8 -> 3 (shared head); 113 = 0..8 with 2 -> 0, 6 -> 2, 6 -> 3, 8 -> 6 (shared tail
6); 235 = 0..9 with 6 -> 0, 6 -> 4, 8 -> 6, 9 -> 6 (vertex 6 shared as tail and as head). So a proof of completeness
through cutwidth 2 must use the shared-endpoint
pairs, i.e. the vertices of backward out-degree or in-degree two, and the crossing chains alone (rank 14, 998 of
1024) are not a complete candidate; the full grammar (rank 17) remains the smallest complete one found. The minimal
rank of a complete candidate therefore lies between 5 and 17, with 14 excluded for these two natural sub-languages.

## 6. Three-smooth ladders and the two-radix sums

The repository's Collatz thread proved that balanced cycles are 3-smooth ladders (THM-4716) and that certain ladders
are complete sets of 3-smooth numbers in a window (THM-4719); S79 observed that every complete grade of the spine
grammars carries both radices (powers of 2 from single runs, powers of 3 from two interleaved runs), and that the
cutwidth-2 counts with few arcs are sums of 3-smooth terms (THM-4748's crossing formula 1 + 2^(p+r+1) + 2^(q+r+1) +
2^(p+q+1) 3^r). The Diophantine lens makes the conjecture of HYP-9470 (C) a representation statement: every odd
number above 21 is a structured sum of 3-smooth numbers of the kind produced by cutwidth-2 words. The classical
result of this type (every positive integer is a sum of 3-smooth numbers none dividing another, Blecksmith, McCallum
and Selfridge 1998, cited from memory and not re-verified here) is proved by a greedy descent; the window
completeness of THM-4719 is the ladder form of the same descent. The proposal recorded in HYP-9477 (C): compute the
exact closed form of the k-arc crossing chain (mode 2 above: rank 14, so a 14-term transfer) as a sum over 3-smooth
terms indexed by the gaps, and test whether its value set is described by a greedy rule on 3-smooth representations;
if so, the completeness of crossing chains would reduce to a representation theorem in the style of Blecksmith,
McCallum and Selfridge.

## 7. What a proof needs, in the paper's assembly form

Transcribing Lemma 1.2 (Assembly) to the spectrum: Chan-Pak follows from (i) G: a Z-indexed family with explicit
counts (the one-apex orbit; KNOWN), (ii) M: multiplication (ordered joins; KNOWN; reduces the conjecture to prime
atoms, 7p, 49, 63, 343), (iii) H: height (floors and the counting lemma; KNOWN), (iv) local coverage (KNOWN), and
(v) a total numeration step on realizations above a bound. The repository's three candidate steps:
- the C2 digit steps h -> 3h + 2, 3h + 4 (klein THM-4862/4865, Session 194): a ternary numeration with digits {2, 4}
  if both children always existed, which Session 194's module barrier shows they do not when the parent is kept as a
  module; the two-branch tree is complete only along explicit lineages;
- the one-apex Z-step (D, E, U) -> (D - E + U, 1, 2U - 1): total, binary, but two-coordinate, which is why it
  programs only signed responses exactly (Session 198, THM-4882);
- the cutwidth-2 arcs as digits (S79): complete to 8191 with at most five arcs, with no proof mechanism yet; the
  three-smooth representation route of section 6 is the candidate mechanism.
The Korselt and Carmichael lenses supply invariants that (v) must respect: the congruence of section 3 constrains
how many labelled realizations a value has at prime orders, and the multiplicity census of section 4 says that above
63 every strong value below 625 has at least six strong realizations, so a numeration step may be chosen among several
realizations of its input - a freedom that Session 194's module barrier shows is necessary.

## 8. Computations and their status

- `spectrum_polynomial_congruence_20261105.cpp`: P_n and F_{n,p} for n <= 7 (FINITE-EXACT; the congruence is
  PROVED in general).
- the multiplicity census: read from THM-4744's committed outputs (no new enumeration).
- `spine_grammar_sublanguages_20261105.cpp` + `spine_grammar_ranks_mod_20261104.py --empty1`: automata, ranks and
  censuses of the three sub-languages (ranks FINITE-EXACT; censuses to 2047 FINITE-EXACT when finished).
- Pre-registration outcomes: Q1 partly confirmed (41 and the multiplicity-10 clause refuted), Q2 rank clause
  "mode 2 < mode 1" refuted, Q3 a check.

## 9. Audit

PENDING.
