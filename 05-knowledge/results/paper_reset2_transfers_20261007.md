# A source-preserving phase compiler for Mersenne reset-two branches

2026-10-07. **PROVED:** the elementary phase compiler, exact finite guard
ledger, and applications of the explicitly named inherited joins.
**FINITE-EXACT:** the declared controls. The four supplied paper headlines
are **TRUSTED PREMISES**, not reverified here. Universal reset-two coverage
and Collatz remain **OPEN**. A smaller dependency is not a ROOT certificate.

[Program](../../04-computation/experiments/paper_reset2_transfers_20261007.py)
and [exact output](paper_reset2_transfers_20261007.out).

## 1. Inheritance and the four paper interfaces

The closest operations are the eight authenticated half-child joins in
[debt-word search, sections2--4](debt_word_join_search_20261004.md), the
uncapped initial-run reset, and the single Mersenne sporadic class in
[terminal lifts, sections3--4](collatz_terminal_lifts_20261007.md).
The canonical hostile is the supplied source2^1459-1, which missed those
last two rules. The corrected near miss is extrapolating a finite word cap
or an unrefined dyadic sign obstruction to every mixed or variable-depth
certificate. The least-used sidecar is the original exponent K, retained
through the nonlinear reset state instead of replacing the source by that state.

Board: **original source / exponent / exact tail / ordered carry / shared
phase / paid child / unresolved branch**. The present contribution is a
compiler for these interfaces and a precisely named restricted-family
coverage increment, not discovery of the eight old identities.

| Supplied source | Read interface (one-based PDF pages) | Useful construction and transfer boundary |
|---|---|---|
| [Maximum number of MUBs in dimension six](C:/Users/Eliott/Downloads/The-maximum-number-of-mutually-unbiased-bases-in-dimension-six-September-24-2026.pdf),76pp | §1.2 pp4--5; shared-coordinate tests §2.3 pp9--10; necessary-condition certificates and unresolved recursion §3.7 pp20--21; bootstrap/atlas §4.6--4.7 pp31--33 | A heuristic may propose a representative or certificate, but acceptance preserves every exact feasible point. Shared parameters cannot be chosen separately. A capped live cell remains unresolved. Figure1 p5 was inspected visually. Our analogue is the exact original-K guard and its explicit complement ledger; the MUB execution is not a Collatz cover. |
| [Threshold parallel repetition for finite-dimensional entangled games](C:/Users/Eliott/Downloads/paper-16.pdf),22pp | §2 pp5--8 | Its sampling lemma works for correlated indicators, but the amplification step obtains independent tensor copies and invokes a separate conditioning bound. Repeating the same source test supplies neither independence nor that bound. We retain exact intersecting residue events. |
| [Weak pinned planar distance theorem](C:/Users/Eliott/Downloads/paper-15.pdf),27pp | outline pp3--4; nested partitions and multiplicity pp9--12; distinct-pair corrections pp15--16; conditional tree transport p21 | A partition, its sampling law, and cell occurrences are separate objects. Giant replacement may duplicate a cell in a multiset; distinct-pair identities remove singleton tails exactly. Figure1 p11 was inspected. For guard coverage we take set union, with authenticated receipts retained separately: duplicate receipts cannot count as additional source mass. No number-field product-formula hypothesis is asserted for our ledger. |
| [Logarithmic Brunn--Minkowski conjecture](C:/Users/Eliott/Downloads/paper-17.pdf),20pp | soft slabs/variance pp3--6; density kernel pp10--12; constant-frame tensor identity pp12--15 | Finite outer constraints preserve an exact nested intersection; measure continuity is a separate conclusion. Affine obstructions are explicitly identified before density is used, and a pointwise constant frame is not differentiated as a moving frame. Figure1 p6 was inspected. Shrinking unresolved measure in our parameter space does not imply an empty integer remainder. No convex-volume inequality is transferred to Collatz. |

The papers supply methods to retain an obligation, not arithmetic premises
for resolving it. Everything used about the following integer maps is
proved directly below or explicitly linked to its earlier certificate.

## 2. The exact reset chart and its invertible finite precision

Write M_K=2^K-1, with K odd and K>=3, and let U be accelerated3n+1 on
positive odd integers. The first K-1 valuations are1, followed by2. Indeed

    U^j(M_K)=3^j 2^(K-j)-1,     0<=j<=K-1,
    z_K=U^K(M_K)=(3^K-1)/2.                            (1)

These initial states grow; they do not reach ROOT. In particular z_K=1mod4.
The next valuation is exactly

    v2(3z_K+1)=1+v2(K+1).                             (2)

For the last equality use v2(3^h-1)=2+v2(h) for even h,
with h=K+1, and subtract one for the denominator2. The elementary
valuation formula follows by repeated squaring from9=1+8: for every
nonzero integer d, v2(9^d-1)=3+v2(d), with the unit denominator understood
when d<0.

There is a stronger coordinate statement. Put K=2j+1 and z_K=4t+1. Then

    t=F(j)=3(9^j-1)/8,
    v2(F(j)-F(l))=v2(j-l).                            (3)

Consequently F induces a permutation modulo2^b for every b>=0, compatible
under reduction. This defines a bijective isometry of Z_2. Its use here is
entirely concrete: it converts a finite native guard on z into a finite
native guard on the *same* K. The archimedean sourceM_K, the reset statez_K,
and the parameterj are distinct marked coordinates. Equation(3) is not a
conjugacy of their Collatz dynamics.

The actual next-letter cells in(2) are

    K=2^(a-1)-1 mod2^a,   a>=2.                       (4)

They have relative natural density2^(1-a) among odd K. A finite split
retains its unsplit tail K=-1 modulo the next power of2. The intersection
of those particular tails is the2-adic value-1, not a positive exponent.
This proves that parsing the *next valuation* terminates on each supplied
positive K. It does not prove that every resulting cell has a paid join.

## 3. Compiling any fixed tail without expanding M_K

Let u be a finite positive valuation word, with carrier

    F_u(x)=(Px+B)/Q,   P=3^length(u), Q=2^A, A=sum(u).

For a nonempty u, substituting(1) shows that its final state is odd exactly
when

    3^K = c mod2^(A+2),
    c=P^(-1)(P-2B+2Q) mod2^(A+2).                    (5)

The inverse is modular. The exact-word cylinder lemma proves that this
final oddness condition gives every positive prescribed intermediate
valuation, if ROOT padding is permitted in the formal identity. It is
proved by reducing the first numerator modulo2^(a1+1), dividing by2^a1,
and repeating. We separately remove ROOT padding below.

**Phase theorem.** Equation(5), with K odd, is solvable iff c=3mod8.
When solvable, it gives one odd residue K0 modulo2^A. Starting with K=1
mod2, determine successive bits: replacing K by K+2^j changes3^K by a
number of exact valuation j+2 for j>=1. Exactly one choice therefore
matches the next bit. There are A-1 lift decisions. This is an arithmetic
operation count, not a constant bit-time claim.

The criterion also follows from(3): the powers3^(2j+1) cover exactly the
residues3mod8 at every compatible precision. The empty tail simply accepts
every odd K>=3.

For an actual first-hit tail, every proper prefix i must exceed1:

    P_i 3^K > 2Q_i+P_i-2B_i.                         (6)

Each is a strict lower bound on K. `compile_phase` retains the least point
of the residue class satisfying all these bounds and K>=3. Thus its output
is an **iff** for the stated actual first-hit tail. It is not merely a
necessary phase filter. A tail may end at1; the next-letter function rejects
that case rather than adding the root self-edge.

If the endpoint after u is greater than1, its next exponent is

    a=v2(3P 3^K+6B-3P+2Q)-(A+1).                    (7)

`last_valuation` evaluates the numerator modulo successively larger powers
of2. It never expands2^K or3^K. Its termination follows because the actual
positive numerator is nonzero. The guard is checked first, and an exact
power-of3 comparison rejects an endpoint equal to1. No uniform constant
precision or ROOT search is claimed.

## 4. Seven old portals become exact Mersenne classes

Use the eight rows(s,e) from the inherited debt-word search. Their full
source/child words at an arbitrary initial run r>=1 are

    source: 1^r,2^e,u,a,
    child:  1^(r-1),2e+s,v,a+delta,                  (8)

where u,v here omit their displayed variable final exponents. They join
the smaller child(n-1)/2, with all native guards and strict ROOT conventions
proved in that note. One can independently check the base carries:
P_source=P_child, Q_source=2Q_child,
2B_child-B_source=P_source. Prepending another1 to each side preserves
the identity because x↦(x-1)/2 commutes with x↦(3x+1)/2.

For M_K take r=K-1. The tail after its first reset2 is2^(e-1),u, where
2^(e-1) in this sentence means repeated valuation letters. Applying(5):

| (s,e) | Tail after first reset2, before variable last a | Exact K class |
|---|---|---|
| (1,1) | (1,6) | empty |
| (1,2) | (2,1,6,4) | 5589 mod8192 |
| (1,3) | (2,2,1,2,9) | 32945 mod65536 |
| (1,4) | (2,2,2,1,14) | 1035585 mod2097152 |
| (2,1) | (6) | 31 mod64 |
| (2,2) | (2,10) | 1289 mod4096 |
| (2,3) | (2,2,10) | 1889 mod16384 |
| (2,4) | (2,2,2,3,1,11) | 1522049 mod2097152 |

The first empty row is structural: z_K=1mod4 cannot have next valuation1.
Every other row's least positive representative already exceeds its strict
cut. All listed preterminal endpoints exceed1 (for example the uniform
bound3^31-1>2^22 dominates the largest tail cost21), so(7) exports a legal
last edge. The inherited child last exponent a+delta>=3 prevents a prior
child root visit; earlier root padding would force all later exponents2.

Each accepted M_K has the marked smaller dependencyM_(K-1). The latter
has even exponent, hence the inherited ordinary reset supplies the further
dependencyM_(K-2). This is a one- or two-rule strict decrease in K, not a
completed certificate without the resulting child's proof.

These seven classes are disjoint: the number of successive2 letters and
then the next letter separate the original rows. The class31mod64 is the
already named sporadic family. The other six are applications of older
identities to this previously restricted routing portfolio, not new identities
or claims of novelty against every stored Collatz rule.

## 5. The separately proved direct two-step exponent reduction

The parent's collision search supplies source reduced prefix

    (2,3,2,2,2,2,1,5,1,1,2), followed by actual c,

and the longer partner

    (2,2,2,1,2,2,2,1,1,3,1,1,1,c+2).

With two extra leading ones on the source, these are identical affine maps
at sources4x+3 and x. The script checks that carry identity for c1..7 and
the actual whole join at K1459. The carries omit the final exponent, and
both denominators acquire the same factor when c is increased, so the
identity at c=1 proves it for every c>=1. General validity is also direct: the source
native cylinder makes the final endpoint odd, the affine equality transfers
that integrality/oddness to the partner, and its final exponent>=3 prevents
earlier ROOT. Source preterminal ROOT is excluded by(6) including its full
preterminal tail.

The tail after the first reset2 is

    u=(3,2,2,2,2,1,5,1,1,2),   sum(u)=21.

Its phase is exactly **K=1459 mod2097152**. For every member, prependK-1
ones on the source andK-3 ones on the partner. Their common future gives
the smaller dependencyM_(K-2) directly. This class is disjoint from the
seven above: its first tail letter is3, while theirs start2 or6.

The source-specific advance is that2^1459-1 now has a checked reduction;
it is not declared ROOT-grounded solely from that reduction. This package
compiles and accounts for the parent's identity; it does not claim to have
discovered it or minimized its word length.

## 6. A finite exact ledger and its unresolved remainder

Normalize j=(K-1)/2. A phase K=r mod2^A is one dyadic cell
j=(r-1)/2 mod2^(A-1), with relative natural density2^(1-A) among odd K.
Our ledger starts from the single j-cell and subtracts a guard by splitting
only the cells containing it. At every split the two children are disjoint
and cover the parent. This proves exact union and complement by finite
induction, including nested or duplicated input guards. Original receipt
labels remain available even when counting ignores duplicates.

Here the ledger counts **dyadic cylinders and their asymptotic masses**.
For a generic compiled phase it does not encode the finite head removed by
`phase.least`; pointwise membership always also uses `accepts`. A generic
zero-cell remainder would therefore still leave those finite head obligations
to check. For the eight named nonempty rows below, the least accepted K is
already the positive residue representative, so there is no such additional
head within the domain odd K>=3. Their displayed complement is pointwise exact.

The seven inherited applicable portals have mass

    16849/524288;
    increment beyond old31mod64 =465/524288.

Including the separate D=2 row gives

    paid named union =33699/1048576,
    unassigned complement =1014877/1048576.             (9)

The complement is represented by86 disjoint dyadic cells, largest normalized
modulus2^20. No2^21-element census is needed. These are densities of exponent
parameters, **not** densities among all integer sources, nor probabilities
that one supplied source converges.

The same relative fractions hold inside the preceding named frozen-policy
residual K=1mod1458, K>1024. Writing K=1+1458j, each phase is a single
j-class modulo2^(A-1), because729 is odd. The exact least exponents and
periods are saved in the output. In particular the new(1,2) family becomes

    K=4470229+5971968t,   t>=0,

and the direct D=2 row becomes K=1459+1528823808t. The old sporadic row is
K=10207+46656t, exactly as inherited. Original exponent and integer source
identities are retained; no giant Mersenne integer is needed for membership.

Equation(9) concerns this **named dyadic union only**. The corrected incoming
three-adic inverse-cycle routing is outside this ledger. For example, the
incoming rule on n=4mod9 applies to M_K when K=5mod6. On the named residual
K=1mod1458, the dependenciesM_(K-2) above have exponent5mod6 and can therefore
enter that additional rule. This observation supports composition of domains;
it is not permission to extend the retracted all-refinements dyadic sign barrier.

### What the paper methods require at the remaining boundary

The MUB finite-recursion interface transfers literally: a processed node
must retain all its unproved children, and a complete finite cover needs
zero unresolved leaves **and no unpaid finite head exceptions**. A discarded
or timed-out child is not paid.
The exact phase compiler makes this condition decidable for any finite
portfolio; it does not make the current86 leaves empty.

Nor does residual measure tending to zero suffice. The nested cells
K=1459mod2^b have masses tending to zero and still contain the same positive
exponent1459 at every stage. This is a synthetic unresolved-ledger hostile,
not a claim that1459 remains uncovered by the new row. The special tail
intersection-1 in(4) has a proved arithmetic exclusion; an arbitrary nested
residual branch has no such automatic exclusion.

Repeated testing does not create parallel repetition. If a single parity
event has probability1/2 and is tested d times on the same parameter, its
intersection still has probability1/2, not2^-d. One must prove a conditional
bound on the actual remaining parameter law or introduce genuinely independent
samples with an appropriate target theorem. Neither establishes that every
individual integer is paid. The retained finite set union avoids this mistake.

A possible finite-base induction would require a source-preserving rule
on every remaining positive exponent, with each dependency exponent smaller,
plus checked base certificates. Current(9) does not meet that premise.

## 7. Reproduction and exact scope

Run the script with `python -B -X utf8` and with `python -B -O -X utf8`.
All checks use exact integers/Fractions and explicit exceptions. The finite
universe contains: every target modulo2^b for b3..12; all words of length0..3
over1..4 against all odd K3..129; finite normalized-chart permutations;
all eight inherited phase rows; four complete moderate-size Mersenne joins;
the actual D=2 join at1459; CRT lifts including parameter10^30; exact ledger
identities and its first4096 positive odd exponents; and malformed, ROOT,
duplicate, correlated-event and shrinking-singleton hostiles. It does not
simulate enormous members or discover ROOT routes.

The production input domains retain exact integer types (bool/float aliases
are rejected), ordered tails and canonical phase fields. A final valuation
is computed only after its source guard and non-ROOT endpoint are checked.
Normal/optimized output must agree. No accepted paper headline, global
Collatz assertion, independent source floor, or new terminal root is inferred
from these finite controls.
