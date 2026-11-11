# Endpoint constraints and infinite seed constructions

**PROVED IN NOTE, with FINITE-EXACT certificates where indicated.**
Same-agent mathematical review and an independent physical-graph counting
implementation; no claim of independent human review or global priority.
Chan–Pak and seeds for every odd coefficient remain **OPEN**.

Write H(B) for the number of directed Hamiltonian paths. Let Q_ab count
ordered pairs of nonempty disjoint paths covering B, the first ending at a
and the second starting at b. Put W_ab=Q_ab+Q_ba. Let U_b count ordered
two-path covers whose second path starts at b, allowing the first path to
be empty. For a nonempty selected set O define

    delta(B,O) = H(B) - sum_{unordered {a,b} subset O} W_ab,
    U(B,O) = sum_{b in O} U_b.

The inherited tail construction adds B -> TT_t -> r, with r -> O and
B\O -> r. For t>=1 its count is

    H(tail(B,O,t)) = delta(B,O) + 2^t U(B,O).                 (0)

This is the KNOWN tail law of
[THM-5126 / unit-offset rays, A](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md),
credited there to Session216. This packet studies its endpoint-preserving
substitutions. A unit seed has delta=1. If U=c*2^j with c odd, it realizes
every 1+c*2^s with s>=j+1.

## E1 The difference lattice also terminates

Let E_k(T) be all tournaments extending a fixed n-vertex T by k named
vertices, with every new incident pair free. Define

    D_k(T) = gcd {H(V)-H(T): V in E_k(T)},

using gcd{0}=0. Adding transitive sources realizes H(T), so this is also
the gcd of all differences between two extension values. Then

    C_k(T)=gcd(H(T),D_k(T)),
    D_(k+1)(T) divides D_k(T),
    D_k(T)=2 for n>=2 and k>=n-1.                             (1)

Here 2 divides 0 in the usual divisibility convention. The first assertion
is the elementary gcd identity; the second follows by padding each child
with a universal sink. All differences are even by Hamilton-path parity.

For the last assertion use the ordinary multilinear extension polynomial
F, not its overlapping positive-cylinder representation. Its nonconstant
coefficients have gcd equal to the gcd of F(x)-F(0): Boolean finite
differences prove both divisibilities. Since H(T) is one of its values,
this is D_k(T). At k=n-1, the alternating old-new path in
[P13 / EXTENSION_CONTENT.md](../modular_centers_collatz_pell_20261111/EXTENSION_CONTENT.md)
has a top monomial with coefficient +2 or -2. Its two orders have the same
sign because their length 2n-2 is even. Hence D_(n-1) divides 2, and parity
gives equality. Sink padding proves all later cases. The n=1 boundary is
D_0=D_1=0, D_k=2 for k>=2; do not apply (1) at n=1. QED.

For fixed k, coefficients of degree at most 2k suffice, so D_k is
computable from the same polynomial-size sample of small supports as C_k.
This strengthens the arithmetic conclusion of P13. It still says nothing
about whether a specified point of the affine lattice actually occurs.

For one added vertex with outgoing set S, endpoint counting gives

    F(S)=H(B)+sum_{i in S}(U_i-E_i)-sum_{{i,j} subset S} W_ij, (2)

where E_i is the number of paths ending at i. Thus D_1 is the gcd of the
linear coefficients U_i-E_i and all W_ij. This response lattice is already
proved in [THM-4118 / ear-response-lattice-and-stateful-unit-component-intervals](../../../01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md).
The new content here is (1)'s general-k terminal statement.

**FINITE-EXACT.** Ordinary bases of orders 1..7 give the following (C_1,D_1)
distribution. All 4,054 children of bases through order six are also
counted as literal tournaments by `endpoint_checks.py`.

| n | Distribution (C_1,D_1): number of ordinary bases |
|---:|---|
|1|(1,0):1|
|2|(1,2):1|
|3|(1,2):2|
|4|(1,2):4|
|5|(1,2):12|
|6|(1,2):55; (3,6):1|
|7|(1,2):452; (3,6):3; (3,78):1|

In this census C_1=1 always accompanies D_1=2; no example with C_1=1 and
D_1>2 is asserted. Even D_1=2 does not force solid support: THM-4118's
strong order-five code8 has H=9 and nonconstant-cut image
{15,17,19,23,25,27,29,33,37,41}, with holes inside its range. The small-code
convention in that theorem need not match this packet's canonical codes.

## E2 A finite and exact test for an infinite seed family

Fix (B,O,v), with n=|B|, v outside O, and O nonempty. Replace v by a
tournament M, preserving every external arc uniformly. Let P_j(M) count
ordered partitions of V(M) into j nonempty directed paths. There are
integers d_j and nonnegative integers u_j, depending only on (B,O,v), with

    delta(B[M at v],O) = sum_{j=1}^n d_j P_j(M),
    U(B[M at v],O)     = sum_{j=1}^n u_j P_j(M).              (3)

**Proof.** In a Hamiltonian path or a two-path cover, contract each maximal
run lying in M to the symbol v. Adjacent v symbols in the same path are
forbidden, since they were one run. Every other vertex occurs once. If
there are j occurrences of v, allocating M's vertices to these ordered
runs is exactly a choice counted by P_j(M). The external arc and endpoint
requirements depend only on the contracted word. Count the Hamiltonian
words as h_j and the endpoint-constrained two-path words as w_j. Then
d_j=h_j-w_j. Count the U-words as u_j.

A single path has at most n runs. Two unrestricted paths could have n+1,
but U fixes the second start outside M; this consumes one of the old
vertices needed to separate runs. Thus U has at most n runs. W fixes the
first end and second start outside M and obeys the same bound. Empty first
paths are included only for U. This proves (3), including its finite bound.
It is the endpoint refinement of the known substitution law in
[THM-3121 / path-cover-walk-content-substitution-kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).

For a transitive block of size a,

    P_j(TT_a)=j! S(a,j),                                    (4)

because each nonempty subset has one increasing path. The matrix of these
values for a,j=1..n is triangular with diagonal j!, hence invertible.
Consequently the following are equivalent:

1. Replacing v by TT_a gives offset 1 for every a>=1.
2. It gives offset 1 for each a=1,...,n.
3. (d_1,d_2,...,d_n)=(1,0,...,0).

This is a complete criterion for this specified operation, rather than an
extrapolation from initial sizes. Call such an unselected vertex a **stable
inflation port**. Equation (3) then gives the stronger identity

    delta(B[M at v],O)=H(M) for every nonempty tournament M. (5)

In particular arbitrary tournament substitution preserves unit offset
exactly when M is transitive. A transitive tournament has H=1; a
nontransitive tournament contains a directed triangle, whose three paths
extend injectively as other vertices are inserted, so its H is at least 3.
This is the precise obstruction to recursively plugging arbitrary new
seeds into a stable port. The independent checks replace each named port
by C_3 and obtain offset 3. QED.

## E3 Certified infinite coefficient families

The following coefficients in (3) are exact. The finite certificates were
computed by quotient-word dynamic programming and independently recovered
from actual substituted graphs at all a=1,...,n using (2) and triangular
inversion. The tuple is (order, canonical code, selected bitmask, port).
Adjacency rows are in `inflation_bank_7.json` and emitted by `construct.py`.
Every displayed row has defect coefficients (1,0,...,0).

| Family | Base tuple | Nonzero u_j, in order | U after a-vertex inflation |
|---|---|---|---|
|seven|(5,258,3,3)|20,7|7*2^a+6|
|three|(6,321,12,4)|64,24|24*2^a+16|
|twentyone|(6,12866,17,5)|88,42|42*2^a+4|
|binary29|(7,524673,5,5)|93,32|2^(a+5)+29|
|ternary30|(6,385,12,0)|89,90,30|30*3^a-1|
|ternary50|(7,1281,12,0)|149,150,50|50*3^a-1|
|ternary78|(7,1537,20,0)|233,234,78|78*3^a-1|
|ternary90|(7,1537,24,0)|269,270,90|90*3^a-1|
|mixed|(5,258,9,1)|32,21,4|4*3^a+9*2^a+2|
|mixed11|(5,258,24,1)|44,36,8|8*3^a+12*2^a-4|

The formulas use P_1=1, P_2=2^a-2 and P_3=3^a-3*2^a+3.
Combining them with (0) proves actual Hamiltonian-count coverage:

| Odd c | Parameter range | All 1+c*2^s realized for |
|---|---|---|
|7*2^r+3|r>=1|s>=2|
|3*2^r+1|r>=1|s>=5|
|21*2^r+1|r>=1|s>=3|
|2^(a+5)+29|a>=1|s>=1|
|b*3^a-1, b in {30,50,78,90}|a>=1|s>=1|
|6*3^r+9*2^r+1|r>=1|s>=2|
|2*3^a+3*2^a-1|a>=1|s>=3|

For the r-families use a=r+1. The factors of two in U are respectively
2,16,4,1,1,2,4 in the table's row order; the odd factors are as displayed.
The four ternary families and binary29 therefore give complete lines,
including s=1. Other rows give tails; their unbounded below-tail problems
are not silently declared solved. These finitely many exponential
coefficient families have only O(log X) members c<=X and do not exhaust
odd coefficients.

Example: ternary30 with a=3 has base order 8, U=809 and delta=1.
Its t=1 closure has order 10 and H=1619. With a=20 it gives the exact
coefficient c=30*3^20-1 and the whole corresponding line, without a large
graph search. Each base adds a-1 vertices; a tail of length t adds t+1.

## E4 What the complete small census classifies

The universe is one canonical representative of every ordinary tournament
on n=1,...,7 vertices, all nonempty O, and every v outside O for which
delta(B,O)=1 before inflation. There are 1,407 such presentations; 634
have stable ports. Among all presentations, 125 have a negative higher
defect coefficient. No positivity premise may be imposed on d_j.

Define P_n=sum_[B] #{(O,v): delta=1 and v is a stable port}. Define Q_n
similarly using an unordered pair of distinct stable unselected ports.
These sums count all marks on an ordinary representative; they do not
quotient the marks by Aut(B). They are independent of its chosen labels.

| n | A000568 ordinary bases | P_n | Q_n |
|---:|---:|---:|---:|
|1|1|0|0|
|2|1|2|0|
|3|2|9|3|
|4|4|26|19|
|5|12|64|66|
|6|56|152|186|
|7|456|381|489|

All 763 separately stable pairs in this census are jointly stable for
**every pair of nonempty transitive block sizes**. This finite-base,
infinite-size statement is stronger than the initial (2,2) test. Its
certificate is the full two-variable analogue of (3), with only the
(1,1) defect coefficient equal to one. `joint_ports.py` computes it exactly.
For two marked vertices the quotient on the marks is acyclic; between
successive occurrences of either mark on a path one needs an unmarked
vertex. Therefore each coordinate is bounded by n for two-path covers,
and the coefficient computation is finite. Product triangular inversion
gives the equivalent finite grid test. No all-order assertion that separate
stability implies joint stability is made.

**Two exact failures of short testing.** These are unit seeds at a=1:

- Base order7, code9281, O-mask17, port2 has d=(1,0,-4,0,...).
  It is a unit seed exactly at a=1,2, since delta=1-4P_3(TT_a).
  Its first defects are 1,1,-23,-143,-599,-2159.
- Base order7, code131713, O-mask68, port5 has d=(1,2,-2,0,...).
  It is a unit seed exactly at a=1,3. Its defect is
  1-2(3^a-4*2^a+5), with initial values 1,5,1,-43,-239,-955.
  The parenthesized expression vanishes at 1 and 3, equals -2 at 2,
  and strictly increases for a>=2 because its next difference is
  2*3^a-4*2^a>0. This proves the complete size classification.

The first failure is invisible to doubling alone. The second shows that
admissible block sizes need not form an initial segment. These statements
classify the actual endpoint obstruction, not just an observed mismatch.

## E5 Finite control of several inflated modules

For a fixed base and a set of marked modules, retain the subset of ordinary
vertices already used, the last symbol, the path-break flag, and the
required endpoint memberships. Reading a marked symbol changes no ordinary
subset; reading an ordinary vertex consumes it. Every transition also
tests the original directed arc. Disallow adjacent copies of the same
marked symbol within one path. This is a finite weighted automaton for
the contracted words. Separate its one-path, W-cover and U-cover accepts.
The weight of a word is the monomial recording visits to each module.

Its multivariate series is rational, by the finite transition matrix
(I-A)^{-1}, interpreted formally. Every repeated-symbol cycle consumes at
least one marked weight, so the formal expansion is well defined. If the
marked subtournament is acyclic, a path can traverse a given mark at most
once between ordinary vertices, and the series is a polynomial with an
explicit finite bound. A marked directed cycle permits arbitrarily many
visits in the contracted language, so no general finite-degree bound can
be inferred merely from the number of marks. Actual finite modules still
truncate the series because P_j(M)=0 for j>|M|.

Apply the linear functional x_i^j -> P_j(M_i) independently in each
coordinate to recover actual path counts. This is not ordinary evaluation
at a scalar x_i. The products of ordered path-cover counts allocate each
module's vertices uniquely to its visits. Thus it gives a rigorous finite
description of the infinite substitution language and identifies exactly
which endpoint data a scalar projection discards. It does not establish
surjectivity of its evaluated counts.

## E6 Endpoint gates remaining after gcd one

For this seed route the necessary and sufficient equations are

    delta(B,O)=1,       U(B,O)=c*2^j for some j>=0.           (6)

Gcd one does not solve either equation. The first equation is a weighted
clique-sum constraint with positive even weights W_ab. For H(B)>1,
every solution O has at least two vertices, and distinct solutions form
an inclusion antichain: adding a vertex to a nonempty set strictly
increases its W-sum. This supports exact branch-and-bound pruning, not a
claim that a numerical subset sum automatically has a graph realization.
For H=1 the nonempty solutions are the singletons. Empty O has zero slope
and is excluded from a positive-c seed.

There are four separate obligations: the extension support must contain
the target; an actual selected endpoint set must satisfy (6); any further
operation must retain its required port/cover state; and the points below
the chosen ray start must also be covered. Stable inflation solves the
second obligation for the E3 families. E1 removes the affine congruence
obstruction after finitely many extensions. Neither result supplies a
uniform legal descent for all remaining odd c.

A sufficient completion would be a finite set of base seed states plus
source-valid inverse rules that take **every** other odd c to a smaller
rank, with retained endpoint state and all guards verified. This proves
termination into seeds. Chan–Pak additionally needs all points below
their ray starts, except 7 and 21, or a construction whose starts are
already 1. The named complete-line families satisfy this additional
obligation individually. The universal rules remain OPEN.
