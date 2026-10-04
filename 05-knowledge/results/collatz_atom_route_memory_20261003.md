# Atoms, signed lifts, and tournaments carrying certified routes

2026-10-03 (America/Denver). Owner-directed continuation of
[the modular tiling/atom note](tiling_modular_atoms_20261003.md).

**Status:** PROVED elementary identities and certificate constructions;
FINITE-EXACT for the stated computations; OPEN for universal coverage of
positive integers. This is a synthesis and proposed storage architecture,
not a novelty claim or a proof of Collatz.

## Inheritance, portfolio, and concept board

Closest proved mechanisms:

- [THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md):
  common-tail combs, finite-prefix realization, and the loss on taking closures.
- [THM-4512, corrected coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md):
  exact valuation words require one more binary bit than integrality alone.
- [Order laws](collatz_procgen_20260922_order_laws.md), raw 4n+1 section,
  and [affine blueprint](collatz_blueprint_20260921_affine.md): signed
  sibling lifts already exist; the atom identification below reconnects them.
- [Trunk chart](collatz_procgen_20260923_trunk_rh.md), Theorem 1.1:
  (4^j-1)/3 is an isometry on the 3-adic integers.
- [Carry composition](collatz_pentagon_operator_coherence_20260930.md),
  Propositions 3-4: ordered affine summaries compose exactly.
- [Recursive entry certificates](entry_20260927_board.md): three rule nodes
  can certify arbitrarily long guarded excursions; node count is not bit cost.
- Concurrent incoming work, integrated before publication:
  [paired prime rows](paired_prime_collatz_20261003.md), section 3, gives
  the four signed steps on 4N+/-1 and the consecutive-pair shortcut laws.
  It reinforces retaining the source coordinate as well as the sign.

Canonical hostile: on the minus sheet 5 -> 7 -> 5, whereas on the plus
sheet 5 -> 1. Corrected near miss: THM-4512's lost final oddness bit and
MISTAKE-554's confusion between two related maps and one map's conjugacy.
Least-used sidecars here: designated root, legal input cylinder, first-hit
convention, ordered carry, and clone ancestry.

| Lane / live object | Operation / retained fact | Loss / decisive probe |
|---|---|---|
| Anchor: atom central complement | Factor z-y=2^k m | Keep sheet; compare plus/minus 5 |
| Modular wedge | Pair moduli 2x and 2x+1 | Same shape is not same arithmetic table |
| Binary parent | ceil(x/2), digit deletion | Parent is not Collatz; retain actual edge equation |
| Niche: tournament route grammar | SCC order, +2 source growth, cloning | Internal changes can violate mod-3 guards |
| Wildcard: route-family splice | Arbitrary head + one shared certified tail | Finite-head realization does not imply all-source coverage |
| Summary/cache | Compose (A,p,S), root and domain | Same clock can have different carries |

The connection contracts are explicit in the sections below. The new
signed chart joins the first three lanes; SCC storage supplies a faithful
operation carrier for a certified basin; the splice supplies actual reuse
across an infinite family. Neither erases the designated-root coverage problem.
Cards used from META-PATTERNS: correct the object before sharpening the
technique; keep local support separate from bounded-height coverage; retain
observer, recurrence class, and finite head. No new meta-pattern promotion.

## 1. Put all three coordinates on one signed atom

Write a family member as (N,epsilon), epsilon=-1 for B and +1 for A:

    x = 2N+(epsilon-1)/2,   y=N,   z=4N+epsilon.
    z-y = 3N+epsilon.                                      (1)

For each F_x, the even modulus 2x and odd modulus z=2x+1 have the SAME
folded domain D_x={(a,b):1<=a<=b<=x}. Thus even/odd neighboring tables
share addresses but not entries. Atom N pairs the odd moduli 4N-1 and
4N+1 around 4N; even-modulus neighbors likewise have an odd midpoint.
Adding 2 to z alternates B_N -> A_N -> B_(N+1). Adding 2 to x preserves
the family letter and increments y by 1; doubling x sends F_x to A_x.
These are different operations and their coordinate must always be named.

The intrinsic parent x -> ceil(x/2) has exact distance bit_length(x-1)
to 1. B/A gives the discarded digit 0/1 of x-1. This is an unconditional
binary route home, but has not thereby become a Collatz route.

For odd positive N, define U_epsilon(N)=oddpart(3N+epsilon). Equation (1)
shows that the central complement stores its actual signed next odd value:

    z-y=2^k m,  k>=1, m=U_epsilon(N), odd m>0, 3 does not divide m.

Conversely every such (m,k) yields one signed edge. Choose epsilon in
{-1,+1} congruent to 2^k m modulo 3, and set

    N=(2^k m-epsilon)/3.

This is positive odd and has exactly the prescribed valuation. It is a
bijection between signed positive odd edges and these target/exponent pairs.
The restriction to odd N is needed for a Syracuse edge; (1) itself holds
for every positive N.
The sign here belongs to the literal numerator at input N. At different
inputs the same central number can have a different shortcut interpretation:
T_+(2N-1)=3N-1, whereas T_-(2N+1)=3N+1, with shortcut
T_sign(n)=(3n+sign)/2 on odd n. Thus a family letter without its source
coordinate is not a universal assignment of a Collatz sign.

## 2. Doubling and the mirrors form one exact ladder

Let J(N,epsilon)=(2N+epsilon,-epsilon). Then

    3(2N+epsilon)-epsilon = 2(3N+epsilon),
    U_(-epsilon)(2N+epsilon)=U_epsilon(N).

In the edge chart J is simply (m,k)->(m,k+1). Applying it twice gives

    J^2(N,epsilon)=(4N+epsilon,epsilon)=(z,epsilon),
    U_epsilon(4N+epsilon)=U_epsilon(N).                     (2)

Thus the atom z-map is the SAME-SHEET common-target lift. The central
complement doubles under J, quadruples under J^2; its odd part is invariant.
For h>=0, repeated same-sheet lifting is

    N_h=4^h N+epsilon*(4^h-1)/3.

The opposite-sheet use is decisively different:
U_epsilon(4N-epsilon)=6N-epsilon, with valuation exactly 1.
The invariance is a conjugacy for the predecessor-lift J, not for the
forward Collatz iteration. A single J changes the first-edge sheet; a
certificate must retain that change. J^2 preserves a fixed-sheet tail.

The earlier central blocks have exact integer self-similarities:

    C_A(4N+1)=4C_A(N)+[[0,1],[1,0]],
    C_B(4N-1)=4C_B(N)-I.

These equalities relate the explicitly reduced central entries at different
moduli; they are not an unproved global table conjugacy.

## 3. A tournament that intrinsically stores a certified odd route

Use U=U_+. For odd target m with 3 not dividing m, set

    k0(m)=2 if m=1 mod3, and 1 if m=2 mod3.
    pred_j(m)=(2^(k0(m)+2j)*m-1)/3,  j>=0.                (3)

These are exactly its positive odd predecessors. Their exponent parity
is forced, and pred_(j+1)(m)=4 pred_j(m)+1. A multiple of 3 can be an
odd source but cannot have an odd predecessor.

Define a guarded word w=(j1,...,jr) by backward decoding from terminal 1
using (3), requiring every inverse step to exist and every reconstructed
source to exceed 1. The latter rule excludes padding by the 1 -> 1
accelerated cycle. The empty word is home 1. Every accepted word is an
actual first-hit route to 1; every odd member of the convergent basin has
one such word, because its forward route and valuations are unique.

Let R_(2j+3) be the regular cyclic tournament on 2j+3 vertices: each vertex
beats its next j+1 cyclic neighbors. Make

    T_w = R_(2j1+3) -> ... -> R_(2jr+3) -> K1,             (4)

with every interblock arc pointing forward. Its strongly connected
components (SCCs) recover these blocks in their unique total order. Their
sizes recover each ji, and the terminal singleton supplies home 1. No
original vertex names or numerical orbit labels are needed. Backward
arithmetic reconstructs the source and all intermediate values.

This is an injective encoding into UNLABELED tournament isomorphism classes
of the odd convergent basin, with a decidable finite certificate grammar.
It is not a claim that arbitrary tournaments or arbitrary words are valid.
An untrusted graph verifier must check the component shapes as well as sizes,
or retain a construction witness; SCC sizes alone do not recognize the grammar.

The arithmetic predicate is preserved because inverse decoding proves each
edge equation. The graph keeps ordered inverse choices and the terminal;
discarding order, root, or arithmetic guards loses that predicate.

### A literal +2 operation

Increasing ONLY THE SOURCE block from R_(2j+3) to R_(2j+5) changes n to
4n+1 and preserves the entire accelerated tail. Example:

    R3 -> R5 -> K1 encodes 3 -> 5 -> 1,
    R5 -> R5 -> K1 encodes 13 -> 5 -> 1.

This can preserve every old edge: embed R_(2a+1) into R_(2a+3) by
i->i for i<=a, and i->i+1 for i>a. The new vertices are a+1 and 2a+2.
Every old winning cyclic interval crosses at most one of these inserted
gaps, so its new length is at most a+1 and it remains a winning interval.

Internal growth is not automatically lawful. Replacing the second R5 by
R7 in the first example would make its target 21; the preceding inverse
step then fails because 3 divides 21. The +2 rule is source-local unless
all earlier guards are reverified. The empty home code has a separate root
rule: its first nontrivial direct odd predecessor is 5, with source block R5.

### A literal doubling operation

Let L(T)=T[TT2], replacing every vertex by an ordered pair, and every arc
by the four parallel cross-fiber arcs. This is transitive cloning; it is
specified separately from the repo's double-round-robin / skew doubling.
For n=2^e u with u odd and certified word w, define

    C(n)=L^e(T_w).

Its strong components have sizes 2^e(2ji+3), followed by exactly 2^e
singleton components. The number of terminal singletons recovers e; the
preceding component sizes recover w. Thus C is injective on the entire
convergent basin and C(2n)=L(C(n)).

For e>=1, halving is intrinsic: in each cloned cyclic component, connect two vertices
if they have identical orientations to every third vertex. Connected
components of this auxiliary relation are the transitive clone fibers.
Within a fiber only adjacent vertices are so related; different cyclic
base vertices cannot be twins, since a regular tournament cannot have
twins with their forced degree difference 1. The terminal transitive block
is recovered in the same way. Pair consecutive positions in each fiber's
unique transitive order and contract. This recovers L^(e-1)(T_w).

For the ORDINARY map C_raw(n)=n/2 if even and 3n+1 if odd, stopped at 1,
the graph rewrite is exact:

- even: contract the clone pairs;
- odd n>1: delete the source SCC, obtaining the code S of odd target m;
  decode k=k0(m)+2j from the deleted block and tail, then output L^k(S).

The latter graph represents 2^k m=3n+1. For the shortcut map it would be
L^(k-1)(S). A recovered rank is e+sum_i(k_i+1), the ordinary remaining
step count. It decreases by one under this rewrite. This rank is defined
on certified codes; universal existence of codes remains Collatz itself.

## 4. Compress a whole family, not just one known trajectory

For a certified positive odd b, all

    h_j=(4^j(3b+1)-1)/3, j>=0

share U(h_j)=U(b). Choose ANY positive valuation word w=(a1,...,ap), p>=1,
and let A=sum ai and S satisfy

    Q_w(n)=(3^p n+S)/2^A,
    S=sum_(i=0)^(p-1) 3^(p-1-i)*2^(a1+...+ai).

Exactly one j0 in {0,...,3^p-1} makes 2^A h_j-S divisible by 3^p. All members

    n_k=(2^A h_(j0+3^p k)-S)/3^p, k>=0,                  (5)

are positive odd, follow EXACTLY word w, and then join the certified tail.

Proof of the phase: binomial induction gives
v3(4^d-1)=1+v3(d), hence v3(h_j-h_l)=v3(j-l) for j!=l.
Thus j mod3^p -> h_j mod3^p is a permutation. Proof of legality: reduce
the integrality condition modulo 3 to reverse the last edge; its value
(2^ap h_j-1)/3 is a positive odd integer. Substitute and repeat through
all p letters. Each reconstructed odd target forces exactly its indicated
2-adic valuation. This avoids confusing a coarse cylinder with an exact word.

The entire family has one affine recurrence:

    n_(k+1)=R n_k+K,
    R=4^(3^p),  K=(R-1)(2^A+3S)/3^(p+1).                 (6)

K is integral because R-1 is divisible by 3^(p+1). Example b=1, w=(1,1):

    j0=4 mod9,
    n_k=(4^(9k+6)-19)/27,
    n0=151 -> 227 -> 341 -> 1,
    n_(k+1)=262144 n_k+184471.

This explicit splice is a repackaging of THM-4501's finite-prefix
realization and the inherited trunk isometry, not a novelty claim.

### Growth of an internal block while retaining the stored head

The splice gives a useful theorem directly on tournament certificates.
Choose a block with p earlier blocks, and increase its index j by 3^p.
Its order increases by 2*3^p. Keep every other block fixed. The result is
again a valid first-hit certificate, with the first p valuation exponents
UNCHANGED. If the original source is n and its first p steps have summary
(A,p,S), the new source is exactly R*n+K from (6). This includes p=0,
where A=S=0, the growth is two vertices and n changes to 4n+1.

Proof: the selected block source h changes to
h'=R*h+(R-1)/3, while its target/tail is fixed. Since
v3(R-1)=p+1, h'-h is divisible by 3^p. Backward replay of the p original
valuation exponents therefore remains integral. All differences in earlier
values are positive; the intermediate target residues modulo 3 remain
unchanged, so their k0 and stored j values also remain unchanged. Thus no
earlier 1 is introduced. Conversely a positive index increment d preserves
this same valuation head only if 3^p divides d, by
v3(h_(j+d)-h_j)=v3(d). This is minimality for the specified valuation head,
not for arbitrary valid routes after changing the head.
For example (0,1) -> (0,3) changes 3 -> 5 -> 1 into 113 -> 85 -> 1.
It uses an index increment 2 at depth 1, but changes the first valuation
from 1 to 2. This is a hostile control against the stronger, false claim
that index increment 3 is necessary merely to keep that compressed prefix.

For the explicit family (5), the complete valuation word is
(1,1,18k+10), its stored word is (0,0,9k+4), and its tournament is

    R3 -> R3 -> R_(18k+11) -> K1.

Increasing the internal final strong block by 18 vertices preserves both
earlier odd steps and corresponds to n->262144*n+184471. This provides a
lawful internal macro where an arbitrary two-vertex internal edit can fail.
Every member reaches 1 in exactly three accelerated odd steps, or
18k+15 ordinary steps. Its tournament has 18k+18 vertices. The proof schema
stores these unbounded return times using one integer parameter and a
shared terminal, without replaying each expanded trajectory.
The factor 3^p is the exact residue precision needed to retain p earlier
inverse steps. It measures preservation of this history, not intrinsic
hardness of Collatz or a universal storage lower bound.

The splice is a legal prefix certificate, not automatically a reduced
first-hit certificate. For b=1 and w=(2), j0=0 gives the prefix 1 -> 1.
Before converting a splice into the tournament grammar, truncate its route
at the FIRST occurrence of 1. This preserves convergence and avoids adding
a forbidden home-cycle block to the canonical encoding.

Hostile control: the same mechanism works on the minus sheet around the
different cycle 5 -> 7 -> 5. Its family
(56*4^(9k+1)+19)/27 starts at 9 -> 13 -> 19 -> 7 -> 5 -> 7.
Therefore attaching every finite head to a root does not show that every
positive integer belongs to that root's basin. Keep sheet AND root.

## 5. Store guards and reusable suffixes alongside algebra

Use three separate records:

1. A rule: signed edge/word, inverse guard, terminal/root identity, and
   proof of the arithmetic identity. Parametric rules retain their bounds.
2. A block summary: (A,p,S) meaning Q(n)=(3^p n+S)/2^A, exact cylinder, and a pointer to a proven
   suffix. Exact plus words use n=(2^A-S)*3^(-p) mod2^(A+1).
3. A graph blueprint: ordered cyclic component parameters j and a binary
   clone exponent e. Expand the graph only when graph observables are needed.

For word u followed by v the summaries compose as

    (A,p,S) * (B,q,R) = (A+B,p+q,3^q S+2^A R).

The order matters: words (1,2) and (2,1) have the same clock (A,p)=(3,2)
but carries 5 and 7. Equal summaries are safe to share only with their
legal domains and endpoint interpretation retained. Proof-certified suffix
pointers and symbolic repeats produce actual reuse; relabeling arbitrary
parity streams alone gives no general compression.
For signed or mixed-sheet words, S is the signed ordered carry; the sign
word and its lawful input/endpoint domains remain part of the rule record.

The expanded tournament has 2^e*(1+sum_i(2ji+3)) vertices; it is not an
efficient byte encoding. Its modular-decomposition blueprint is compact.
For a search implementation, memoize first proven descent below a retained
source, not only arrival at 1; recursively linked smaller certificates then
prove termination by strong induction. Keep OPEN sources separate from
verified terminal families. A library's coverage remains an independent goal.

An attractive but unproductive alternative is orienting word pairs by their
affine commutation defect D_u S_v-D_v S_u, D_w=2^A-3^p. Since S_w>0 for
nonempty plus words, this only orders D_w/S_w. It is a total preorder,
with genuine ties; quotienting ties gives a transitive tournament. That
construction stores a scalar order, not the missing return structure.

The classical parity-coordinate precedent is Bernstein--Lagarias,
[*The 3x+1 Conjugacy Map*](https://websites.umich.edu/~lagarias/doc/bernstein.pdf),
equation (1.5): the entire parity itinerary is a 2-adic coordinate and
Collatz becomes a shift. This motivates retaining itineraries, but a
shift representation alone does not prove positive-integer convergence.

## Verification and next obligation

Reproduce with:

    python3 04-computation/experiments/collatz_atom_memory_20261003.py
    python3 04-computation/experiments/collatz_route_tournaments_20261003.py

The first script independently replays 10000 signed edges, 5010 inverse
chart pairs, 10000 binary-parent chains, and 2040 exact splice instances
(head lengths1..4, letters1..4, seeds1/5/27, parameters0/1). It tests carry
composition, exact cylinders and residue permutations as separate paths.
Seeds are verified by bounded explicit iteration. There are no external
tables or assumed convergent arbitrary inputs. The second script checks
the guarded tournament language, graph/SCC recovery and rewrite controls;
its companion output records its exact finite universe.
In that script the full j-alphabet 0..4 at lengths 1..5 gives 694 valid
and 3211 rejected words. All 3174 valid word/position pairs pass the
internal-growth law; 65134 bounded phase controls check its minimality
for the fixed exponent head. The graph tests use 18 carriers, 12 intrinsic
halvings, 30 ordinary odd-step rewrites and relabeling controls. The special
family is replayed at k=0..3, including three induced internal +18 extensions.

The next substantive search is a SMALL GUARDED RULE GRAMMAR whose covered
integer sets strictly increase under useful operations while leaving a
precisely described residual. The present construction supplies lawful
storage and infinite certified families. It does not prove universal coverage,
nor does it assert that all inverse words or all internal graph mutations work.
