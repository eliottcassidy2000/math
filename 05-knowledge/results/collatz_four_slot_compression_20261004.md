# Four-slot compression: divisor balance, trace recursion, and paid Collatz exits

2026-10-04. **PROVED:** the explicit profile maps and their loss boundaries,
the ordered-carry identities, six all-height paid word families, their
all-depth permuted extension with constructed root certificates at every
depth, and a new source cell adding1/8 of the
specified power-three exponent domain.
**INHERITED PROVED:** divisor balance, the tournament constructions, the
proper Collatz rank, and the preceding coverage banks. **FINITE-EXACT:**
the recorded audits and the larger minimum-exit probe. **OPEN:** the
all-parameter sharper exit rule, universal positive coverage, and complete
negative-basin classification. No novelty or literature-priority claim.

The connection has a practical consequence: the smallest switched word121a
has profile(2,1,1). Retaining its order-sensitive carry proves that its safe
exit starts at a=3, improving the previous sufficient threshold4. This adds
the entire family `187+256t ->119+162t` and the exponent class `27+64t`
for sources3^k. The permuted word1123 also adds `7+256t ->5+162t`.
It extends to every balanced profile P^hQR by the paid words1^h2a below.
The previous specified exponent-domain coverage30.2206 percent
becomes42.7206 percent. These are paid dependencies, not supplied root proofs
for every resulting child.

## 1. Inheritance and the shared object

Closest mechanism: the [paid-controller theorem](collatz_paid_portrait_controllers_20261004.md)
and the [divisor balance classifier](seam_prime_balance_20260925.md).
Canonical hostile: words12 and21 have equal multiplier but different carry;
the new exit1212 grows from123 to157. Corrected near miss: equality
`F=S+U` is a count identity with U contained in S, not a disjoint partition
of divisors into the original S and U. Least-used sidecar: the exact defect
created by changing the order of two operations.

Anchor: pay a previously weak Collatz exit. Niche: an operation-aware bridge
from prime multiplicity to a four-event controller. Wildcard: the same
repeated slot in the reciprocal trace square and the strong four-tournament.
The board is **multiplicity / order / carry / guard / quotient / original rank**.

Relevant inherited constructions:

- [Divisor balance](seam_prime_balance_20260925.md), sections1--3: the shapes
  `p,p^3,p^2qr`, and the need to retain prime names under multiplication.
- [Four-color decoder](duck_decoder_20260925.md), section2: a different,
  already proved ten-letter symmetric-square model. Its transported color
  action is not a divisor-poset automorphism. We keep that distinction.
- [THM-1960, tournament substitution](../../01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md):
  the strong four-tournament has a two-vertex module; scalar path counts
  lose information under cyclic substitution.
- [Intrinsic halving decoder](decoder_halving_20260925.md), sections1--4:
  nested regular odd tournaments grow by two, while raw substitution and
  pair-halving have an explicit compatibility obstruction.
- [Quadratic operation atlas](quadratic_escape_rank_atlas_20261004.md):
  word doubling squares a multiplier and folds to the Chebyshev map.

The reusable local object is a **four-slot multiset with one repeated
species**: `{P,P,Q,R}`. This specifies multiplicities; it does not yet
specify arithmetic values, chronological order, arcs, or legal inputs.

While this work was in progress, incoming commit61f6e81861 supplied
[four-channel kernels and exact tournament modes](recursive_small_structures_20261004.md).
Its reciprocal four-product expansion independently gives the same trace
mechanism below; its Mode B recovery specifies source/sink enclosure.
Its [guarded dependency kernel](collatz_recursive_dependency_kernel_20261004.md)
also constructs completed families with retained terminal proofs. The
principal additions here are the divisor-profile map, the carry extremum,
and the resulting new paid domains, rather than a priority claim for
the shared four-product algebra.

## 2. The controller parameters have an exact divisor-profile image

Keep the controller integers q,r separate from prime labels P,Q,R.
Use the positive odd map U(n)=(3n+1)/2^v, with v=v2(3n+1);
a word records these exact valuations in chronological order.
For q,r>=1 and a>=3, its word is

    w(q,r,a)=(12)^q 1^r a.

Send valuation1 to formal letter P, valuation2 to Q, and the marked exit
a to R. Forget word order but retain the exit label. Then

    ab(w)=P^(q+r) Q^q R,                                 (1)
    length(w)=2q+r+1,    cost(w)=3q+r+a.                  (2)

This is an actual map from words to a commutative monomial. Evaluating P,Q,R
at any three distinct primes gives a divisor lattice with the same exponent
profile. It is a code for the operation pattern, not the prime factorization
of the source n. The value of a must stay attached to R if the word is to be
reconstructed. At a=2 the distinct-species hypothesis fails.

At q=r=1, (1) is precisely P^2QR. More generally, put k=q+r. The divisor
lattice is

    C_(k+1) x C_(q+1) x C_2.

For the encoded integer N=P^k Q^q R, let F count proper nontrivial divisors,
S_k count those which are k-free, and U count those which are prime. Here
k-free means that no prime exponent is at least k. Because k>q,

    F=2(k+1)(q+1)-2,
    S_k=2k(q+1)-1,   U=3,
    F-S_k-U=2q-2.                                       (3)

**Proof.** Each divisor is one exponent vector in the displayed box.
Remove its bottom and top to get F. A k-free divisor omits the top P-face;
there are2k(q+1) such vectors including the unit, giving S_k. Each of the
three prime axes supplies one proper prime divisor. Subtraction proves(3).

Thus q=1 is exactly the balanced lane of this generalized statistic;
increasing r deepens its distinguished axis while preserving balance.
Increasing q also widens the second axis and incurs defect2(q-1).
At q=r=1, k=2 is ordinary squarefreeness and (3) is10=7+3.
This gives q and r different structural roles. It does not make the
divisor defect a decreasing Collatz rank.

## 3. The strong four-tournament and the twelve-state lattice

Take vertices P0,P1,Q,R with arcs

    P0 -> P1,   {P0,P1} -> Q -> R -> {P0,P1}.             (4)

This is `C3[TT2,1,1]`, the strong four-tournament. Its unique pair module
is{P0,P1}, and its score sequence is(1,1,2,2). The three module sizes
are(2,1,1), the same multiplicity object as in section2.

There is a precise lattice map. Consider only the event dependency
`P0 precedes P1`, with Q,R independent. An order ideal chooses zero,
one, or two P-events, and independently chooses Q,R. Send it to

    (i,j,l) -> P^i Q^j R^l.

This is a lattice isomorphism from its12 ideals to the divisors of P^2QR:
intersection/union of ideals become gcd/lcm. Removing bottom and top
leaves10. Seven remaining ideals have i<=1. The other three form the
proper portion of the extra P-face and have the same count as the prime
axes. This is the mechanism behind the count balance.

Two information losses matter. The dependency poset retains only the
within-module order from(4); it does not retain the cyclic intermodule
arcs. Also, the16 arbitrary subsets are not the12 ideals. Their occupancy
map has four fibres of size2 and eight of size1. It does not preserve
intersection: {P0} and{P1} have the same occupancy but disjoint intersection.
Choosing prefix ideals repairs this, with the P-order retained.

Nor does the controller word121a become a Hamiltonian path of(4) under
these color names. A path P-Q-P would require Q to distinguish the two
members of a module, which is impossible. The actual correspondence is
the multiplicity and ideal lattice, not the full directed dynamics.
The four-core has five Hamiltonian paths, independently enumerated.

## 4. Why the same four slots give J^2-2, and where +2 belongs

Start with two weights u,v. Their Cartesian square has four ordered slots

    u^2,   uv,   vu,   v^2.

Commutativity identifies the two mixed weights. Before numerical
specialization their multiplicity pattern is again(2,1,1). Set z=u+v
and d=uv. Then the exact two-coordinate recursion is

    (z,d) -> (z^2-2d,d^2).                              (5)

The term2d is the sum of the two mixed slots. It is not an arbitrary
adjustment. At (u,v)=(lambda,lambda^(-1)), d=1, so

    J(lambda)=lambda+lambda^(-1),
    J(lambda^2)=J(lambda)^2-2.                          (6)

This is the classical Chebyshev recurrence in the normalization
`J(lambda^m)=2T_m(J(lambda)/2)`; see the
[NIST recurrence table](https://dlmf.nist.gov/18.9.T1).
The preceding repo atlas already applied(6) to affine word doubling.
The present four-slot expansion identifies exactly which two terms
were compressed. At degenerate weights further numerical collisions
can occur; formal slot names still distinguish their roles.

For the tournament recursion, two actual operations recovered from the
repo have size shadows:

    substitution square T[T]: X -> X^2,
    Mode B source -> T -> sink: X -> X+2.

The first squares the vertex set. The second adjoins a source and sink,
preserving the Hamiltonian-path count by adjoining those endpoints to
each path. The incoming note repairs the arm count in
[THM-291, Mode B](../../01-canon/theorems/THM-291-mode-b-multilinear-recursion.md):
the two arms and apex contain2n-1 new free tiles. The regular two-vertex
extension from the halving-decoder note is a separate +2 construction,
with different crossing arcs. For (6),
`J(lambda)^2=J(lambda^2)+2` has the same square-and-two arithmetic because
the two mixed products are both1. These are explicit operation formulas,
not a conjugacy of both tournament operations to both Collatz operations.

Indeed the exact plus-two transport under J is

    J(lambda+2)=J(lambda)+2-2/[lambda(lambda+2)],          (7)

when the denominators are nonzero. Equation(7) exposes the missing
correction in a supposed common scalar model. Another loss is
`J(lambda)=J(lambda^(-1))`: the trace alone cannot tell expansion from
contraction. Retain the selected multiplier, affine carry, and legal word.
This is how the small identity can recur at many scales without claiming
that the different native dynamics are interchangeable.

There is a second exact tournament square with a different observable:
the order join of two copies of T has `H(T -> T)=H(T)^2`, while its order
is2X. Every Hamiltonian path must finish the first block before entering
the second, so concatenation is a bijection of path sets. In contrast,
substitution squares the order but need not square H: the inherited
`H(C3[C3])=3159` differs from `H(C3)^2=9`.
Thus the transferable primitive is **compose two copies and specify the
observable**. Saying only "square" drops that specification. The +2
extension above refers to vertex order, not a universal increment of H.

## 5. Repair the commutative profile by its exact order defect

Write f_a(x)=(3x+1)/2^a, and apply words from left to right. For words u,v,
compare the two words u,a,b,v and u,b,a,v. Both have the same length L
and total cost A. Write their formal actions as `(3^L*x+B)/2^A`.
Then

    B_(u,a,b,v)-B_(u,b,a,v)
      =3^len(v) 2^cost(u) (2^a-2^b).                    (8)

**Proof.** The middle two maps have equal slope, with constant difference
`(2^a-2^b)/2^(a+b)`. Precomposition by u changes neither that constant
difference nor its independence of x. Postcomposition by v multiplies it
by `3^len(v)/2^cost(v)`. Clear the common denominator2^A.

Adjacent swaps therefore reconstruct every order defect. In particular,
among all permutations of a fixed valuation multiset, increasing order
minimizes B and decreasing order maximizes B. With a fixed prefix, apply
the same assertion to the remaining suffix. This comparison concerns
formal affine actions. It does not authorize swapping a trajectory's
actual valuations: the exact input cylinders change with B.

There is a second precise compression loss for our three registers.
The clock map `(q,r,a) -> (L,A)` in(2) has integer kernel generated by
`(1,-2,-1)`. For r>=3 and a>=4, replacing the registers by
`(q+1,r-2,a-1)` preserves both clocks, but changes B by

    2^(3q+2) [3^(r-1)-2^(r-1)] > 0.                    (9)

This follows either from(8) or the closed carry

    B(q,r,a)=5*3^(2q+r+1)-12*3^r*2^(3q)-2^(3q+r+1).

For example121115 and121214 both have `(3^L,2^A)=(729,2048)`;
their carries are925 and1085. Thus even these structured words cannot be
compressed to their two clocks without losing their guards and endpoints.

## 6. A four-slot payment theorem and the new coverage

**PROVED.** Let a>=3. Every actual positive odd Collatz word which starts
with1 and is a permutation of the multiset `{1,1,2,a}` reaches a smaller
positive odd integer in four steps. It also strictly lowers the inherited
root rank. There are six distinct word orders for each a.

**Proof.** Each has multiplier81/2^(a+4). By(8), fixing its first1,
the largest carry is attained at `(1,a,2,1)` and is

    B_max=45+14*2^a.

The source has n=3 mod4. It cannot be3: the first four valuations at3
are1,4,2,2, which do not have the specified multiset. Therefore n>=7.
For every a>=3,

    (2^(a+4)-81)*7-B_max=98*2^a-612>0.

So every permitted carry gives `81n+B<2^(a+4)n`, proving strict descent.
The orbit is positive and its final oddness comes from the exact guard.
For n=3 mod4 the inherited energy is3(n-1)^2/4; every smaller positive
odd m has energy at most3(m-1)^2/4, with root1 handled separately.
Hence the original lexicographic rank decreases. No convergence premise
is used.

This four-step arithmetic statement uses U(1)=1. If a word reaches1
earlier, a first-hit certificate stops there:151 has valuations1,1,10
to1, so the fourth letter2 in(1,1,10,2) would be root padding. The
specific new1213 and1^h2a families have no such intermediate hit: their
initial runs grow, and a valuation2 can reach1 only from1 itself.

Each word has a single exact residue class modulo2^(a+5). All these
four-step cylinders are disjoint, since they specify different actual
prefixes. Their union has odd-relative density

    6 sum_(a>=3) 2^(-(a+4))=3/32.                       (10)

For natural density, truncate the exit size: the tail lies in a finite
union of cylinders demanding a high valuation at one of positions2--4,
whose densities tend to zero. This justifies the infinite sum. Much of
this union was already handled by earlier rules;3/32 is not all new gain.

The specific new word is1213. Its entire source and endpoint are

    n=187+256t ->119+162t=m,   t>=0.                     (11)

The intermediate states are `281+384t`, `211+288t`, and `317+432t`;
their next exact valuations are2,1,3. The source's first valuation is1.
Thus all four guards hold for every parameter, and n-m=68+94t>0.
The boundary is sharp in a: word1212 sends123 to157 and does not pay.

**Disjoint addition to the preceding named banks.** The new source cell
has initial word1213. The rational single-anchor bank with h=1,q>=2
starts1212; all its other h, and the -17 bank, disagree earlier. The
switched bank either has q>=2 (1212), or q=1,r>=2 (1211), or q=r=1 with
final valuation at least4. The checkpoint155 mod2048 starts1211. Finally
the old binary16 bank has debt index e=1 here and requires suffix16 or6
after the first reset2, both excluded by the actual suffix13.
So(11) adds exactly1/128 odd-relative binary coverage.

The permutation1123 supplies a second whole cell:

    7+256t ->5+162t,   t>=0.                             (11b)

Its intermediate states are `11+384t`, `17+576t`, and `13+432t`;
the next valuations are1,2,3. It pays since the source minus child is
2+94t. It is disjoint from(11). The only potentially matching rational
anchor bank has h=2 and starts1121, while(11b) starts1123. The other
banks disagree earlier, and the same debt-one suffix test excludes the
binary16 rows. Thus the two cells together add1/64 binary coverage.
In particular the old bank's hostile7 now has a paid route to5.

**All-depth permutation transfer.** The profile argument yields a larger
bank, not just these two cells. For every h>=2 consider

    v(h,a)=1^h 2 a,
    L_2=3,
    L_h=min{a>=3:2^(h+2+a)>2*3^(h+2)},  h>=3.           (11c)

Require the actual last valuation to be at least L_h. At h=2 the
four-slot theorem applies. At h>=3 the endpoint is

    y=alpha*(n+1)-2^(-a-1),
    alpha=3^(h+2)/2^(h+2+a)<1/2,

so every actual positive source n>1 has0<y<n and pays its original rank.
This is a uniform proof over all h, a, and source heights.

The word has exactly the balanced k-free profile P^hQR at k=h. It is a
permutation of the earlier switched word `12 1^(h-1) a`. Both have the
same two clocks; the earlier word's carry exceeds the new one by

    12*(3^(h-1)-2^(h-1)).

Thus the profile leads to a lawful family of reordered candidates; the
carry difference and the independently rebuilt input guards decide what
can be transferred. No trajectory is reordered in place.

Each h contributes a coarse cylinder of odd-relative density
`2^(-(h+L_h+1))`. Initial run length distinguishes h. The other added
banks have incompatible initial runs or demand a1 after the first reset2,
whereas this bank demands at least3. Against binary16, the only overlap
has next valuation exactly6. It occurs precisely at h=2,...,6; remove
its exact cylinder, of odd-relative density `2^(-(h+8))`, for each h.
For example h=2 removes391 mod2048. Therefore the disjoint new bank has

    delta_perm=sum_(h>=2)2^(-(h+L_h+1))
               -sum_(h=2..6)2^(-(h+8))
              approximately0.019033158994196.          (11d)

The h>H numerical tail is at most2^(-H-4), using L_h>=3. Its sources
also lie in the deep initial-ones cylinder `n=-1 mod2^(H+2)`; this
supplies natural-density convergence. All sources are7 mod8, so the
entire bank is disjoint from1213 and from every power of3. The cell(11b)
is already included and must not be counted twice.

We have `3^27=187 mod256`, and3 has order64 modulo256. Hence the entire
new power-three family is

    k=27+64t,  source n=3^k,  t>=0.                     (12)

Inside the specified exponent domain k=3 mod8, (12) occupies exactly1/8.
The second cell is7 modulo8 and contains no power of3, whose residues
modulo8 are1 and3.
It is disjoint from the preceding additions by the source-prefix proof.
Adding1/8 to that frozen domain's previous coverage gives

    0.427206155834282... .                              (13)

The combined binary addition is the preceding binary addition plus
`delta_perm+1/128`, an increment of approximately0.026845658994196.
With binary16 and the strengthened original-source ternary bank, the
specified residual odd-relative density becomes

    0.109660841806587... .                              (14)

CRT gives the increment `(delta_perm+1/128)*(1-delta_ternary)` outside that bank;
the inherited controlled tails justify its natural density. Exact rational
enclosures and the prior-artifact hash are recorded in the JSON.
Other checkpoint positions and other repo rules are not exhausted by this
comparison. These numbers describe paid dependency guards, not densities
of hypothetical counterexamples or completed root proofs.

## 7. A reusable method, and the next test it produces

The discovery procedure used here is concrete:

1. Expand one local composition into distinguishable slots.
2. Record which slots coincide and form its multiplicity object.
3. Map that object into another setting using the native operation there.
4. Compute the exact defect of the information being forgotten.
5. Use a bound on that defect to prove the target predicate on legal inputs.

Here the four-slot pattern led from divisor balance to the controller's
valuation multiset. The defect was the affine carry. Its extremal bound
proved payment at the previously omitted exit3. This is a proof transfer
through an explicit intermediate object, not a transfer from the word
"balance" or a matching cardinality.

The small construction can be nested. Repeated blocks remain compressed
only when their interfaces carry the data needed by the next operation:
prime-support overlap, tournament ports/path covers, or affine carry and
valuation guards. Replacing an object by a smaller number of containers
alone is not a well-founded decrease.

There are also **completed root certificates at every depth** in the
new permuted family. For every h>=2 and t>=0 set

    a=3^(h+1)*(2t+1)-1,
    n=2^(h+1)*(2^(a+1)+1)/3^(h+2)-1.                   (15a)

Then n is a positive odd integer and its exact first-hit word is1^h2a,
ending at1. To prove integrality, induction by cubing gives
`2^(3^j)=-1 mod3^(j+1)`; raising to the odd power2t+1 preserves this
congruence. Here a is even and at least26. Write
`n+1=2^(h+1)M` with M positive odd. The first h valuations are1.
Their endpoint is `(2^(a+2)-7)/9`; its next step has exact valuation2
and reaches `c=(2^a-1)/3`. Finally3c+1=2^a, so the next valuation is a
and its endpoint is1. The preceding states exceed1, giving first hit.
For h=2,t=0 this is `13256071 --(1,1,2,26)-->1`.
These exits exceed the thresholds in(11c).

These are explicitly constructed completed subfamilies, not a conclusion
that every child in the density computations is rooted. They implement
the same terminal-interface principle as the concurrent dependency
kernel through a different word and an elementary closed congruence.

A sharper open target is now testable. For every q,r>=1 let a0 be the
least a>=2 for which `2^(3q+r+a)>3^(2q+r+1)`. Let P,D,B be the resulting
formal coefficients, and let c be the least positive member of its
coarse source cylinder `c=-B/P modD`. The sufficient all-height inequality

    (D-P)*c>B                                           (15)

would pay that entire cylinder, including actual larger final valuations.
The previous general theorem uses a factor-two safety margin. We proved
the sharper rule at q=r=1. The exact probe finds (15) throughout
q,r=1..128, but this is **FINITE-EXACT**, not an all-parameter theorem.
The remaining question is an arithmetic bound on the least guarded residue;
no global Collatz conclusion follows from this finite success.

## 8. Reproduction and scope

Run:

    python3 -B 04-computation/experiments/collatz_four_slot_compression_20261004.py
    python3 -O -B 04-computation/experiments/collatz_four_slot_compression_20261004.py

The [script](../../04-computation/experiments/collatz_four_slot_compression_20261004.py)
writes matching [JSON](collatz_four_slot_compression_20261004.json) and
[output](collatz_four_slot_compression_20261004.out). It uses explicit
exceptions, exact integers/rationals, no convergence oracle, and no imported
prior implementation. The prior density artifact is a named input with a
recorded SHA256; its strengthened ternary interval is independently rebuilt.

Declared universes: all16 four-event subsets and all12 ideals; all24
four-core vertex orders; profiles q,r<=16;1600 rational weight pairs;
576 adjacent-swap contexts;2016 same-clock register moves; all six paid
orders for a=3..64 with symbolic progression guards and five literal
parameter controls each;4096 members of each new cell; all64 exponent phases plus101
progression controls and eight expanded powers; the saved prior binary
cylinders; the all-depth permuted bank through h=64 with all exact-six
overlaps removed and a rational tail bound; and16384 minimum-exit parameter
pairs; and21 literal completed-family controls at h=2..8, t in{0,1,3}.
Negative controls retain
the raw-subset intersection failure, the impossible P-Q-P tournament path,
the same-clock distinct carries, the unpaid1212 exit, and the root-padding
boundary at151.

The proofs supply the infinite quantifiers in sections2--6. The finite
minimum-exit probe in section7 has no such promotion.
