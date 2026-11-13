# Arithmetic clocks, physical births, and a cylinder obstruction

**PROVED IN NOTE:** the scoped two-adic slope obstruction, the physical
lift of the inherited signed inverse clock, and the hereditary-spectrum
operation rules. **FINITE-EXACT:** the explicitly declared controls in
`check.json`. **INHERITED:** the signed inverse-fibre classification, the
binary-shell stencil, and the supplied404 ray. **REFUTED:** an all-k
unit-ready405 return in the exact indegree-one interface, by the incoming
complete k17 exclusion reviewed in [K17_INCOMING_AUDIT.md](K17_INCOMING_AUDIT.md).
**OPEN:** broader interfaces, complete Chan--Pak coverage, and Collatz convergence.
No external novelty claim is made.

## 1. Inheritance and the useful distinction

The closest arithmetic mechanism is
[thirtysix_digits, section3](../../thirtysix_digits_20260925.md): its signed
half-step, every iterate, and full inverse-germ bijection already give the
formulas below. The closest physical mechanisms are
[fixed balanced interfaces](../../changing_core_returns_20261113/returns/FIXED_BALANCED.md)
and the literal404 source in
[birth budget, section4](../../changing_core_returns_20261113/returns/BIRTH_BUDGET.md).
The canonical hostile is a valid physical birth whose arithmetic source
becomes divisible by3. The corrected near miss is using an address map as
a scalar Hamilton-count map. The least-used sidecar is the actual induced
subtournament count spectrum.

Our board is actual source / deletion parameter / address language /
halving clock / hereditary count spectrum. Each map below says which
target it preserves: an odd arithmetic target, a Hamilton count, and a
deletion count are different objects.

## 2. Recovered signed clock, including all guards

For positive odd n and sign e in {+1,-1}, put

    U_e(n)=oddpart(3n+e),    a=v2(3n+e)>=1.

The inherited bijection sends (n,e) to (m,a), where m is positive odd and
3 does not divide m. Conversely, e is the centered residue of 2^a m
modulo3, and n=(2^a m-e)/3. The exact half-step is

    D(n,e)=(2n+e,-e),
    D^j(n,e)=(2^j n+J_j e,(-1)^j e),
    J_j=(2^j-(-1)^j)/3.                              (1)

It fixes m and increases a by1. In particular D^2=(4n+e,e),
D^3=(8n+3e,-e), and D^6=(64n+21e,e). This is an inverse-fibre clock,
not the forward time of either fixed-sign Collatz map. The inverse
((n+e)/2,-e) is legal as an odd source exactly when n+e=2 mod4,
equivalently a>=2. At a=1 the inverse stops, retaining its unresolved
target m.

The further requirement 3 does not divide n has an exact guard:

    2^a m != e mod9.                                 (2)

For e=+1 the allowed residues are4,7 mod9; for e=-1 they are2,5.
For each fixed sign there are two allowed exponent classes modulo6.
The complete six-state clock (n mod3,e) is

    (1,+) -> (0,-) -> (2,+) -> (2,-)
          -> (0,+) -> (1,-) -> (1,+).

On deleting the two zero-residue states, the first-return cycle is

| source | next admissible state | physical exponent gap | source update |
|---|---|---:|---|
| (1,+) | (2,+) | 2 | 4n+1 |
| (2,+) | (2,-) | 1 | 2n+1 |
| (2,-) | (1,-) | 2 | 4n-1 |
| (1,-) | (1,+) | 1 | 2n-1 |

This table follows directly from D modulo3. Its four admissible transitions
consume six physical half-steps. A four-state quotient without those edge
lengths loses the physical clock. It only controls one fixed inverse
fibre: recursive inverse generations require additional ternary digits,
as proved in
[reverse-tree pieces, section1](../../collatz_mod6_20260917_reverse_tree_pieces.md)
and illustrated by targets7,25 in
[predecessor sections](../../inverse_predecessor_sections_20261005.md).

The user's sequence2,4,5,7,8,10,... is the positive3-free integers with1
omitted. Positive3-free integers are uniquely2^v m with m odd and3-free;
the valuation v is necessary to recover the original integer. They are
not the odd-source domain alone. Ordered from1, their j-th element is
floor((3j-1)/2).

Two cheap hostiles prevent changing clocks silently. Starting at(n,e)=(1,+),
the correct D^3 gives(11,-), still with odd target1. The alternatives
(8n+1,-)=(9,-) and(8n-1,-)=(7,-) have odd targets13 and5. Also, a literal
forward odd rule n->cn+e for even c>=4 remains odd forever and is strictly
increasing on positive integers. It is not the target-preserving sibling
operation4n+e in(1), whose forward problem still uses multiplier3.

## 3. The exact physical lift, and what it does not construct

Use the inherited rooted state(h,s,q,F), with h the full Hamilton count,
s the number of Hamilton paths starting at the root, q the root-deletion
count, and F the ordered two-path-cover count whose first nonempty path
ends at the root (the second path may be empty). Suppose an actual source
has balanced state(h,t,t,h). Add a new vertex beaten by the retained root
and beating every other old vertex, retaining the old root. The inherited
birth identity gives

    R:(h,t,t,h) -> (2h-t,t,t,2h-t).                    (3)

Thus A=h-t doubles, while t is fixed. If A=2^a m with m positive odd and
3-free and a>=1, encode e=A mod3 in{+1,-1}, n=(A-e)/3. Equations(1),(3)
give a commuting square: the actual birth R maps to D. It preserves the
odd arithmetic target m and the physical deletion count t. It changes h.
The graph and root are retained as essential sidecars; an integer germ
does not supply a graph realizing its proposed(h,t).

The lift extends to any nonzero positive A after retaining b=v3(A) and
encoding A/3^b. Both a and b are fixed/advanced exactly as stated: R adds1
to a and fixes b and the odd part after removal of3^b. This rescaling
does not create a physical source. In particular the desired amplitude
405*2^k has b=4 and normalized odd part5. Its simple signed arithmetic
clock is compatible with the previously proved absence of a balanced
source(h,71,71,h) on that ray.

The actual404 family supplies a decisive positive/hostile control:

| physical h | t | amplitude A | encoded(n,e) | n mod3 |
|---:|---:|---:|---|---:|
|879|71|808|(269,+)|2|
|1687|71|1616|(539,-)|2|
|3303|71|3232|(1077,+)|0|
|6535|71|6464|(2155,-)|1|

All four literal graphs have balanced state(h,71,71,h), and each is the
physical R successor of the previous one. Pruning n divisible by3 deletes
the legitimate third graph. Consequently the four-state first-return
model must retain the skipped birth count and intermediate graph if it
is used as a physical compiler. Neither the existence of this arithmetic
model nor its exact finite controller proves growing-t target coverage.

## 4. A two-adic slope obstruction even with changing cores

We use the actual two-flank compiler and its four-point detector from
[central compiler C1--C2](../../central_modular_chanpak_20261110/compiler/PROOF.md).
For a named core K and cut, d=H(K), the exact count has the form

    Phi(A,B)=d(g(A)+g(B)-2)+qAB+eA+sB+c
             +u(A mod2)+v(B mod2),                    (4)
    g(0)=g(1)=1,
    g(2z+1)=g(z),  g(2z)=g(z)+z+(z mod2).

All coefficients belong to the actual core and cut. On every fixed binary
cylinder A=2^r z+a, the coefficient of g(z) remains exactly d. Hence,
holding the other free address fixed, the stencil at z=4j,...,4j+3 is

    Phi_0-Phi_1-Phi_2+Phi_3=-2d.                      (5)

**Corollary (two-adic slope rigidity).** Suppose a proposed exact identity
between two actual named-core compilers holds for every pair of free
addresses in a binary cylinder, or for all sufficiently large such
addresses. Both source and destination addresses are fixed suffix appends
of these free variables, possibly with their roles interchanged. If

    Phi_new=lambda Phi_old+c,                        (6)

with constant rational lambda, then lambda=d_new/d_old. In particular
v2(lambda)=0, because every tournament Hamilton count is positive odd.

**Proof.** Apply(5) in one free variable to(6). It gives
d_new=lambda d_old. Oddness gives the valuation conclusion. The same
argument works beyond any fixed finite address threshold. QED.

Consequences, with their exact boundaries:

* No such rule implements2Phi-71 on even one full cylinder. Allowing an
  unbounded bank of changing cores cannot repair this obstruction.
* A policy depending on finitely many low address bits has a cylinder on
  which its rule is fixed. Therefore it cannot implement that recurrence
  everywhere by these rules. A fixed finite composition has the same
  limitation, since suffix maps pull low-bit guards back to low-bit guards.
* A uniform exact-H repair in this language must keep d_new=d_old. This
  is only necessary: equality of d does not establish(6).

This is a new consequence of the inherited detector, not a new detector.
It does not prohibit a single numerical target hit, a thin address
language, unbounded-height guards, arbitrary full-address transducers,
selective nonuniform module edits, or graph operations outside this
compiler. The actual R ray in section3 has one constrained amplitude;
it is not a full two-address cylinder.

The inherited
[square-graceful rebase](../../h13_paired_descent_20261111/square_graceful/PROOF.md)
is a concrete escape: deletion counts13,9,5 yield the same actual counts
77*2^(k-1)-31 only on B=2^k-2. That address language has density zero.
Its continuation Z retains the full equality; an intervening O then Z
reveals the discarded response. A pointwise repair and a repeatable
uniform repair therefore have different proof obligations.

Another constructive escape is the inherited
[offset transport](../../multiplicand_saturation_chanpak_20261111/chanpak/OFFSET_TRANSPORT.md):
actual additions to distinct transitive modules commute as

    L(D,U,h)=(2D-1,2U,2h-1),
    R(D,U,h)=(D,U,2h-D),    h=D+2^a U.

This genuinely grows an offset, but its D is a marked-interface offset,
not automatically the balanced deletion parameter t. The stronger
arbitrary-module readout H(M)+(alpha+beta*2^a)Q2(M) couples H(M) and Q2(M);
choosing those statistics independently would lose the actual source.

## 5. A complete obstruction for the recovered module carriers

The preceding arbitrary-module interface can be tested at every module
size, not just on its transitive examples. Write P2(M)=Q2(M) for ordered
spanning two-path covers allowing either bin empty. For every tournament
M with H(M)=71, the inherited complete71 census gives

    oddpart(P2(M)) in {253,67,247}.                    (7)

Here is the all-order reduction. Hamilton counts multiply over ordered
SCCs. Since71 is prime, M has exactly one nontrivial strong component S
of count71, with transitive prefix and suffix. P2 likewise multiplies
over ordinal sums, by independently restricting each of its two bins;
P2(TT_a)=2^a. Thus P2(M)=2^a P2(S). The complete strong71 head in
[birth budget](../../changing_core_returns_20261113/returns/BIRTH_BUDGET.md)
has P2(S) in{1012,1072,1976}, whose odd parts are respectively253,67,247.
Its inherited order bound and complete census, not this note's small
control, justify that these are all possibilities.

**Corollary.** No readout71+c P2(M), with H(M)=71 and positive integer c,
equals405*2^k+71 for any k>=0. Otherwise oddpart(P2(M)) would divide405;
none of the three possibilities does.

In particular all three physical carriers recovered above fail on the
entire405 ray if the inserted module has H71, even allowing arbitrary
module order and arbitrary positive tail length. Their coefficients
c=alpha+beta*2^a are positive integers. This is a new operation-specific
consequence of the inherited71 census, not a new census. It does not
exclude these carriers with a changing module Hamilton count, other
interfaces, or the presence of the target counts themselves. The useful
next change is therefore to the actual module's count or response
interface, rather than its size alone.

## 6. Hereditary count spectra as an operation obstruction

For a tournament T let

    I(T)={H(T[S]): S subset V(T)},   H(empty)=1.

These are induced-subtournament counts, not the global Hamilton spectrum.
If T embeds induced in W, then I(T) is a subset of I(W). Every vertex-only
birth word retains the initial tournament induced, irrespective of root
changes, new-vertex cuts, and the number of births. Therefore if a complete
target family F satisfies h not in I(W) for all W in F, no source with
Hamilton count h can reach F by vertex-only births. Arbitrary initial
same-h preparation does not help: the prepared graph itself is still an
induced h-count witness in the final graph.

For ordinal sum T->U there is an exact product rule

    I(T->U)={ab:a in I(T), b in I(U)}.                  (8)

Indeed every induced subset splits into two ordinal factors, and a
Hamilton path must traverse each nonempty factor in order, giving the
product of their counts. Conversely every pair of subsets can be chosen.
Since I(TT_a)={1}, dominated transitive padding leaves I unchanged.
This invariant is not determined by H, even at fixed order: at order6,
C3[TT3,1,1]->TT1 has H9 and I={1,3,5,9}, whereas C3->C3 has H9 and
I={1,3,9}. Thus a same-count repair can change a later birth obstruction.
Thus a complete strong-core target classification plus this padding rule
can turn a finite induced-subset calculation into an all-order birth-only
obstruction. Completeness and the exact padding normal form are separate
dependencies: checking a convenient sample of targets is insufficient.

This is the appropriate invariant for the parent's881-to1691 audit. It
does not itself certify the census or assert nonexistence of1691-count
tournaments. It identifies why changing old arcs after a birth or deleting
vertices can escape a birth-only obstruction: the initial induced graph
need no longer survive.

## 7. Reproduction and claim boundary

Run `python -B .../geometry/check.py` and the same command with `-O`.
The checker verifies the signed inverse/guard tables over the declared
finite universe, the stencil on binary cylinders, the four literal404
graphs by independent endpoint subset DP, their actual birth relation,
all eleven inherited strong71 core P2 values, and the hereditary product
identity on small literal tournaments. These
are controls for the proofs; finite ranges do not replace their quantified
arguments. `source_pins.json` hashes every directly used inherited file.
No source family, conjecture, or global metagraph is declared classified.
