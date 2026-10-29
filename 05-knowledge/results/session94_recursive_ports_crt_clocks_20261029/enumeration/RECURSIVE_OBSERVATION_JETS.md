# Recursive observation jets: SCC ports and prime-root child slots

**Status: PROVED scoped recurrences; FINITE-EXACT implementations.**
The inherited A000568/Burnside and prime-composition identities are not new.
The contribution here is a small scalar-jet compiler for particular joint
observations: occupied strong components, named prime-root child slots, and
selected root-cycle weights. It avoids constructing the full aggregate prime
cycle index when those observations, rather than that index, are the target.
It does not generate prime representatives, individual automorphism groups,
root-switch rows, or transporters. No new A000568 record or measured speedup
against the strongest scalar counter is claimed.

Source pin: `88c969273f2b` (2026-10-29). Output provenance: Session94 enumeration
lane. The independent audit lane reviewed both derivations and separately
checked ordinary marked SCC counts through five vertices by labelled Burnside.

## Inheritance and scope board

* Session86 `A000568_COMPOSITION.md`: all-order Burnside, aggregate prime
  plethystic recovery, canonical strong-module species, and the boundary
  between aggregate indices and class/group inventories.
* Session87 `FAST_COUNTS.md`: degree-local scalar counter and distinct fixed
  vertex formula; Session89 `counting/README.md`: fused positive channels.
* `singleton_pair_reconstruction_20261028/README.md`: ordinary and strong
  count prefixes already reach350; this session does not claim a new record.
* THM-4816, `prime-cycle-index-tower-order-eleven-twelve-thirteen-types-from-non-rigid-primes-and-the-ledger.md`:
  the corrected, audited finite 11/12 type inventories used below.
* Session93 fixed-point laws distinguish joint port actions from independent
  marginals. The present problem concerns all tournaments, not just W_k/C_k.

Board: Burnside channels; canonical SCC order; prime-root child actions;
formal derivatives; divisor/Frobenius substitutions; faithful-carrier losses.
Canonical hostile remains the order12 same-cyclic-triples 81/27 pair. Our
smallest new type-check hostile is C3[TT2,1,1]: three root children but four
vertices. These quantities must not be conflated.

## 1. From existing channels to all ordered port moments

For r>=1 tournament layers and s>=0 binary vertex attributes, aligned under
one permutation, retain the inherited channel

    F_n(r,s) = sum_(odd lambda |- n) 2^(r E(lambda)+s length(lambda))/z_lambda.

Let R_j(n) be the sum over unlabelled aligned structures of the number of
orbits of ordered j-tuples of vertices under their actual joint automorphism
group. Repetitions are allowed. R_0(n)=F_n(r,s), including the empty object.
Using Stirling numbers of the second kind gives the all-order identity

    R_j(n) = sum_(i=0)^min(j,n) {j brace i}
             2^(r binom(i,2)+s i) F_(n-i)(r,s+r i).                 (1)

Proof: a fixed j-tuple has a partition into i equality classes. The i ordered
distinct named vertices are individually fixed by an automorphism. Their
internal orientations contribute r binom(i,2) bits, their attributes s i,
and their incidences with each remaining vertex contribute r i additional
binary attributes. Equivalently, remove i singleton cycles in Burnside.
This extends the inherited distinct-point formula; it is not a replacement
for the actual joint action on a particular tournament.

Only J+1 positive channels are needed for every j<=J. They can use Session89's
shared traversal. The checker deliberately uses independent odd-partition
sums at a small bound, not the fastest existing engine.

## 2. Canonical ordered components give a marked reciprocal

Form the union digraph of all r layers, retaining both arcs if layers disagree.
Its SCC quotient is a transitive tournament: between distinct SCCs all layers
must agree, and the condensation is acyclic and complete. Thus there is one
intrinsic ordered sequence of components. For r=1 this is the familiar
strong-component decomposition. For r>1, **union-strong is not the same as
strong in every layer**; two opposite arcs on two vertices are the minimal
counterexample.

Set A(x)=sum_n F_n(r,s)x^n and R_j(x)=sum_n R_j(n)x^n. Let I_j(x) count the
same ordered j-port decoration restricted to a single nonempty union-strong
component. Form exponential marker series in the named-port labels:

    A(x,t)=sum_j R_j(x)t^j/j!,   S(x,t)=sum_j I_j(x)t^j/j!.
    A(x,t)=1/(1-S(x,t)).                                          (2)

The port labels are distributed among the ordered components; no component
permutations remain. Multiplication of exponential marker series is exactly
the binomial label distribution, including possible collisions within one
component. Formal differentiation of (2) yields

    I_0 = 1-A^-1,
    I_j = A^-2 R_j
          - A^-1 sum_(i=1)^(j-1) binom(j,i) R_i I_(j-i), j>=1.     (3)

All arithmetic is in formal power series, with A(0)=1.

### Exactly b occupied components

Let G_b,j count structures with j ordered ports visiting exactly b SCCs.
Boundary conditions are important:

    G_0,0=A;   G_0,j=0 for j>0;
    G_b,j=A sum_(i=1)^j binom(j,i) I_i G_(b-1,j-i), b>=1.           (4)

The first A is the arbitrary unvisited gap preceding the first occupied
component, and i>=1 prevents counting an unoccupied component as occupied.
Equivalently,

    sum_(b,j) G_b,j u^b t^j/j!
      = A / (1-u A(S(x,t)-S(x,0))).                               (5)

At u=1 this is A(x,t), so sum_b G_b,j=R_j. In particular

    same-SCC ordered pairs = R_2 - 2 R_1^2/A,
    different-SCC ordered pairs = 2 R_1^2/A.                      (6)

This is an aggregate joint observation computed without a class bank. It is
not obtainable by squaring one-point counts at the same degree, nor does it
restore individual joint port order.

Ordinary tournament examples:

|n|classes|strong classes|ordered pairs|same SCC|different SCC|
|--:|--:|--:|--:|--:|--:|
|2|1|0|4|2|2|
|11|903753248|884464590|109287476096|108891621298|395854798|
|12|154108311168|152310149735|22187481894592|22147346586994|40135307598|

## 3. Prime-root child slots from scalar jets, without the full prime index

Now return to one ordinary tournament layer. Use symmetric-function cycle
indices B=Z_T (including empty), F=B-1, S=1-B^-1, and the inherited aggregate
prime index P (prime quotients have at least3 slots). Canonical strong-module
composition gives

    K := S-p_1 = P composed with F.                               (7)

Let delta denote evaluation p_l=x^l and Psi_i the Adams substitution
p_l -> p_(i l). Critically, partial differentiation with respect to p_1
changes F but **does not change Psi_i F for i>1**. Define

    A_j = delta(partial_1^j B),
    J_j = delta(partial_1^j K),
    H_j = delta((partial_1^j P) composed with F).

The A_j are directly available from the inherited colored channels:

    A_j(x) = 2^binom(j,2) sum_(n>=0) F_n(1,j)x^n.                 (8)

Here the derivative has removed j units of degree; this is NOT R_j(x) from
section1. Compute derivatives of S by (3), using A_j rather than R_j, then
subtract x from J_0 and1 from J_1. For j>=2 no singleton subtraction remains.
The ordinary one-variable Faà di Bruno chain rule now gives

    J_j = sum_(h=1)^j H_h Bell_(j,h)(A_1,...,A_(j-h+1)).           (9)

Since Bell_(j,j)=A_1^j and A_1(0)=1, this is triangular:

    H_j = [J_j-sum_(h=1)^(j-1) H_h Bell_(j,h)] / A_1^j.           (10)

Finally,

    L_j(x)=(A_0(x)-1)^j H_j(x)                                  (11)

counts ordinary strong tournaments with j **ordered DISTINCT named children
of their canonical prime root**, modulo the actual automorphism group.
A derivative of P names and removes a root slot; F replaces the decoration
in that slot. Canonical decomposition guarantees that these are actual child
blocks, not arbitrary modules. L_0=S-x counts all strong nonsingletons,
including inflated prime roots; it does not count prime classes only.
Repeated child marks can subsequently be recovered by a Stirling transform.

Finite results (columns j=0,1,2,3,4):

    n3:   1,1,2,2,0
    n4:   1,3,6,6,0
    n11:  884464590,9630536734,95353397744,849834467648,6733390029112
    n12:  152310149735,1817710219292,19886726641176,
          197801301532516,1770733526329944

The n4 strong tournament is C3[TT2,1,1], with trivial automorphism group:
three fixed root children, four fixed vertices. This catches the most likely
misreading of the formula.

### Computational benefit and exact limitation

For a fixed observation order J, obtain J+1 scalar channels through N and
perform truncated univariate convolutions. The displayed straightforward
Bell implementation uses O(J^2 N) stored coefficients and O(J^3 N^2)
coefficient operations after those channels, rather than constructing every
odd-partition monomial of the full prime index through N. These are arithmetic
operation/storage bounds, not bit-complexity or measured end-to-end speedups.
Big-integer costs and the inherited channel engine remain. No improvement of
the state quotient or of bare A000568 counting is claimed.

This is an observer-specific compiler, not universal recovery from a few
moments. For example, regular C9 and C3^2 actions on9 points have identical
fixed-point moments of every order (1 identity fixes9;8 others fix0), but their
cycle indices differ and their binary-color orbit counts are60 versus64.
This is an action-level loss control, not a claim that these particular
regular actions are full groups of two specified prime tournaments. The
repository's81/27 and prime81 group-shape controls remain stronger warnings
against promoting aggregate data to concrete group inventories.

## 4. Divisor-triangular cycle jets: retaining the normalization and clocks

There is a broader first-derivative compiler, with additional weighted
channels rather than only the F_n(1,s) channels. For odd d define

    D_d = delta(partial_(p_d) B),
    h_d = delta((partial_(p_d) P) composed with F).

The exact chain rule is

    partial_(p_d)(Psi_i F) =
        Psi_i(partial_(p_(d/i)) B) if i divides d, and0 otherwise.

There is no extra i factor: p_d is an independent power-sum variable. Thus

    h_d = [D_d/A^2 - 1_(d=1)
           - sum_(i|d, i<d) h_i D_(d/i)(x^i)] / D_1(x^d),        (12)

where A=delta B. The denominator has constant1. Even derivatives vanish for
both tournament and prime-tournament indices. Formula(12) works over Q[[x]];
these are **rational aggregate cycle weights, not class counts**.

Removing one d-cycle from an odd partition gives an explicit source:

    D_d(x) = 2^((d-1)/2)/d *
      sum_(odd mu) [2^E(mu)/z_mu]
         product_(l in mu) (2^gcd(d,l) x^l).                      (13)

The factor1/d is essential. For d=1 this is the binary-color channel. For
other d the extra per-cycle weight is 2^gcd(d,l), not a constant palette size.
The existing future-GCD response state can still carry it: adding k cycles
of length l adds k gcd(d,l) to the exponent, without changing the future
responses. A shared positive-channel implementation may reuse traversal
metadata but must retain independent arithmetic weights. This packet does
not benchmark or replace that production engine.

The Adams substitutions are a precise arithmetic connection: only divisors
of a cycle length can feed its derivative. This is a formal composition
clock, not a Lucas prime-distribution or Collatz drift theorem.

### The actual 11 and12 ledger check

The calculation gives

    h_11(0)=32/11 = 3*(10/11) + 10/55.

The THM-4816 order11 inventory has three C11 prime classes and one Paley
prime with automorphism order55. Each contributes ten11-cycles normalized
by its group order, exactly explaining this value.

It also gives [x^3]h_9=16/3. The inherited order12 census contains eight
prime classes with C9 automorphism group, six order9 elements of cycle type
3^1 9^1 each. Their contribution is8*(6/9)=16/3. The competing cycle type9^1 1^3
is excluded for a prime tournament by the inherited isolated-cycle module
constraint, so it is the3^1 9^1 coordinate being read here. The full small prime index
check confirms this coefficient. This verifies the ledger's aggregate
cycle weight; without the inherited group-type census, the number16/3
alone would not identify eight classes or their groups.

## 5. Reproduction and bounded evidence

Run from this directory:

    python3 check_marked_components.py
    python3 check_prime_slot_jets.py
    python3 check_divisor_jets.py

Each was also run with `python3 -O`, with identical JSON outputs; all validity
checks use explicit exceptions, not removable asserts.

* `check_marked_components.py`: series through20 and four ports, four channels
  (r,s)=(1,0),(2,0),(2,1),(1,1);480 literal coefficient comparisons. Literal
  class/automorphism controls through orders5,4,3,4 respectively, with up to
  three ports. Checks every occupied-component count, not only totals.
* `check_prime_slot_jets.py`: scalar jets through20 and four distinct child
  slots;215 literal prime-root slot counts over every strong class through6.
  Modules are found independently by subset tests and overlap exclusion,
  then actual automorphisms act on their maximal proper strong blocks.
* `check_divisor_jets.py`: d=1,3,5,7,9,11, compared with the separately coded
  inherited full prime-index inversion through total weight15;60 rational
  coefficient comparisons in their valid truncated ranges. Checks h_1 against
  the separate Bell compiler. The reference reconstruction is a consistency
  round trip, not an independent discovery of the Burnside count formula.

Receipts: `receipt.json`, `slot_receipt.json`, `divisor_receipt.json`.
The corresponding `_optimized.json` files record the identical -O runs.

## Loss ledger / stopping boundary

Source: exact Burnside-weighted joint actions and canonical decomposition.
Target: declared finite-order aggregate port or cycle observations.
Preserved: component occupancy, ordered distinct child marks, selected
normalized cycle weights, and every formal source contribution.
Lost: representative classes, individual groups, coherent transports,
orientations of visited ports, and higher cycle/slot correlations outside the
chosen observation menu. Sidecar for those tasks remains the actual action,
child-isotype stabilizer and named transports. The81/27 pair must not be
identified by a cyclic-triple-only or moment-only carrier.

This closes the present bounded enumeration lane: three proved extraction
recurrences, literal hostiles, and cross-checks against the current family
machinery. It does not claim unrestricted prime classification, metagraph
connectivity, or any major open conjecture has been settled.
