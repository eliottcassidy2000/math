# Five, eight, and nine: exact mechanisms and a guard-transfer boundary

2026-10-04. **PROVED:** the elementary transfer identities and the all-depth
clock-versus-contraction obstruction below. **RECOVERED PROVED:** the named
repository mechanisms, with their original scopes. **FINITE-EXACT:** the
declared independent controls. No assertion concerns universal Collatz
coverage, or a numerical pattern shared by unrelated objects.

Artifacts: [script](../../04-computation/experiments/five_eight_nine_transfer_20261005.py)
and [saved output](five_eight_nine_transfer_20261005.out).

## 1. Four useful connections, with their types retained

| Mechanism and source | Exact map/predicate | Loss and sidecar | Coverage consequence |
|---|---|---|---|
| Discriminant5, field9, clock8: [golden prime clocks, sections1–3](golden_prime_clocks_20261003.md) | Multiplication by phi permutes the8 nonzero elements of F9 | Modulo reduction loses exact value and rational-integer membership; retain the exact pair and radix | A finite register can check a known guard; the clock does not find a paid route |
| Signed anchor−5 and the factor9/8: [negative-cycle shadow](collatz_negative_cycle_shadow_20261004.md), with [THM-4527, Syracuse owner labels](../../01-canon/theorems/THM-4527-syracuse-in-owner-labels-microcosm-and-difference-spectrum.md) | Actual word12 sends n+5 to(9/8)(n+5) on its exact binary cylinder | Endpoint integrality loses the last oddness bit; retain v2(n+5), actual source and original comparison source | At27 mod64 there is exactly one12 block, followed by a forced two-one corridor |
| Golden period8 and the arithmetic carry9n+5: [THM-4528, golden beta-map and holonomy](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md), [arithmetic realization gate](denominator_arithmetic_filter_20261004.md) | The phase(1+phi)/3 has parity10100000; its rational Collatz return is(9n+5)/64 | A valid golden cycle loses the arithmetic integrality predicate; retain halving cost and rational denominator | The fixed point is1/11, so a complete golden phase cycle supplies no integer-cycle or source-coverage certificate |
| Eight presentations, a five-element fiber, and path product9: [four-state tournament codec, sections1–2](tournament_recursive_four_state_20261004.md) | Eight fixed-path masks form four classes; the strong class has five masks and H=5; joining the two diamonds gives H=9 | Class labels lose named-flip response; H loses ordered SCC data | The graph can store an ordered proof word, but its scalar count supplies no arithmetic guard |

The closest proved mechanism is the exact affine carry with its source
guard. The canonical hostile is a legal modular or golden phase with a
noninteger arithmetic realization. The corrected near misses are a unit
clock confused with a contracting forward map, and a smaller predecessor
of the current state confused with a smaller child of the original source.
The useful sidecars are halving cost and the immutable source. Our live
board is **quadratic ring / clock / ordered itinerary / signed anchor /
binary fuel / ternary address / original-source debt**.

These are four different roles for small numbers, rather than four
instances of a single untyped operation. The proofs below identify which
ones actually share an affine numerator and which ones do not.

## 2. A genuine5–9–8 object: the golden field

Let O=Z[phi], phi²=phi+1. The polynomial's discriminant is5, which is2
modulo3 and is not a square. Consequently O/(3) is the field F9.
In the basis(1,phi), multiplication by phi is

    M(a,b)=(b,a+b), det(M)=−1.

The eight powers of phi are

    (1,0),(0,1),(1,1),(1,2),(2,0),(0,2),(2,2),(2,1).

They are exactly the nonzero field elements. In particular phi^4=−1,
phi^8=1. The norm a²+ab−b² alternates1,2 along this clock, because the
norm of phi is−1. Here5 is a discriminant,9 is a field cardinality, and8
is a multiplicative period. This is one explicit object carrying all three.

The role of5 changes at modulus5: eta=phi−3 is nonzero but eta²=0.
That ring has25 elements and the phi clock has order20. It is not F25.
Thus even a two-coordinate register needs its ring and unit predicate.
The inherited exact order at ternary precision a is

    ord_(O/(3^a))(phi)=8*3^(a−1), a>=1.

The script verifies the base field directly and every pair modulo3^a for
a=1,...,4. These register clocks do not preserve the rational-integer line:
even1 shifts to phi. The [literal golden digit reader](golden_digit_carry_20261003.md)
therefore retains the exact coefficient pair and radix position, not only
a clock phase. Literal source digits, Zeckendorf digits and an orbit's
parity itinerary remain different input types.

## 3. The same5/9 carry has different dynamics at budgets8 and64

For U(n)=oddpart(3n+1), the two-letter valuation word(1,a) has affine map

    F_(1,a)(n)=(9n+5)/2^(a+1).                         (1)

The first numerator3n+1 creates the carry1; the second creates
3*1+2=5. The last exponent changes the denominator but not this carry.
Its exact odd source cylinder is

    9n+5 = 2^(a+1) mod2^(a+2).                        (2)

The fixed point is5/(2^(a+1)−9). In particular:

| Valuation word | Ordinary cyclic parity | Affine return | Fixed point |
|---|---|---|---|
|12|10100|(9n+5)/8|−5|
|15|10100000|(9n+5)/64|1/11|

These parity strings use the ordinary branches n/2 and3n+1, so their
lengths are5 and8. The negative−5 cycle is an actual signed integer cycle.
The positive1/11 cycle is an actual rational-parity cycle with odd
denominator11, not a positive-integer cycle. Word15 is instead a perfectly
valid finite positive-integer route on35 mod128:35->53->5. That guard
excludes27. Retaining only the numerator9n+5 loses exactly this distinction.

The golden connection is also literal. Put theta=(1+phi)/3. Under the
upper golden map x->phi*x−1_{phi*x>1}, it has exact period8 and itinerary
10100000. This can be checked by eight exact quadratic-pair operations.
Its phase in(O/3) follows the8-cycle from section2. Replacing each parity
digit by its ordinary arithmetic branch gives the return in(1) with a=5,
whose fixed point is1/11. Its arithmetic orbit is

    1/11,14/11,7/11,32/11,16/11,8/11,4/11,2/11.

This independently reproduces the q=3 hostile in the arithmetic-realization
note: golden period completeness does not imply integer realization. The
general arithmetic gate retains the ordered carry K and requires
(2^e−3^s)|K. Here55 does not divide5.

THM-4528 concerns the complete **future parity itinerary** Theta(n), not
the positional embedding n->(n,0) in O/(3). Its convergence equivalence is
a tail predicate. None of the finite phase calculations decides that tail.

## 4. A new all-depth transfer obstruction, and a useful fuel law

Let G=F_12. Shift z=n+5. Then

    G(n)+5=(9/8)(n+5).                                 (3)

On the ternary ring R_a=O/(3^a), multiplication C by9/8 is well-defined,
because8 is a unit, and

    C^ceil(a/2)=0.                                    (4)

In contrast M, multiplication by phi, is a permutation of R_a.

**Clock-transfer obstruction.** The only function f:R_a->R_a satisfying
f(Mx)=C f(x) for every x is the constant zero function. More generally the
same conclusion holds for any permutation on a finite source set and any
nilpotent target endomorphism. Iterate the alleged identity through the
nilpotence exponent: f(M^q x)=0. Since M^q is onto, f is identically zero.
No linearity assumption is needed. In unshifted coordinates, the constant
is−5. Replacing M by any power, including the five ordinary steps of12,
does not repair the obstruction: it remains a permutation.

Thus a rotating golden digit register cannot be turned into this forward
ternary observer by a nonconstant fixed state quotient. This does not
forbid using a golden reader to store known congruence guards. It also does
not contradict the itinerary semiconjugacy in section3, whose source,
target and retained future information are different.

There is nevertheless an exact useful two-prime law for every legal G
step on positive odd integers:

    v2(G(n)+5)=v2(n+5)−3,
    v3(G(n)+5)=v3(n+5)+2.                              (5)

Consequently the number of complete actual12 blocks is exactly

    floor((v2(n+5)−1)/3).                             (6)

The extra oddness bit matters: at n=3 the formal expression G(3)=4 is
integral but even; the actual valuations are(1,4), not(1,2). Binary fuel
is consumed while ternary divisibility accumulates. Neither valuation
change by itself is a numerical payment against the original source.

## 5. What this says about27 mod64

For n=27+64t, t>=0, v2(n+5)=5. Thus exactly one12 block is available,
and

    G(n)=31+72t =7 mod8.

It cannot immediately enter H's155 mod2048 guard, whose members are3
mod8. Its next two actual valuations are both one, so the entire source
cell has the exact four-letter prefix

    (1,2,1,1),  endpoint=(81n+85)/32.                  (7)

This is a forced exit corridor, not a repeatable12 loop. It is a concrete
place to search for a new guarded smaller-child interface; the existing
H/G credit account alone does not supply one for every member of this cell.

There is also a tempting circularity. Every legal G endpoint x satisfies
x=4 mod9. The inherited smaller-predecessor rule

    m=(8x−5)/9

then gives a positive odd m<x with common future x. But when x=G(n), this
formula gives m=n exactly. For27->31 it returns27. It proves descent from
the *current*31 while paying nothing against the immutable source27.
The rule is recovered in
[the connectivity note, Proposition3](collatz_connectivity_from_rigidity_20261001.md)
and remains useful on unrelated sources; only this composition is an identity.
For q repeated G blocks, the complete inverse similarly returns n:
8^q(x+5)/9^q−5=n. The new source-aware selector must preserve this equality
boundary rather than count it as newly paid coverage.

There is a genuine positive composition on the219 mod256 quarter-child
subcell, developed in the concurrent
[paid-guard budget note](paid_guard_budget_20261005.md). Put n=219+256t.
The checkpoint in(7) is4K+1, where

    K=(81n+53)/128=139+162t<n.

That checkpoint is5 mod8, so U(checkpoint)=U(K). This K is4 mod9, and now the same
inverse-G rule produces a *new* smaller original-source child:

    L=(8K−5)/9=(9n−3)/16=123+144t<K<n,   G(L)=K.        (8)

The ternary sidecar makes the exact boundary visible:

    K+5=18(8+9t), v3(K+5)=2;
    L+5=16(8+9t), v3(L+5)=0.

Exactly one inverse G is available here. Unlike27->31->27, the quarter-child
switch changes the common-future branch before undoing G. This is genuine
source-relative progress, conditional on the displayed binary subcell.

The same equation interfaces with the inherited credit potential:

    L+5=(9/16)(n+5)+2,
    ((L+5)/(n+5))(9/8)^4 <=6561/7168<1,                (9)

with equality at n=219. The ratio decreases thereafter. K earns three
of these credits and L earns four, with the same bound because G(L)=K.
The paid-guard note owns the resulting adaptive payment theorem; this
package independently checks its algebra and common-future interface.
The subcell occupies one quarter of27 mod64. It does not cover the other
three quarters or itself provide a rooted terminal certificate.

## 6. The tournament interpretation stores different data

Fix the path0->1->2->3 and let a,b,c flip02,13,03. The8 masks have
isomorphism fibers

    {0}, {a}, {b}, {ab,c,ac,bc,abc}.

Their Hamiltonian-path counts are respectively1,3,3,5. Thus the strong
class has five presentations and five Hamiltonian paths; these equal
counts use its trivial automorphism group. The four-class projection is
not predictive for arbitrary named flips: ab and c are isomorphic, but
flipping a yields b and ac, which are not.

Let + and− be the two diamond classes. Both order joins +->− and−->+
have eight vertices and Hamiltonian count9=3*3. Their ordered SCC size
lists are(3,1,1,3) and(1,3,3,1), so the tournaments retain the order that
the count loses. This is an exact8/9 storage boundary. It supplies no
intrinsic arithmetic orientation on integers and no Collatz guard. The
sidecar is the ordered component decomposition plus the chosen arithmetic
serialization. The script independently computes all eight masks and
the two eight-vertex path counts rather than assuming those values.

## 7. Reproduction and scope

    python -B 04-computation/experiments/five_eight_nine_transfer_20261005.py
    python -B -O 04-computation/experiments/five_eight_nine_transfer_20261005.py

Finite universes: all7380 quadratic-ring pairs modulo3^a for a=1,...,4;
all81 linear maps F3²->F3² as a finite control of the stronger arbitrary-map
obstruction; the8 golden phases and8 rational cycle states; words(1,a)
for a=1,...,8 with three exact integer sources each; all odd n<=1023 and
all their585 legal G-repeat prefixes, including empty prefixes; the first64
members of27 mod64 and64 positive quarter-child/inverse-G compositions;
all8 fixed-path masks with all24 relabellings each;
and exact dynamic-programming Hamiltonian counts for the two8-vertex joins.

The all-depth proofs do not rely on these finite cutoffs. The concrete
consequence for coverage is an exact forced corridor, a test that rejects
a circular predecessor payment, and the verified positive interface(8).
The field and tournament correspondences
remain useful storage/verification mechanisms, not a proof of a paid exit
for every supplied source. No literature-priority claim is made.
