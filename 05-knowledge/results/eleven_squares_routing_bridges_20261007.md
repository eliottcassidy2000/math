# Labelled routing, marked packing, and exact ordered Collatz packets

2026-10-07. **CITED / ACCEPTED PREMISES:** the supplied *The Fibonacci
Geometry of Eleven Squares* and the conditional multiplication interfaces
below; their proofs are not re-audited here. **PROVED:** the elementary
guard/cut/credit composition and permission distinctions below.
**FINITE-EXACT:** the declared controls. **OPEN:** universal Collatz
completion and an exhaustive cover by grounded guards. This package does
not claim new guard coverage or a new multiplication complexity theorem.

Artifacts: [exact program](../../04-computation/experiments/eleven_squares_routing_bridges_20261007.py)
and [output](eleven_squares_routing_bridges_20261007.out).

## 1. Inheritance and the useful permission

The closest existing mechanism is the native progression fibre product in
[carry interfaces, section 3](collatz_carry_interfaces_20261004.md).
That note already proves associative composition of the paid letters'
native guards, including empty intersections. The present extension adds
a strict original-input threshold and a minimum initial credit, then uses
the combined summary in balanced ordered storage with local updates.
The affine carry is inherited from
[reciprocal four-channel kernels](reciprocal_four_channel_kernel_20261004.md).
No priority claim is made for affine composition, generalized CRT, or
segment trees.

The canonical hostile is `(1,2)` versus `(2,1)`: equal multipliers, unequal
carries and source classes. The corrected near miss is replacing a native
guard by output integrality after cancellation. The least-used coordinate
here is the greatest original-input threshold needed to prevent an earlier
ROOT visit. Our board is **named records / ordered carry / native guard /
strict ROOT boundary / funding / supplied source**.

Douglas Colkitt's [routing audit](https://github.com/CrocSwap/integer-mult-bounds/blob/5fa9e9772e67414b8c712401379d2c2b4e4550ab/docs/nonadjacent-axis-audit.md)
uses direct nonadjacent field interchanges already allowed by the accepted
upstream interface. Axis reversal uses `floor(d/2)` swaps; exposing and
restoring all target axes uses at most `2(d-1)`. Axis names, validity masks,
and coefficient payloads travel together. The numerical processing order
is unchanged; numerical stages are not assumed to commute. These schedules
replace the counted quadratic number of adjacent moves with a linear
number of permitted full-field moves. The claimed multiplication result
retains the source's conditional status.

The snapshot read was commit `5fa9e9772e67414b8c712401379d2c2b4e4550ab`.
The [repository README](https://github.com/CrocSwap/integer-mult-bounds/blob/5fa9e9772e67414b8c712401379d2c2b4e4550ab/README.md)
states the resulting exponent saving `2^-78`; we do not derive or test that
algorithmic bound. Our small tuple controls test only the elementary
labelled schedules used to describe the transfer.

This separates an implementation restriction from a mathematical one.
Restricting already-permitted swaps to adjacent positions creates an
avoidable routing cost. Permuting chronological Collatz instructions would
change the map and its domain, so that is not an available primitive.

## 2. What Eleven Squares actually preserves

The supplied local source is
`C:/Users/Eliott/Downloads/eleven-squares.pdf`, PINGYOU LTD,
*The Fibonacci Geometry of Eleven Squares*, 35 pages. Its Figure 1 was
inspected visually, and the contact/mark statements below were read in
sections 2--4. We use the source as an accepted premise.

For the specified two-parameter Trump contact family, the source reduces
feasibility to two signed clearances `g_A,g_B >= 0`. Its Theorem 3.5 and
Corollary 3.6 give a semialgebraic bijection between that feasible family
and marked tori, preserving both clearances and the enclosing side `T`.
The period is

    M_w C_- Z^2,   M_w=I+wQ,   C_-=4I-Q,
    Q=[[0,1],[1,1]],   Q^2=Q+I.

The canonical lift `z` of the mark retains the contact data through

    M_w z = xi + lambda e,
    g_A = Lz-d,   g_B=lambda/s^2,   T=3+chi*z.

For fixed `w`, different admissible `T` have the same unmarked torus but
different marks. Thus the unmarked torus is a decisive quotient-loss
example, rather than an alternative carrier for the full packing problem.
The source explicitly does not identify each individual planar square
with one torus cell. Its minimum is the minimum in this stated contact
family; this use does not enlarge that domain to all eleven-square
packings. The square torus at `w=1/3` is outside the physical parameter
window, another reason to retain the admissibility data.

The quarter-turn in Theorem 4.1 transports the frame, cells, labels,
canonical lift and clearance functionals together. It then preserves the
two inequalities and `T`. This is the same *type of lossless operation*
used below: change a representation and pull back every predicate.
There is no map here from square configurations to Collatz sources and no
deduction of arithmetic coverage from a packing density.

## 3. Exact partial-affine packets

A packet carries

    F(x)=(P*x+B)/Q,       P,Q positive integers,
    x = r mod m,         m positive,
    x > ell,
    credit change delta, minimum input credit c.

The finite lower cut `ell` is rational; `ell=-infinity` is also permitted.
Inputs and credits are integers, and credits are nonnegative. The native
progression has integral outputs:

    Q divides P*r+B,     Q divides P*m.                  (1)

Use a canonical rational affine triple by dividing `P,Q,B` by their common
gcd, **without cancelling the native modulus**. A packet is an arithmetic
interface, not an assertion that an arbitrary affine map is a Collatz
step. Each leaf retains its separately proved operation type and proof.
For a token summary require `c>=max(0,-delta)`.

**PROVED composition.** For a first packet `u` and a second packet `v`,
pull the second guard to the original input:

    P_u*x = Q_u*r_v-B_u mod (Q_u*m_v).                  (2)

Let `g=gcd(P_u,Q_u*m_v)`. If `g` does not divide the right side, the
composite has empty domain. Otherwise divide (2) by `g`, invert the
remaining coefficient and intersect its residue class with
`x=r_u mod m_u` by generalized CRT. Incompatible residue classes again
give EMPTY. For a compatible result retain its exact residue and modulus.
The other summaries are

    P=P_v*P_u, Q=Q_v*Q_u, B=P_v*B_u+B_v*Q_u,
    ell=max(ell_u,(Q_u*ell_v-B_u)/P_u),
    delta=delta_u+delta_v,
    c=max(c_u,c_v-delta_u).                            (3)

An absent cut contributes `-infinity`. The formula for `c` is the usual
maximum prefix deficit; it includes all internal requirements carried by
the child packets.

**Proof.** At an admitted input, (1) makes `F_u(x)` an integer. Condition
(2) is exactly `F_u(x)=r_v mod m_v`, not a weakened divisibility test.
Since `P_u/Q_u>0`, the second lower inequality pulls back to the one in
(3). The second stage sees credit `credit+delta_u`, giving the final
formula. Substituting the two affine maps gives the first line. These are
necessary and sufficient conditions for the ordered composition.
Every compatible residue class contains arbitrarily large integers, so a
finite lower cut introduces no further emptiness test. Associativity
follows either from these three-stage conditions or from substitution in
(2)--(3). EMPTY is absorbing, and the unrestricted identity packet is a
two-sided identity. Signed carries and nonunit coefficients are allowed;
positive slopes are essential to the one-sided cut formula. Upper bounds
or negative slopes require a more general interval carrier.
Identity and associativity are semantic statements; canonical generated
packets also have exact dataclass equality. All well-typed `empty=True`
records denote the same empty partial map, and composition canonicalizes it.

This extends, rather than replaces, the inherited aligned-progression
formula. It is particularly useful when the proof includes an immutable
source bound or a no-earlier-ROOT condition.

### Strict first-hit boundary for ordinary valuation words

For a positive integer `a`, the ordinary odd step

    F_a(x)=(3x+1)/2^a

has exact native guard

    x=(2^a-1)*3^(-1) mod 2^(a+1),    x>1.              (4)

The cut excludes starting a further step at ROOT. Compose the leaves in
chronological order. If the word has length `r`, total valuation `A`, and
prefix affine data `(P_i,Q_i,B_i)`, with the empty prefix `(1,1,0)`, its
packet has full source guard and cut

    x=(2^A-B_r)*(3^r)^(-1) mod 2^(A+1),
    x > max_{0<=i<r} (Q_i-B_i)/P_i.                   (5)

The congruence includes the final oddness bit. Exact composition proves
that it is equivalent to all the valuation guards. The cut in (5) is
equivalent to every preterminal value being greater than 1. A nonempty
admitted positive trajectory therefore avoids ROOT before its endpoint;
the endpoint can equal 1. The packet itself does not promise that it does.

For example `(4,2)` has source residue `5 mod128` and strict cut `5`.
It rejects `5 -> 1 -> 1` but accepts `133 -> 25 -> 19`. Omitting the
cut silently turns a first-hit certificate into a padded one.

## 4. Balanced compilation is lawful; reordering is a different theorem

Store the ordered leaves in a balanced binary tree and put their packet
composition at each internal node. Associativity gives the same root
packet as sequential compilation. Replacing one leaf changes only its
ancestors: at most `ceil(log2 d)` packet combinations for `d` leaves.
Construction uses `O(d)` combinations and storage. Empty leaves are not
deleted from the labelled proof history merely because their summary is
EMPTY. Leaf labels and proof provenance are retained by the tree.

These are counts of exact algebraic operations and tree nodes. The
integers, moduli and rational thresholds grow with the expression; no
fixed-bit memory, Turing-machine speedup, bounded-time orbit algorithm,
or shorter decoded proof is claimed. A recomputed source guard must be
tested at the **same supplied integer**, not at a new representative of
that guard.

Two minimal controls distinguish the permitted operations:

    word 12: (9x+5)/8, source11 mod16;
    word 21: (9x+7)/8, source 9 mod16.                 (6)

Their slopes and valuation multisets agree; their maps and domains do
not. In contrast, moving a *labelled storage record* to another physical
slot, or changing parentheses while keeping chronological labels, does
not replace `12` by `21`.

The native paid library provides a second control. Inherited L is
`(9x-3)/16` on `219 mod256`, with credit4. Odd output alone holds on the
larger class `27 mod32`. The source `27` gives the formal output15 but
does not satisfy the native dependency proof. Our carrier retains256.
It computes `LB=EMPTY` and

    LG=(81x+53)/128 on219 mod256,
    net credit3, minimum starting credit0.

The credit awards are inherited local theorems. This compilation does not
prove a new award, supply a ROOT suffix, or infer realizability from a
balanced formal word. An unpaid source or an empty native intersection
remains rejected regardless of the stored affine expression.

### A whole-block change can be legitimate without free swaps

The concurrent parent found the two words

    w=111112223113,    v=121313121121.

They have equal length12 and cost19, and their carries satisfy
`B_v-B_w=4*3^12`. The exact packet guards give

    n=257727 mod1048576,
    F_w(n)=F_v(n-4).                                  (7)

This program independently checks the displayed identity and the first16
positive lifts. Both actual words have no earlier ROOT, and every source
prefix in these controls grows above `n`. The all-height family and its
coverage belong to the parent's separately proved replacement; no finite
control is promoted to that theorem here. Equation (7) illustrates the
available permission: supply a checked whole-block relation and transport
its source guard. Equality of a multiset, slope, or geometric shape is not
that permission. A smaller child's ROOT certificate is still a premise
when the relation is used for conditional induction.

## 5. The repelling-point chart retains a translation coordinate

The earlier [two-orbits synthesis, section 3](oai3_two_orbits_twos_and_threes_20261007.md)
routes to THM-4568 and the accepted matrix-multiplication growth profile.
Only the following elementary fixed-point identity is needed here:

    c_a=-1/(3-2^a),
    F_a(x)-c_a=(3/2^a)*(x-c_a),   a>=0.              (8)

For `a=0` this is unhalved tripling, center `-1/2`. For `a=1` its center
is `-1`, and two such letters give `(9x+5)/4` on `7 mod8`, with slope
`9/4`. The two centers must not be identified. When the next exponent is
different, changing from its predecessor's fixed-point chart introduces
a translation. Dropping that transition drops the ordered carry in (6).
There is no deduction here equating a Collatz multiplier and the exponent
of a tensor-complexity bound.

The previous Ellison control locates the exceptional **exponents**
`13,14,16,19,27` in an inequality comparing powers of2 and3; they are not
automatically the same domain as Collatz starting integers. Likewise,
the source's eleven labels and affine group do not identify a Platonic
object or make every chronological Collatz permutation lawful. Those
suggestions need a source-preserving map before they can support a proof.

## 6. Consequence and exact verification boundary

The source-to-target map in this package sends a labelled sequential
receipt to a balanced tree carrying the same partial affine operation.
It preserves the exact admitted integer sources, endpoint, native
congruences, strict ROOT boundary and credit requirements. The tree retains
the chronology/proof labels; its scalar root summary alone need not retain
every local proof type. The strongest current benefit is cheaper lawful
recompilation and exact rejection of broken interfaces.

For a fixed family of leaf programs, reassociation preserves the union of
their admitted source sets. It therefore cannot turn partial coverage into
global coverage. New coverage must come from a newly proved primitive or
whole-block relation, a genuinely enlarged local guard with proof, or a
newly certified endpoint. The eleven-square marking and direct-axis
routing suggest how to retain the necessary information during those
changes; they supply none of those arithmetic premises themselves.

Reproduction:

```text
python -B -X utf8 04-computation/experiments/eleven_squares_routing_bridges_20261007.py
python -B -O -X utf8 04-computation/experiments/eleven_squares_routing_bridges_20261007.py
```

The explicit universe is axis dimensions0..32; every pair of residue
classes for moduli1..9; 1,000 triples from ten generic signed-cut packets;
every valuation word on1..4 of length0..5, with literal replay of selected
positive lifts/off-guards and both end-leaf edits; and every H/G/A/B/L
program of length0..4 with native replay and three input-credit controls.
Additional controls cover nonunit congruence pullback, earlier ROOT,
cancelled native guards, the fixed-point charts and the stated16 block
lifts. Thirteen malformed exact-input cases are rejected, including bool
aliases, floating cuts, nonintegral progressions and unfunded leaf summaries.
The checks are explicit exceptions and remain active under `-O`.
