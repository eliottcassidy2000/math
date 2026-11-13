# A three arc module mirror and its quadratic norm output

**PROVED:** the universal construction and identities below, including the
precise primitivity criterion. **FINITE-EXACT:** the accompanying literal
controls. **OPEN:** prescribed-target coverage and Chan–Pak. The one-apex
reflection and its transitive specializations are inherited objects; the
new interface retains a named root through a constant three-arc change
with an arbitrary module and its converse.

## 1. The construction and the information it retains

For a nonempty tournament M put

    h = H(M),
    p = P_M(2) = sum_S H(M[S]) H(M[V(M)\S]),
    k = p−2h,

with H(empty)=1. Thus k counts ordered covers by two nonempty paths.
Write M^op for the converse. Its h,p,k are identical, by reversing every
path. Form the ordinal chain of seven blocks

    r -> s -> M -> t -> M^op -> u -> v,

where the five named letters are singleton vertices. Add a new apex a.
In T, a beats precisely the blocks r,t,M^op. In U, it beats precisely
r,s,M^op,u. These rules specify all arcs; each graph has 2|M|+6 vertices.

Exactly the three arcs a-s, a-t and a-u reverse. All other arcs, including
the complete star of the named root r and its sole predecessor a, remain
fixed. Reflection of the seven blocks together with converse identifies
T with U^op, so H(T)=H(U). This uses the specified converse pairing of
the two modules; it does not assume M is self-converse.

Both tournaments are strong. Choose Hamilton paths through M and M^op.
Together with the singleton blocks, their concatenation is a spanning
path from r to v; a→r and v→a close a spanning directed cycle.

Use the rooted state (H,s_r,q_r,F_r), where s_r counts paths starting at r,
q_r=H(T−r), and F_r counts ordered two-slot covers whose first path is
nonempty and ends at r. Then, writing Z=12p²+8hp−h²,

    T: (Z, 6hp−h²,    6hp−h², Z),
    U: (Z, 4p²+4hp, (2p+h)², Z+h²).                   (1)

The preserved scalar is H; the repaired coordinate is q_r−s_r:

    0 -> h².                                         (2)

Thus the arbitrary-module construction is unit-ready exactly when h=1,
equivalently when M is transitive. It must not be described as a universal
unit-ready repair. For M=TT_j, j≥1, set B=2^(j+1); then (1) is exactly
the constant-three-arc family in [APEX_REARM §3](history/APEX_REARM.md),
with Z=3B²+4B−1. In particular M=TT1 gives the 63 repair.

## 2. A direct boundary proof of the formulas

A Hamilton path splits uniquely into its part before a and its part
after a. Each of these parts visits the ordinal blocks in increasing
order. Each module is therefore used in at most two runs. A module used
entirely before or entirely after contributes h; one used on both sides
contributes k. Allowing either side, including an empty side, contributes
p. Every singleton contributes 1 on its chosen side.

An assignment of the seven blocks to the two sides is legal precisely
when the last nonempty prefix block loses to a and the first nonempty
suffix block is beaten by a. Empty prefixes and suffixes are allowed.
These two boundary conditions are sufficient, since all intervening
inter-block arcs point in the required direction.

For clarity, there is an equivalent sum by the two boundary blocks.
If every vertex precedes or follows a, test just the last or first block.
A prefix ending immediately before a suffix starts gives an adjacent
L,W cut. Otherwise the first suffix block i lies strictly before the
last prefix block j, with colors W,L. A module at either endpoint
contributes h+k=p−h; an internal module contributes p; a module outside
the interval contributes h. An internal singleton contributes 2 and an
endpoint or external singleton contributes 1. There can be no i=j,
because that block cannot both beat and lose to a.

For T the colors are W,L,L,W,W,L,L. The two endpoint placements and
the adjacent L,W cut contribute 3h². The remaining W,L intervals give

    h² + 2h(p−h) + 4p² + 8p²
        + hp + 2hp + h(p−h) + 2h(p−h).

Their sum is Z=12p²+8hp−h². Converse-reflection gives the same H for U.
After deleting r, the colors for T are L,L,W,W,L,L. Endpoint and
adjacent-cut contributions total 2h²; its W,L intervals contribute
3hp+3h(p−h). Hence q_T=6hp−h².

For U−r the colors are W,L,L,W,W,L. Endpoint and adjacent-cut
contributions total 3h². Its W,L intervals contribute

    h(p−h) + hp + 4p² + 2h(p−h) + h²,

so q_U=4p²+4hp+h².

Since r beats every remaining vertex except a, the HPs of T−r that do
not extend to a path beginning at r are exactly those beginning at a.
Their number is zero in T−r, whose first carrier block s beats a, and
h² in U−r, where a beats s. This proves both s coordinates and (2).
The same concatenation argument as [APEX_REARM §2](history/APEX_REARM.md)
gives F=H+q−s for any root of indegree one. All four coordinates in (1)
follow.

The two-run coefficient matrix gives a second finite proof interface:

| graph and observer | h² | h k | k h | k² |
|---|---:|---:|---:|---:|
| T, H | 63 | 30 | 26 | 12 |
| T, s | 11 | 6 | 0 | 0 |
| T, q | 11 | 6 | 0 | 0 |
| U, H | 63 | 26 | 30 | 12 |
| U, s | 24 | 10 | 10 | 4 |
| U, q | 25 | 10 | 10 | 4 |

Here the middle columns distinguish which of the two modules is split.
The companion script independently enumerates quotient run words and
prefix/suffix assignments, compares these exact integer coefficients,
and substitutes k=p−2h symbolically. No fitted polynomial is used.

## 3. The exact quadratic norm and the primitive exit

The inherited N/R root-child operations add x with root→x→every other
old vertex; N retains x and R retains the old root. Their full state maps
are stated in [primitive §3](primitive/PROOF.md) and proved in
[ROOTED_RESETS R1](../two_port_return_induction_20261103/ROOTED_RESETS.md).
For a root of indegree one, with F=H+δ and q−s=δ, their composition NR
gives the balanced state

    (3H+2δ, H, H, 3H+2δ).

Apply it to U. The output is therefore

    (Y,Z,Z,Y),
    Y = 3Z+2h² = (6p+2h)²−5h².                       (3)

In particular this is the explicit product of conjugates

    Y = (6p+2h+h sqrt(5))(6p+2h−h sqrt(5)).

The connection with quadratic norms is an exact algebraic identity for
an actual graph constructor. Its input remains the physical pair (h,p)
of a tournament. A representation Y=x²−5y² does not supply this graph
unless y=H(M) and (x−2y)/6=P_M(2) for a suitable M. In particular neither
Pell recurrences nor units of the quadratic ring automatically preserve
the physical input set.

The balanced output is primitive exactly when

    gcd(h,p)=1 and 3 does not divide h.                 (4)

Indeed its content is gcd(Y,Z)=gcd(h²,Z), since Z is odd. A common prime
of h and p divides both h² and Z; and 3|h implies 3|Z. Conversely, if
both conditions in (4) hold, every prime ell dividing h satisfies
Z≡12p²!=0 modulo ell. Therefore gcd(h²,Z)=1. This proves both directions.

The distinction between unit readiness and primitive continuation is
productive: h>1 loses the former but can retain the latter. For the strong
four-vertex module Q=(2,4,9,3), h=5,p=32, so Z=13543 and Y=40679; the
defect is 25 and the final balanced content is 1. For M=C3, h=3,p=12,
Z=2007 and Y=6039, the final content is 9. These are positive and hostile
controls for (4), respectively.

More generally take r ordinal copies of Q followed by TT_a, with
4r+a≥1. Slotwise multiplication gives h=5^r and p=2^(5r+a). Every such
module satisfies (4), hence gives a primitive balanced norm output.
The module size is unbounded and can be nontransitive. This two-parameter
family is still sparse; the claim is a source of actual primitive states,
not coverage of every prescribed odd integer.

There are exact congruence limits. Since p is even and h is odd,
Z≡7 (mod 8) and Y≡7 (mod 8). Under the primitive guard 3 does not
divide h, equation (3) also gives Y≡2 (mod 3), hence

    primitive Y ≡ 23 (mod 24).                         (5)

Furthermore Y is a square modulo 5. In particular this chart cannot
cover any target congruent to 2 or 3 modulo 5. These are actual image
restrictions, not failures of a bounded search. The quadratic norm is
useful as one chart in a larger construction, not as a universal formula.

Nevertheless exact membership in this chart is a finite source problem.
Cutting a Hamilton path at each of its |M|+1 positions gives distinct
ordered two-slot covers, so p≥(|M|+1)h. In particular Z≥63h² and
Y≥191h². For a prescribed target Y, enumerate the positive odd integers
h≤sqrt(Y/191). A candidate must have an integer

    x=sqrt(Y+5h²),  x≡2h (mod 12),
    p=(x−2h)/6≥2h,

and an actual module of order at most floor(p/h)−1 with response (h,p).
Conversely any such module supplies the target by the literal construction.
For primitive outputs also impose (4). The arithmetic candidate test and
the bounded physical-module question are both necessary; dropping the
latter changes the problem. This is a finite inverse description of
this chart, not a claim that every target passes it.

## 4. The question this changes

The original same-H move at 63 comes from reflection, a familiar global
symmetry. Its useful coordinate is local: the same root star survives
while a deletion count changes by one. Arbitrary module inflation exposes
exactly what is lost: that one becomes h². The next physical operation
turns this loss into a quadratic norm and, under (4), a primitive balanced
source. Thus the fixed-source problem can be replaced by a question about
the realizable two-slot response set

    R = {(H(M),P_M(2)): M a nonempty tournament}.

Two concrete remaining questions are **OPEN**:

1. Can a target equation such as (3), or the independent linear wrapper
   Y=117p+29h, select a physical point of R for every remaining prime or
   seven-prime target after a bounded collection of additional wrappers?
2. Can a constructor move between these source charts while preserving
   a primitive marked interface and decreasing an explicit target height?

The exact feedback law in [APEX_REARM §6](history/APEX_REARM.md) says why
one chart alone cannot be assumed closed: NR changes the feedback rank
to tau(T−r)+1. The deleted-atom models in
[semigroup §5](semigroup/PROOF.md) say why localization or square/cube
attainment cannot replace the physical response-set question.

## 5. Reproduction and scope

Run from the repository root:

    python3 -B 05-knowledge/results/sixty_three_bridge_20261113/mirror_module.py

The script checks both exact symbolic compilers and six named literal
modules through order 14, including transitive, cyclic and strong
four-vertex inputs. A subset DP independently counts all four rooted
coordinates; the three changed arcs and the unchanged root star are
checked on the actual adjacency rows. The NR conclusion uses the proved
state map, not a purported large-order census. The receipt records the
universe and every literal pair. All-parameter conclusions follow from
the boundary proof above.
