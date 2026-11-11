# Completing the entire 35-grid with offset tails and a product replacement

**PROVED:** every number `N(35,r,s)=1+2^s+35·2^(r+s)`, r,s≥1, is the
Hamilton-path count of an explicit tournament. The new constructions close
all four formerly unresolved levels s=2,3,4,5; levels s=1 and s≥6 are
inherited. **FINITE-EXACT:** the literal controls and specified exploration
menus below. **OPEN:** the other residual grids in general, and Chan–Pak.

## 1. The inherited interface and the changed operation

The [seed-recursion packet](../../seed_recursion_births_20261111/PROOFS.md),
R10, constructs

    N(35,r,s)=1+(35·2^r+1)2^s, r≥1,

at s=1 and every s≥6. It leaves the four uniform levels s=2,3,4,5 open.
The same packet's R3 completes the entire d=3 grid, including the coefficient
7 ray `1+7·2^t` for every t≥1.

The nearest failed model is precise: a fixed one-port transitive-inflation
kernel starting at H285 cannot produce the full s=2 line in R10's exact
inventory decision. That exclusion concerns a uniform graph architecture,
not the values in the line. An actual serial composition can change the
architecture and preserve the target.

The live concepts are actual parent graphs / ordinal product / finite heads
versus repeatable tails / endpoint capacities / balanced binary shells.
The useful underused coordinate is the factorization of the **whole target**,
including its constant term, rather than only its odd slope.

## 2. The three offset seeds close levels 3, 4 and 5

A selected-set seed need not have offset one. For any actual tournament B
and selected set O, define H,U,W as in the next section. The inherited
nonempty-tail construction has the exact count

    H(T_t(B,O)) = D+2^t U,     D=H(B)-W_B(O), t≥1.    (1a)

Its graph is B→TT_t→z, with z→O and B-O→z. Retaining D instead of
requiring D=1 gives the following three small literal sources.

| s | B outgoing rows | O | (H,U,W) | D | tail length t |
|---|---|---|---|---|---|
| 3 | (26,44,57,48,34,1) | {1,2} | (31,140,22) | 9 | r+1 |
| 4 | (62,124,120,112,32,64,17) | {4,6} | (27,70,10) | 17 | r+3 |
| 5 | (30,124,24,112,96,69,5) | {2,3} | (95,280,62) | 33 | r+2 |

All three bases are strong, and O is nonempty. The resulting tail graph
is strong: from its last vertex one reaches O and then all of B, while
B reaches the tail and then that last vertex. Every r≥1 gives a nonempty
tail, so there is no absent-module exception. Substituting the displayed
statistics into (1a) gives, respectively,

    9+140·2^(r+1) = 1+2^3+35·2^(r+3),
    17+70·2^(r+3) = 1+2^4+35·2^(r+4),
    33+280·2^(r+2)= 1+2^5+35·2^(r+5).               (1b)

This proves all three uniform levels. The finite source statistics are
exact counts of paths and two-path covers on six or seven vertices. The
all-r conclusion is the already proved tail identity, not numerical fitting.

The useful change is structural: these targets have offsets 9,17,33.
Searching only unit-offset seeds discards the three sources above before
their tail operation is even used. No physical halving or rearming is
needed at these levels. The offset and selected set remain part of the state.

## 3. The product replacement theorem for level 2

For every r≥1,

    N(35,r,2)=5+140·2^r=5(1+7·2^(r+2)).              (1)

If A and B are tournaments and every cross arc points from A to B, every
Hamilton path of their ordinal join traverses A first and B second. Restriction
and concatenation are inverse bijections, so

    H(A→B)=H(A)H(B).                                 (2)

The four-vertex tournament with outgoing rows

    A=(2,4,9,3)

has H(A)=5. Apply (2) to A and the inherited coefficient-seven constructor
at exponent r+2. This proves every value in (1), with no hypothetical seed
halving or same-count rearming.

Here is a completely explicit version, with no constructor lookup needed.
At r=1 use the H57 graph

    (2,64,67,39,15,23,57),

and at r=2 use the H113 graph

    (34,48,11,67,13,28,55).

These literal heads come from the completed d=3 grid's `three_s3` and
`three_s4` prefix records. For r≥3, inflate vertex4 of

    (2,40,3,5,15,29)

to a transitive two-vertex module, and keep vertices2,3 selected. The resulting
seven-vertex marked base B has

    H(B)=33,    U_B(O)=112,    W_B(O)=32.             (3)

Here U counts ordered spanning two-path covers with second path nonempty
and starting in O, the first path allowed empty. W requires both paths
nonempty, the first ending in O and the second starting in O.

Append a transitive tail of length r−2 after B, then a vertex z. Direct
B→tail→z; let z beat precisely O in B and let the other B vertices beat z.
The inherited exact selected-set tail identity gives

    H=H(B)−W_B(O)+2^(r−2)U_B(O)
     =1+112·2^(r−2)=1+7·2^(r+2).                    (4)

The tail is nonempty because r≥3. Joining A to this graph gives (1).
The seed in (3) is the size-two instance of R3's certified `three` family;
its three finite statistics can also be counted directly on seven vertices.

The emitted tournament has eleven vertices at r=1,2 and r+10 vertices
for r≥3. It is an ordinal product, so strong connectivity is not claimed.
The Hamilton-count problem allows this graph. Every vertex and arc is
specified by the construction.

## 4. Why the fixed-kernel exclusion is respected

This realization has two finite heads and a different repeatable tail.
It does not claim one fixed marked H285 graph whose single port can be
inflated to produce every subsequent member. The repeatable template starts
at r=3 with H1125. Its lower two values are supplied separately.

This is exactly the distinction between a finite constructive family and
an unnecessarily fixed physical genealogy. The new saturation theorem for
strict integer addition does not automatically transfer to tournaments:
its binary operation still needs an actual graph implementation. Here
ordinal multiplication already has such an implementation, so the arithmetic
factorization transports without missing port conditions.

A general safe transfer is immediate. If `2^s+1` divides d, then

    1+2^s+d·2^(r+s)
      =(2^s+1) [1+(d/(2^s+1))·2^(r+s)].              (5)

Whenever the factor count `2^s+1` and the full required inner ray have
actual tournaments, ordinal join fills that entire residual-grid level.
This is a constructor lemma with two explicit realization hypotheses,
not an assertion that arbitrary scalar factors are already available.
For d35,s2 those factors are 5 and the proved coefficient-seven ray.

## 5. Complete constructor and reproduction

`realize35.py` emits a literal tournament for every positive r,s. Its
construction choices and vertex orders are:

| s | construction | order |
|---|---|---|
| 1 | inherited H=3+70·2^r one-port kernel | r+8 |
| 2 | ordinal product from section 3 | 11 for r=1,2; r+10 thereafter |
| 3 | offset-9 tail | r+8 |
| 4 | offset-17 tail | r+11 |
| 5 | offset-33 tail | r+10 |
| ≥6 | inherited unit seed U=35·2^(r+5)+32 and tail s−5 | r+s+7 |

For the last row, the base is `(4,1,18,7,75,31,47)`, vertex5 is inflated
to TT_(r+5), and vertices0,3 are selected. These are exactly incoming R10
rows and marks. The first row uses its quoted nine-vertex kernel with port4.
The level-2 construction need not be strong; all the other listed outputs
are strong. Strong connectivity is not required for the count theorem.

Run the complete checker or emit any individual instance:

    python3 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/realize35.py
    python3 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/realize35.py --r 5 --s 4

The checker derives (H,U,W) directly for all three new offset seeds and
counts the full Hamilton paths in all 21 outputs with r=1,2,3 and s=1,...,7,
through order17. Literal rows and strongness are retained in
`realize35_receipt.json`. The independent level-2 constructor is also retained.

## 6. Product controls and the nearest tested birth lane

Run

    python3 05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/realize35_s2.py

The checker counts the finite seed statistics and literal whole-graph paths
for r=1,...,6, giving285,565,1125,2245,4485,8965. Its receipt includes every
literal row. The all-r statement follows from (2) and (4), not extrapolation.
To emit any member, pass `--r` followed by a positive integer.

Before the product replacement, `head_lift_probe.py` tested every newborn
cut of the existing s1 source: inflate vertex4 of
`(8,5,9,322,15,31,55,127,247)` to TT_r, whose count is `3+70·2^r`, and
ask for a single newborn taking it to twice its count minus one. For r=1,...,7
there are no such cuts. This is only a finite observation about those exact
physical sources. The positive replacement above explains why that local
failure need not obstruct the target family.

## 7. An actual operation transporting the offset

The subsequent [offset-transport theorem](OFFSET_TRANSPORT.md) supplies
marked carriers for3,21,35 whose unselected transitive-module growth realizes
`(D,U)→(2D−1,2U)` at every positive length. Together with tail growth these
are two commuting physical graph operations. Their all-module extension
constructs `H(M)+(a+b·2^t)Q_2(M)` for three explicit pairs(a,b), preserving
the path-cover sidecar rather than claiming independent scalar addition.
This adds a reusable construction interface; the complete35-grid coverage
proved above remains unchanged.
