# 1024 connects the 2024 graph, but cannot give a Hamilton path

Scale (8,15,17) by 3, giving the sum alphabet 576,2025,2601
on 1…2024. As in the [reflection orbit theorem](../inverse_fences_landing_20261116/FENCE_AND_SQUARE.md),
the two reflections modulo 2025 are x↦−x and x↦576−x; their
product is translation by 576, with gcd(2025,576)=9. The five
components are the ±r residue classes modulo 9 for r=0,…,4:
one 224-vertex path with ends 288 and 576, and four 450-cycles.
The largest label is 2024; 1024 is not an extremal label or a
component size.

The user's 1024=10³+3·2³=32² nevertheless has a physical role.
Adding the edges of sum 1024 connects the five components.
Their quotient (ignoring internal edges) is 0–2–4–3–1, as
1024≡7 mod 9. Example edges, one for each seam, are
(9,1015),(2,1022),(4,1020),(6,1018).

**But all these added edges together still cannot give a Hamilton
path.** Suppose one exists. Color old edges not used red and new
edges used blue. The blue edges form a matching. Pair red and
blue incidences at each vertex as far as possible. This gives
alternating trails, with unpaired blue incidences only possible
at old endpoints 288 and 576, and unpaired red incidences only
at new endpoints. Total red and blue edges are equal, since both
the original disjoint graph and a Hamilton path have 2023 edges.

An alternating cycle would have red sums 576 or 2025 (its vertices
all ≤1023 because they meet a blue edge). In an alternating cycle,
one application of a blue reflection followed by a red reflection
translates by −448 or +1001. If the red counts are a,b then
448a=1001b, forcing b≥64 if the cycle is nonempty. There are only
11 old 2025 edges with both endpoints at most 1023. So no alternating
cycle exists. An unpaired blue end at 288 goes to 736 then can
only continue with red to 1289 or 1865, out of blue range.
From 576 the blue goes to 448, then red to 1577 (out of blue
range) or 128; in the latter case blue then goes to 896, whose
old neighbors 1129 and 1705 are again out of blue range.
Thus there is no blue-to-blue trail, so count balance precludes
red-to-red trails too. Every blue edge lies on one of the
described trails, at most three in total. That cannot join five
old components.

## A smaller square is a better repair

The least new square sum that allows a Hamilton path with
exactly four cross-component new edges and no other new edges
is 529=23². For s below 529, the endpoint 576 cannot move,
so start there, exit the old path at 288, and visit each
cycle once (the four bridges contract to a tree and the
Hamilton traversal visits its endpoints successively).
The old cut inside every intermediate cycle has both ends
below s, hence can only have sum 576. Arrivals are
y_j=s−288+(j−1)(s−576), j=1,…,4, so 4s−2016>0 is necessary.
The next square is 529. The four actual bridges are
(288,241),(335,194),(382,147),(429,100); the internal cycle
cuts are (241,335),(194,382),(147,429). At 100 one may
leave out either incident cycle edge to obtain the final path.
This minimality has the stated four-edge scope.

## Where 24 and Ramanujan actually enter

2024=C(24,3). In the binary Golay code on 24 coordinates, 759
octads have 56 triples each, and every triple is in 21 octads
(by the Steiner property); this gives the same 2024. Prescribing
two coordinate bits leaves 2^(12−2)=1024 codewords, since the
dual has no weight-one or weight-two word. These are independent
Leech/Golay interpretations, not an identification of the
five graph pieces with code orbits.

There is an exact arithmetic resonance in Ramanujan's tau:
τ(2)=−24, 576=τ(2)², and 2024=2¹¹+τ(2). A two-bridge
interchange using sum s and two old cut sums A,B requires
A+B=2s, since the four endpoints are the same on both sides.
At s=1024 and A=576, the required partner is
2¹¹−τ(2)²=−τ(4)=1472, by the prime-power recurrence. It is
not an available old sum. Using 2025 instead requires the
nonsquare 23. The direct repair uses **23²**, not the missing
color 23. This is an exact missing-partner calculation, not a
modular-form construction of the Hamilton path.
For the repo's actual Leech–Ramanujan theorem (a theta-series
identity, not just the repeated numeral), see
[THM-4758](../../../01-canon/theorems/THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md).

## The small-scale success has an exact cutoff in this template

Scale by any odd k≥3 instead, with A=64k², B=225k²,
d=k² and r=(d−1)/2. The same reflections give one path of length
224 in vertices (ends A/2 and A) and r cycles of 450 vertices.
For an added square color s with 2s<B, a minimum-budget repair
using exactly r new edges exists if and only if gcd(s,d)=1 and

$$
A-\frac{A}{2r}<s<A+\frac{A}{r-1}.
$$

At that budget the contracted walk is a single pass through the
old pieces. Every intervening old cut has both ends below s,
so its sum must be A. The first exit x is either A/2 or A;
successive arrivals are y_j=j(s−A)+A−x. Below A only the first
exit is possible, and y_r>0 is the lower bound above. Above A
the best exit is A, and y_(r−1)<A is needed for the last internal
cut, giving the upper bound above. These bounds also suffice for
all the intermediate edges by the explicit arrival rule. Residues
are js mod d, and they visit all ± residue classes precisely
when s is a unit modulo d (if not, a nonzero common divisor
blocks any escape to those other classes). Use the old sum-A
cut at each intermediate component and either neighboring cut
on the final cycle.

For odd k≥9 even the neighboring squares (8k−1)² and (8k+1)²
are outside this interval: their distances 16k−1 and 16k+1
exceed the respective bounds 64k²/(k²−1) and 128k²/(k²−3).
Thus the minimum-bridge scheme genuinely stops in this restricted
low-color range; it says nothing against larger square colors
or multi-pass repairs.
