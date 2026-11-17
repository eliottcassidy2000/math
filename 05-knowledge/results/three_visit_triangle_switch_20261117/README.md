# The triangle switch: an exact +2 with three visits and nonuniform continuation

**Proved.** One edge reversal produces an exact gap of twice the Hamiltonian count of a prescribed quotient. When the quotient is transitive the gap is exactly two, with no length bound and independent of how individual exterior vertices contact the two ports. There is also a literal three-count window with a block that can be visited three times. This does **not** prove Chan–Pak.

## 1. A two-port triangle switch

Let `Q` be any tournament with distinguished vertex `d`. Replace `d` by a triangle `b→c→y→b`, giving every exterior vertex of `X=Q−d` the same orientation to all three triangle vertices that it had to `d`. Add ports `x,a`, with
```
x→b, a→b, c→x, a→c, x→y, y→a.
```
For each individual exterior vertex choose, arbitrarily and independently, whether it dominates **both** ports or loses to **both** ports. Let minus have `x→a` and plus have `a→x`; the rest is unchanged. Then
```
H(plus) − H(minus) = 2 H(Q).
```

**Proof by contraction.** Paths not using the flipped arc are identical. In every remaining path contract the consecutive ports to one letter `q`. Its incoming arcs read the first port of the original pair, and outgoing arcs read the last. This is an exact bijection to Hamiltonian paths in the associated directed graph (which may have double or missing arcs). The two contractions agree everywhere except on the triangle. In minus, `q` can be entered from `c` and can exit to `b,c`; in plus it can be entered from `y` and can exit to `b,y`.

* If neither immediate neighbor of q belongs to the triangle, use the identity. If exactly one belongs, choose a cyclic rotation mapping the required singleton incoming set or two-element outgoing set onto the analogous set in the other graph. This rotation preserves all other constraints, because the triangle is a module relative to X.
* If both neighbors belong to the triangle, the only possible local strings are `c q b` in minus and `y q b` in plus. When the third triangle vertex touches neither end of this string, swap the names c and y; no triangle–triangle adjacency is used, so this gives a bijection. Otherwise, the third triangle vertex is immediately before or after the string. Neither placement works in minus (`y→c` and `b→y` are false); both work in plus (`c→y` and `b→c` are true).
* Each of the two extra four-letter strings in plus starts and ends in a triangle vertex, with the same exterior contacts as d. Compress the whole string to d. Each string therefore contributes exactly H(Q).

In particular, choose a transitive Q and add exterior vertices one at a time, preserving its order, with any sequence of shared port contacts. The same flip remains an exact +2 under *every* such continuation. Indeed arbitrary modifications of the exterior or the position of d are also allowed if Q stays transitive. The two ports must see each exterior vertex in the same direction, but the direction may vary across vertices. This is not the one-fixed-root operation, and the proof does not require that either port's individual continuation profile be invariant.

The calculation is a finite-difference rule under composition. For an ordinal quotient `M→d→V` it reads exactly `2H(M)H(V)`, allowing the scale to change with independently substituted M and V. If V is transitive, *any* external contact word for V is allowed, not just a uniform run. This gives a physical way to carry a unit through a nonuniform seam, rather than bounding a scalar error at one isolated size.

## 2. An explicit three-visit window, and what each contact changes

Use a transitive spine with apices x→y. A contact digit has bit 0 set when x beats that spine vertex, and bit 1 set when y beats it. Compare words
```
2 1^m 0 1        3 1^m 3 0        2 1^m 3 0 .
```
Let B occupy the repeated slot instead of a transitive block, with all external orientations unchanged. Let W_j(B) count ordered decompositions of B into j nonempty directed paths, with W_1(B)=H(B). The three Hamilton counts are exactly
```
27 W_1 + 21 W_2 + 4 W_3
29 W_1 + 21 W_2 + 4 W_3
31 W_1 + 21 W_2 + 4 W_3.
```
For transitive B of size m≥1, its ordered j covers number j! S(m,j), the surjection count. The window is consequently
```
4·3^m + 9·2^m − 3,  −1,  +1         (same center)
```
where only the last two variants differ by a single edge. For a direct three-visit witness at m=3, number the six spine vertices 0,...,5 and use x=6,y=7. In both last variants,
`1,4,5,6,2,7,0,3` is Hamiltonian, and the block {1,2,3} is visited in the three singleton runs. In the switch proof these last two variants have Q=M→d with M the repeated slot, and both ports beat all M. A new transitive exterior V with contacts as in §1 preserves the same two-endpoint difference, in particular when m and |V| grow at unrelated rates.

**Deriving the coefficients.** A path is split by x,y into three increasing spine runs (possibly empty). Fix apex order and assign the three fixed spine positions to these runs. For each nonempty subset of runs receiving vertices of B, checking the resulting path uses only the fixed contacts and that subset. Assigning paths of a W_j cover to those j visits lifts the compressed trace bijectively. There are at most three visits. Thus the coefficient vector is determined exactly by the three transitive tests at sizes 1,2,3, whose profile matrix has rows (1,0,0), (1,2,0), (1,6,6). The counts are respectively (27,69,177), (29,71,179) and (31,73,181). The supplied verifier also compares with independent full vertex-permutation or subset counts, including nontransitive B. Neither a recurrence fit nor extrapolation is needed.

### A growing continuation needs more than three reads

Fix the middle gadget G (word `3 1^m 3 0` before replacing B by its contents). For an ordered cover with k paths there are at most k+2 visits to B: leaving B and returning in the same path uses at least one apex, since without either apex the fixed spine is ordered around B. The transfer W_k(G[B]) is an exact nonnegative integer linear combination of W_j(B) for j≤k+2. The top coefficient is positive for every k: take the exhibited three-visit compressed trace as one cover member and k−1 B-only singleton-slot traces as the rest. A lift requires sufficient block size, already encoded by W_{k+2}(B).

Consequently if G is substituted d times, then m↦H(G^d[TT_m]) has leading exponential base 2d+1, with strictly positive coefficient. Those reads, considered as functions of positive m, are linearly independent by taking the largest leading base. **No fixed finite-dimensional linear profile** on these transitive inputs can recover every iterated read by linear functionals. This is not an obstruction to nonlinear encodings, or to truncation for a fixed number of wraps. It explains why a three-coordinate exact reader for one seam is not automatically a fixed-size transfer system for arbitrary nesting.

## 3. Exact scope and a target-facing question

For homogeneous V of size k where every vertex of V beats both ports, the two-module read has a positive 3^m 3^k term. The coefficient array on i^m j^k (i,j=1,2,3, m,k≥1) for minus is
```
  71   -80    22
-116   187   -76
  42  -104    66
```
and plus changes only 71→73. This homogeneous subfamily is sparse (O(log² X) distinct counts up to X): for m,k≥3 the positive 66 W_3(TT_m)W_3(TT_k) trace contribution bounds m+k by O(log X); when either width is 1 or 2 the remaining variable has a strictly positive leading coefficient at base 3. For k=0 the §2 formula applies. The bound is **not** asserted for arbitrary contact words. In particular these homogeneous rails do not yield all odd fibres or a formula for the Mandelbrot target (2^p+1)²+2.

### Unequal growth really does help placement, but not exact placement

One nonuniform continuation is r ordinary sinks followed by one new vertex beating both ports and beaten by every other old vertex. For m,r≥1 the two counts are

    L(m,r) = (28·3^m + 29·2^m − 15)·2^r + (4·3^m + 3·2^m − 5)
    L(m,r) + 2.

The compressed-trace coefficient array, on i^m j^r (i=1,2,3 and j=1,2), has rows (−5,−15), (3,29), (4,28); the + branch changes just −5 to −3. The surjection expansion of §2 gives the formula.

For each fixed ε>0 and every sufficiently large real X, there are such physical *pairs* just below and just above X within εX. Write A_m=28·3^m+29·2^m−15. The phases log₂ A_m mod 1 are dense on every tail of m, since log₂ A_m−m log₂3 approaches log₂28 and log₂3 is irrational. Choose a fixed finite bank giving a fine net in both directions and choose r so A_m 2^r is inside the required multiplicative interval with some margin. The additive terms and +2 are bounded on that bank, so for sufficiently large X both endpoints fit. This is relative placement, not bounded additive error.

In fact neither branch ever attains Q(p)=(2^p+1)²+2 for integer p≥1. For r≥3, m=1 gives 5 or 7 mod 8, and m≥3 gives 7 or 1 mod 8. For m=2, only minus can have phase 3 mod 8: its read is 353·2^r+43, and for p≥3 the 2-adic valuation of Q(p)−43 is exactly 3, forcing r=3. That read is 2867, between Q(5) and Q(6). The small targets Q(1),Q(2)=11,27 are below both branches. For r=1 or 2 with m≥6 the counts exceed Q(4); for p≥5, Q(p)≡3 mod 64. When r=1 only plus survives mod 8 and its phase mod 64 would require 3^m≡7 mod 16; when r=2 only minus survives and would require 3^m≡5 mod 16. Neither is a power of 3. Finally for m=1,...,5 the minus reads at r=1,2 are (267,521), (749,1455), (2073,4019), (5801,11235), (16497,31931); none of these or their +2 partners equals a Q(p). Thus stable repair, two independent exponential scales and relative placement can still fail absolute placement.

The new physical inverse problem is to place the **center** of the switch pair. For a transitive quotient the center is H(minus)+1, and a prescribed odd P is hit precisely when some quotient/order/contact word places that center at P−1 or P+1. This is stronger than placing the residue or producing any pair at distance two, and is not proved here. Changing the address pattern or the quotient is a stronger next test than increasing both widths.

For the original conjecture and attribution, see [Chan–Pak, §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf). The exclusions at 7 and 21 are not new claims here. See also the [homogeneous-width obstruction](../homogeneous_mandelbrot_width_20261117/README.md) and the [curried birth diagrams](../curried_births_tooth_diagram_20261117/README.md).

## Reproduction

Run `python replay.py` and `python two_slot.py` from this directory. The former checks the spine counter against literal permutations, and the arbitrary-quotient formula and substitution responses against independent subset path counts. The latter counts two-slot compressed traces and independently replays both arrays on literal tournaments.
