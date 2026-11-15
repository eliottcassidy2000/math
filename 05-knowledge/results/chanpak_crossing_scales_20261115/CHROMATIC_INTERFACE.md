# The exact bounded-checksum interface, and two persistent unit-step switches

**Status:** the finite compiler and its asymptotic exact interpretation for the gap-separated grammar are all-order results. The palette census for checksum at most 99 is finite-exact. Two persistent +2 switches follow. The exclusion certificates are for a specified construction and checksum budget, not for the tournament spectrum.

## Why keep the full interface?

OpenAI's [Euclidean no-five-coloring proof and Lean statement](https://github.com/openai/math/blob/main/lean/docs/158.md) cover arbitrary colorings, without assuming regularity of the input color classes; the result leaves six or seven colors. The related but distinct [linear List Hadwiger result](https://github.com/openai/math/blob/main/lean/docs/157.md) bounds list chromatic number by a universal constant times clique-minor size even though the ordinary Hadwiger conjecture fails. In [Shukla's November preprint](https://arxiv.org/abs/2611.06269), three exact finite unit-distance witnesses improve the size records **within specified small disks** for five chromaticity, with coordinates, 5-colorings and DRAT obstructions to four colors. These scoped finite certificates do not decide six versus seven colors.

The repo's [six-color interface profile](../six_colour_interface_profile_20261015f.md) already distinguishes the full boundary relation from its marginal projections and notes that cross-color compatibility needs the actual labeled support. We use the same design principle here: the two sides of a glue must report the same **normalized absolute leading coefficient** and satisfy the low checksum, not merely hit the high and low projections separately. This is a list-compatibility analogy and an exact arithmetic compiler, not a reduction to unit-distance graph coloring.

## The exact two-coordinate signature

Use the [two-core anchor glue](../chanpak_quadratic_anchor_atlas_20261115/README.md) and [binary anchor decoder](../chanpak_heisenberg_target_20261115/ANCHOR_DECODER.md). The anchor has a distinguished apex over a transitive spine; its word is 1 0^g v with g≥1, and v has length k≥1. The bit 1 means that the anchor apex beats that spine vertex. The rest of the glue has R→W→U, U→p and R−p→U, with p the anchor apex. It is strong, since the anchor apex reaches the first spine vertex, every vertex of R reaches W, and W returns through U.

Write u=h(U), h=h(W), x=Q(U)−u, q=Q(W), d=xq, t=uh, and G=x(q−h)+2t. Let b be the value of v and r be the value of its k-bit reversal. Denote by A(v) the path count of the apex tournament with this suffix alone. Then the exact path count is
```
H = 2^g L + c
L = d·2^k + t(2^k−r)
c = G + db + t(A(v)−1).
```
Indeed B(anchor)=1+2^(g+k)+b and H(anchor)=2^g(2^k−r)+A(v). The second identity follows by counting the odd cycles using the first 1 or from the seam identity; there is no seam adjacency correction since g≥1. Inserting both counts into the two-core formula gives the stated signature.

The actual leading color is λ=oddpart(L), and its clock is j=g+v₂(L). Two such signatures have the same set of allowable sufficiently large clocks when they have the same λ. Equality of raw L is unnecessarily restrictive; the 59 example below uses coefficients 118 and 59. Padding v on the **left** by zeros leaves b, A(v), c unchanged and multiplies L by a power of two: the bit reversal shifts to the left by exactly the padding width. Thus retain only the ordinary binary expansion of b (take v=0 when b=0). No information used by the low checksum or the normalized color has been lost.

## The finite-list theorem at all heights

For an odd low checksum c let Λ(c) be the set of normalized colors of all physical signatures with that c. The following finite enumeration computes it exactly:

1. Enumerate U,W with positive orders a,b' and total order at most 1+floor(log₂(c−1)); for c≤1 the set is empty.
2. For each pair, enumerate suffix values 0≤b≤floor((c−G)/d) when G≤c. Use the shortest representative as above. Compute A(v), c and λ; retain exact matches.

To prove step 1, note Q(X)=sum over vertex subsets S of h(X[S])h(X[V−S]), with empty path count one. The terms at the two ends contribute 2h(X); all terms are positive by Rédei. Hence Q(X)−h(X)≥2^|X|−1. At total order a+b', the product (2^a−1)(2^b'−1) is minimized at a=1 or b'=1; for example this follows from (2^(a−1)−1)(2^(b'−1)−1)≥0. Therefore
```
c≥G≥(2^a−1)(2^b'−1)+2≥2^(a+b'−1)+1.
```
Step 2 follows by A(v)≥1. Conversely every retained signature constructs physical tournaments for every g≥1, with positive order on both sides of each port.

Here is an explicit uniform interpretation at separate heights. For 1≤c≤Y, set M=ceil(log₂(Y+1)) and J=M+ceil(log₂(6Y+1))+1. At every j≥J, a target N=λ2^j+c (λ positive odd) has a witness within this grammar with **some** checksum at most Y if and only if λ belongs to Λ(c). For this, d≤2G≤2Y, t≤d/2. If b>0, 2^k≤2b and L≤3db≤3Y; if b=0 then L=2(d+t)≤6Y and λ=d+t≤3Y. A witness at low checksum ≤Y has normalized color ≤3Y, so its clock cannot be ≤M when N≥2^J. Reducing the two expressions modulo 2^M shows equality of the low checksums; then uniqueness of factoring out twos shows equality of normalized colors and clocks. Conversely use the shortest representative for a listed color and set g=j−v₂(L); by L≤6Y this g is positive and 2^g>Y, so the low part has no carry. This proves the certificate is complete for that budget, not an extrapolation from a finite range of target counts.

Consequently exact adjacent odd counts at any sufficiently high common clock use the interface intersection Λ(c)∩Λ(c+2). A list of high colors and a separate list of obtainable low checksums cannot certify adjacency; the joint interface is necessary.

## Two literal +2 switches forever

Here TT_m denotes the transitive tournament of order m, C3 a directed 3-cycle, and S4 the strong tournament with bit rows (8,1,3,6). Its profile (h,Q) is (5,32). Use the signatures below, with column k determined by the displayed suffix.

| U | W | suffix v | L | c |
| --- | --- | --- | ---: | ---: |
| S4 | TT1 | 0 | 118 | 37 |
| TT3 | TT1 | 10 | 59 | 39 |
| TT3 | TT1 | 110 | 117 | 97 |
| C3 | C3 | 0 | 234 | 99 |

For every j≥2, the first pair realizes both 59·2^j+37 and 59·2^j+39 (choose g=j−1 and g=j, respectively). The second pair similarly realizes 117·2^j+97 and 117·2^j+99 (choose g=j and g=j−1). For j≥8 the low checksum of every displayed arm is strictly less than its clock modulus. These locks genuinely place a unit step on the odd number line; a +32 and a −30 on different windows cannot be subtracted to assert this. None of the four underlying readouts is just an arbitrary algebraic coefficient: the cores and words are explicit.

For Y=99 the theorem leaves each core of order at most six and the two core orders sum to at most seven. A literal subset-DP census of every labelled tournament through six vertices, followed by the exact compiler, gives 34 nonempty checksum lists and 95 distinct incidences (c,λ). The only nonempty adjacent-checksum intersections are exactly Λ(37)∩Λ(39)={59} and Λ(97)∩Λ(99)={117}. In particular there is no three-entry consecutive-odd run in this bounded-checksum grammar at a common sufficiently high clock. This does **not** rule out larger checksum budgets, other anchors, other glues, or unbounded-low representations of the same integer.

**Replay:** `python replay_chromatic_interface.py` uses literal DP for core and eight small whole-glue witnesses, then the exact finite list compiler. The outside results motivate retaining labeled compatibility and explicit forbidden regimes. They are not used as assumptions in the arithmetic proof, and no statement about a new chromatic bound is made.
