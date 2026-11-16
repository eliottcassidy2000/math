# Keeping the absolute target at a real boundary

## New exact handoffs

At height $s$ the old +4 child asks for the three counts
$L_s,L_s+2,L_s+4$, with $L_s=756\cdot2^s-33$. Changing the deletion layer
reopens three more heights:

| Height | Exact three counts | One shared $(E,O,g)$ |
|---:|---|---|
| 6 | 48351, 48353, 48355 | (46567,46576,9) |
| 7 | 96735, 96737, 96739 | (73697,73734,37) |
| 8 | 193503, 193505, 193507 | (183047,183056,9) |

The displayed witnesses, and a second height-seven continuation (77999,78036,37), are strong, with outgoing row masks in
[WITNESSES.txt](WITNESSES.txt), last vertex the root and vertex 0 the low bit.
They meet the prescribed absolute count. They are re-encodings, not induced
extensions of a chosen original graph. None of these displayed
continuations is three times an integral continuation: the height-six
and height-eight even coordinates are nonmultiples of three, and height
seven has $g=37$. Together with the earlier height-four handoff at $g=45$
and height-five at $g=15$, this gives exact full-continuation windows at
heights 1, 2, 4, 5, 6, 7, 8. Height 3 remains open: a separate exact scan rules out every
shared continuation with $g\le73$, including signed-zero and
nonstrong hosts. The [earlier ordinal-factor and signed-zero replay
method](../target_permuted/ROW_SCAN.md) works for any of these smaller
deletion counts by using all strong component factors through 73.
Replay **c++ -std=c++17 -O3 height3.cpp -o /tmp/height3** and
**/tmp/height3**; no `triple` or `signzero` lines are produced.
The lower bound 75 on odd $g$ is not a global impossibility of height 3.
A universal selector remains open.

Replay the literal rows with
```sh
c++ -std=c++17 -O3 literal.cpp -o /tmp/verify-boundary
awk '$2==6||$2==8' WITNESSES.txt | /tmp/verify-boundary 9
awk '$2==7' WITNESSES.txt | /tmp/verify-boundary 37
```

### An actual discrete intermediate-value path

The three height-eight witnesses live on the **same** labeled vertex
set. Turn the first into the second by reversing edge 2–4. Turn the
second into the third by reversing the directed triangle
$2\to4\to18\to2$, where 18 is the root. The two moves add exactly
2 each and leave all three continuation coordinates unchanged. Neither
move touches vertex 17, whose sole outward neighbor is the root. The
[terminal-slot law](../target_permuted/README.md) then preserves the
same sequence at every transitive inflation length of vertex 17: the
Hamilton outputs are $183047\cdot2^b-172591+d$ for $d=0,2,4$, with
common $(E,O,g)=(183047\cdot2^{b-1},183047\cdot2^{b-1}+9,9)$.
These new inflated windows do not hit another old-row height: the
left endpoint is 1 modulo 4 at $b\ge2$, whereas every old left endpoint
is 3 modulo 4.

There is a simple discrete intermediate-value criterion. In a graph of
literal edits preserving a specified continuation, require every edge
on a path to change the odd count by at most 2 in absolute value. If
two endpoints bracket an odd target, the first step reaching or crossing
it must hit it. This is a proof about the integer lattice, not a
continuity assumption on orientation space. The height-eight chain
is a literal instance. The complete flip graph does not have its
small-jump or same-continuation hypotheses.

## What real boundary limits lose

The motivating [1969 paper](https://www.ams.org/tran/1969-141-00/S0002-9947-1969-0243078-8/S0002-9947-1969-0243078-8.pdf)
studies limits obtained along an approach arc, not just a point's
location on the boundary. The [linked Lean blueprint's theorem statements](https://github.com/micah541/kaczynski-mathematical-works/blob/main/TJK/MainTheorems.lean)
also distinguish the existence of a good approach from a general limit,
and point out how the descriptive complexity can change when the domain
changes from disk to cube. These theorems are inspiration, not hypotheses
imported into the tournament proof.

Here the exact real chart is straightforward. For a terminal signature
$(F,D,g)$ put $t=2^{-b}$. The count ray is the line
$(t,tH)=(t,F+Dt)$. All intercepts at fixed slope approach the same
boundary point $(0,F)$. Knowing that a curve reaches that point does
not decide whether it attains a prescribed integer fiber. The
transverse coordinate $(tH-F)/t$ is exactly **$D$**, and admissible
points also require a dyadic address. Thus the chart that retains the
first-order approach, rather than just its endpoint, is necessary even
in this simplest straight-line family.

There is a genuinely quantitative theorem behind this distinction.
[COMPACTIFICATION.md](COMPACTIFICATION.md) proves that normalized
rooted readers with $g\le G$ and bounded relative defect reduce to
finitely many terminal charts, with at most
$$
\exp\left(\frac{9\log2}{2(\log5)^2}(\log G)^2+O(\log(G+2))\right)
$$
core presentations for fixed relative-defect bound. For fixed $G$
they leave a gap of positive **relative** length in every sufficiently
high octave, although all nonconstant rays tend to the same normalized
profile. In particular such roots with $g\le\exp(c\sqrt{\ln X})$
hit only $o(X)$ counts below $X$ for every fixed $0\le c<0.911$ (the
exact threshold is in the proof). This does not exclude larger
deletion counts or unbounded relative defects.

The lesson transfers to other exact inverse problems: a scaled
approximation or a measurable set of good directions gives no
uniform control of the **unscaled** residual. A nested compactness
argument would need compatible physical choices and a tight final
error, for instance $|H-N|<2$ for odd $H,N$. Merely having infinitely
many approximate approaches at each successive scale is weaker.
The graph path criterion supplies this missing control when it
applies; it is not implied by density of projected limits.

## Exact exclusion for one natural slope-preserving strategy

A tempting normalized strategy on the old row is to keep
$2E=756\cdot2^s$ and look for the three intercepts
$D=H-2E=-33,-31,-29$. This strategy **cannot** work with $g\le100$,
for any height. In particular the middle intercept $-31$ is absent
at this mantissa in the stated deletion range.

Here is a finite search certificate with no bound on suffix length.
In a strong normalized root deletion write $A$ for old outcontacts,
$S$ for weighted paths starting in $A$, and $K$ for weighted even paths
starting and ending in $A$. Then $S=2E$ and $D=g-2K$; these identities
are from the [terminal decoder](../target_permuted/NEW_READERS.md).
Strip every singleton deletion component after the last contact
component. Its optional position doubles $S$ but leaves $g,K$
unchanged, even if followed by other SCCs. This preserves the odd
part of $S$, which is **189** for the desired slope. Among $J$ singleton
components strictly between first and last contact components, fixing
one contact at either end and choosing an even subset between them
gives $K\ge2^{J-1}$ when $J\ge1$. Therefore at most
$3+\lfloor\log_2 K\rfloor$ singleton components survive. If the
last remaining deletion component is a singleton it must be an
outcontact, otherwise it too would strip off. This yields a bounded
search in $(g,K)$ before the dyadic suffix is programmed.

The compiler generates the strongly connected deletion components
of counts at most 100 by adding one vertex to each orientation up to
isomorphism and pruning induced counts greater than 100. Degree-ordered
canonical relabeling allows all equal-degree permutations, so does
not discard a type. Busch's bound excludes orders ten and higher.
It builds all ordered component lists and all allowed singleton slots:
53,269 normalized deletion layouts. Exact parity convolution assembles
the start weights and even endpoint weights. The root contact search
prunes partial $K$ only by its nonnegative contributions, and partial
$S$ only by the sum of remaining weights. This exhausts all candidates
at the required mantissa. [ALIGNED.txt](ALIGNED.txt) gives the complete
$(g,D)$ support in the nine tested offsets, with one smallest $S$ and
a literal reduced rooted row for each. None has $D=-31$.

For nonstrong hosts at this missing middle intercept, normalization
forces a terminal root. Every prefix count divides $2E$ and $H-2E=-31$,
and is odd. Since $\gcd(189,31)=1$, the prefix multiplier is one;
singleton padding cannot add a solution. A terminal singleton root has
$E=0$ and is irrelevant. Hence the exclusion is not just for strong
hosts. It rules out a finite deletion band **with the old exact slope**,
not other continuations at the same target. The successful new rows
explicitly use other slopes.

Replay with **c++ -O3 -std=c++20 aligned.cpp -o /tmp/align-sign**
then **/tmp/align-sign 3 100 189**. The recorded layouts and reduced
profiles are finite; the unbounded-height conclusion follows from the
stripping argument, not a census of many heights.
