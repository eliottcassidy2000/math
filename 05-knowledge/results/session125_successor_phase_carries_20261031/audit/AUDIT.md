# Independent audit: phase obstruction and radix-49 successor carries

Session 125, 2026-10-31. **PASS.** This audit writes only its own scratch
directory. It imports no producer module, does not edit repository sources or
remotes, and never compiles or executes the external C source. No permission
denial or unresolved execution blocker occurred.

## 1. Independent finite input

`independent_audit.py` constructs square-sum edges by enumerating squares and
extends oriented paths from vertex 1, using only necessary residual-degree and
connectivity tests. This differs from the sibling independent onset solver's
degree-subset enumeration. It independently recovers:

- Order 41: 49 eligible paths, 38 positional-parity masks.
- Order 42: 275 eligible paths, 132 masks.
- Exactly three common masks, with nice-pair multiplicities 1,4,1:
  `1187208845135`, `1453217412461`, `1487448155501`.
- The first two masks admit order 43 continuations but none at 44; the third
  admits no order 43 continuation. The receipt includes independently generated
  positive witnesses, not paths copied from external tables.

The negative certificates are independently rechecked. For the first mask,
degree-two vertices 18,42 force both edges at 7, and degree-two vertices 16,44
force both edges at 20. Vertex 29's allowed neighbors are only 7,20,35, so its
required internal degree 2 is impossible. The second mask gives internal
vertex 44 only neighbor 37. The third gives internal vertex 43 only neighbor 38.

The exhaustive searches at 43 and 44 are redundant checks; the short degree
certificates already establish nonexistence. All search pruning is necessary,
and every completed path is counted once because orientation starts at 1.

## 2. Finite phase defects: scope and proof check

If a common-position-parity original-port family exists at every order
`n>=n0`, then each shared label's phase becomes permanent. Counting odd slots
when a label arrives forces every label greater than `n0` to have its numerical
parity. Thus the defect set is finite, contained in `[1,n0]`; comparing the
odd-slot and odd-label counts gives equally many odd and even defects.
Endpoints at sufficiently large odd and even orders force ports 1,3,8 to be
nondefects.

Nondefect-to-nondefect edges have odd-square sums, hence residue sum 1 modulo 8.
Their four residue components are `{0,1}`, `{2,7}`, `{3,6}`, `{4,5}`. A
defect-to-nondefect edge has even-square sum, so a single defect can attach
only to one of the two supergroups `A∪D8` and `B∪C`.

At a sufficiently large even order, every residue component has nondefect
vertices. A path beginning at 1 and ending at 8 begins and ends in A. It needs
at least five surviving nondefect runs, hence four intervening defect runs.
At least two runs cross between the supergroups and therefore each uses at
least two defects. Consequently there are at least six defects. An independent
finite weighted-walk check returns 6 for A-to-A and 4 for A-to-C.

**Quantifier boundary:** the six-defect statement is a necessary condition,
not an existence theorem. The exact 41..44 obstruction excludes a coherent
family beginning at 41, not every later-start family or every cycle family.

## 3. Exact within-pair edits and complete finite audit

For the pinned inert integer arrays, the audit validates all 1,992 base nice
pairs at orders 41..2032, including coverage, square edges, endpoints and shared
position parity. It checks all 196 transition words, all 98 member-change states,
the forced threshold rule, and an actual expanded nice pair for every state.
Only one full block changes member in each state. The common input phase
keeps all shared labels in the same block, so old internal edges and connector
edges cannot exchange roles. Therefore

`d_pair(n)=d_pair(parent(n))+K(digit,parent parity)`.

The complete audited constant range is `9<=K<=44`; maximum base path deletion
count is 1699, and maximum base cycle count 1700. This yields logarithmic
within-pair deletions but does not identify the longer member with the next
pair's shorter member.

As an additional check, any potential table connector/internal collision with
canceling signed offsets decodes only to parent pair `{1,3}` or loop `{8,8}`.
The former is forbidden by the original pair phase and the latter by simplicity.

## 4. Audit strengthening: the carry recurrence is exact

Write `H(r)=|E(L_(r-1))\E(S_r)|`. At a full carry `t=49r+24`, both output
words consist of all 49 offsets of one order-r parent path, plus 24 singleton
atoms. All internal transformed labels are at least 25. The decoder

`center(z)=floor((z+24)/49)`

recovers the unique parent label. The union of lifts of an undirected parent
edge `{x,y}` is

`{{49x+c,49y-c}: -24<=c<=24}`,

independently of that edge's input phase. Different parent edges have disjoint
lifts, so internal deletions multiply by 49 even when parent phases differ.

There is a stronger separation than the initial `49H±72` estimate. Both
parents have endpoints 1 and `q`, where q is 3 for odd r and 8 for even r.
A singleton-incident connector has a label<=24 and cannot be internal. Every
other connector decodes to two labels from `{1,q}`. If equal, these would be
a parent loop; if distinct, they would be the closing edge `{1,q}`. Neither
can be an internal edge of a Hamiltonian path of order r>2 with endpoints 1,q.
Thus **all carry connectors are disjoint from both full internal lifts**.

Consequently

`H(49r+24)=49H(r)+K_carry`, with `0<=K_carry<=72`.

For the original fixed words, the complete connector comparison gives
`K_carry(even)=40`, `K_carry(odd)=52`. Independently expanded samples are

| Seed | Depth 0 gap | Depth 1 gap | Depth 2 gap |
|---:|---:|---:|---:|
|42|18|922 at 2082|45218 at 102042|
|43|9|493 at 2131|24209 at 104443|

The common closing edge is absent internally in both paths, so closing two
same-order representatives leaves H unchanged.

## 5. Canonical and selector quantifiers

If H(42)=H(43)=0, the fixed starting vertex 1 turns undirected edge equality
into equality of oriented paths. The three base pairs at 41,42,43 would then
glue to the forbidden coherent family 41..44. Hence at least one of H(42),
H(43) is nonzero. Same-endpoint Hamiltonian paths cannot differ by only one
old-edge deletion: equal vertex degrees force the replacement edge to have
the same endpoints. Therefore one seed gap is at least 2.

Along `t_k=49^k(t0+1/2)-1/2`, exact carry separation gives
`H(t_k)>=2*49^k`. If `B(n)` is the within-pair cycle deletion count and
`D(n)` the canonical short-member cycle deletion count, directed edge-set
triangle inequalities give precisely

`H(n+1)-B(n)-1 <= D(n) <= H(n+1)+B(n)`.

The minus 1 is essential: reversing a comparison between an n-cycle and an
(n+1)-cycle adds one deletion. Since B is O(log n), the universal canonical
bound is **limsup D(n)/n >=4/87**, improving the initial 1/87 bound.

For the pinned even seed 42,

`H(t_k)=(113*49^k-5)/6`, `t_k=(85*49^k-1)/2`.

The same two-sided triangle estimate gives the **actual ray limit 113/255**,
not merely a lower bound, for canonical deletion cost divided by old order.

The local two-step selector corollary also passes. At a carry r, put
`A=close(S_(r-1))` and `Z=close(L_r)`. Then

`delta(A,Z)>=H(r)-B(r-1)-B(r)-2`.

Any allowed selector at r-1 differs from A by at most H(r-1)<=72; any allowed
selector at r+1 differs from Z by at most H(r+1)<=72, because both neighbors
are noncarry orders. Their two-step gap is therefore at least

`H(r)-B(r-1)-B(r)-146`.

Any middle cycle, even one outside the two-representative bank, must split
that gap over two successor steps. At least one step is half this bound.
Thus allowing an arbitrary middle cycle, while restricting its two neighbors,
still gives universal limsup at least 2/87, and the pinned seed 42 yields at
least 113/510. This does **not** exclude cheap successor families free to use
other cycles at the neighboring orders too.

For a genuine two-representative bank selector, the stronger case split keeps
the **full**, not half, lower bound. Set `U=close(L_(r-1))`, `V=close(S_r)`.
If the chosen middle is V, its incoming step has cost at least
`H(r)-B(r-1)-73`: compare the selected predecessor to U through
`close(S_(r-1))`. If the chosen middle is U, its outgoing step has cost at
least `H(r)-B(r)-73`: compare the selected successor to V through
`close(L_r)`. The constants 73 are exactly 72 for the neighboring representation
gap plus 1 for reversing an adjacent-order deletion comparison.

Therefore **every two-representative selector** satisfies the full universal
limsup lower bound 4/87, and the pinned 42-ray full lower bound 113/255. Merely
switching which stored representative names the current integer cannot
remove the large carry; it can move that carry to the other adjacent step.

## 6. Exact pinned global limsup: finite maximum, not just one ray

The independent audit also completes a bounded finite maximization. Every
order t>=42 has a unique backwards carry reduction until either:

1. its root s lies in `[42,2033]`; or
2. s>=2034 is a noncarry order.

The exceptional root 2033 must be included. Its predecessor uses a stored
base pair while its own pair is recursive, so the noncarry 72 bound cannot
be applied there. The audited value is `H(2033)=1861`.

Let kappa be 40 or 52 according to the root parity. Along a rooted carry ray,

`H(t_k)=49^k(H(s)+kappa/48)-kappa/48`,
`t_k=49^k(s+1/2)-1/2`.

Its limiting ratio is

`R(s)=(H(s)+kappa/48)/(s+1/2)`.

For every root s>=2034, `H(s)<=72`, hence
`R(s)<=(72+52/48)/(2034+1/2)=877/24414`.
This is below the already established 42-ray ratio 113/255. The maximum must
therefore occur among the 1,992 roots 42..2033. Exhaustive exact rational
comparison gives

`max R(s)=22273/24042`, attained at `s=2003`, `H(2003)=1855`.

Here the odd carry constant is 52. The first descendant is 98171 with gap 90947.
The exact finite-k ratio differs from R(s) by
`(R(s)/2-kappa/48)/t_k`, uniformly O(1/t_k). The uniform logarithmic
within-pair bound and canonical triangle squeeze show
that this is the **exact global canonical limsup**, not a finite-size maximum
or an extrapolated empirical rate.

It is also the exact global limsup for **every** two-representative selector.
For the upper bound, an arbitrary bank step at n costs at most
`H(n)+B(n)+H(n+1)`. Outside the finite boundary region, at most one of n,n+1
is a carry, while the other representation gap is at most 72. Thus its
normalized limsup is at most the same M. The preceding middle-choice case
split gives limsup at least M along the maximizing carry ray. Hence switching
between the two representatives does not even change this asymptotic worst
case for the original fixed constructor.

## 7. Reproduction and source policy

From this packet's root, with the separately obtained pinned C source:

```sh
python3 audit/independent_audit.py /path/to/squares.c --output independent.normal.json
python3 -O audit/independent_audit.py /path/to/squares.c --output independent.optimized.json
cmp independent.normal.json independent.optimized.json
```

All correctness checks use explicit exceptions, not Python assertions.
The script accepts only inert integer-initializer syntax after verifying the
whole-file SHA256:

`602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e`.

The 17,321,973-byte external source and its arrays are not redistributed.
The separate source URLs, pin, and no-redistribution policy are recorded in
[../SOURCE.json](../SOURCE.json), inherited from the Session 124 Gerbicz reconstruction. Reproducing the numeric table checks requires that
separately acquired pinned source; the mathematical carry separation and
universal lower bounds do not rely on any particular external base choice.

No result here identifies centered radix 49 with Fibonacci/Zeckendorf
numeration, nor claims a global no-go theorem for unrestricted graph families.
