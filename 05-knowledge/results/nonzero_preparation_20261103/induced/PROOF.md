# A prepared source cannot take its next +2 by any induced extension

2026-11-03. **PROVED from FINITE-EXACT computation and a CITED order bound — [independent audit](AUDIT.md) passed.** The all-order conclusion combines a finite exact
induced-extension certificate with a CITED strong-order minimum. The positive
all-clock clone continuation is proved separately below. No reusable
nonzero two-branch class is asserted.

## 1. Actual source and the two conclusions

Take the original-clock t=3 output of the proved
[prepared-four continuation](../../multiport_grammar_20261103/component/PREPARED_FOUR.md).
Its eleven old labels, in order, are

```text
p0, D0,D1,D2,D3,D4,D5, u, r, x, y.
```

The literal out-rows are

```text
S=(1790,1668,1672,1778,1670,1942,1718,1280,1631,128,512),
H(S)=1429, F_r(S)=1563, r=8.
```

In particular the old triangle has `u -> y -> x -> u`. Every vertex of D
beats the triangle, r beats x and y, and the old D port is D4, label5.
An **induced extension** means any tournament containing these eleven named
vertices with every old arc unchanged. There is no restriction on the
incidences or number of the new vertices. The finite search classifies strong
final cores; an SCC argument then covers weak endpoints too. Intermediate
induced graphs may always be weak.

The consequences of the complete certificate are:

1. No induced extension of S, strong or weak, has Hamilton count
   `4289=3*1429+2`, at any order. This does not even require an F constraint.
2. Every strong induced extension of S having Hamilton count
   `4291=3*1429+4` is one of the following, up to relabelling only the new
   vertices: a single new vertex with incoming old set `{p0,D4,r}`, or two
   new vertices enlarging the old singleton y to a transitive three-vertex
   block with the same outside incidences. The first has F_r=4677, while
   each of the six labelled two-vertex placements has F_r=4689=3*1563.
3. Every weak induced extension with count4291 is precisely a transitive
   prefix, one of those strong cores containing all of S, and a transitive
   suffix, in that order. All prefix/suffix vertices are new. In particular
   dropping final strongness does not restore the missing +2.

The all-order scope is specific to this actual eleven-vertex source. It is
not an impossibility theorem for every tournament of H1429, every larger
clock, or all operations. It rules out a uniformly applicable induced +2
recipe whose domain includes this source. If every old vertex must remain, an actual +2 construction must change
at least one old arc, even if a weak final graph is allowed. Temporary helpers and deletions of newly introduced
vertices do not escape an induced-final-graph obstruction.

The positive two-birth placement extends to every original clock t>=3,
with both intermediate births strong and every old arc/name retained:

```text
(189*2^t-83, 189*2^t+51)
   -> (567*2^t-245, 567*2^t+153).
```

Thus the actual scalar chain at t=3 continues
`157 -> 475 -> 1429 -> 4291`. The last root defect is398.

## 2. Inheritance, hostiles and preserved information

The inherited C2 operation, exact joint response, all-auxiliary stopping
boundary, and transitive clone-block formulas are proved in
[multiport defect grammar, construction/RESULT.md](../../multiport_defect_grammar_20261103/construction/RESULT.md),
Sections1–5. In particular its H17-to-H53 example already exposes the
one-time clone continuation. We apply that existing identity to the new
physical prepared family; the clone identity itself is not claimed new.

The old all-auxiliary obstruction is only about the C2 recipe. This search
allows every induced incidence, including all dual cuts, and never treats
that old theorem as an exhaustive extension classifier. The corrected near
miss is also concrete: the one-birth endpoint below hits H4291 at one clock
but fails the all-clock identity and the retained-root target. The least-used
sidecar is the full ordered two-path joint endpoint matrix K, which makes
an arbitrary birth a quadratic function of its entire actual cut.

The search preserves all old labels and the complete old relation. It does
not identify graphs with equal H, equal H/F, or equal seven-coordinate
assembly signature. No isomorphism quotient, new-label quotient, or port
uniformity filter is used. The count monotonicity used for pruning is proved
below; final strongness is checked only after the count target is reached.

## 3. Exact arbitrary-birth law and safe pruning

For any current tournament T, let e_i and s_i count Hamilton paths ending
and starting at i. Let K_ij count ordered covers of T by two nonempty paths,
the first ending i and the second starting j. For a new vertex v with
incoming old set I, deleting v from a Hamilton path gives the exact formula

```text
H(T+v)=sum_(i in I)e_i + sum_(j outside I)s_j
       +sum_(i in I,j outside I)K_ij.                 (3)
```

The first two terms are v-last and v-first. The third is v-internal, with
its two adjacent old endpoints retained. These cases are disjoint and their
inverses insert v at precisely the indicated join. Diagonal K_ii is zero.
This proof makes no strongness assumption.

With cut bits b_i, put

```text
A_i=e_i-s_i+sum_j K_ij, B_ij=K_ij+K_ji.
H(I)=H(T)+sum_i A_i b_i-sum_(i<j)B_ij b_i b_j.       (4)
```

The producer evaluates all cuts by removing the least set bit i:
`H(I)=H(I-{i})+A_i-sum_(j in I-{i})B_ij`. It independently checks (3),
(4), and direct forward path counting on every one-birth extension of every
labelled nonempty tournament through order4: 1,098 full cases.

For monotonicity, insert a new vertex into each old Hamilton path by a fixed
rule: insert before the first vertex it beats, or append if it beats none.
If insertion is not at the first position, its preceding vertex beats the
newcomer. Deleting the new vertex recovers the old path. Thus this is an
injection and `H(T+v)>=H(T)` for every tournament and every cut.

It follows that a partial birth word with H>4291 has no descendant attaining
4287,4289 or4291. Pruning that partial word loses no requested endpoint,
including endpoints whose intermediate graphs are weak. Conversely, every
labelled induced extension by d vertices has a unique birth-cut word in the
fixed order of those d new labels. Every prefix of a target endpoint passes
the monotone bound and is therefore explored.

## 4. The complete finite universe and all-order reduction

At depth d, the unfiltered labelled universe has
`2^(11d+d(d-1)/2)` extensions. The source includes every incoming subset at
every retained parent. Neither early strongness nor a selected port palette
is imposed. Exact search levels are:

| New vertices d | Raw universe | Tested cuts after safe prefix pruning | Retained H<=4291 |
|---:|---:|---:|---:|
|1|2^11|2,048|37|
|2|2^23|151,552|164|
|3|2^36|1,343,488|766|
|4|2^50|12,550,144|4,176|
|5|2^65|136,839,168|26,520|

At each of these depths the count4289 occurs zero times even before imposing
strongness. For count4291, the strong counts are respectively `1,6,0,0,0`.
The six depth2 objects are exactly the six total orders on the block
`{old y,new11,new12}`. Four depth2 strong endpoints with H4287 provide zero
controls; they are retained in the receipt but are not additional nonzero
returns. Depth0 is excluded directly by H(S)=1429.

**CITED order bound.** Arthur H. Busch, *A Note on the Number of Hamiltonian
Paths in Strong Tournaments*, EJC13 (2006), N3,
[Theorem1 p3 and Corollary1 p4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/),
proves the exact strong-order minimum

```text
m(3q)=3*5^(q-1), m(3q+1)=5^q, m(3q+2)=9*5^(q-1), q>=1.
```

I independently recovered the primary statement. Its values are strictly
increasing, and `m(17)=9*5^4=5625>4291`. Hence every strong endpoint of either
requested count has order at most16, so contains at most five vertices beyond
S. The complete finite universe therefore proves the stated all-order
strong-core consequences. The repository also records this cited bound
in [c4 regeneration](../../c4_port_returns_20261103/regeneration/PROOF.md);
no new strong-order extremal theorem is asserted here.

### Removing the final-strongness hypothesis

The following reduction is independent of the search and applies to any
strong source B with odd Hamilton count h>=3. Suppose an induced extension T
has `H(T)=3h+d` for d in {2,4}. Every old vertex lies in a single strong
component C of T, since all old directed paths remain. Hamilton paths
traverse the linearly ordered SCCs in their forced order, so
`H(T)=H(C)*product_(J!=C)H(J)`. The insertion injection gives H(C)>=h.

Every tournament has odd Hamilton count (the inherited Redei parity theorem,
also recalled on p1 of Busch). If the outside product exceeded1, it could
not equal3 because `3h+d` is not divisible by3. Thus it would be at least5,
forcing `H(C)<h`, since `(3h+d)/5<h` for h>=3. This is impossible. The outside
product is1, every outside strong component is a singleton, and H(C)=H(T).
Thus C is one of the strong cores already bounded and classified above.
The other vertices form an arbitrary transitive prefix and suffix around C.

For the present +4 classification this also predicts the complete number
of labelled endpoints at any new-vertex depth d>=1: it is
`d!*(4d-3)`. The one-new-vertex core contributes `d!*d`; the six labelled
two-new-vertex cores contribute `3*d!*(d-1)`. The exact depth1–5 totals
`1,10,54,312,2040` in the search agree with these independent SCC assembly
counts. There is no analogous +2 core, hence no +2 induced extension at
any order, strong or weak.

If the prefix has a vertices, the retained-root response is `2^a F_r(C)`;
a transitive suffix has no effect on it. Therefore the +4 extensions that
also triple F_r use a clone core with no prefix, and may have any transitive
suffix. Requiring final strongness removes that suffix as well. This final
statement is an endpoint classification, not a legal-operation claim about
all possible ways of constructing it.

## 5. The uniform positive clone continuation

Let B be any old tournament, with retained root r and auxiliary a. Use the
inherited quantities h,F,S,G for this actual marked source. Add x beating
only a, and a nonempty transitive block Y_m with

```text
B -> Y_m -> x -> a,  every vertex of B-{a} -> x.
```

The inherited exact clone formulas are

```text
h_m=h+2^m S, F_m=F+2^m G.                            (5)
```

The proof retains ordered slots: if x is terminal all Y occurs in one
forced block and contributes h or F; otherwise x is immediately followed
by a and the Y labels can be assigned independently to the two allowed
slots, contributing the factor2^m. These are actual vertex assignments,
not replacement by a scalar-equivalent graph.

At m1 to m3, subtraction gives

```text
(h_3,F_3)-3(h_1,F_1)=(2(S-h),2(G-F)).                 (6)
```

Thus enlarging the first C2 image's singleton y to a transitive block of
size3 repeats its original response defect once. In the actual prepared
family, B has `S=h+2,G=F`, so (6) is exactly `(4,0)` at every clock. One may
append each new vertex after the existing y-block, with incoming set all
previous vertices except x. Each newborn then has incoming vertices in the
old strong graph and an outgoing arc to x, so both births are strong.
The other five label orders describe the same abstract induced extension
relative to the old names.

This is the inherited clone mechanism applied to the actual source. It does
not repeat indefinitely. Asking another transitive clone enlargement to
repeat the same defect `(2(S-h),2(G-F))` requires, already from H,

```text
2^(m')=3*2^m+2.
```

For positive integer m this has only `m=1,m'=3`: when m>=2 the right side
has2-adic valuation1, while it exceeds2. Therefore the isolated clone step
is not a reusable same-letter channel.

There is also a direct all-auxiliary stopping statement for m>=2 when B is
nontrivial strong. Every old vertex has at least three outgoing neighbours.
Every nonfinal Y vertex either has at least three, or has exactly two and
the residual graph contains B. The dual root-gap formula from READY §1 gives
start gap at least6 in either case. The final Y vertex has its sole
out-neighbour x, and its start gap is exactly H(B)>=3 by deleting it and
counting paths ending x. Finally x has start gap0 because deleting x leaves
`B -> Y_m`, so no Hamilton path can end at a in B. Thus no auxiliary has
start gap1 or2, and no C2 digit2/4 is prepared by the clone continuation.
The nontrivial-strong hypothesis matters: for singleton B the final clone's
gap is1. This statement concerns C2 readiness, not all possible next recipes.

## 6. The one-birth hostile to single-clock inference

The unique depth1 strong H4291 endpoint has incoming old cut289, namely
`{p0,D4,r}`. It adjoins a vertex to D that is beaten only by D4 inside the
core, while the prefix and root beat it and it beats the old triangle.
The new seven-vertex core and its port4 have literal signature

```text
D'=(2,4,57,3,75,27,47),
(h,p,f,g,j,e,s)=(41,712,47,128,0,2,9).
```

At arbitrary prefix parameter z=2^k this gives

```text
H'=2277z-263, F_r'=2277z+123.
H'-(3H+4)=9z-18, F_r'-3F_r=9z-30.
```

At k1 it therefore has H4291 but F4677, below the prescribed tripled F by12.
At every larger k it fails the H identity as well. It must not be promoted
to an all-clock recipe. In development, the initial scout printed target
hits only at depths>=2; the explicit depth1 recovery exposed this additional
H-only point. The final certificate includes every depth1 target explicitly.
The initial literal-source transcription also failed its tournament gate
and was repaired before any frozen receipt. Neither issue affected the
constructively generated scout's pruning or its depth2–5 counts.

## 7. Exact source contract and replay

`induced_complete.cpp` is standalone. It reconstructs the endpoint and joint
cover matrices at every retained parent, evaluates every cut, and recurses
without quotienting. Direct forward counting independently verifies every
reported strong target graph and its retained-root F. Every such graph is
checked to be a tournament with exactly the eleven old induced rows. A
separate small exhaustive birth universe validates the cut polynomial.

All arithmetic is integral. Parents have order at most15 and children at
most16. Path counts are bounded by16!, and a two-slot count by17*16!, both
below2^63. At a parent of order n, the sum of all K entries counts two
nonempty slots and is at most (n-1)n!. Each signed cut update starts from
a nonnegative child count plus A_i and then subtracts nonnegative B terms;
its absolute intermediates are bounded, for example, by 4*17*16!<2^63. Unsigned path counts
and signed cut updates therefore incur no overflow. The program uses
always-active exceptions, not assertions removed by optimization.

Reproduce from the worktree root:

```bash
c++ -std=c++17 -O3 04-computation/nonzero_preparation_20261103/induced_complete.cpp \
  -o /tmp/induced_complete_O3
/tmp/induced_complete_O3 > /tmp/induced_complete_O3.json
c++ -std=c++17 -O0 04-computation/nonzero_preparation_20261103/induced_complete.cpp \
  -o /tmp/induced_complete_O0
/tmp/induced_complete_O0 > /tmp/induced_complete_O0.json
```

Optimized and unoptimized C++ replays both pass **150,939,321 always-active
gates** and produce the identical **2,615-byte** receipt. The retained-row
stream FNV-1a64 digest is `8281660521594346102`; it is a deterministic check
of this traversal, not an isomorphism hash.

- Source `induced_complete.cpp`: 8,894 bytes, SHA256
  `3e2647525cc79a33788031e5b2e11dad425eb6744fa7cefb93215bc0ba3d83f1`.
- Receipt `induced_complete.json`: 2,615 bytes, SHA256
  `9247ce871dde82aaca8b441e23fb2a07bd69c0a018f13e619f7a33af6c551cfa`.

Separate literal-permutation and actual-family controls are recorded by
`induced_controls.py`. Normal and optimized Python replays pass **2,261
gates** with identical 1,272-byte receipts. They independently recount every
reported strong endpoint, validate all 1,098 small cut cases by literal
permutations, recover the one-birth core signature, and check the clone
response and every-auxiliary stopping boundary at k1,2,3. The singleton
source hostile correctly fails that stopping boundary. Reproduce with:

```bash
python3 -B 04-computation/nonzero_preparation_20261103/induced_controls.py
python3 -B -O 04-computation/nonzero_preparation_20261103/induced_controls.py \
  --output /tmp/induced_controls_optimized.json
```

These are bounded independent controls; the all-order negative conclusion
uses the complete search and its proved/cited reduction. Full frozen pins
are recorded in [MANIFEST.json](MANIFEST.json).
The independent analytic/source/replay audit passed. The precise stopping
object is an induced-source obstruction at H1429; other old-arc surgeries,
other source representatives, or a preparation reached at another clock
remain available. The finite arbitrary-cut result is stronger than the
initial small port-uniform menu but is not a universal birth obstruction.
