# Ninja means weak B: head-cut radix-25 successor carries

**Status:** PROVED structural carry mechanism; FINITE-EXACT weak-B onset input;
VERIFIED owned A/B certificates and bounded examples. This is a separate
extension, not a renaming of the original-port Gerbicz Q49 theorem.

## 1. Correct object and scope

A weak-B pair `Pair(n)=(S_n,L_n)` consists of square-sum Hamiltonian paths
on `[1,n]` and `[1,n+1]`. Both start **1,3**, end at **8 for even order**
and **10 for odd order**, and agree on the position parity of every common
label. This is the current “ninja pair” convention.

The Session 124 owned operators A and B both preserve weak B. They use
radix 25, all offsets -12 through 12, singleton labels 1 through 12, and
may split a transformed input path only at the head edge 1–3. A preserves
an additional prescribed prefix based on its hook 12; B uses hook 8 and a
longer prescribed prefix. These extra prefix properties are not needed for
the carry theorem.

**The objects here are paths, not directly closed cycles.** In particular,
`1+10=11` is not square, so an odd weak-B path cannot simply be closed back
to 1. The owned separate cycle-consumer certificates are a different map;
no successor-cost theorem for their output cycles is claimed here. The
original-port Q49 cycle-bank theorem stays separate and unchanged.

The base interval is 42 through 1061, with first recursive pair at
`1062=25*42+12`. At larger n put `m=floor((n-12)/25)` and `d=n-25m`, so
`12<=d<=36`. The homogeneous A and B constructors used for the numerical
examples select their corresponding owned paired word at every recursion
level. No external Gerbicz C source is needed.

## 2. A finite fixed-phase obstruction for weak B

The independent weak-B onset catalogue gives exactly six common odd-position
masks for the first pair at 42/43. Bit v−1 records label v. If these phases
were kept unchanged in a family of consecutive paths, every new label v>42
would have to acquire natural parity `v mod 2`: subtract the odd-position
counts `ceil(n/2)` at consecutive orders.

The following elementary degree certificates kill all six possibilities:

| Common mask at 42/43 | Obstructed order | New internal vertex's only phase-compatible neighbor |
|---:|---:|---:|
| 1585793692427 | 46 | 18 |
| 1599483269979 | 46 | 35 |
| 1602705117979 | 47 | 34 |
| 3372531545867 | 46 | 18 |
| 3389445068555 | 47 | 34 |
| 3526867243787 | 47 | 34 |

These are exact neighborhoods in the phase-filtered square-sum graph:
46 has arithmetic neighbors 3,18,35; 47 has arithmetic neighbors 2,17,34.
After applying the fixed phase, only the displayed neighbor remains. The
new vertex is not a prescribed endpoint (8 at order 46, 10 at order 47),
so it would need path degree two. Therefore no weak-B phase-coherent family
exists across the whole interval **42 through 47**. We do not claim that
each displayed order is its phase's first failure; no larger search was
needed.

The catalogue dependency is the separately checked weak-B onset packet,
[the published weak-B certificate](../../session125_origin32_numeration_20261031/pair_taxonomy/weak_b_certificate.json), whose
complete degree-subset and endpoint-growth searches agree. This checker
revalidates all six positive 42/43 pair witnesses and the displayed negative
neighborhoods. It does not infer completeness merely from six examples.

For any proposed globally fixed weak-B phase, the finite defect set would
necessarily include labels **3 and 10** (phase even and odd respectively),
unlike the original-port case. The original-port six-defect theorem is not
being transferred to this changed port convention.

## 3. Separate the common head edge from the moving body

Let `E(P)` denote a path's undirected edge set and
`delta(P,Q)=|E(P) minus E(Q)|`. Every weak-B path contains the head edge
`e0={1,3}`. Put `R(P)=E(P) minus {e0}` and let
`T={1,3,8,10}`.

No edge of R(P) has both endpoints in T. Among distinct labels of T, the
only square-sum pairs are 1–3 and 1–8. The first was removed; the second
cannot occur because vertex 1 already has its sole path edge 1–3. This
argument is uniform in the path's order and terminal parity.

The lift at offset c sends the label x in position j (starting at 1) to

```
25*x + (-1)^j*c.
```

A full block contains all lifted R(P) edges and one lifted 1–3 edge.
A split head/tail pair contains exactly the same lifted R(P) edges, but
the 1–3 carrier has been cut. Reversing either piece does not change its
undirected internal edges. Thus the moving body R(P) is unaffected by all
permitted head splits and orientations.

Use nearest-center decoding `nu(z)=floor((z+12)/25)`. It recovers x from
each transformed label. A boundary edge is either a connector between atoms
or an uncut 1–3 carrier. Its endpoints decode only into
`{0,1,3,8,10}`: 0 denotes a singleton, and the other possibilities are head,
tail-start, or terminal ports. A lifted R(P) edge cannot have both decoded
endpoints in this set. Hence **every lifted body edge is disjoint from every
boundary edge**, even across different input paths, phases, cut sets, or
choices of operator A versus B.

There are exactly **61 boundary edges** in every valid output word. With h
split offsets there are `25+12+h` atoms, hence `36+h` connectors, and
`25-h` uncut head carriers. Their total is 61. Equivalently, at output
order N the body has N−62 edges and the Hamiltonian path has N−1 edges.
This count is independent of the cut pattern.

## 4. Exact body/boundary edit laws

Write `d(n)=delta(S_n,L_n)` and `H(t)=delta(L_(t-1),S_t)`.

**Within a pair.** The new input label m+1 has phase m+1. Its offset-c
image is `25(m+1)+(-1)^(m+1)c`; consequently the longer input member is
used exactly when

```
25+(-1)^(m+1)*c <= d+ext,  ext=0 or1.
```

Precisely one offset changes member between the two output paths. On their
common labels its map agrees, because the input pair agrees on phase.
All other body-edge sets coincide. Since 1–3 is common in the parent pair,
removing it does not change the parent deletion distance. The body/boundary
disjointness therefore gives an exact recurrence

```
d(n)=d(m)+K_pair,  0<=K_pair<=61.
```

The finite base interval implies `d(n)=O(log n)`. For the pinned owned bases
the maximum base value is 813 at n=1054. The actual A/B word constants are
respectively in 4..22 and 3..22.

**Noncarry bank boundary.** In the recursive range, if t is not congruent
to 12 modulo 25, the two representations `L_(t-1)` and `S_t` use the same
parent pair and the same member at each offset. Their bodies agree, so
`H(t)<=61`. The actual noncarry maxima over the owned words are 37 for A
and 36 for B.

**Full carry.** At `t=25r+12`, r>=43, the previous pair's longer output uses
`L_(r-1)` at all 25 offsets; the new pair's shorter output uses `S_r` at all
25 offsets. These parent paths may have different phases. Nevertheless each
parent body edge {x,y} contributes the phase-independent set

```
{{25*x+c,25*y-c}: -12<=c<=12}.
```

Changing the endpoint phase replaces c by −c and leaves this union unchanged.
Distinct parent edges have disjoint lifted sets by nearest-center decoding.
Thus body-edge differences multiply by exactly 25. No boundary edge can
cancel them. This proves the exact carry law

```
H(25r+12)=25*H(r)+K_carry,  0<=K_carry<=61.
```

Different cut sets can only change the nonnegative boundary correction;
they do not reduce the amplified parent disagreement. For the owned
homogeneous A and B tables, direct boundary-set comparison gives the stronger
identity **K_carry=0 in both parent parities**.

## 5. A universal obstruction for the weak-B constructor bank

Suppose arbitrary valid weak-B bases are chosen on 42..1061, together with
valid radix-25 head-cut words having the standard coverage rule above.
If all of `H(43),H(44),H(45),H(46)` were zero, the shared-start path equalities
would glue Pair(42) through Pair(46) into a phase-coherent family at every
order 42 through 47, contradicting Section 2. Thus one seed
`t0 in {43,44,45,46}` has H(t0)>0. Distinct Hamiltonian paths with the same
endpoints cannot differ by only one deleted edge: equal degree vectors force
the replacement edge to equal the removed edge. Therefore H(t0)>=2.

Along `t_k=25^k(t0+1/2)-1/2`, the carry law gives
`H(t_k)>=25^k H(t0)`. For the canonical path choice S_n, directed
set-difference triangles give

```
H(n+1)-d(n)-1 <= delta(S_n,S_(n+1)) <= H(n+1)+d(n).
```

More generally choose at each order n either bank representative `S_n` or
`L_(n-1)`. At a large carry r, neighboring bank gaps `H(r-1),H(r+1)` are
at most 61. If the middle choice is S_r, its incoming edit costs at least
`H(r)-d(r-1)-62`. If it is L_(r-1), its outgoing edit costs at least
`H(r)-d(r)-62`. One of these two steps is large regardless of the choice.
Consequently **every such weak-B path-bank selector** satisfies

```
limsup_(n->infinity) delta(P_n,P_(n+1))/n >= 4/93.
```

The constant follows from `2/(46+1/2)=4/93`. This is universal for the stated
head-cut Q25 constructor architecture and base interval, not for all weak-B
path families, not for arbitrary square-sum cycles, and not an impossibility
of a different successor algorithm.

## 6. Bounded owned examples, with no global maximum search

The owned small bases give `H(43)=0, H(44)=9, H(45)=14, H(46)=14`.
For both homogeneous A and B, exact carry multiplication gives:

| Seed | Carried orders | H values |
|---:|---|---|
| 44 | 44,1112,27812 | 9,225,5625 |
| 45 | 45,1137,28437 | 14,350,8750 |

On the seed-45 ray, `t_k=(91*25^k-1)/2` and `H(t_k)=14*25^k`.
Thus the canonical predecessor-step deletion ratio tends **28/91**.
For any two-bank selector, the maximum of the two steps adjacent to t_k,
divided by t_k, has the same limit 28/91: the preceding case split gives the
lower bound, and `delta(P_n,P_(n+1))<=H(n)+d(n)+H(n+1)` gives the upper.
This is a particular carry-ray limit, **not a claimed global limsup value**.

Concrete predecessor-step deletions at order 28436 to 28437 are 8787 for A
and 8769 for B; their within-pair costs at 28436 are only 42 and 24.
The bounded controls stop here: no global maximizing seed scan, new
constructor search, or unbounded path census is part of this extension.

## 7. Reproduction and dependencies

Both checkers default to the owned sibling result directories relative to
their source location, including `pair_taxonomy/weak_b_certificate.json` in
the Session 125 origin packet. Every input path can be overridden explicitly.
The checker uses Python's standard library, reads only owned repository JSON
and gzip-root base data, and never needs or executes external C. It checks
all 1020 base pairs, all 100 A/B paired states by actual expansion, exact
body/boundary recurrences, both parity carry rays through depth two, the
bank-selector inequalities, the six phase obstructions, and hostile controls.
All checks use explicit exceptions and remain active under `python -O`.

```
python check_weakB_carry.py \
  --operator-dir /path/to/session124_mixed_radix_square_rays_20261031 \
  --base-dir /path/to/square_sum_B_base_interval_20261019 \
  --small-bases /path/to/session124_square_prefix_threshold42_20261031/SMALL_BASES.json \
  --weak-b-catalogue /path/to/weak_b_certificate.json \
  --output normal.json
```

Repeat with `python -O` and a different output file, then compare receipts.
Compressed base chunks, their decoded data, the small bases, and every used
operator certificate are verified against the existing owned packet's pins.
The weak-B catalogue is a separately verified finite dependency; its six
positive witnesses and the new negative degree certificates are rechecked here.
