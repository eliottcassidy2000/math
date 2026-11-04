# The nested source, an exact one-exit compiler, and a new-port obstruction

2026-11-03. **PROVED** context identities and finite arithmetic reduction;
**FINITE-EXACT** coefficients, hostile examples and audit receipts.
Reusable two-branch closure and the unrestricted `{7,21}` spectrum claim remain
**OPEN**. No theorem ID is allocated here.

The useful progress is a complete description of a much larger failed route:
one may replace the five-vertex bottom core of the literal 170807 source by
**any tournament of any order**, and independently choose **any positive outer
transitive-prefix length**; this entire family still misses both 512423 and
512425. The successful literal switches genuinely change the exterior
incidence. A five-vertex hostile identifies precisely why a second physically
named port must then be retained.

## Inheritance and scope

Closest proved mechanism: the single-port zero-ready identity in
[ZERO_READY.md](../../triangle_return_170807_20261103/ZERO_READY.md), together
with the maximal-block-run proof of **THM-3121**,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The latter is an August 2 mechanism recovered here, not a new priority claim.
The least-used relevant sidecar is the deletion count `q=H(C-w)`.

The inherited hostile is the actual 81/27 pair, whose equal H does not give
equal birth response. The corrected near miss is the previously classified
two untouched-module enrichment: neither its scalar agreement nor clock
resynchronization proves a reusable carry. This packet keeps the source
adjacency and names instead of substituting an equal-H tournament.

Live concepts: nested zero-ready layers; named exit ports; ordered path-cover
slots; deletion monotonicity and parity; affine return equations. Their
connection here is exact: the slot interface compiles the layers, deletion
monotonicity constrains its integral values, and new exits force extra named
boundary information.

## 1. A universal one-exit context identity

Let C be any nonempty tournament with a specified vertex w, and E any disjoint
tournament, possibly empty. Every E vertex dominates `C\{w}`. Let `L⊆V(E)` be
the vertices beaten by w; all remaining E vertices beat w. No order or strong
connectivity assumption is needed.

Use `H(empty)=1`. For an induced E-subtournament X, let `e_L(X)` and `s_L(X)`
count its Hamilton paths ending and starting in L respectively, and set both
to zero on empty X. Put

```
h = H(C),       s = start_w(C),       q = H(C-w),
F = E_w(C,2),
a = H(E),       b = e_L(E),
p = sum_(X subset E) H(X) s_L(E\X),
d = sum_(X subset E) e_L(X) s_L(E\X).
```

Here F counts ordered two-slot spanning path covers whose first slot is
nonempty and ends at w; the second may be empty. Then

```
H(C inside E) = a h - b s - d q + p F.                  (S1)
```

Proof for `|C|≥2`. A Hamilton path has at most two maximal C-runs, because
only w can leave C and w is used once. One C-run in final position contributes
`a h-b s`: entering C from a predecessor in L forbids paths starting at w.
A single nonfinal C-run ends at w and contributes `p end_w(C)`. Two C-runs
contribute `p(F-end_w(C))-d q`: the subtracted cases are exactly those where
the first C-run is the singleton w, but its preceding E vertex belongs to L.
The end-count terms cancel. For singleton C, the final-run contribution is
`a-b` and the nonfinal contribution is `p-d`, giving the same formula since
`h=s=q=F=1`. Empty E gives H(C) directly.

This is a four-dimensional linear observer interface, not a complete graph
encoding. Its linear dimension cannot generally be lowered. The checker
exhibits four independent actual profiles `(h,s,q,F)` and four independent
actual context coefficient vectors. For example the profile basis is

```
(1,1,1,1), (1,0,1,2), (1,0,1,4), (3,1,1,3),
```

and the coefficient basis is

```
(1,0,0,0), (1,-1,0,1), (1,0,0,2), (1,-1,-2,3).
```

This dimension statement concerns linear coordinates and this observer family;
it makes no claim about arbitrary nonlinear encodings.

## 2. The literal source is five nested zero-ready layers

The bottom core, with its actual physical names, is

```
C=[6,7,8,14,31], rows=[14,16,10,18,5], w=31.
```

Equivalently `C=C3(TT2[6,8],TT2[14,7],{31})`, with block cycle
`[6,8] -> [14,7] -> 31 -> [6,8]`. Its profile is
`(h,s,q,F)=(11,1,1,13)`.

Starting at C, apply the following wrappers in order. In each wrapper the
listed prefix dominates the existing strong core, and the new root is beaten
only by the listed in-port.

| Prefix in order | New root | Sole in-port | Result H |
|---|---:|---:|---:|
| 4→2 | 22 | 31 | 63 |
| 54 | 33 | 4 | 211 |
| 1 | 21 | 33 | 633 |
| 40 | 32 | 33 | 2195 |
| 50→51 | 39 | 31 | 170807 |

Every internal vertex except 31 loses to every exterior vertex. Vertex31
beats just exterior22 and39. Therefore S1 applies even if the whole bottom C
is replaced by an arbitrary marked tournament. The exact exterior compiler
gives, in coordinate order `(h,s,q,F)`,

```
H(B)    = (53,0,0,124) · (h,s,q,F) = 53h+124F,
E31(B,2)= (0,0,-487,3280) · (h,s,q,F) = 3280F-487q.
```

Here B is the 14-vertex suffix when C has its actual five vertices. Thus an
outer transitive prefix of any size `s0≥1` gives

```
H_s0(C) = 53h + (124+3280·2^s0)F - 487·2^s0 q.       (S2)
```

The current source has `s0=2`, hence coefficients `(53,0,-1948,13244)`.
The coefficient computation is exact subset DP on the fixed exterior, not
a fit to example cores. S1 establishes its validity for all cores. The B
endpoint row also follows identically by subtracting H(B) from the outer
two-prefix expression and dividing by four, using the inherited zero-ready
identity. No dependence on the particular five-vertex C is introduced.

Replacing C arbitrarily can make B non-strong. The counting identity still
holds: for P=`TT_s0 -> B`, adjoining the root beaten only by w gives
`H(Z)=H(P)+E_w(P,2)-start_w(P)`. The nonempty prefix forces
`start_w(P)=0`, while `H(P)=H(B)` and the ordered-slot factorization gives
`E_w(P,2)=2^s0 E_w(B,2)`. Thus S2 is valid for arbitrary B. Strongness is
needed for the inherited canonical zero-ready classification, not for this
algebraic count extension.

Even `(h,F)` is insufficient in this context. A rooted directed triangle has
`(h,s,q,F)=(3,1,1,3)` and gives H=37943. A dominating root above a directed
triangle has `(3,3,3,3)` and gives H=34047. Both have `(h,F)=(3,3)`. This is
a hostile test of the quotient, not a proposed replacement of the literal
source. The latter core is not strong; S1 deliberately allows it.

## 3. Neither desired child occurs anywhere in the S2 family

For every nonempty marked tournament C,

```
1 <= q <= h <= F,          h and q are odd.             (S3)
```

For `q≤h`, insert w into each Hamilton path of C-w: prepend if w beats the
first vertex; otherwise insert just before the first vertex beaten by w, or
append if there is none. Deleting w recovers the original path, so this is
injective. For `h≤F`, split a Hamilton path immediately after w. Its prefix
ends at w and its suffix may be empty, and concatenation recovers the path.
Positivity and oddness are the classical Rédei statement recorded in
**THM-001**, `01-canon/theorems/THM-001-redei.md`; empty deletion has q=1.

Because `487·2^s0>53`, equations S2 and S3 imply

```
H_s0(C) >= (177+2793·2^s0)F.                            (S4)
```

For either N=512423 or N=512425, this proves `1≤s0≤7` and bounds F as follows:

| s0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| maximum F | 88 | 45 | 22 | 11 | 5 | 2 | 1 |

For every pair `1≤h≤F` in these boxes, the equation uniquely determines

```
q = [53h+(124+3280·2^s0)F-N] / (487·2^s0).
```

There are 5289 pairs per N. The exact arithmetic receipt checks all of them,
using divisibility before comparison. At N=512423 none yields an integer
`1≤q≤h`. At N=512425 the only such quadruple is

```
(s0,h,F,q)=(1,73,79,20),
```

which violates the necessary oddness of q. This is a bounded arithmetic
certificate after an all-order reduction; **no bound on |C| was assumed**.

There is a shorter fixed-s0=2 check. The bounds force `39≤F≤45`, and

```
h = 1801(N-13244F) mod1948,       53·1801 = 1 mod1948.
```

The least nonnegative residues are respectively
`1687,555,1371,239,1055,1871,739` and
`1393,261,1077,1893,761,1577,445`. Every residue exceeds its F, so S3 excludes
both targets without using parity in this fixed-prefix case.

This exclusion is limited to S2's fixed four inner wrappers and their fixed
incidences. It does not exclude changes to those wrappers, nontransitive
outer prefixes, new exit edges, or the actual successful literal switches.

## 4. The new exit has a measurable missing port

The first literal +2 switch, branch0 in `literal/realized_depth3.json`, changes
named edges `{4,31},{6,22},{39,31},{39,4}` and uses outer prefix
`50→51→55→56`. Its bottom-core exits are now

```
6 -> {22},       31 -> {4,22}.
```

Reverse just the bottom arc `6→8` in the source. This is a root31-preserving
relabeling of the core, so **every one-port observer at31**, not just S1's four
coordinates, sees the same rooted core isomorphism type. Both complete sources
have H=170807 and `(h,s,q,F)=(11,1,1,13)`.

Apply the same four named raw edge changes and the same prefix extension to
both. Independent vertex-level DP gives

```
original source -> 512423,
6/8 hostile     -> 654021.
```

The first additional named joint field already changes:

```
J_C(31,6) = 7 versus 3,       J_C(6,31) = 1 versus 1,
```

where J counts ordered disjoint spanning two-slot covers with first slot
ending at the first name and second starting at the second name. The physical
port6 profile `(h,start6,q6,F6)` changes from `(11,4,5,12)` to `(11,3,5,17)`.
Thus forgetting port6 is a demonstrated failure, rather than a conjectural
need for more data. This comparison concerns raw named edge and prefix
operations. The old native Hamilton-path guards are **not** transplanted to
the hostile source; such guards would have to be compiled again.

## 5. The exact multi-exit successor interface

Let D be the set of core vertices that beat at least one exterior vertex,
with all exterior vertices dominating C\D. Every nonfinal maximal C-run in
a Hamilton path ends at a different vertex of D. Consequently the number
of C-runs is at most `|D|+1`.

There is an exact finite interface for a fixed D: keep ordered spanning
path covers of C with up to `|D|+1` nonempty slots, recording each slot's
first and last vertex as its physical name in D or the symbol `*` for C\D.
Let `R_C(sigma)` count such covers for an ordered boundary signature sigma.
Let `K_E(sigma)` count ordered exterior path slots before, between and after
those core slots, with nonempty slots between successive core slots, and
the actual allowed named entry/exit edges. Initial and final exterior slots
may be empty. Then

```
H(C inside E) = sum_sigma K_E(sigma) R_C(sigma).         (S5)
```

Proof: contract each maximal C-run and each intervening E-run of an actual
Hamilton path. The remaining records are uniquely an ordered core cover,
an ordered exterior cover, and their boundary signature. Conversely those
records concatenate uniquely whenever their boundary edges are allowed.
No core permutation factorial is inserted: the slots are already ordered.
This is THM-3121's maximal-run mechanism with physical port decorations;
uniform module substitution alone would lose these decorations.

Equivalently enumerate an exterior permutation and choose gaps, including
the two boundary gaps, for the nonempty C slots. At an interior gap the
slot's start must be beaten by the preceding E vertex and its end must beat
the following E vertex. Exterior adjacencies are required exactly in the
unselected gaps. Summing the corresponding ordered core-cover responses
gives S5 without assuming a Hamilton order on E.

The literal bounds are sharp. Verified full named Hamilton paths in the
receipt have the following maximal C-runs:

| Branch | Core exit ports D | Consecutive C-run lists |
|---|---|---|
| +2 branch0 | 6,31 | [6]; [7,31]; [8,14] |
| +2 branch1, no permanent birth | 6,7,31 | [6]; [7]; [14,31]; [8] |
| +4 branch4 | 14,31 | [6,7,31]; [14]; [8] |

These examples certify actual three- and four-run contributions. Alone they
would not prove that every possible lower-dimensional identity fails; the
separate named-port hostile proves a concrete failure of one-port closure.

The incoming exact endpoint matrix `G=(I-ZA)^(-1)Z` at
`05-knowledge/results/session180_170807_returns_20261104/ports/PROOF.md`
(read from origin/main commit 838b01cec4) supplies every S5 joint field by
squarefree products, and updates it under the actual named edge changes.
It is a sufficient growing sidecar, not a bounded-size compression. The
literal/phase packets own native guards, F49/F343 recomputation and the
81/27 control; S1 does not identify or transport their separate phases.

## 6. A testable replacement for scalar return matching

A candidate family should specify its named core operation, exterior state,
and joint response map. If a returned response satisfies `R'=M R+v`, and its
exterior compiler is K', the reusable carry condition is

```
K'M = 3K,          K'v = c,          c in {2,4}.         (S6)
```

These are sufficient coefficient identities. They are necessary on a
full affine span of admissible profiles; on a restricted family they need
hold only modulo the affine relations of that family. The returned named
port language and legal guards remain separate requirements. Testing S6
on profiles spanning that affine hull gives a decisive finite linear test;
testing it at the actual single `(11,1,1,13)` does not.

The strongest next target is an exterior-changing transform that controls
the new joint fields and returns to a specified multi-exit language. The
present packet proves why arbitrary bottom-core growth inside the old
exterior cannot substitute for that task.

## Reproduction

```
python 05-knowledge/results/zero_ready_strong_switch_20261103/structure/check_structure.py
python -O 05-knowledge/results/zero_ready_strong_switch_20261103/structure/check_structure.py --output 05-knowledge/results/zero_ready_strong_switch_20261103/structure/verification_optimized.json
```

The checker uses explicit failures, not Python assertions. It checks all
2175 combinations of labeled cores and exteriors through order3, every
core root and exterior exit subset, including singleton and empty boundaries;
scans 5405 marked cores through order5; independently counts the literal
source, the same-profile hostile sources, and the two diverging target
graphs; checks the nested components and all named maximal-run paths; and
performs the complete finite arithmetic reduction. All universes and inherited
filters are stated above. JSON receipts pin the actual input and checker
bytes. No orbit enumeration, isomorphism filter, or unproved spectrum
exclusion is used.
