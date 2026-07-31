---
id: THM-2825
title: "Nearest half-step common/right collar and semantic parity boundary"
status: >
  PROVED + VERIFIED-EXACT; CORE COLLAR INDEPENDENTLY HOSTILE-AUDITED.  On the
  complete 567-cell full-semantic common/right bank, exactly 193 cells are
  nonempty.  Their unscaled equal-copy relation is a 193-block disjoint union
  of complete bipartite graphs with 63,308 common pieces, 587 right pieces,
  and 195,517 edges.  Circular half-step distance selects a unique cellwise
  collision-free right-to-common collar at +h; delayed coefficient equality
  is exactly even displacement parity, so +h reverses all 587 semantic values
  while +2h is the unique nearest same-value collar.  The parity symbol is a
  norm-one Schur tripotent, but the metric collars are absent from the coarse
  cell/colour algebra.  The rooted half-step forest factors +2h as a
  non-idempotent transverse/tangent operator and yields an exact
  M_3 tensor I_587 linking algebra with a punctured V_4 grading.  Native
  factors and the source carrier still fail uniformly and endpoint data do
  not repair the physical typing.  No common allocated atom, row exclusion,
  or LRC(14) conclusion follows.
source: root/nearest-half-step-common-right-collar-2026-07-28
depends_on:
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2818-right-cofiber-positive-copy-stratification-and-alternating-half-step-chains
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
related:
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2814-projective-allocation-square-holonomy-and-idempotent-provenance-no-go
  - THM-2819-nonwrapping-endpoint-conjugacy-odometer-wrap-death-and-sharp-target-eleven-face
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - HYP-9046-signed-idempotent-decomposition-for-the-lrc-wall
script: 04-computation/lrc14_nearest_half_step_common_right_collar_thm2825.py
output: 05-knowledge/results/lrc14_nearest_half_step_common_right_collar_thm2825.out
script_sha256: bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a
output_sha256: c4a31e5ee0aa5af69faa3efbe315d0968a85cba49c2d77c0ca93a229bc39fa0c
secondary_script: 04-computation/lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.py
secondary_output: 05-knowledge/results/lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.out
secondary_script_sha256: af924103c30a45eeab88520d6b569c00a5a8c9cafbc8052b3585dd86e07b5dac
secondary_output_sha256: e610c5bffb720e074662f2222a50b2ce461c1b1293e946aa260faf898c4347b7
independent_script: 04-computation/lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.py
independent_output: 05-knowledge/results/lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.out
independent_script_sha256: 5a74cfa0b9a5f30e15460e7e1bee17214193d1f307e79f1eeb68b8781b7d23d6
independent_output_sha256: b9846e35c615b6ff0183f6df47ba8e8f99796a195afb6480f922254e35cbe82d
schur_script: 04-computation/lrc14_nearest_half_step_common_right_collar_schur_audit_thm2825.py
schur_output: 05-knowledge/results/lrc14_nearest_half_step_common_right_collar_schur_audit_thm2825.out
schur_script_sha256: 90728bcd82a0552742f29fb34c9754c1a454c4b3ebc177121679ae2d7ec5b83a
schur_output_sha256: 6f90f04c1efa0a6be5f19bdb5e9a353dcda5e7ebb70cbabb8a34f9d62fc5b1d4
path_script: 04-computation/lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.py
path_output: 05-knowledge/results/lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.out
path_script_sha256: 7f0780e70161cbfafc4499d02a7d1c5aa40366b6dfa9935b00dc652e3d54c8e0
path_output_sha256: 6cff846944c59d8e243e74afe8829046feb288f0cf23da1c998037dcf70411b2
hash_basis: LF-normalized bytes
---

# THM-2825 -- scale selects a two-coloured common/right boundary collar

**PROVED + VERIFIED-EXACT; CORE COLLAR INDEPENDENTLY HOSTILE-AUDITED.**

THM-2818 decomposes the nonzero right cofiber into positive equal interval
copies and identifies a half-step parity chain at target twelve.  That result
uses one source-label section.  It does not compare the right pieces with the
common part throughout the full THM-2806 semantic bank.

The full comparison has two sharply different answers.  If scale is forgotten,
every common piece in a cell is tied with every right piece.  If circular
half-step distance is retained, every right piece has a unique nearest common
mate.  The delayed coefficient is not preserved by this first neighbor: it is
the parity character of the displacement.  The second neighbor is the first
semantic-preserving one.

This produces a canonical thin collar and a two-colour rectangular
decomposition.  It also records why the collar is not yet the typed physical
cospan required at the LRC(14) allocation wall.

## 1. The full semantic bank

Retain THM-2806's seven clocks and take

```text
S={0,1,2,3,8,9,10,11,12},
U={3,4,5,6,7,8,9,10,11}.                              (1)
```

For every labelled cell

```text
c=(e,s,u) in {0,...,6} x S x U                         (2)
```

intersect the inherited weighted source and target carriers with the full
semantic section.  Pull the target carrier back by `-SHIFT`, and write

```text
M_c=source intersect pulled-target,
R_c=pulled-target minus M_c.                            (3)
```

Decompose both objects into maximal weighted half-open intervals.  The theorem
is about these **full-semantic weighted** pieces.  Before this semantic cut,
the raw carrier has 63,895 right pieces and no equal-length common mate.  Thus
the relation below is created by the retained semantic refinement, not by bare
endpoint coincidence.

Of the 567 cells in `(2)`, exactly 193 are nonempty and 374 are empty.
Moreover,

```text
sum_c |M_c|=63308,          sum_c |R_c|=587.            (4)
```

Every nonempty cell consists of equal pieces of length

```text
L=26444880=6h/91,                                      (5)

h=T/(2*13^5)=401080680.                                (6)
```

The weight is `27581135604` at clock one and `27580222516` at clocks two
and three.  No other clock is nonempty.  The complete cell-type census
`(clock,|M_c|,|R_c|; number of cells)` is

| clock | common | right | cells |
|---:|---:|---:|---:|
| 1 | 189 | 3 | 54 |
| 1 | 213 | 2 | 18 |
| 1 | 239 | 2 | 9 |
| 2 | 190 | 3 | 7 |
| 2 | 426 | 4 | 28 |
| 2 | 474 | 2 | 14 |
| 2 | 526 | 2 | 7 |
| 3 | 382 | 5 | 7 |
| 3 | 404 | 4 | 28 |
| 3 | 409 | 3 | 7 |
| 3 | 452 | 2 | 7 |
| 3 | 504 | 2 | 7 |

This table sums independently to `(4)` and to 193 cells.

## 2. Forgetting scale gives complete bipartite ties

For `r in R_c` and `m in M_c`, declare

```text
r E_c m

iff length(r)=length(m), weight(r)=weight(m),
    and left(r)-left(m) lies in h Z.                    (7)
```

Then in every nonempty cell

```text
E_c=R_c x M_c.                                         (8)
```

Hence the global labelled relation is a disjoint union of 193 complete
bipartite graphs.  It has

```text
sum_c |R_c||M_c|=195517                                (9)
```

edges.  The possible degrees on the right are

```text
189,190,213,239,382,404,409,426,452,474,504,526.       (10)
```

Their right-vertex multiplicities are respectively

```text
162,21,36,18,35,112,21,112,14,28,14,14.               (11)
```

Thus the unscaled relation contains no unique mate at all.  It is not an
intrinsic tournament: every pair in a cell is tied by the retained
observable, and orienting those ties would erase the geometry needed below.

## 3. Circular scale selects the one-step collar

The half-step circle has

```text
T/h=2*13^5=742586                                      (12)
```

positions.  For a relation edge define the centred signed displacement

```text
delta(r,m)=(left(r)-left(m))/h in Z/(2*13^5),           (13)
```

represented in the centred interval.  Every right piece has a unique common
piece at minimum `|delta|`, namely

```text
j_1(r)=[left(r)+h,right(r)+h;weight(r)],               (14)

delta(r,j_1(r))=-1.                                    (15)
```

Within each labelled cell, `j_1` is injective.  Consequently its image is a
587-piece thin boundary collar inside the 63,308 labelled common pieces.  The
direction is load-bearing:

```text
j_1:R -> M,                                            (16)
```

not a surjection from the abundant common bank onto the right cofiber.

The collision-free statement is cellwise, or equivalently on the disjoint
labelled union.  It does not claim that erasing `(e,s,u)` leaves 587 distinct
physical intervals.

The direction in `(16)` is sharp.  If instead one asks for a nearest right
piece for each of the 63,308 common pieces, exactly 63,066 have a unique
nearest right piece and 242 have an exact two-way tie.  Thus the reverse
nearest rule is not a function on the full common bank.

## 4. The delayed coefficient is exactly displacement parity

For a piece `I` at clock `e`, let

```text
v_e(I)=(source delayed carry 12,
        translated-target delayed carry 6).             (17)
```

Throughout the bank,

```text
v_e(I) in {(C,C),(0,0)},       C=103478815440.          (18)
```

For all 195,517 relation edges,

```text
v_e(r)=v_e(m)       iff       delta(r,m) is even.       (19)
```

The exact edge census is

```text
even/equal     97661,
odd/opposite   97856,                                  (20)
```

with no even/opposite or odd/equal exception.  In particular, `(14)` reverses
the value on every one of the 587 right pieces:

```text
(C,C) -> (0,0) on 573 pieces,
(0,0) -> (C,C) on 14 pieces.                           (21)
```

The fourteen zero right pieces are seven target-labelled copies of one
physical interval at `(e,s)=(2,3)` and seven target-labelled copies of one
physical interval at `(e,s)=(3,2)`.

The unique nearest common piece with the **same** semantic value is

```text
j_2(r)=[left(r)+2h,right(r)+2h;weight(r)],             (22)

delta(r,j_2(r))=-2.                                    (23)
```

It is again cellwise injective.  Thus

```text
2h=T/13^5=802161360                                    (24)
```

is the semantic period, while `h` is its nontrivial two-colour square root.
The metric collar is therefore a `Z/2`-graded boundary object, not merely a
nearest-neighbor matching.

This is the depth-five analogue of THM-2698's central half-odometer grammar:
translation by `h` represents the nontrivial class of
`C_(2*13^5)/C_(13^5)`, and `(19)` identifies the delayed value with its
nontrivial character.  But `j_1` is only a partial arrow from `R` into `M`,
not an action on a closed free two-state fibre.  It therefore identifies the
missing `C_2` sidecar without constructing THM-2698's still-required semantic
bibundle.

## 5. Rectangular idempotents and the exact Schur reframe

Regard `(8)` as a Boolean matrix with rows the labelled right pieces and
columns the labelled common pieces.  Each cell is one all-one rectangular
block, and different cells share neither a row nor a column.  Therefore the
whole relation is a contractive idempotent Schur multiplier.

Partition its edges into

```text
E_equal={(r,m):v(r)=v(m)},
E_opposite={(r,m):v(r)!=v(m)}.                          (25)
```

Within a cell, each of these is a union of the rectangles obtained by
partitioning the rows and columns by the two values in `(18)`.  These
rectangles again share neither row nor column.  Hence both indicators in
`(25)` are contractive idempotent Schur multipliers, and the parity sign is
the exact two-term decomposition

```text
(-1)^delta = 1_(E_equal)-1_(E_opposite).                (26)
```

The nearest collars `j_1` and `j_2` are matching sub-idempotents of
`E_opposite` and `E_equal`, respectively.

There is a sharper norm statement.  Give every row and column the sign of its
value in `(18)`.  On the relation, `(19)` factorizes the parity symbol as

```text
S(r,m)=(-1)^delta=epsilon_R(r) epsilon_M(m).             (27)
```

Consequently its Schur action is the rectangular projection conjugated by
diagonal sign unitaries:

```text
M_S(X)=D_R M_E(X) D_M.
```

It has Schur norm exactly one, not merely the norm-at-most-two estimate from
`(26)`.  Pointwise,

```text
S^2=1_E,                  S^3=S.                         (28)
```

This is a proved finite instance of the blocky signed-decomposition language
suggested by HYP-9046 and the Beke--Goh--Hatami--Jaffe--Naylor theorem.  It is
not a proof of HYP-9046's uniform Schur-norm obligation for the large-diameter
LRC wall.  Here the norm-one block structure is already visible; what it
forgets is precisely the metric scale that distinguishes `(14)` from the
other 195,517 ties.

### 5a. Exact offsets, the coarse-algebra no-go, and the common-side torsor

Resolve every relation edge by the forward offset

```text
k=(left(m)-left(r))/h=-delta(r,m).
```

There are exactly 1,043 nonempty offset masks: 461 negative and 582
positive, with range `-560,...,582`.  Every offset mask is a partial matching,
so each is a norm-one Schur idempotent, and distinct offsets are
Schur-orthogonal.  Exactly two offsets have all 587 rows in their domain:

```text
k=1 gives j_1,              k=2 gives j_2.
```

Every other offset supports at most 527 rows.  Thus, once the offset
resolution is retained, domain saturation and parity select the two collars
without a nearest-distance comparison.

This does not make scale recoverable from the coarse Schur algebra.  The
algebra generated by labelled cells and the two semantic colours has 414
rectangular atoms.  Each collar strictly cuts all 207 atoms of its parity.
Therefore no Schur polynomial in the cell/colour masks can equal either
collar.  Contractivity and semantic parity are exactly blind to the metric
coordinate that selects them.

The images of `j_1` and `j_2` are disjoint 587-piece common collars.  Pairing
them over their common right root gives a fixed-point-free involution on
1,174 common pieces; it swaps adjacent pieces and reverses semantic value.
Equivalently, if `J_i` is the row-to-column matching matrix, then ordinary
matrix composition gives

```text
V_ij=J_i^* J_j,            V_ij V_kl=delta_(j,k) V_il.
```

The four `V_ij` generate `M_2 tensor I_587` on this common-side collar.  Its
Pauli generators satisfy `X^2=Z^2=Q` and `XZ=-ZX`.  This is an abstract free
`C_2` torsor only **after** the metric collars have been supplied.  It lives
entirely on the common side and is not the missing physical source-carrier
bibundle.

### 5b. The rooted path forest and a transverse/tangent factorization

In each labelled cell put an arrow `x -> x+h` whenever `x` is a right or
common piece and `x+h` is common.  The resulting graph is a disjoint union
of 685 finite paths:

```text
587 cofiber-rooted paths cover 54,754 common atoms;
 98 common-only paths cover the remaining 8,554 atoms.
```

The rooted common lengths range from two to 288; exactly 60 have length two.
The common-only paths have lengths

```text
13,39,65,91,117,143
```

with multiplicities `14,14,14,14,14,28`.

Let `V` be the free coefficient module on all labelled right and common
pieces, let `N` be the arrow operator, and let `P` project onto the common
pieces.  This is a module on the labelled interval decomposition, not the
original physical allocation module.  Write

```text
d=P N (1-P),                 a=P N P.
```

Then

```text
N=d+a,       rank(d)=587,       rank(a)=62,623,
d^2=d a=0.
```

If `J` is the sign of the delayed semantic value, exact alternation on all
685 paths gives

```text
Jd=-dJ,                       Ja=-aJ.
```

The same-content collar factors rather than teleports:

```text
S=a d=P N^2(1-P).
```

It has rank 587, satisfies `S^2=0`, commutes with `J`, and obeys
`[P,S]=S`.  On the matched right/image summand,
`S^*S=Q_R` and `SS^*=Q_2`; hence `S+S^*` is an entrywise-nonnegative
self-adjoint involution there.  It is not a positive-semidefinite operator
and not an idempotent.

The right root, odd collar, and even collar together give a still cleaner
linking algebra.  On the 587-dimensional right-root space set

```text
W_0=the inclusion into R,       W_1=d,       W_2=S.
```

Their images are pairwise disjoint and

```text
W_i^* W_j=delta_(i,j) I_R.
```

Hence `E_ij=W_i W_j^*` are matrix units generating
`M_3 tensor I_587` on the 1,761-dimensional initial ladder.  Its compressed
half-step operator is

```text
B=E_10+E_21,              B^2=E_20=S,              B^3=0.
```

Relative to each right root, semantic parity is the grading
`diag(1,-1,1)`, while source-carrier absence/presence is
`diag(-1,1,1)`.  Thus the absent source carrier is represented exactly as a
third graded summand, not silently identified with a common atom.  This
`M_3` algebra still depends on the metric-selected labelled collars and does
not descend to a physical allocation action.

Together the two gradings form a punctured `V_4` square:

```text
R :  (semantic,carrier)=(0,0),
M1:  (semantic,carrier)=(1,1),
M2:  (semantic,carrier)=(0,1),
missing corner          =(1,0).
```

The arrows have the three nonzero degrees

```text
deg(d)=(1,1),            deg(a)=(1,0),            deg(ad)=(0,1),
```

and their `V_4` addition is exactly the factorization `ad=S`.  The fourth
corner cannot be manufactured by relabelling the present right bank: its
semantic populations are 573 live and 14 zero, so no semantic-flipping
permutation of the 587 right roots exists.  Completing the square therefore
requires an enlarged cofiber or signed multiplicity even before endpoint and
native-factor typing are imposed.

The path forest is sharply nilpotent:

```text
N^289=0,                     a^288=0.
```

Moreover `a^(k-1)d` is translation by `kh` exactly on roots whose path has
length at least `k`.  In particular,

```text
rank(a^13 d)=527,
```

while the other 60 roots die after the `+2h` collar.  This is precisely the
`527/60` direct `+14h` boundary: on the long paths the odometer displacement
factors as the first transverse `+2h` collar followed by six
common-tangent `+2h` moves, whereas the short paths terminate.  It is not a
composition of seven typed `j_2:R->M` arrows.

The 98 common-only paths carry an additional signed prime-power unit.  For a
path of length `13m`, with odd `m in {1,3,5,7,9,11}`, reduce its residue word
in `F_13[C_13]`.  With

```text
N_13=sum_(r=0)^12 X^r,       A_13=sum_(r=0)^12 (-1)^r X^r,
```

the raw mask `R` and live-colour mask `L` satisfy

```text
R=m N_13,                    2L-R=A_13,
(1+X)A_13=2,
L^(-1)=(1+X)(1-m N_13).
```

Thus semantic alternation turns every rootless augmentation-zero norm block
into a full-spectrum signed unit.  The inverse is not a positive physical
transport.

## 6. Prime-power unit mass gives full signed ancestry on each colour

Choose the least left endpoint in a cell as origin.  The half-step lattice
then identifies with

```text
C_(2*13^5)=C_2 x C_(13^5),                              (29)
```

where parity is the `C_2` coordinate and `+2h` generates the `13^5` factor.
The exact parity-mass types `(M_even,M_odd;R_even,R_odd)` are:

| clock | common masses | right masses | cells |
|---:|---:|---:|---:|
| 1 | `(94,95)` | `(3,0)` | 54 |
| 1 | `(106,107)` | `(2,0)` | 18 |
| 1 | `(119,120)` | `(2,0)` | 9 |
| 2 | `(95,95)` | `(1,2)` | 7 |
| 2 | `(213,213)` | `(4,0)` | 28 |
| 2 | `(237,237)` | `(2,0)` | 14 |
| 2 | `(263,263)` | `(2,0)` | 7 |
| 3 | `(191,191)` | `(1,4)` | 7 |
| 3 | `(202,202)` | `(4,0)` | 28 |
| 3 | `(205,204)` | `(3,0)` | 7 |
| 3 | `(226,226)` | `(2,0)` | 7 |
| 3 | `(252,252)` | `(2,0)` | 7 |

Every displayed nonzero mass is prime to thirteen.  The piece weights are
also units modulo thirteen: `11` at clock one and `1` at clocks two and three.
THM-2839 therefore applies separately to every nonempty parity mask:

```text
each common colour has all 13^5 ancestry characters;

the sole right colour has all 13^5 characters in 179 cells;

both right colours have all 13^5 characters in the remaining 14 cells.
                                                               (30)
```

Equivalently, every nonempty colour mask has translate rank `13^5=371293`
over `Q` and `F_13`.  This is a signed reconstruction theorem.  Since every
mask here has nonsingleton support except the one-piece colour classes, its
group-algebra inverse is not a nonnegative physical multiplier.

The `C_2` factor records an additional sharp boundary.  A mask
`A_0+sA_1` in `F_13[C_2 x C_(13^5)]` is a unit exactly when both character
components `A_0+A_1` and `A_0-A_1` have nonzero augmentation.  Exact counts
give

```text
R_c is a full doubled-group unit in all 193 cells;

M_c is a full doubled-group unit in 88 cells and fails the sign
component in 105 cells;

M_c union R_c is a full doubled-group unit in all 193 cells.     (31)
```

Thus the right cofiber and pulled target retain full signed spectrum even
where the common mask loses the parity character.  The missing datum is not
ancestry rank; it is positive multiplier access and physical source support.

## 7. Native-factor and carrier boundaries are uniform

Let the six native factors be

```text
(E3,clock,q1,q2,c2,c3).                                (32)
```

Both `j_1(r)` and `j_2(r)` lie in all six factors in all four
native/pulled frames, for every `r`.  No right piece does.  The right pieces
have exactly four hole types:

| hole in the first/fourth frame | pieces |
|---|---:|
| `E3,c2` | 37 |
| `E3` | 319 |
| `q1` | 14 |
| `c2` | 217 |

The two middle frames are complete.  Thus metric translation repairs bare
factor membership only by moving from the right cofiber into the common
bank; it does not prove that the two pieces have one native origin.

The thirteen-twist carrier masks are even sharper.  Write `delta_0` for
presence only at twist zero.  Uniformly over all 587 triples,

```text
R:   (source,target)=(empty,delta_0),
M1:  (source,target)=(delta_0,delta_0),
M2:  (source,target)=(delta_0,delta_0).                 (33)
```

Hence neither collar is a co-supported physical translation of the right
piece.  This is exactly the missing hypothesis isolated by THM-2820's affine
Hasse-jet boundary: `(14)` supplies a numerical displacement, but its source
carrier is absent.

The literal THM-2791 ancestry chambers, by contrast, are preserved by both
`j_1` and `j_2` on all 587 pieces.  This positive sidecar sharpens the
failure: ancestry provenance survives, while native factor and source-carrier
typing do not.

## 8. Endpoint data do not recover the semantic period

The independent endpoint audit checks all 587 triples against the full
`F_13^2` endpoint-address bank.  The target endpoint masks of `r` and either
`j_i(r)` are identical.  On the source endpoint masks the Hamming-distance
census, for **both** `i=1` and `i=2`, is

```text
distance 0: 74,
distance 9: 187,
distance 10:245,
distance 81:81.                                        (34)
```

Only the 74 equal cases admit an address translation, namely zero; the
remaining 513 admit none.  Thus the endpoint-pair census is identical for the
semantic-reversing `h` collar and the semantic-preserving `2h` collar.  A
coarse endpoint signature cannot recover the parity sidecar.

The exact endpoint-profile censuses of `j_1` and `j_2` also agree.  This is
an aggregate equality unless the audit separately proves pairwise equality;
no pairwise claim is made here.

## 9. Relation to the odometer wrap

The physical target shift satisfies

```text
13*SHIFT=5615129520=14h=7(2h).                         (35)
```

So THM-2819's full odometer wrap is seven semantic collar periods.  This
arithmetic identity explains why the support-level wrap and the coefficient
parity can coexist.  But direct translation by `14h` lands back in `M_c` for
only 527 of the 587 right pieces; 60 have no such common landing.  All 527
landings preserve the delayed value, as parity predicts.

Therefore `(35)` does **not** compose seven typed `j_2` arrows, even at the
finite-bank support level.  Equation `(33)`, the native-factor holes, and
`(34)` additionally show that the required intermediate physical maps are
absent.  THM-2819 independently records the wrap death at the typed
endpoint/factor level.

## 10. Exact verification and scope

The primary exact companion reconstructs all 567 cells, verifies `(4)--(24)`,
and uses explicit failures rather than Python assertions.
It pins the THM-2818 dependency by LF-normalized SHA-256 and passes ordinary
and optimized Python.

The physical companion independently retains all native-factor, twist-carrier,
and endpoint masks in `(32)--(34)`.  The hostile audit reconstructs the bank
through a separate lower-level path and checks the relation, metric, semantic,
direction, ancestry, prime-power masses `(29)--(31)`, and wrap landing without
importing the primary script.

The Schur companion verifies all 1,043 exact-offset matchings, the 414 coarse
semantic atoms, the norm-one factorization, the collar-algebra no-go, and the
common-side `M_2` composition.  The path companion verifies all 685 path
components, the transverse/tangent operator identities, the
`M_3 tensor I_587` and punctured-`V_4` laws, sharp nilpotence, every rank
plateau, and the 98 rootless `C_13` inverses.  Both pass ordinary and
optimized Python with byte-identical stored output and no assertion nodes.

What is proved is:

1. the complete bipartite unscaled relation;
2. the unique scale-selected `+h` collar;
3. the global semantic parity law and unique `+2h` same-value collar;
4. the norm-one Schur tripotent, exact-offset resolution, and coarse-algebra
   no-go;
5. the rooted path factorization, `M_3` linking algebra, punctured `V_4`
   obstruction, and rootless prime-power unit;
6. the full signed prime-power ancestry spectrum on every nonempty colour;
   and
7. the exact factor, carrier, endpoint, reverse-direction, and wrap
   boundaries.

What is **not** proved is a common allocated atom, a common endpoint origin,
a lawful equivariant translation, a row exclusion, the missing THM-2772
incidence, or LRC(14).  Section 5b does construct the requested
non-idempotent, gauge-transverse **labelled coefficient** sidecar.  The
remaining obligation is its physical descendant retaining the source
endpoint and native factors.  The cheapest next test is to decorate the 587
rooted paths by their source endpoint-defect classes
`0:74,9:187,10:245,81:81`, ask whether those classes are lawful path
coboundaries, and compare any surviving potential with THM-2847's physical
q3/q11 transverse target edge rather than with another untyped interval
translation.
