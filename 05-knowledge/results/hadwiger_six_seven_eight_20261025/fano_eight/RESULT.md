# Eight frames retain what the seven-point completion forgets

2026-10-25. **PROVED in the stated scopes; [independent audit PASS](../audits/PALETTE_AUDIT.md).**
The affine-code and graph-cohomology ingredients are classical;
the claims below are elementary deductions in the explicitly marked palette
setting. **FINITE-EXACT** controls supplement the proofs. No unit-distance
realization, six-color decision, or Hadwiger–Nelson lower bound is asserted.

## 1. Inheritance and the missing operation

[Session61's palette-deck packet](../../session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md)
already proves the six-label pairing / fifteen Lagrangians / two Fano
completions correspondence. Its section 3 proves that fixed ordered-pair
transport loses `2 b1` obstruction bits when only completion parity is
kept, and explicitly leaves pair permutations to a local-coefficient
extension. During this session the independently audited concurrent
[Session62 twisted-pair packet](../../session62_agent_theme_bridges_20261025/holonomy/TWISTED_PAIR_LIFTS.md)
arrived at `origin/main` in checkpoints `9807ae1678` / `9ef7014b8f`. It
already proves the exact local-coefficient extension, its section counts,
and the endpoint/conic cube comparison. Sections 4 and the corresponding
cube interpretation below are therefore **recovered concurrent results**,
not new claims of this packet. Our additional explicit structures are
the `C2 x S4` tetrahedral reduction, the affine-design/basis sidecar of
index 28, and the common-fixed-point witness with complete finite controls.
The [older Fano-code packet](../../arithmetic_braids2_20260917_fano_code.md)
already constructs the extended Hamming code as affine evaluations; the
code name and construction are recovered, not new. The
[root-deletion cube](../../six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md)
and [free simplex/conic action](../../six_seven_geometry_20261025/fano_codes/RESULT.md#6a-the-same-subgroup-acts-regularly-on-the-eight-central-k3-conics)
provide a separate geometric realization of the eight-state translation
action. We retain that distinction.

Closest proved mechanism: graph gain equations plus the marked palette
pairing. Canonical hostile: every scalar parity passes while even frame
transport does not close. Corrected near miss: each generating loop may
have a fixed frame but the loops may have no *common* fixed frame. The
least-used sidecar is pair-permutation monodromy. It can remove some
obstructions rather than merely add complexity.

Anchor: make the six-palette lift faithful. Niche: the eight affine-frame
geometry. Wildcard: a three-cycle of pair labels makes the scalar parity
test complete on a single cycle. The live board is:

| Object | Operation | Preserved / missing data |
|---|---|---|
| Six paired colors | choose one endpoint per pair | eight endpoint frames; a pairing is still required |
| Two Fano completions | split frames by total parity | loses the four-state even component |
| Eight affine points | retain all affine hyperplanes | seven direction classes, but no chosen cube basis |
| Unordered pair seams | allow `S3` channel permutations | actual monodromy modifies obstruction dimension |
| Full frame gluing | solve one common affine fixed-point problem | separate loop solvability is weaker |
| K3 conic cube | compare its marked translations | no palette proper-coloring predicate transported |

## 2. Six facets, seven directions, and eight endpoint frames

Fix three pairs of named colors `P_j={a_j,b_j}`, `j=1,2,3`. A frame chooses
one endpoint of each pair, so the frame set is the affine cube `A=F2^3`:
the j-th bit is zero for `a_j` and one for `b_j`. Each individual palette
color identifies one of the six coordinate facets `{x:x_j=b}` of this
eight-point cube. This is a literal incidence map: a frame lies on the
three facets for its chosen colors.

Add a seventh named label `*`. The two Fano completions of the fixed
pairing have three common lines `(*,a_j,b_j)` and respectively the four
transversal lines whose endpoint frames satisfy

```
p(x)=x_1+x_2+x_3=eta,       eta=0 or 1.                    (F1)
```

Thus the eight frame points are exactly the eight off-central transversal
triples distributed between the two completions. The completions and
the frame cube have different vertex sets; equation (F1) is the actual
correspondence between them.

For every nonzero linear functional `n` on `A`, the pair of affine
hyperplanes `n(x)=0,1` is one parallel class. There are seven such classes
and fourteen four-point blocks. The sixteen affine binary functions

```
{ x -> n(x)+b : n in (F2^3)^*, b in F2 }
```

form the length-eight extended Hamming / `RM(1,3)` code with weights
`1+14 z^4+z^8`. This is the inherited affine-evaluation construction.

The fourteen nonconstant supports also form a Steiner system `S(3,4,8)`.
Indeed three distinct affine binary points `u,v,w` determine the unique
four-point plane `{u,v,w,u+v+w}`. Any affine line has only two points,
so the three are affinely independent. This proves both existence and
uniqueness of the block through each triple. Choose an origin, and the
seven blocks through it, with that origin removed, form a Fano plane.
The seven parallel classes are instead the *dual* Fano direction set.
An identification of point and direction Fano planes requires a chosen
duality; it is not built into the design.

The six colors single out three independent direction classes from these
seven. There are `168/6=28` unordered bases of the direction space.
Therefore the affine design is not a faithful carrier of the palette
unless those three directions are retained. Its full automorphism group
is `AGL(3,2)`, of order `8*168=1344`. The stabilizer of the chosen three
directions is the palette group

```
W = F2^3 semidirect S3,        |W|=48,                    (F2)
```

acting by independent endpoint swaps and permutations of the three pairs.
To verify the full automorphism assertion without a citation: after
composing with a translation, a design automorphism fixes zero. The
unique block through `0,u,v` has fourth point `u+v`, so such an
automorphism preserves binary addition and is linear. The converse is
immediate. The checker independently exhausts all `8!` point permutations.

**Hostile to forgetting the basis.** The shear
`(x_1,x_2,x_3)->(x_1+x_2,x_2,x_3)` preserves all fourteen blocks, but sends
the cube edge `000--010` to `000--110`, a nonedge. Hence an affine-design
identification alone does not preserve cube adjacency, its K3 component
intersection matrix, or the designated six palette facets.

## 3. The hidden four-state object is a tetrahedron

Let `c=111` and let `E=ker p={000,011,101,110}`. Since there are three
coordinates, `p(c)=1`, and there is an `S3`-equivariant direct sum

```
A = F2 c direct-sum E,
x = p(x)c + (x+p(x)c).                                   (F3)
```

For an affine palette transport `g(x)=pi(x)+a`, decompose
`a=p(a)c+a_E`. Then its action is

```
(eta,e) -> (eta+p(a), pi(e)+a_E).                         (F4)
```

Coordinate permutations act faithfully as all `GL(2,2)` transformations
on the two-dimensional `E`: they permute its three nonzero points.
Translations by `E` give every affine permutation of its four points.
There are `4*6=24`, exactly `S4`. Consequently (F4) is an explicit group
isomorphism

```
W = C2 x S4.                                             (F5)
```

The first factor acts on the two Fano completions. The second acts on
the four choices of full frame *within* a specified completion. Thus
the missing coordinate is a marked tetrahedron, not just an unspecified
pair of binary bits. This description remains valid when pair labels
are transported nontrivially around cycles.

No new Euclidean tetrahedron is claimed: these are four endpoint frames.
An injective unit-distance `K4` in the Euclidean plane is not supplied.
The tetrahedron's role here is the permutation action on four objects.

## 4. Unordered-pair transport: the exact local-coefficient theorem

Let `Gamma` be a finite connected gluing graph, considered as a
one-dimensional complex. Choose an orientation for each edge. A seam
`e:u->v` carries `pi_e in S3` and `a_e in A`, with frame equation

```
x_v = pi_e x_u + a_e.                                   (F6)
```

Reversing the edge uses `pi_e^-1` and `pi_e^-1 a_e`. The permutations
define a rank-three local coefficient system `A_rho`, whose rank-two
even subsystem is `E_rho`. There are no two-cells and hence no extra
flatness requirement on the assigned edge permutations. A vertex frame
change adds a local-coefficient coboundary. Complete compatible endpoint
frames exist precisely when

```
[a]=0 in H^1(Gamma; A_rho).                              (F7)
```

This asks for a section of the bundle of endpoint choices. It does not
ask for a global ordering or identification of all six palette labels:
that stronger trivialization target requires the prescribed permutation
holonomy itself to vanish. A nontrivial three-cycle can have compatible
endpoint choices while preventing such a label trivialization, as the
incoming Session62 packet also emphasizes.

The equivariant split (F3) induces the split exact sequence

```
0 -> H^1(Gamma; E_rho) -> H^1(Gamma; A_rho)
  -> H^1(Gamma; F2) -> 0.                               (F8)
```

The last map is total parity, and its vanishing is exactly compatibility
of the completion choices. Put `b1=|edges|-|vertices|+1` and
`epsilon=dim(E^M)`, where `M<=S3` is the group of permutation holonomies
around loops, based at one vertex. Then

```
dim H^1(Gamma; E_rho) = 2(b1-1)+epsilon.                 (F9)
```

**Proof.** For any rank-r local system on a graph, the cochain dimensions
are `r|V|` and `r|E|`. There is no `C^2`. Rank-nullity gives
`dim H^1=r(|E|-|V|)+dim H^0`; parallel sections are exactly the vectors
fixed by the monodromy group. Apply this to rank two. The sequence (F8)
is split because the underlying coefficient sequence already splits
equivariantly by (F3). This proves the theorem for arbitrary graph size,
not just the bounded receipt. QED.

The possible invariant dimensions are elementary:

| Pair-permutation holonomy M | epsilon | Hidden obstruction dimension |
|---|---:|---:|
| trivial | 2 | `2 b1` |
| generated by one transposition | 1 | `2 b1-1` |
| contains a three-cycle (`C3` or `S3`) | 0 | `2 b1-2` |

The latter rows require enough cycles to generate that group. A tree has
trivial monodromy, `b1=0`, and no obstruction. When a compatible full
frame exists, its set is a torsor for `H^0(Gamma;A_rho)` and contains
exactly `2^(1+epsilon)` frames.

This recovers Session61's `2 b1` formula when the pair labels do not
permute, and restates the incoming Session62 unordered-pair theorem. Pair
permutations can *reduce* the hidden obstruction: on a connected
unicyclic graph with a three-cycle as permutation holonomy, (F9) is zero.
There the completion-parity condition is necessary and sufficient for
a full endpoint-frame lift. For identity or transposition holonomy,
two or one hidden bits remain, respectively.

## 5. Individually soluble loops need not share a frame

Choose a spanning tree, propagate a root frame along it, and compose
the affine transports around the remaining `b1` loops. Equations (F6)
are soluble iff these affine loop maps have a *common* fixed point
in the eight-frame cube. Once completion parity closes, (F5) reduces
the remaining question to a common fixed point of permutations of
the four-point tetrahedron.

Here is a small hostile with two loops. Use integer masks for bits and
let `pi` cycle bit positions `0->1->2->0`. Take

```
g_1(x)=pi(x),            g_2(x)=pi(x)+101.
Fix(g_1)={000,111},      Fix(g_2)={011,100}.              (F10)
```

Both translations have even parity. Both affine maps have order three
and individually admit two frames, one of each completion parity.
Their fixed sets are disjoint. Hence both completion tests and both
separate loop-solvability tests pass while the full gluing fails.
The generated affine group contains the nonzero even translation
`g_2 g_1^-1`, which is fixed-point-free. The missing test is common
root-frame compatibility, not another scalar parity.

The bouquet of two loops is merely a compact presentation. Subdivide
each loop into a triangle, assigning identity transport on the first
two edges and its stated gain on the third, to obtain the same witness
on a simple graph. This is an exact frame-constraint example, not an
assertion that arbitrary such seams occur in a unit-distance assembly.

## 6. What reaches coloring, Fano geometry, and the K3 cube

The six-color compiler in Session61 first fixes an existential ternary
assignment and its same-label graph. Our theorem solves the more general
*specified frame-transport system* once an actual interface construction
has produced equations (F6). It does not prove that every graph-coloring
boundary relation is affine, or replace that existential ternary choice.

The eight-frame construction supplies a concrete common carrier:

| Source -> target | Exact map / preserved predicate | Loss and needed sidecar |
|---|---|---|
| Six paired colors -> eight frames | frame selects one color from each pair | keep pairing; frame is not a vertex color |
| Eight frames -> two Fano completions | parity (F1) | loses tetrahedron coordinate; no equivariant full-frame selector |
| Palette frames -> affine `S(3,4,8)` | all affine hyperplane supports | full affine group loses three selected directions |
| Pair seams -> `C2 x S4` | equivariant split (F4) | parity alone loses common tetrahedron fixed point |
| Unordered seams -> graph cohomology | local coboundary equation (F6) | retain linear monodromy, not just cycle rank |
| K3 conic cube -> frame cube | mark its three sign directions | preserves translations and Walsh characters; no coloring transfer |

The prior K3 packet has seven nonconstant Walsh modes, with its three
singleton modes actual divisor relations and its other four negative
directions. The selected three cube directions distinguish those modes.
The affine-design shear above destroys this distinction while preserving
the fourteen blocks. Hence the larger eight-point symmetry is useful
but needs the marked basis to preserve the geometric intersection form.

There is no contradiction between the two warnings: the full affine
design may forget *too much basis data*, while the parity completion may
forget *too much affine translation data*. The eight frames retain both
once the palette pairing and its three directions are specified.

## 7. Reproduction and exact scope

From the repository root:

```
python -B 05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/check.py 05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/receipt.json
python -B -O 05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/check.py 05-knowledge/results/hadwiger_six_seven_eight_20261025/fano_eight/receipt.optimized.json
```

The outputs agree byte for byte. The source imports no producer/repository
implementation and uses explicit exceptions. It checks all 48 palette
transports on both Fano completions and the tetrahedron; all 40,320
permutations of the eight points against the design and cube; all 1,344
constructed affine maps; all 28 unordered direction bases; all
`1+48+48^2=2353` normalized bouquet transports through cycle rank two;
and all 2,304 transport assignments on a two-edge tree. The linear proof
of (F9) is not extrapolated from this finite universe. The receipt gives
the complete blocks/codewords, section-count census, explicit hostile
maps, and source SHA-256.

Remaining task for Hadwiger–Nelson: obtain such frame constraints from
an exact unit-distance assembly while retaining its full boundary
color relation. Nothing in this packet supplies that assembly or changes
the plane's chromatic bounds.
