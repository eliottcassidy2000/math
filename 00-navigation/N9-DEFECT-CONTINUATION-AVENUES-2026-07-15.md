# The `n=9` defect continuation frontier

**Status:** working atlas after THM-828/832/834/838/842/843/846/853,
codex-2026-07-15-S13.  The exact static join, its linear Cech core, one
centered-CF continuation, its first closed lift/face composites and complete
fixed-size identity-adjoined return monoid, and the
collision-fibre-to-node placement are complete.  This file records the
object dictionary, preservation ledger, and ranked pull list for continuation
and LRC-facing agents.  It is deliberately explicit about which statements
concern tilings, tournament nodes, complement lines, defect sectors, or LRC
metric packets.

## The object we actually found

The raw lower codec on the `2^27` canonical apex-zero size-nine tilings has 58
disjoint doubletons and no larger fibres.  Every doubleton is
`{u,sigma u}`.  Its difference

```text
D=(1+sigma)u
```

passes through four successive ambient objects:

```text
F_2^28
  -> raw-S2 syndrome code, dimension 14
  -> reflection-fixed syndrome code, dimension 6
  -> span of realized sectors, dimension 4
  -> 11 occupied nonzero points carrying 58 literal pair fibres.
```

The four-dimensional space has basis

```text
0192486, 08c2c0c, 11b4600, 4483414.
```

Four of its fifteen nonzero points have empty final fibres.  One is killed by
the size-eight face-`H` support and three pass through B compatibility but fail
literal raw-S2 realization.  The coordinate graph on the eleven occupied
points is connected with 14 edges and diameter five.  These edges are
syndrome-coordinate flips, **not** merged-metagraph edges.

Every face, every pair overlap, and the ten-cell triple overlap retains all
four coordinates.  The associated constant Cech complex has cohomology

```text
(dim H0,dim H1,dim H2)=(4,0,0).
```

Thus the four-dimensional object is a holonomy-free linear core decorated by
nonlinear witness fibres.  The linear core glues perfectly; the `H_8` kernel
and literal histogram decide which bases exist.

## Three sorts, then two more

The merged metagraph is not one graph.  At fixed size it contains at least
three mathematical sorts:

```text
tilings X_n  --q-->  converse-merged tournament nodes M_n
    |
    v
complement lines L_n = X_n/<kappa>  --boundary--> Sym^2(M_n).
```

THM-828/832 add two further sorts:

```text
raw-address cells K_n = X_n/~_Lambda,
decorated defect sectors (D,F_D),  F_D=set of literal equal-address pairs.
```

Here `u ~_Lambda v` means `Lambda_n(u)=Lambda_n(v)`; this is a nonlinear
kernel-pair quotient, not the kernel of a linear map.

Nodes classify tournaments and forget the observer Hamiltonian path.  Lines
retain two tiling endpoints and colour.  Raw-address cells compare recursive
face observations.  A sector groups several cells by a common reflection
difference.  None can silently substitute for another.

## Exact maps currently available

1. **Tiling to node.**  Build the fixed-path tournament, canonically label its
   ordinary isomorphism class, and merge its converse class.  THM-834 gives
   the exact subset-DP function on the 116 size-nine collision endpoints;
   THM-828's unsigned-16 atlas supplies it on all `2^21` size-eight tilings.
2. **Node to tilings.**  Generate directed Hamiltonian paths modulo
   automorphisms for each ordinary class in a merged node, or scan an exact
   atlas.  This is a genuine inverse fibre, not one canonical representative.
3. **Tiling to raw address.**  Evaluate the three face observations
   `H_(n-1)`, upper colour, and raw reflection-layer histograms.
4. **Collision pair to defect sector.**  Compute `D=u xor v`; read four pivot
   bits in any A/B/C intersection; apply THM-832's four-bit decoder.
5. **Sector to collision fibre.**  Map the four-bit coordinate back to `D`
   and filter THM-828's exact witness table.  Sector masses are
   `2,2,2,2,2,4,4,4,4,6,26` in sorted order.
6. **Collision endpoint to chirality.**  On the 116 endpoints only, take the
   sign of the first nonzero layerwise skew moment `J_tau`.  This one bit
   separates all 58 doubletons; define a default before using it elsewhere.
7. **Collision pair to a sufficient node coordinate.**  The marginal merged
   node has 53 values and its complement-partner marginal has 54.  Their
   ordered coordinate `P^->=(nu(u),nu(u xor FULL))`, and also its unordered
   version `bar P`, has 58 singleton values.
8. **Centered-CF copy.**  THM-838 derives the literal coordinate copy
   `Phi:X_9->X_10`.  It preserves all 58 `Q` and `bar P` cells, retains defect
   rank four, and exposes exactly the `span(b0,b1)` half to target raw-S2.
9. **B3 continuation.**  THM-842 proves that B gap contraction maps the 58
   cells injectively to 58 lower `Q` and 58 lower `bar P` cells.  The unordered
   A/C endpoint deck has 29 doubletons; adjoining B restores 58 singletons.
10. **Full n8 node-axis shell.**  THM-843 expands the 155 face nodes in the
    exact blue/black complement-line graph.  Black radius one has 3,308/3,528
    nodes; blue re-enters through exactly two mixed portals.  Monotone
    unrestricted and `blue* black*` reach coincide on 3,069 nodes and contain
    every one of the 50 balanced nodes.
11. **First closed CF lift/face words.**  THM-846 computes
    `F_R=R_10 Phi_(9->10)` literally.  A/C jointly omit only the apex; B is an
    idempotent projection omitting the seam `x_64 xor x_73`.  On the 58-cell
    bank, ordered A/C is injective while unordered A/C is exactly the 29
    theta orbits.  Bare node marginals fail different THM-840 operation-
    kernel tests; coupled `bar P` and literal Q pass for B.
12. **Closed-return monoid.**  THM-853 enumerates all 46,126 empty-or-nonempty
    word endomorphisms generated by the three closed returns; 46,125 come
    from nonempty words.  Their Cayley graph has six
    rank-one coordinate sinks and only 144 present Q equality partitions, but
    20,419 future-Q classes.  B is a Q congruence at every state; A/C each
    fail at exactly 76 mirrored states.  The return closure contains 8,940
    masks after complements; its complete node classification remains open.

The missing map is now multiword continuation, not static lookup: given a
decorated sector or node fibre, what decorated fibres can lift, delete, flip,
or undergo an arbitrary centered-CF word?

## Preservation ledger

| carrier | preserves | destroys / cannot decide |
|---|---|---|
| literal tiling | all fixed-path arc bits and its tournament class | LRC metric/owner/wall data and future continuation not supplied |
| ordinary class | the unmarked tournament, distinct from its converse unless self-converse | vertex labels, observer path, complement partner |
| merged node | isomorphism up to converse, fibre mass if weighted | literal path, line endpoint coupling, deletion continuation |
| complement line | two endpoint tilings modulo simultaneous swap, colour | a chosen endpoint orientation |
| raw `Lambda_9` cell | three lower node-pairs, `UABC`, raw S2 | chirality on exactly 58 cells |
| bare `D` sector | reflection difference, four Cech coordinates, all restrictions | base tiling, node tuple, `H_8` truth, S2 realization |
| decorated `(D,F_D)` | sector plus every literal collision witness; exact A/B/C action on the n9 bank | general lift/deletion/CF semigroup beyond audited words |
| curvature tuple | longitudinal C3 flow, `(q0,q1,S,epsilon)` | mirror position; all 58 pairs lie at `epsilon=0` |
| literal `Q=L/<sigma>` orbit | reflection-safe line action | projected node/edge descent can split (THM-813) |
| Bezout-owner stalk | exact `GL_2(Z)` inverse-owner transport | metric image of the LRC witness interval |
| LRC packet | gap metric, owners, wall schedule if retained | no safe quotient is currently proved to the tournament carriers |

## How the old viewpoints now fit

### Kernel pairs and Cech descent

THM-818 replaced node labels by the equality relation
`R_8={(x,y):H_8(x)=H_8(y)}` with literal overlap projections.  THM-828 reduced
the full three-way join to a difference-indexed A/C join; THM-832 shows the
surviving linear carrier is already constant on the entire B3 nerve.  The
lesson is recursive: quotient labels locate equality, kernel rows retain the
witness needed to compose it.

THM-830 makes the three roles intrinsic at
the marked-path level: a chord `(a,b)` survives vertex deletion with role
multiplicities `n-a`, `a-b-1`, and `b-1`, summing to `n-2`.  Its exact replay
through `n=14` reports binary deck rank nine at `n=9` and proves the gap
tournament, half-defect carrier, and deletion/seam calculus used below.

THM-842 applies those roles to the actual defect stalk.  A/C are endpoint
deletions and B is gap contraction.  B is the unique literal-Q-pure role and
is injective on all 58 cells.  The unordered A/C deck has 29 doubletons,
generated by the apex-relative antipode
`theta=sigma(nonapex-complement)`.  Theta preserves `D` and the full skew word,
so this is not chirality loss.  Endpoint deletion converts theta to ordinary
`sigma kappa`, already quotiented by Q; B transports theta itself and retains
the sheet.  This is the first exact fibre-valued meaning of the difference
between the B term and the A/C terms in the seven-term recursion.

### `A+B+C-D-E-F+G`, `A+B-C`, and the odd seven-term mode

The full recursion is the B3 face nerve.  On the defect carrier it reads
exactly `4+4+4-4-4-4+4=4`; the face-exclusive corners are zero and the common
intersection retains everything.  This does not make Venn ranks a nonlinear
codec: identical exact-region weight profiles carry fibre masses 2, 4, and
26.  The even/odd half-tiling recursions should next be applied as operations
on decorated fibres, not scalar identities.

### Black-edge curvature disintegration

THM-811/814/827 show that `q` is the correct conditional flow coordinate.
THM-828 supplies its exact boundary: every false palindrome preserves the full
`(q0,q1,S,epsilon)` tuple and lies at `epsilon=0`.  Curvature controls movement
but cannot orient the two mirror sheets.  Positional chirality is a transverse
coordinate missing precisely on the Smith balance wall.

### Transitivity flow and the wiggly dictionary

THM-810 now gives exact organizing coordinates for ordinary classes and
edges:

- ordinary-class stationary mass is `H/(|Aut| 2^m)`;
- its weighted wiggly degree is `m H/|Aut|`;
- a flip changes the score axis by `4(d_loser-d_winner)+8`;
- margin-two flips are exactly level edges;
- one unit of electrical current crosses every transitivity-axis cut.

Incoming THM-833 sharpens “flow” to an exact OU law.  If
`x=sum_i(2d_i-(n-1))^2`, then a uniformly chosen arc flip satisfies

```text
E[Delta x | T] = -8/(n(n-1)) (x-n(n-1)).
```

The dynamical centre is the random-tournament value `x_*=n(n-1)`, not the
maximally regular score class.  For the size-eight face support in THM-834,
only the first two of eleven sector barycentres lie on the transitive side of
`x_*=56`; nine lie on the overdistributed side, where mean drift reverses.

For a non-self-converse merged node, the corresponding fibre mass and degree
have an extra factor two; a self-converse node does not.  This supplies an
objective node order: score-axis value, score sequence,
blue/mixed/pure-black type, ordinary `H/|Aut|`, then canonical word.  It
also supplies an exact local edge orientation.  THM-834 places the 58 n9
pairs and their 348 distinct face tilings: they occupy 155 n8 nodes and 87
distinct exclusively-black projected complement-line node pairs.  The
selected defect bank therefore has no blue rail; its black `|Delta x|`
histogram is `(0:32,8:52,16:72,24:4,32:13,48:1)`.

THM-843 shows precisely where that local statement stops.  The full black
radius-one shell contains 3,308 nodes; the two mixed seed nodes are the only
blue portals, supplying 25 external blue neighbours.  The coloured union has
3,310 nodes and contains the entire `x<56` distributed tail.  Endpoint
curvature refines the 47,872 incident black projected supports to 58,937
cells, but every fibre is still even under reflection.  Thus the selected
bank is internally black yet globally almost spanning, and curvature still
does not replace mirror position.

### Continued fractions and owner transport

THM-778 says a continued-fraction digit is safe only with its token action.
THM-812 constructs a literal coordinate-copy action; THM-813 proves that
projected node/edge cells cease to be functorial while literal reflection
orbits remain safe.  THM-829 closes inverse-owner transport by
`v'=Av,b'=bA^-1`.  THM-838 carries the rank-four core through the next
centered consecutive word exactly: rank `4->4`, raw-S2 kernel rank two, all
58 literal `Q` orbits and coupled `P` cells distinct.  Static chirality is
not functorial—eight signs reverse—so the next fibre product must retain the
layerwise moment word and owner/Bezout sidecar.

THM-846 closes that lift with A/B/C.  The resulting three kernels are not one
phenomenon: endpoint A/C forgets the apex on the full cube; B forgets the seam
`x_(6,4) xor x_(7,3)`; unordered A/C forgets the theta sheet on the finite
bank.  The seam repairs the literal B map and the direct-node source
congruence, but it neither repairs the partner-node congruence nor refines the
target node marginal.  Coupled `bar P` is the operation-safe metagraph address
for this named word.  This is the exact example the operation-indexed language
was designed to distinguish.

THM-853 iterates those fixed-size returns without pretending that they are
the increasing-size CF cocycle.  The exact algebra has 46,126 maps, 31,075
Cayley SCCs, and six rank-one horizons on the central tiles
`(5,3),(6,3),(6,4),(7,3),(7,4),(7,5)`.  Reflection pairs the outer four and
fixes the two coordinates of the B seam.  Static Q equality takes 144 values;
future equality under all return words takes 20,419.  This is the first exact
Myhill--Nerode-sized warning in the defect programme: small present fibres do
not imply a small continuation state.

Incoming THM-840 makes the general logical boundary operation-indexed: a
kernel can be a congruence for monotone insertion while failing deletion,
replacement, or exact evaluation.  THM-842/843 give the tournament-side
instance.  Bare node identity, endpoint deletion, B gap contraction,
blue/black line expansion, and a CF lift must be audited as five different
relations rather than called one recursive address.

### Positional moments and additive reconstruction

At `n=9` the 58 residual full lower-address collisions are exactly
reflection-chirality pairs.  One skew sign repairs these actual fibres.
THM-825's stronger per-state `(M0,M1,M2)` codec is
unconditionally literal through n15 and first fails at n16 on
`{0,4,5}` versus `{1,2,6}`.  This is the clean additive-combinatorics ladder:
actual minimal separator, robust moment separator, then literal word.

### Coding, matroid, and sheaf language

The 14D syndrome is a binary parity-check code; reflection cuts it to 6D; the
realized sectors span a 4D subcode.  The 28 cell-coordinate functionals
define a representable rank-four binary matroid; its seven Venn blocks have
ranks `(0,0,0,3,2,3,4)`.  Treating the seven blocks themselves as elements
instead gives a rank polymatroid, not this matroid.
The punctures are not linear-code holes: three have many compatible bases and
die nonlinearly at literal S2.  This makes “linear code plus nonlinear stalk”
the most economical current description.

THM-839 now resolves the nonlinear stalk locally.  On all 636 post-B bases,
every active difference layer toggles two ordered states and their reflected
coordinates.  Raw-S2 equality is exactly the complement-pair condition
`chi_tau=1`.  The dominant missing point has 504 bases but one uniform
`tau=7` parity obligation kills all of them.  Reflection acts with 388
endpoint-swap fixed bases and 124 two-cycles; only 58 fixed bases also satisfy
every balance bit.  Thus the final 58 are a balanced fixed locus, not merely a
reflection-fixed locus.

## Critical next computations, in order

1. **Return-closure node map.**  Canonicalize THM-853's 8,940 masks with the
   exact size-nine subset-DP isomorphism key, then minimize direct node,
   partner node, and coupled `bar P` as three separate partition-output
   automata.  Compare the result—not just its node count—with THM-851's
   constant-composition refinement.
2. **n10 codec preflight.**  Carry the rank-four sectors through the next
   relation join before building a full n9 atlas.  Retain second moments as
   the unconditional repair and test whether the residual kernel grows by
   new position modes or by node-fibre continuation.  Share the resulting
   node-coloured composition table with THM-851's independently reserved
   coherent-refinement calculation; that thread asks for constant composition
   numbers, while this one asks for preservation under a named lift.
3. **General CF-word action.**  Factor the literal coordinate action of two
   consecutive continued-fraction digits, determine its semigroup on the
   four-bit core, and classify which words preserve the coupled `P` descent.
4. **Metric fibre product.**  Couple literal `Q` orbit, skew chirality,
   THM-829 Bezout row, owner/root labels, and actual witness intervals.  Prove
   a preservation lemma before claiming any implication for LRC(14).

## Past niche threads worth pulling

- **THM-796 phase square:** coherent descent phase survives identical coarse
  node/curvature data; compare its square with the four defect coordinates.
- **THM-814 fixed-layer swap:** the n8 black codec fails on one positional
  kernel; test whether it is the first shadow of the n9 skew sectors.
- **THM-830 deletion-deck mirror current:** THM-842 applies its A/B/C roles and
  isolates the theta sheet.  THM-846 answers the first closed word: endpoint
  order retains the sheet, B transports it, and the B seam is a distinct
  coordinate.  Pull this thread next for two-word compositions.
- **THM-549/550 half-tiling parity:** blue lines are the reflection-fixed
  half-tiling; use even/odd recursion on the defect basis, not only cell count.
- **HYP-3809 parity skeleton / even black graph:** test whether sector mass or
  cube degree refines black even-degree pairing.
- **THM-849 black self-line obstruction (reserved):** if its n8 obstruction
  survives audit, intersect its Klein-four realizer weights with THM-843's
  218-node omitted cap rather than treating self-converse symmetry as a node
  scalar.
- **THM-852/854 and S14 H-drift corrections:** the proposed self-line/SC
  equality is refuted at n8 (`404!=176`), and THM-854 replaces the invalid
  involution corollary by an exact witness-order/path-XOR parity law.  Exact H
  drift requires the OCF endpoint coordinate `K`, not H alone.  Use this as a
  triple warning for the return atlas: converse-merged nodes lose same-versus-
  converse, unmarked symmetry loses witness action, and one scalar level loses
  future drift.  Test OCF `K` and witness-path XOR on the 8,940 masks only after
  node keys are available.
- **THM-850 chi-seven B3 charge (reserved):** test the proposed face weights
  `(1,2,4) mod 7` on the eleven defect sectors and on the closed A/B/C words.
  A residue law for face roles would be a new sidecar candidate for the
  `A+B+C-D-E-F+G` recursion, but no link is assumed before its proof lands.
- **THM-851 coherent node-coloured defect closure (reserved):** compare its
  coarsest constant-composition refinement with the operation kernels of
  THM-846.  Equitable composition and deterministic continuation are distinct
  tests and may demand different minimal states.
- **THM-805/826 Tutte-Farey profiles:** both succeeded by retaining labelled
  atoms before summing; use them as models for a fibre-valued generating
  function, not as a claim of direct equivalence.
- **THM-499/505 OCF correlation tower:** the spectrum sees totals while the
  OCF sees disjointness correlations; analogously the syndrome sees parity
  while the stalk sees base realization.
- **FKN transitive-corner threads:** view the four generators as higher-order
  perturbations from the transitive vacuum and compute their Walsh support.
- **Smith/electrical current:** use the `epsilon=0` localization as a genuine
  disintegration boundary and seek an antisymmetric conjugate potential.
- **Matroid circuits and cocircuits:** compute the restriction matroid,
  puncture indicator, and chirality as a nonlinear orientation of its topes.
- **Automata/Myhill-Nerode:** treat lift/delete/CF words as operations and
  minimize decorated fibres by future truth, not present equality.
- **Owner-labelled wall movies:** test whether first skew layer aligns with
  first changed owner or wall crossing in the ten-wall CF schedule.

## Guardrails

- On the canonical apex-zero bank `X_9^0`, `Lambda_9(u)=Lambda_9(v)` implies
  `v in {u,sigma u}`; this is not an all-size claim, and the raw address does
  not globally descend to the reflection quotient.
- A sector is not an isomorphism-class node, and cube adjacency is not a
  blue/black metagraph edge.
- Full rank on an overlap proves linear recovery, not nonlinear fibre purity.
- Reflection-converse merging does not authorize independent endpoint swaps;
  line and deletion couplings swap simultaneously.
- Static injectivity does not imply lift/deletion/CF continuation.
- Bare `D->d_B` is not a valid continuation stalk: all four missing sectors
  have B footprints before their A/C or S2 death certificates (THM-842).
- No tournament carrier above currently preserves the LRC loneliness metric.
  The fourteen-runner case remains open.
