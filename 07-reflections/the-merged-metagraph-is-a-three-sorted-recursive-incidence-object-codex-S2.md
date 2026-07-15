# The merged metagraph is a three-sorted recursive incidence object

**codex-2026-07-14-S2/S11. Status:** supporting synthesis for THM-796:
general elementary identities, exact `n=3..7` atlas audit, and a preservation
prescription.  General statements are separated below from finite census
evidence.

## The graph drawing suppresses the primary objects

Fix the Hamiltonian path `n -> n-1 -> ... -> 1`, and let

```text
E_n = {(x,y): x-y >= 2},             m_n = |E_n| = C(n-1,2),
X_n = F_2^{E_n},                     fixed-path tilings,
G_n = tournament isomorphism classes,
M_n = G_n / (T ~ T^op),              converse-merged nodes.
```

There are really three sorts, not one graph:

```text
                    q_n
        X_n ---------------------> M_n
         |
         | quotient by tau_n
         v
        L_n = X_n / <tau_n>.
```

Here `tau_n(t)=t+1` flips every non-path tile.  A line `ell=[t]` in
`L_n` has node boundary

```text
partial_n(ell) = {q_n(t), q_n(tau_n t)} in Sym^2(M_n).
```

The two entries may agree.  A self-loop is therefore not an exceptional
graph convention: it is a genuine complement orbit whose two tiling
presentations lie in the same merged node.  Parallel lines are also genuine.
The simple node support graph forgets both facts.

Let

```text
sigma_n(x,y) = (n-y+1,n-x+1)
```

be anti-diagonal grid reflection.  It commutes with `tau_n`, and it sends the
tournament represented by a tiling to the converse class.  Consequently
`q_n sigma_n=q_n`.  A line is blue if `sigma_n(t)=t`, and black otherwise.
This is a color on **line instances**, not merely on unordered node pairs: at
`n=6,7`, one node pair can support both colors.

The intrinsic inverse-fibre theorem supplies another description of the top
sort.  A tiling presentation of `[T]` is a directed Hamiltonian path of `T`
modulo `Aut(T)`.  Thus `X_n` is the finite space of tournament objects together
with an observer path; `M_n` forgets the observer and then folds converse.
This is why the tiling sort is mathematical data rather than a visualization
artifact.

The complete fixed-size node shadow is still useful when its weights are not
discarded.  Put `F_u=q_n^{-1}(u)` and define the coloured half-edge kernel

```text
A_n^chi(u,v)
  = #{t in F_u:q_n(tau_n t)=v and color([t])=chi}.
```

It is symmetric, its row sum is `|F_u|`, and a diagonal entry is twice the
number of loop lines.  Hence row normalization is a reversible Markov kernel
with stationary mass `|F_u|/2^m_n`.  This reconstructs fibre volume and
coloured incidence, but it does not reconstruct which parallel line was used,
which class sheet a merged loop traversed, or how either endpoint continues
under deletion.  “Weighted” is therefore necessary but not synonymous with
“recursively sufficient.”

## Endpoint deletion is the honest recursion

There are two natural face maps on tilings:

```text
d_top    : X_n -> X_{n-1},  delete the first path vertex n (drop x=n),
d_bottom : X_n -> X_{n-1},  delete the last path vertex 1 (drop y=1, relabel).
```

Each drops exactly `n-2` free coordinates, so every lower tiling has
`2^(n-2)` extensions through either face.  Direct coordinate calculation gives

```text
d_top tau_n       = tau_(n-1) d_top,
d_bottom tau_n    = tau_(n-1) d_bottom,
d_top sigma_n     = sigma_(n-1) d_bottom,
d_top d_bottom    = d_bottom d_top.
```

The last equality is understood where both composites are defined.  Repeated
endpoint deletion is therefore a genuine two-faced recursive system on
presentations.  Reflection exchanges the faces and complementation is natural
for both.

Because deletion commutes with `tau`, it induces honest maps

```text
bar d_top, bar d_bottom : L_n -> L_{n-1}.
```

Every lower line has exactly `2^(n-2)` high-line lifts: its two tilings have
`2 * 2^(n-2)` extensions, paired freely by high complementation.  This is a
regular recursion of line **instances**.

There is generally no map `M_n -> M_{n-1}`.  Different tilings in one merged
node can delete to different lower nodes.  The exact node datum is the weighted
span

```text
M_n  <-q_n-  X_n  -q_(n-1)d_top->  M_{n-1}
```

or its matrix

```text
P_n(a,b) = #{t in X_n: q_n(t)=a and q_(n-1)(d_top t)=b}.
```

It has exact marginals

```text
sum_b P_n(a,b) = |q_n^{-1}(a)|,
sum_a P_n(a,b) = 2^(n-2) |q_(n-1)^{-1}(b)|.
```

Substitution `t=sigma_n(u)` and the reflection identity show, for every `n`,
that the top and bottom matrices are equal after converse merging.  This is
the precise sense in which the two endpoint deletions become indistinguishable
at node level.  They are exchanged by a symmetry; they are not the same map on
tilings.

The exact atlas shows how strongly a node fails to have one parent:

| n | support of `P_n` | number of high nodes by number of lower targets |
|---:|---:|---|
| 4 | 5 | `1:1, 2:2` |
| 5 | 19 | `1:3, 2:5, 3:2` |
| 6 | 112 | `1:5, 2:4, 3:11, 4:5, 5:8, 6:1` |
| 7 | 1312 | `1:11, 2:16, 3:31, 4:62, 5:33, 6:70, 7:49` |

## Even the weighted node span is not continuation-complete

Normalize a row of `P_n` to obtain the distribution of the lower node after
deleting a uniformly chosen tiling in a high node fibre.  It is tempting to
multiply these node kernels for repeated deletion.  That is wrong: conditioned
on the first high node, the reached tilings are not uniformly distributed
inside their middle node fibre.

The exact two-step comparison gives:

| n | unequal high-to-low entries | high rows with an error | maximum error |
|---:|---:|---:|---:|
| 5 | 16 | `8/10` | `1/2` |
| 6 | 86 | `32/34` | `1/2` |
| 7 | 1778 | `270/272` | `1/2` |

At `n=5`, for example, one address-ranked high node reaches a particular
`n=3` node with actual probability `1`, while the product of the two normalized
node matrices predicts `1/2`.  The ranks are atlas coordinates, but the
inequality is invariant content.

Thus the tiling reached after one step is recursively necessary.  A span can
be composed exactly only while retaining its witness in `X_{n-1}`; replacing
that witness by its node forms the wrong fibre product.  In automata language,
node identity is not a continuation equivalence.  In probability language,
the node projection is not a sufficient statistic for the deletion process.

## The exact recursive edge object is endpoint-coupled

For a high line, choose one representative `t` only long enough to write

```text
(a_0,a_1;b_0,b_1)
  = (q_n(t),q_n(tau t); q_(n-1)(dt),q_(n-1)(d tau t)).
```

Changing the representative to `tau t` performs the **simultaneous** swap

```text
(a_0,a_1;b_0,b_1) ~ (a_1,a_0;b_1,b_0).
```

Let `tilde E_n^{alpha,beta}` count these simultaneous-swap orbits.  This is
the exact endpoint-coupled tensor: it knows which high half-edge descends to
which low half-edge.  Sorting `(a_0,a_1)` and `(b_0,b_1)` independently gives
a smaller tensor `E_n`; that operation additionally identifies straight and
crossed couplings and is only a marginal.  The loss is already visible in the
support counts:

| n | coupled support | independently sorted support |
|---:|---:|---:|
| 4 | 4 | 3 |
| 5 | 31 | 30 |
| 6 | 464 | 455 |
| 7 | 15112 | 15074 |

Diagonal pairs retain self-loops in both tensors.  The projected tensor has
the exact marginals

```text
sum_{b,b',beta} E_n^{alpha,beta}(a,a';b,b')
    = W_n^alpha(a,a'),

sum_{a,a',alpha} E_n^{alpha,beta}(a,a';b,b')
    = 2^(n-2) W_(n-1)^beta(b,b'),
```

where `W` counts unoriented line instances, including each loop once.  The
half-edge kernel above agrees with `W` off the diagonal and equals `2W` on the
diagonal.  The tensor is the
direct answer to “what happens to an edge as `n` changes?”, while `tilde E`
is the answer to “what happens to each of its two half-edges?”  A high
edge can have several low images, a low edge can receive several high node
supports, and loop status can appear or disappear.  The verifier records all
four loop transitions rather than treating loops as zero-length edges.

The finite census illustrates the information loss in simple support:

| n | merged nodes | blue cross/loop lines | black cross/loop lines | colored/plain supports | loop supports |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | `1/0` | `0/0` | `1/1` | 0 |
| 4 | 3 | `1/1` | `2/0` | `3/3` | 1 |
| 5 | 10 | `8/0` | `20/4` | `20/20` | 2 |
| 6 | 34 | `30/2` | `456/24` | `187/183` | 10 |
| 7 | 272 | `256/0` | `16014/114` | `6126/6076` | 44 |

There is a second exact parity layer.  Reflection `sigma_n` acts on line
instances.  A blue line is fixed pointwise.  A black line cannot be fixed as a
line: `sigma_n(t)=tau_n(t)` would disagree on every reflection-fixed tile.
Hence reflection-line orbits have size one for blue and size two for black.
Because `q_n sigma_n=q_n`, the two black lines in such an orbit project to the
same merged node pair.  It follows that, including diagonal pairs,

```text
every projected black-pair multiplicity is even,
every black half-edge degree is even (a loop contributes two).
```

Blue half-edge degree is the number of grid-symmetric tilings in the node
fibre.  It is odd on every self-converse node and zero on every non-self-
converse merged node.  These statements are asserted directly for every atlas
level, so the line involution, node quotient, and loop convention are checked
together rather than as separate censuses.

## Colour is the zero fibre of a linear defect

There is a more exact description than a blue/black bit.  The line set is the
binary quotient code

```text
L_n = F_2^{m_n}/<1>.
```

If `f_n=floor((n-1)/2)` is the number of reflection-fixed tiles and
`h_n=(m_n+f_n)/2=floor((n-1)^2/4)`, define

```text
delta_n([t]) = t+sigma_n(t).
```

This is complement-independent and gives the exact sequence

```text
0 -> Fix(sigma_n)/<1> -> L_n -> im(1+sigma_n) -> 0,

dim Fix(sigma_n)/<1> = h_n-1,
dim im(1+sigma_n)    = (m_n-f_n)/2.
```

Blue is the zero-defect fibre.  Every nonzero defect labels an equally large
black coset, of size `2^(h_n-1)`.  The impossible equation
`sigma_n(t)=tau_n(t)` is now visible algebraically: an image defect vanishes on
reflection-fixed coordinates, while the all-one word does not.  Reflection
therefore pairs black lines freely even after conditioning on an exact defect
and an exact node pair.  This strengthens “black multiplicity is even” from a
global colour fact into a defect-by-defect law.

## Blue and black form a renormalization channel, not nested subgraphs

Put

```text
h_n = floor((n-1)^2/4),
B_n = 2^(h_n-1),
K_n = 2^(m_n-1)-B_n.
```

For `n>=4`, the exact color transition matrix, with high color indexing rows
and deleted color indexing columns, is

```text
                 low B                              low K
high B       2^(n-3)                         B_n - 2^(n-3)
high K       2^(n-2) B_(n-1) - 2^(n-3)      2^(n-2) K_(n-1) - (B_n-2^(n-3)).
```

Row sums are `(B_n,K_n)` and column sums are
`2^(n-2)(B_(n-1),K_(n-1))`, as required by regular line lifting.  The exact
values are:

| n | BB | BK | KB | KK |
|---:|---:|---:|---:|---:|
| 4 | 2 | 0 | 2 | 0 |
| 5 | 4 | 4 | 12 | 12 |
| 6 | 8 | 24 | 120 | 360 |
| 7 | 16 | 240 | 1008 | 15120 |

The formula has a short structural proof.  On the retained staircase,

```text
sigma_n(x,y)     = (n-y+1,n-x+1),
sigma_(n-1)(x,y) = (n-y,n-x).
```

Their composition translates `(x,y)` to `(x+1,y+1)`.  Therefore a high tiling
and its top deletion are both blue exactly when the mask is constant on every
difference diagonal `x-y=d`.  There are `n-2` such diagonals, so there are
`2^(n-2)` masks and `2^(n-3)` complement-lines.  This proves the BB entry.  The
other entries follow from the row and column sums.

This identifies the recursively stable blue-to-blue core with the old
**difference-striped** or circulant-like carrier.  Among all high blue lines,
the fraction that remains blue after deletion is

```text
2^(n-h_n-2),
```

equal to `1, 1/2, 1/4, 1/16` for `n=4,5,6,7`.  Most blue structure turns black
under size descent.  Conversely, a low blue line has high blue lifts only when
it is difference-striped, and then it has exactly two.  Thus the blue
half-tiling at fixed `n` is real, but blue subgraphs are not a nested tower in
`n`; their common recursive spine is much smaller.

## Symmetric deletion turns the colours into a defect tower

The one-end recursion is the right microscope for half-edge coupling.  The
two-end Mode-B map

```text
p_n=d_top d_bottom:X_n->X_(n-2)
```

is the right microscope for colour ancestry because it commutes with both
`tau` and `sigma`.  It induces exact sequences

```text
0 -> F_2^(2n-5) -> L_n -> L_(n-2) -> 0,
0 -> F_2^(n-2)  -> Blue_n -> Blue_(n-2) -> 0,
0 -> F_2^(n-3)  -> Defect_n -> Defect_(n-2) -> 0.
```

The Mode-B high/low colour channel is consequently triangular:

```text
BB = B_n,                         BK = 0,
KB = (2^(2n-5)-2^(n-2))B_(n-2),  KK = 2^(2n-5)K_(n-2).
```

A nonzero parent defect can never be repaired by adding boundary tiles; it is
inherited blackness.  A blue parent can acquire a nonzero kernel defect; this
is fresh blackness.  The black sea is therefore recursively stratified by the
first level at which its defect appears, together with the later images of
that defect.  That ancestry is invisible in the binary colour.

Loops require one further bit.  When `[t]` is a merged loop, its two endpoints
may lie in the same unmerged tournament class or in the two converse classes
folded into one node.  Call this same-sheet versus converse-switch holonomy.
At `n=6` the 24 black loops split `6+18`; at `n=7` the 114 black loops split
`44+70`.  Under Mode-B `6->4`, loop status itself has transition counts

```text
cross/cross 366, cross/loop 120, loop/cross 18, loop/loop 8.
```

So neither “is a loop” nor “which class sheet it traverses” is hereditary.
The diagonal of a merged adjacency matrix is not enough.

## The recursive object

The most economical faithful description is a graded involutive incidence
system:

1. `X_n` is the presentation stalk, with commuting endpoint faces and the
   involutions `sigma_n,tau_n`.
2. `M_n` is a quotient base, connected across sizes by a weighted
   correspondence rather than a function.
3. `L_n` is the complement-orbit sort, with regular deletion maps.
4. `partial_n` attaches each line instance to a node pair, including the
   diagonal.
5. `delta_n` refines colour into a linear defect and remembers black ancestry.
6. Loop incidence carries same-sheet/converse-switch holonomy.
7. The lower-line pair plus one coherent phase bit is an exact recursive line
   address.
8. `Xi_n` jointly retains the upper endpoints, both face endpoints, and the
   three-colour word.
9. `tilde E_n` is the one-face endpoint-coupled tensor; `E_n` is its
   independently sorted node-edge marginal.

One can call the resulting node picture a colored Bratteli diagram, but its
edge fibres and witnesses must remain attached.  The recursive failure of the
node Markov product proves that a bare Bratteli matrix is too coarse.

This also challenges the default Tournament Analysis vertex choice.  The
natural “vertices” are alternately tournament classes, Hamiltonian-path
presentations, and complement-line events.  Their native relation is symmetric,
weighted, and looped, so forcing it into a binary tournament destroys exactly
the multiplicity and fixed-point data that control recursion.  The preserved
predicate is not merely current adjacency.  It is continuation under every
future deletion, extension, complement, reflection, and quotient step.

The groupoid reading makes this failure unsurprising.  A tiling is a
tournament with a marked Hamiltonian path.  Deleting a terminal path vertex is
an honest operation on marked objects.  Forgetting the path and then asking
for a node parent is not cartesian: different path orbits in one node have
different parents.  Arbitrary vertex deletion also needs the vertex's
automorphism orbit and stabilizer; internal path deletion needs a repaired
Hamiltonian path because its shortcut can face backward.

Deletion is not reversible, so the whole tower is not itself a groupoid.  The
accurate phrase is **groupoid-valued recursion**: complement, reflection, and
isomorphism transport form the reversible core inside a fixed stratum;
deletion and certified quotient/descent are functors or spans between strata.
A decreasing condensation graph must be proved, not inferred from the word
“recursion.”

This suggests a recursive definition rather than another hand-picked tuple.
Take the operation alphabet `{d_top,d_bottom,tau,sigma}` and a declared set of
terminal observations.  Two marked states have the same exact address only
when every legal continuation word produces the same observations.  This is a
finite-state/Nerode description at a bounded level.  Weighted line-WL is one
computable approximation to it, but the non-Markov audit proves that node
identity alone is not stable under continuations.

The first exact minimization makes that definition operational.  On literal
lines in `L_3 disjoint-union ... disjoint-union L_7`, refine a chosen current
observation until it is stable under top deletion, bottom deletion, and
reflection.  At `n=7` the labelled-face cell counts are

```text
node boundary   6076 -> 16359       Xi             16031 -> 16382
+ colour        6126 -> 16359       Xi+sheet+delta 16270 -> 16382
+ loop sheet    6128 -> 16359       phase address  16384 -> 16384
+ exact defect  8091 -> 16359       literal line   16384 -> 16384.
```

The first four stable partitions are equal, as are the two `Xi` partitions.
This is a useful warning about the word “preserve.”  Colour, class sheet, and
defect are mathematically different coordinates with different algebraic
laws, but for this bounded labelled-face observation language they add no
line separation beyond the bounded future tree of node boundaries.  They cannot
be discarded globally: changing the operation alphabet, asking a direct
colour/defect theorem, changing the face-successor semantics, moving to
`n>=8`, or pulling back to metric LRC states changes the equivalence relation.  Necessary
information is therefore not an absolute list of fields.  It is a congruence
relative to `(operations, terminal predicates, gauge)`.

Only two `Xi` collision pairs survive.  Both are phase-torsor deck pairs:

```text
0x12ca/0x12cb over face lines (53,150),
0x146c/0x146d over face lines (150,163).
```

In each pair the exact reflection defect is the same (`0x06a6`) and the
coherent phase bit is opposite.  The deck move is the upper apex flip, and
that coordinate is fixed by reflection, so `1+sigma` cannot see it.  This is
the precise reason the defect tower cannot replace the phase torsor.  The four
lines form a commuting involution square: phase swaps within the displayed
pairs, while reflection swaps the two pairs.  They are parallel black lines
between one non-self-converse node, whose two ordinary class sheets contain
151 tiling presentations each, and a self-converse node with 57 presentations.
Phase changes only the marked-path presentation inside fixed ordinary endpoint
classes; reflection exchanges the two 151-element converse sheets.  The
class-identifying relabelling for each phase pair is unique, fixes both path
endpoints, and cycles all five interior vertices.  The two interior 5-cycles
are reflection-conjugate.  This reframes the next search: classify endpoint-
fixed interior permutations that absorb simultaneous complement of both
faces.  Such a permutation is an algebraic witness for phase blindness and
may be enumerable without constructing the full next-size node atlas.

This suggests a better recursive ontology: node data are observations, lines are
continuation states, tilings are witnesses, and phase is descent data gluing
the two face witnesses.  The object is organized by what operations can be
performed, not by which marks happen to appear in one graph drawing.

If top and bottom successor cells are retained only as an unordered multiset,
the same carriers give `8310,8310,8310,8311,16380,16380,16384,16384` cells.  This
does not quotient current states by reflection, and ordered fields already in
the current `Xi` observation stay ordered; it changes only the refinement of
future branches.  The large loss shows that named future face order is real
transport data unless the target theorem is invariant under exchanging the
two ends.

The S9 carrier Tournament Analysis makes the information comparison explicit.
Its vertices are

```text
fibre_size, colour_degree, line_support_row, line_weighted_row,
lower_face_support, lower_face_weighted, lower_face_normalized, exact_node.
```

The pairwise observable is separated unordered node pairs; total retention
and retention per cell are the two switches; the displayed list is the tie
Hamiltonian path.  At `n=7` both tournaments are transitive with score
histogram `{0:1,...,7:1}`, zero directed triangles, singleton SCCs, and one
Hamiltonian path, but the gauges flip 18 edges.  This is diagnostic telemetry,
not a replacement for the three-sorted relation.

The continuation-carrier tournament uses node boundary, colour, sheet,
defect, `Xi`, `Xi+defect`, phase, and exact line as vertices.  Its pairwise
observable is the sign of the difference in `n=7` line pairs separated after
stable labelled refinement; total separation and separation per cell are the
switches.  Both tournaments are transitive with score histogram
`{0:1,...,7:1}`, no directed
triangles, singleton SCCs, and one Hamiltonian path, but they flip 20 edges.
Total retention orders the unique path
`phase, exact, Xi, Xi+defect, node, colour, sheet, defect`; economy follows
the declared carrier order.  Transitivity is forced by scalar ranking and is
therefore a consistency check.  The switch reversal is still the point: a
nearly exact recursive quotient can be more economical than an exact name,
but only the latter closes the final phase pairs.

## What this changes about the four-coordinate LRC object

The metagraph is not the HYP-6815 suspension

```text
{(u,t,c,lambda):u=ct and Phi_(A,R)(u,t)>=lambda}.
```

It is a constructible atlas placed over phase-order strata of that suspension.
The theorem-facing LRC predicate is nonemptiness of the integer-slope section
at `lambda=1/14`.  The three-sorted object preserves combinatorial order,
observer-cut fibres, complement transport, and recursive defect ancestry.  It
destroys metric widths, owners, endpoint side, scale/residue, wall chronology,
inverse winding, and carry.  Those fields remain a stalk above the atlas until
one proves fibre-purity or a reconstruction/annihilation theorem.

This gives a practical preservation hierarchy:

```text
simple support
 < weighted coloured node incidence
 < named line + defect + loop sheet
 < simultaneously coupled half-edge transport
 < marked deletion witness
 < owner-labelled metric wall state.
```

The hierarchy is recursive: what must be retained at one level is exactly the
information needed to make the next intended operation well-defined.

## Cross-thread convergence: component obligations and three holonomies

The independently developed LRC14 component-obligation synthesis reaches the
same preservation law from circle geometry rather than tournament tilings.  It
proposes a typed packet

```text
(P,C,O,I;A,M,W,K)
```

consisting of a normalized packet/chart, a safe component or endpoint germ,
active obligations, their labelled incidence, and arithmetic, metric,
transport, and proof decorations.  This is not yet a minimality theorem, but
it identifies what the metagraph atlas must sit **over**.  The correspondence
is structural:

```text
merged node M_n          quotient chart/base P,
tiling half-edge X_n     marked observer/witness transport W,
line instance L_n        reversible event/obligation pairing,
line-node incidence      component-obligation incidence I,
defect and phase         arithmetic/transport decorations A,W,
continuation observation proof predicate and residual debt K.
```

The comparison also corrects an ambiguity.  “Holonomy” has at least three
typed meanings:

```text
h_class : same or converse tournament class sheet of a merged loop,
h_red   : reduced token/sheet return modulo the declared gauge,
Delta_M : metric translation of the component and endpoint phases.
```

THM-796 computes `h_class`.  THM-794 shows why the other two cannot be merged:
a full packet can have trivial reduced return while translating its metric
base by a nonzero amount.  A closed loop in a finite collision automaton need
not be a closed loop in the LRC suspension.

This leads to a sharper quotient test.  Before discarding a field, prove one
of the following: the theorem-facing observation is constant on its fibres;
the field is reconstructed by retained transport; a metric/Fourier/Hall/dual
certificate annihilates its effect; every legal continuation has the same
terminal result; or the ambiguity is emitted as a named residual obligation.
This five-way rule unifies the metagraph stalk, controlled forgetting, and the
component-obligation program without pretending that any one finite quotient
is the underlying four-coordinate object.

## Reproducibility

Run:

```bash
python3 04-computation/merged_metagraph_recursive_three_sort_audit_codex_S2.py --stdout
python3 04-computation/three_sorted_metagraph_recursion_codex_S9.py
python3 04-computation/three_sorted_metagraph_continuation_minimization_codex_S11.py --stdout
node 04-computation/test_tournament_tiling_explorer_line_api_codex_S9.js
```

The independent verifiers consume the exact `n=3..7` address atlases, check
every tiling and line instance, verify the general identities and formulas
above, and write the S2, S9, and S11 result artifacts.  The JavaScript test
checks the browser's literal tiling/line/node inverse API through `n=6`.
