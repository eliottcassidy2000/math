---
id: HYP-6880
title: Joined B3-face and folded-B2 sidecars form a recursively useful metagraph address
status: OPEN GENERAL/TRANSPORT CLAIM; exact through n=8, with n9 kernel/punctures, B3 continuation, and one n9->n10 defect-bank transport classified by THM-801/809/828/832/834/838/839/842; exact gap-tournament and half-defect carrier proved by THM-830
source: codex-2026-07-15-S12/S13/S11
depends_on: [THM-553, THM-796, THM-801, THM-809]
related: [THM-805, THM-811, THM-812, THM-813, THM-814, THM-818, THM-828, THM-830, THM-832, THM-833, THM-834, THM-838, THM-839, THM-842, HYP-2685, HYP-3234, HYP-6825, HYP-6870, HYP-6885]
---

# HYP-6880 — joined B3/B2 metagraph address

The full three-subtriangle chart preserves overlap ownership and the missing
gap-contraction face.  The mirror-folded two-coordinate chart preserves the
reflection-orbit phase that face-node projection forgets.  The hypothesis is
that their join is a useful recursive address for complement lines and their
node-pair fibres beyond the current finite atlas.

## Finite result now proved

THM-801 completed the exact audit:

```text
n       lines       Xi cells    Omega cells    Omega+S2 cells
4           4              4              4                 4
5          32             32             32                32
6         512            509            510               512
7      16,384         16,031         16,308            16,384
8   1,048,576               -              -         1,048,576
```

Thus `Omega+S2` is injective on literal lines through `n=8`.  At `n=7`, the
gap face removes 277 of `Xi`'s 353 collision excess and the `S2` vector from
the B2 chart removes the remaining 76.  Inside coloured node-pair fibres,
`S2`, `S3`, and their join leave respectively 172, 1,368, and 16 excess
collisions; the two charts preserve genuinely different information.

The same computation found exact compact companions:

- the gap-face support row resolves five of endpoint deletion's eight support
  twin pairs at `n=7`;
- the boundary-curvature polynomial `K_u(z)` separates 238/272 nodes, or
  249/272 with total `C3`;
- THM-811's richer black orbit polynomial `(C3,H_x)` separates all 272 nodes,
  while node pair plus `(q0,q1,|epsilon_Smith|)` separates only 7,248/8,064
  black edge orbits; positional `(B2,B3)` data remain necessary;
- the common-core node pair improves `Xi_7` to 16,110 cells but is weaker than
  the gap face;
- the joined address remains a static codec, not a proved recursive Markov
  state.

The finite question is exact: disintegrate literal complement lines first by
their `Omega_n` upper/three-face node tensor, then by raw mirror crossing-layer
counts, retaining loop multiplicity only at the endpoint-incidence marginal.
Record every remaining collision rather than replacing a multigraph fibre by
simple adjacency.

THM-809 proves the `n=8` row by the strictly weaker lower-only key

```text
Lambda=((A,A'),(B,B'),(C,C'),UABC,S2),
```

where all six nodes are `n=7` face nodes.  No upper `n=8` node label is used.
Its collision-excess ladder is

```text
B3 lower nodes + colour                 418
+ tau=3                                 252
+ tau=4                                 148
+ tau=5                                  74
+ tau=6                                  52
+ tau=7                                   0
+ fixed tau=8                             0.
```

Every collision before the last step is a double collision.  Thus the new
size-three crossing layer `tau=7` is the first complete separator, while the
fixed layer is redundant after it.  The literal B3 lower-line triple is also
injective on all `2^20` lines, independently checking the Cech descent.

The 418 base doubletons have first-separator counts
`166,104,74,22,52` on `tau=3,...,7`.  Their literal face-difference patterns
are `A+B:4,A+C:44,B+C:4,A+B+C:366`; none has both endpoint faces `A,C`
equal, and none contains one of the four S11 `n=7` phase-square residuals as a
face.  The `n=8` ambiguity is therefore a new genealogy rather than a direct
lift of the two-face phase obstruction.

The all-size claim remains deliberately open.  Even an injective finite codec
need not be a continuation-complete state: THM-796 already proves failure of
strong node lumpability.  A recursive theorem would have to transport the
sidecars and their action on the rooted Hamiltonian-path fibre, not merely
name a static signature.

THM-814 identifies the first positional warning at the same size.  On black
reflection orbits, projected node pair plus orbit-symmetrized raw `(B2,B3)`
has sixteen double collisions at `n=8`.  All swap two positions in the fixed
`ABC` layer while preserving its count.  THM-809's lower-node `Lambda` remains
injective because it retains more face position; a fixed-layer first moment
also repairs all sixteen.  Therefore any recursive minimization of
`Omega+B2` must keep positional data precisely when a layer ceases to be
reconstructible from its count.

## Exact carrier supplied by THM-830

THM-830 identifies the object on which a recursive minimization should now be
performed.  A size-`n` tiling is exactly a labelled tournament on the `n-1`
path gaps.  The `B` face is its cut/window tournament and the discarded
adjacent arcs are a literal seam word.  More globally, if `C_n` is the
half-coordinate cube with fixed-boundary trace, then

```text
X_n ~= C_n x_trace C_n.
```

Writing a tiling as common boundary `z`, first half `u`, and additive defect
`delta`, the merged-node label

```text
c_n(z,u,delta)=pi_n(u,u+delta,z)
```

determines the full node/tiling/edge tensor.  Blue is `delta=0`; every nonzero
`delta` is black and is itself a size-`n-1` blue word.  Thus `S2` is a population
shadow of a literal defect word, while `B3` is a face shadow of a literal gap
tournament.  The next minimizer should refine or disintegrate this tensor,
not search for an unrelated tuple of scalars.

The unprojected recursion already has exact composition numbers.  With
`O_n=Blue_n/<1>` and `D_n=Blue_(n-1)`, complement lines are the arrows of
`O_n semidirect D_n`; exact-defect matrices obey

```text
A_delta A_epsilon=A_(delta+epsilon).
```

Each boundary-trace component carries the binary Hamming scheme on
`h_(n-1)` coordinates.  Its weight intersection numbers and Krawtchouk
spectrum are explicit in THM-830.  Thus the open coherence problem is not in
raw defect composition.  It is exactly whether the nonlinear colouring
`c_n` together with the ordered sidecars can be refined to an equitable
pushforward of that algebra.

The preservation boundary is already exact.  Total weighted deletion roles
collide at `n=5`; even deck totals plus `C3` collide between distinct score
spectra at `n=6`.  Conversely, on the whole one-flip shell the address
`(gap,|mirror depth|)` is a complete merged-node and reflection-orbit line
codec.  THM-814's sixteen fixed-layer collisions are separated by the gap
component of this weighted deck.  Therefore literal root positions, kernel
rows, or an equivalent action on `c_n` remain necessary beyond radius one.

## Exact completion and minimization target

For a literal core `c`, write THM-796's reconstruction variables as
`(p_L,p_H,a)` and define

```text
I_(u,c)(p_L,p_H,a)=1[reconstruct(c,p_L,p_H,a) in F_n(u)].
```

The full anchored Boolean Möbius transform of this function is invertible and
therefore continuation-safe if its substitution action is retained.  The
research question is whether `Omega+S2` identifies a small closed sector of
these coefficients.  Static equality must be tested separately from equality
under future deletion, complement, reflection, and lift operations.

THM-828 now decides the first previously open finite size `n=9`.  Raw counts
do not remove every lower-node collision: they leave exactly 58 disjoint
doubletons, all reflection mates.  A single antisymmetric positional
chirality bit repairs them all.

THM-818 removes the exponential outer scan from that question.  Define the
oriented face observation

```text
H_8(x)=(nu_8(x),nu_8(kappa x),chi_8(x))
```

and retain its kernel relation `R_8`, not just its fibre labels.  Exact census
gives `876,512` fibres and `5,997,416` relation rows, of which `3,900,264` are
off-diagonal.  Ordered pairs of candidate `n=9` tilings with equal three-face
data are exactly the compatible triangles in three copies of `R_8`; both
apexes must then be zero, the remaining upper colour bit and `S2` must agree
(the upper bit plus the three face colours is `UABC`).  The three face
projections are therefore part of the address's recursive structure.  Fibre
cardinalities alone destroy precisely the gluing information needed at the
next size.

THM-828 factors this join through a 14-dimensional raw-S2 difference code and
then through the unique gluing of its `A/C` faces.  It finds 9,540 literal
overlap matches, 636 B-compatible candidates, and 58 final pairs.  Their
eleven difference masks span a punctured four-dimensional reflection-defect
cube; one absent sector is forbidden by face-H support and three pass B
compatibility but fail literal S2.  The raw address therefore determines the
reflection orbit but needs one sheet-orientation bit to determine the literal
tiling.

THM-832 identifies the recursive linear core of those differences.  Each of
the three faces, each pair overlap, and even the triple overlap retains all
four survivor coordinates; explicit decoders recover the global sector from
any one chart.  The resulting constant Cech system has no positive-degree
cohomology.  Nevertheless equal regional support profiles carry different
literal fibre masses, so this holonomy-free linear carrier cannot replace the
nonlinear kernel rows or their base witnesses.

This also identifies the positional boundary exactly.  At
`n=9`, each nonfixed layer has at most three positions and the four-position
fixed layer has only three free after apex orientation, so per-state count and
first moment reconstruct the literal layer.  At `n=10`, `{0,3}` and `{1,2}`
first defeat count-plus-sum.  Static raw-S2 injectivity is now refuted at
`n=9`, while a minimal one-bit skew repair is proved there.  Continuation
purity under lift, deletion, and arbitrary centered-CF words remains open.

THM-838 decides one centered-CF continuation on the entire defect bank:
literal rank stays four, target raw-S2 sees exactly two dimensions, all 58
literal reflection orbits and coupled endpoint-node cells remain distinct,
and eight static chirality signs reverse.  THM-834 gives the exact quotient
placement.  The 116 endpoints occupy 53 converse-merged nodes, five of which
contain two collision pairs.  The bare node quotient is insufficient, as is
its 54-valued complement-partner marginal, but their coupled `bar P` coordinate
has 58 singleton values.  The six face presentations occupy 155 size-eight
nodes and select 87 exclusively black projected complement-line node pairs.
The defect cube is therefore neither a node set nor the metagraph: it is a
decorated carrier whose stalk maps into node and line fibres.

THM-839 explains every puncture in that stalk.  One sector has empty A/C
face-kernel support.  On each of the other three, an active raw-S2 layer is a
two-toggle packet and equality is exactly one complement-pair parity
obligation.  Uniform killer layers are `tau=8,7,7`; in particular one `tau=7`
bit kills all 504 bases of `54a5692`.  The 636 post-B bases decompose under
reflection into 388 endpoint-swap fixed points and 124 two-cycles.  Exactly
58 fixed points also satisfy every ordered `01/10` balance bit.  The nonlinear
stalk is therefore a balanced reflection-fixed locus, not an unexplained
punctured code.

THM-842 decides the first deletion/gap-contraction continuation on that stalk.
B is the unique B3 role descending to literal Q on each source cell and maps
the 58 cells injectively.  The reflection-exchanged A/C endpoint deck has 29
doubletons.  Their involution is the apex-relative antipode
`theta=sigma(nonapex-complement)`, which preserves `D` and the entire upper
skew word; the missing bit is therefore affine sheet context, not chirality.
Endpoint deletion converts theta to ordinary complement-reflection, while B
transports theta itself.  The unpositioned induced-tournament deletion deck
has exactly the same 53-cell partition as the bare upper node, whereas the
position-marked deck is injective.  This proves that recursive usefulness
requires path position and sheet action even after static node completeness.

## Continued-fraction connection

THM-778 supplies the exact analogy and the guardrail.  THM-808 now supplies an
actual transported coefficient: on a degree-one prime-sheet fibre, a centered
Christoffel owner-count block acts on the duplicate root by
`d'=d-sum c_a w_a^(-1)`.  The same block and same least-path mask can still
have two different target masks, so the root/owner lift is genuinely needed.
A continued-fraction
digit is useful only with the induced substitution on its labelled token
fibre; a metagraph address is useful only with the induced action on its
path/core stalk.  HYP-6880 therefore asks for a transported address, not a
larger tuple of static scalars.  MPA-38 proposes the first direct test: apply
centered Christoffel substitutions to the low/high leg variables and their
Möbius coefficients.

THM-812 now performs that direct test for the first nontrivial centered word
`d=(1,2)`.  Its coordinate-copy action `X_5->X_6` is complement/reflection
equivariant.  Bare nodes are not functorial, but every one of the 20 projected
coloured edge cells has a unique, distinct target cell.  On the coefficient
side, the induced general law is subset-image pushforward; degree at most
three loses one fixed-core target-node pair and three explicit quartic
coefficients repair it.  The remaining open transport claim is therefore not
whether a CF digit can act at all, but whether a finite coefficient sector
containing those quartics is closed under arbitrary Euclidean words and the
`Omega+B2` recursion.

THM-813 supplies the first negative transport boundary.  The immediate
`X_6->X_7` centered copy splits 51 projected colour/node-pair edge cells, even
though it remains exact on literal tilings and complement lines.  The safe
recursive carrier is the staircase-reflection orbit `Q_n=E_n/<sigma>`; descent
to a projected edge cell is an additional constancy property, not functorial
by definition.  The two-step node-indicator audit still closes at degree four,
so the live finite-sector target is now: quartic plus the `Q`-orbit sidecar,
with closure tested under coordinate-image saturation.

THM-838 tests the next consecutive word on the first actual nontrivial
kernel.  The literal and `Q` actions remain exact; the special 58-cell coupled
`P` coordinate also descends without collision, even though bare projected
edge cells fail globally by THM-813.  The raw histogram is not stable: only
ten of 58 pairs remain equal, and two defect directions become visible at
target layers `tau=5,7`.  This replaces the vague transport target by a
finite one: classify the semigroup action of CF words on the decorated
four-bit stalk and its endpoint coupling, not just on a node or S2 vector.

Assumption challenge: the useful third face is not induced deletion of an
original tournament vertex.  THM-830 proves that it is an honest tournament
on path gaps/cuts, with adjacent seams retained separately.  It preserves the
gap tournament and exact face-overlap ownership but destroys literal induced-
subtournament ancestry in the original marked tournament; that loss must be
stated whenever the face is pushed to merged nodes.
