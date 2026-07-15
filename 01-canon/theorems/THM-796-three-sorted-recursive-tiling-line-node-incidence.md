---
id: THM-796
title: Three-sorted recursive incidence of tilings, complement lines, and converse-merged tournament nodes
status: PROVED (general pullback/torsor, incidence, defect/parity, one-face and Mode-B recursions, and colour laws) + FINITE-EXACT (node/line coupling, loop holonomy, non-lumpability, bounded continuation minimization, and census through n=7)
source: codex-2026-07-15-S9/S11 (independent S9 and S2 audits, reconciled)
depends_on: [THM-280, THM-345, THM-643, THM-781, THM-793]
related: [THM-477, THM-785, THM-790, THM-801, HYP-6815, HYP-6825, HYP-6865]
verification:
  - 04-computation/three_sorted_metagraph_recursion_codex_S9.py
  - 05-knowledge/results/three_sorted_metagraph_recursion_codex_S9.out
  - 05-knowledge/results/three_sorted_metagraph_recursion_codex_S9.json
  - 04-computation/merged_metagraph_recursive_three_sort_audit_codex_S2.py
  - 05-knowledge/results/merged_metagraph_recursive_three_sort_audit_codex_S2.out
  - 04-computation/three_sorted_metagraph_continuation_minimization_codex_S11.py
  - 05-knowledge/results/three_sorted_metagraph_continuation_minimization_codex_S11.out
  - 04-computation/mobius_cech_metagraph_codec_codex_S12.py
  - 05-knowledge/results/mobius_cech_metagraph_codec_codex_S12.out
  - 03-artifacts/visualizations/tournament-tiling-explorer.html
  - 04-computation/test_tournament_tiling_explorer_line_api_codex_S9.js
---

# THM-796 — The three-sorted recursive metagraph

## 1. The object is a span, not one graph

Fix the directed Hamiltonian path `n -> n-1 -> ... -> 1`.  Put

```text
S_n = {(x,y):1 <= y < x <= n and x-y >= 2},
m_n = |S_n| = binom(n-1,2),
X_n = F_2^{S_n}.                                           (1.1)
```

An element of `X_n` is a fixed-path tiling, equivalently a tournament with the
displayed Hamiltonian observer path.  Define the commuting involutions

```text
a_n(t) = t+1                                      all-tile complement,
r_n(t)(x,y) = t(n-y+1,n-x+1)                     grid reflection. (1.2)
```

By THM-280, `r_n` realizes tournament converse after relabelling; `a_n` is the
antipodal complement operation inside the fixed-path cube.  Let

```text
I_n = tournament isomorphism classes,
M_n = I_n/(T ~ T^op),                            merged nodes,
L_n = X_n/<a_n>,                                 line instances,
pi_n:X_n->M_n,       lambda_n:X_n->L_n.           (1.3)
```

Then `pi_n r_n=pi_n`, and the exact object is

```text
                         pi_n             lambda_n
                 M_n  <-------  X_n  ------------>  L_n.   (1.4)
```

Tilings are its named, oriented half-edges.  A line `ell=[t]` has boundary

```text
partial_n(ell)={pi_n(t),pi_n(a_nt)} in Sym^2(M_n).         (1.5)
```

The two entries may agree.  A loop is therefore a free two-cycle of tilings
inside one node fibre, not a fixed tiling or a drawing convention.  Parallel
lines remain distinct elements of `L_n`.

For `u in M_n`, let `F_n(u)=pi_n^{-1}(u)`.  THM-781 gives the intrinsic inverse

```text
F_n(u)=disjoint union_(c in C(u)) HP(T_c)/Aut(T_c),
|F_n(u)|=sum_(c in C(u)) H(T_c)/|Aut(T_c)|.                (1.6)
```

Thus a node returns a fibre of observer-path orbits, never a fictitious unique
tiling.  In particular `a_n` does not descend to a node involution: tilings in
one node fibre can have complements in several different nodes.

For a colour `c` and unordered node pair `{u,v}`, define the literal line fibre

```text
G_n^c(u,v)={ell in L_n:partial_n(ell)={u,v}, colour(ell)=c}. (1.7)
```

A line is blue when `r_n(t)=t`, and black otherwise.  This is endpoint-
independent because `a_n,r_n` commute.  The browser explorer exposes the maps

```text
tilingToMergedNode, mergedNodeToTilings,
tilingToComplementLine, complementLineToTilings,
complementLineToMergedNodes, mergedNodeToComplementLines,
mergedNodePairToComplementLines.                           (1.8)
```

The executable browser atlas covers `n=3..6`; the same relations through
`n=7` are stored in the theorem JSON.

## 2. Exact two-face pullback of tilings

For `n>=4`, let `d_L:X_n->X_(n-1)` delete path endpoint `1`, and let `d_H`
delete endpoint `n`.  Let `alpha(t)=t_(n,1)` be the apex bit.  Then

```text
t |-> (d_Lt,d_Ht,alpha(t))                                 (2.1)
```

is a bijection

```text
X_n ~= (X_(n-1) x_[X_(n-2)] X_(n-1)) x F_2(apex),         (2.2)
```

where `(x_L,x_H)` is compatible exactly when

```text
d_H(x_L)=d_L(x_H).                                        (2.3)
```

The common value is the interior tiling on vertices `{2,...,n-1}`.  It fills
the overlap; the two faces fill the disjoint legs; and `alpha` fills the only
missing tile `(n,1)`.  This proves both injectivity and surjectivity.

The pullback is equivariant:

```text
a_n(x_L,x_H,alpha)=(a_(n-1)x_L,a_(n-1)x_H,1-alpha),
r_n(x_L,x_H,alpha)=(r_(n-1)x_H,r_(n-1)x_L,alpha).          (2.4)
```

In particular,

```text
d_s a_n=a_(n-1)d_s,
d_L r_n=r_(n-1)d_H,
d_Ld_H=d_Hd_L.                                            (2.5)
```

The common core is exactly THM-793's Mode-B restriction.  Fixing it leaves
`2(n-3)+1=2n-5` free leg/apex bits.  This is a fixed-path-coordinate pullback,
not a pullback of unmarked tournament classes; deleting a path endpoint does
not mean deleting a tournament source or sink.

## 3. The line tower is a phase torsor

Every complement line has a unique endpoint with apex bit zero.  Taking that
endpoint in (2.2) gives the canonical set bijection

```text
L_n ~= X_(n-1) x_[X_(n-2)] X_(n-1).                       (3.1)
```

Each face descends to an honest map

```text
bar d_L,bar d_H:L_n->L_(n-1),                              (3.2)
```

and every lower line has exactly `2^(n-2)` upper lifts through either face.
Equation (3.1) does not make upper lines a fibre product of lower **lines**.
For every `n>=5`,

```text
(bar d_L,bar d_H):L_n ->
L_(n-1) x_[L_(n-2)] L_(n-1)                               (3.3)
```

is uniformly two-to-one: a canonical `C_2` endpoint-phase torsor.  In the
apex-zero model, the deck operation simultaneously complements both tiling
faces.  Equality of their bare core lines admits exactly two coherent endpoint
choices.  At `n=4` the empty core is degenerate and the fibre has size four.

There is an explicit recursive phase address.  Let `c=d_H(x_L)=d_L(x_H)` and
let `epsilon(ell)` be the first fixed coordinate of this common core.  Its
value changes under simultaneous complement, so

```text
ell |-> (bar d_Lell,bar d_Hell,epsilon(ell))                (3.4)
```

is a bijection onto compatible lower-line pairs plus one phase bit.  The deck
mate can equivalently be obtained by flipping only the apex of the apex-zero
upper endpoint and then taking its complement line.  Explicitly, if `s_0(ell)`
is the apex-zero endpoint and `e_(n,1)` is the apex basis vector, then

```text
rho_n(ell)=lambda_n(s_0(ell)+e_(n,1)).                    (3.5)
```

Thus the information
lost by bare-line recursion is exactly a simultaneous endpoint phase, not an
unnamed node scalar.

The exact torsor audit is

| n | upper lines | lower-line-pair support | fibre histogram | failures |
|---:|---:|---:|---:|---:|
| 5 | 32 | 16 | `2:16` | 0 |
| 6 | 512 | 256 | `2:256` | 0 |
| 7 | 16,384 | 8,192 | `2:8192` | 0 |

## 4. Lines are a quotient code and colour is zero defect

The line sort is the binary vector space

```text
L_n=F_2^{m_n}/<1>.                                        (4.1)
```

Let

```text
f_n=floor((n-1)/2)                         fixed tiles,
h_n=(m_n+f_n)/2=floor((n-1)^2/4)           reflection orbits,
B_n=Fix(r_n)/<1>,
D_n=im(1+r_n).                                           (4.2)
```

The reflection defect

```text
delta_n([t])=t+r_n(t)                                     (4.3)
```

is representative-independent and gives the short exact sequence

```text
0 -> B_n -> L_n --delta_n-> D_n -> 0,
dim B_n=h_n-1,          dim D_n=(m_n-f_n)/2.               (4.4)
```

Its kernel statement is the definition of `Fix(r_n)`; surjectivity is by the
definition of the image.  Each reflection two-cycle contributes one image
coordinate and each fixed tile contributes none, proving the dimensions.

Blue lines are exactly the zero-defect fibre.  Every defect fibre has size

```text
b_n=2^(h_n-1),                                            (4.5)
```

and the number of black lines is

```text
k_n=2^(m_n-1)-b_n.                                        (4.6)
```

There is no setwise-fixed black line.  Every image of `1+r_n` vanishes on the
reflection-fixed tile coordinates, whereas the all-one word does not, so
`r_n(t)=a_n(t)` is impossible.  Reflection acts with singleton line orbits
exactly on blue lines and free two-element orbits on black lines.

Since `pi_nr_n=pi_n`, the two lines in a black reflection orbit have the same
node boundary and the same nonzero defect.  Hence, for every `n`,

```text
every black node-pair multiplicity is even, loops included;
every black half-edge degree is even;
both statements remain true after conditioning on an exact defect.         (4.7)
```

A blue tiling represents a self-converse class.  Conversely, reflection acts
on the odd THM-781 fibre `HP(T)/Aut(T)` of a self-converse node and therefore
has an odd number of fixed tilings.  On a non-self-converse merged node it
exchanges the two class fibres and has none.  Thus

```text
u is self-converse      iff its blue half-edge degree is odd,
u is non-self-converse  iff its blue half-edge degree is zero.              (4.8)
```

Pure-black nodes are exactly the non-self-converse nodes; self-converse nodes
are pure-blue or mixed.

## 5. Projected incidence and the reversible node shadow

Define the incidence multiplicity

```text
J_n(u,ell)=#{t in ell:pi_n(t)=u} in {0,1,2}.               (5.1)
```

Then

```text
sum_u J_n(u,ell)=2,
sum_ell J_n(u,ell)=|F_n(u)|.                               (5.2)
```

For colour `c`, define the directed half-edge kernel

```text
A_n^c(u,v)=#{t in F_n(u):pi_n(a_nt)=v,colour([t])=c}.      (5.3)
```

Complementing `t` proves

```text
A_n^c(u,v)=A_n^c(v,u),
sum_(v,c) A_n^c(u,v)=|F_n(u)|,
A_n^c(u,u) is even.                                       (5.4)
```

The literal line multiplicities satisfy

```text
|G_n^c(u,v)|=A_n^c(u,v)       if u!=v,
|G_n^c(u,u)|=A_n^c(u,u)/2.                                (5.5)
```

Therefore node fibre size is coloured multigraph degree with loops counted
twice.  The row-normalized total kernel

```text
P_n(u,v)=sum_c A_n^c(u,v)/|F_n(u)|                         (5.6)
```

is reversible with stationary mass `|F_n(u)|/2^m_n`.  It reconstructs fibre
volume and weighted coloured incidence, but not literal parallel-line
identity, class-sheet holonomy, or recursive ancestry.

Three loop notions must remain distinct: an ordinary class-self line, a loop
only after converse merging, and a simple support loop.  At `n=6` there are 8
ordinary class-self lines but 26 merged-node loops, split into 2 blue and 24
black.

## 6. Faces give a weighted node correspondence, not a parent map

Faces do not descend to functions `M_n->M_(n-1)`.  The exact node object is
the span

```text
M_n <-pi_n- X_n -pi_(n-1)d_s-> M_(n-1)                   (6.1)
```

or its matrix

```text
R_n(u,v)=#{t in F_n(u):pi_(n-1)(d_Lt)=v}.                 (6.2)
```

It obeys

```text
sum_v R_n(u,v)=|F_n(u)|,
sum_u R_n(u,v)=2^(n-2)|F_(n-1)(v)|.                       (6.3)
```

Reflection exchanges the faces, so the `d_L,d_H` matrices agree after
converse merging.  This is equality of weighted correspondences, not equality
of maps on tilings.

Divide a nonzero row by the gcd of its entries, key the entries by the already
ordered lower nodes, and call the primitive vector `Rbar_n(u)`.  This is an
objective recursive face coordinate.  It complements, rather than replaces,
THM-785's horizontal `C3` coordinate.

| n | nodes | support-row cells | weighted-row cells | primitive-row cells |
|---:|---:|---:|---:|---:|
| 4 | 3 | 2 | 3 | 3 |
| 5 | 10 | 7 | 10 | 10 |
| 6 | 34 | 34 | 34 | 34 |
| 7 | 272 | 264 | 272 | 272 |

Thus relative face multiplicities distinguish all audited nodes, though this
is not a general injectivity theorem.  HYP-6865 supplies an independent
horizontal coordinate: harmonic voltage on the unmerged local-flip resistor
network is perfectly concordant with score variance through `n=6` and has
`92502/92634` pairwise concordance in the floating `n=7` audit.  The Smith
graph is not the complement-line graph; voltage, `C3`, and `Rbar_n` are
related coordinates, not an asserted graph identity.

## 7. The universal node/line/face coupling tensor

For an upper line choose its apex-zero endpoint `t` and write

```text
u =pi_n(t),       u'=pi_n(a_nt),
l =pi_(n-1)(d_Lt),l'=pi_(n-1)(a_(n-1)d_Lt),
h =pi_(n-1)(d_Ht),h'=pi_(n-1)(a_(n-1)d_Ht).               (7.1)
```

Let

```text
Xi_n(u,u';l,l';h,h';ULH)                                 (7.2)
```

count upper lines with these ordered endpoint-node pairs and upper/low/high
colour word.  Then:

- summing lower data gives the upper line fibres `G_n^c`;
- adding the two endpoint contributions gives `R_n`, for either face;
- summing by `ULH` gives the colour atoms in Section 9;
- reflection swaps the low/high data and the last two colour letters.

Thus `Xi_n` is one joint count-valued object whose marginals cannot disagree.
It remains a quotient: several literal lines may occupy one tensor cell.

| n | lines | nonempty `Xi` cells | multiplicity histogram | collision cells | max |
|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 4 | `1:4` | 0 | 1 |
| 5 | 32 | 32 | `1:32` | 0 | 1 |
| 6 | 512 | 509 | `1:506,2:3` | 3 | 2 |
| 7 | 16,384 | 16,031 | `1:15704,2:309,3:10,4:8` | 327 | 4 |

All marginals and the reflection law have zero failures.  The actual pair of
lower lines plus the coherent phase bit from (3.4) is exactly invertible.

For a single face, a representative-free coupling retains

```text
(u_0,u_1;v_0,v_1)~(u_1,u_0;v_1,v_0).                    (7.3)
```

This simultaneous swap records which upper half-edge descends to which lower
half-edge.  Independently sorting `(u_0,u_1)` and `(v_0,v_1)` identifies
straight and crossed couplings and is only a marginal:

| n | simultaneous-coupled support | independently sorted support |
|---:|---:|---:|
| 4 | 4 | 3 |
| 5 | 31 | 30 |
| 6 | 464 | 455 |
| 7 | 15,112 | 15,074 |

## 8. Complete static node fingerprints are not recursive states

For `s in F_(n-1)(v)`, define the individual lift count

```text
Q_(u,v)(s)=#{t in F_n(u):d_Lt=s}.                          (8.1)
```

Node-level strong lumpability would require this to be constant over every
lower node fibre.  It fails:

| n | nonuniform nonzero blocks | all nonzero blocks | max lift range |
|---:|---:|---:|---:|
| 4 | 0 | 5 | 0 |
| 5 | 11 | 19 | 3 |
| 6 | 76 | 112 | 6 |
| 7 | 1,206 | 1,312 | 6 |

There is a second, coarser manifestation.  Row-normalizing `R_n` and
multiplying consecutive node kernels assumes that reached tilings are uniform
inside the middle node.  Direct two-step deletion disagrees with this Markov
product as follows:

| n | unequal high/low entries | erroneous high rows | maximum error |
|---:|---:|---:|---:|
| 5 | 16 | `8/10` | `1/2` |
| 6 | 86 | `32/34` | `1/2` |
| 7 | 1,778 | `270/272` | `1/2` |

The lower tiling, its Hamiltonian-path/automorphism orbit, or a genuinely
continuation-equivalent refinement must survive.  A node is a base address,
not the whole recursive stalk.

## 9. Closed one-face and three-face colour recursion

### 9.1 One face

For `n>=4`, count lines by `(upper colour,deleted colour)`.  With high colours
as rows and low colours as columns,

```text
BB=2^(n-3),
BK=b_n-BB,
KB=2^(n-2)b_(n-1)-BB,
KK=2^(n-2)k_(n-1)-BK.                                   (9.1)
```

An upper line and one face are both blue precisely when the mask is constant
on every difference diagonal `x-y=d`.  The two relevant reflections compose
to unit translation along each diagonal.  There are `n-2` diagonal bits and
complement identifies opposite words, proving `BB=2^(n-3)`; the other entries
follow from row and column totals.  The difference-striped/Toeplitz subspace
is the recursive blue-to-blue spine.

### 9.2 Both faces

Let

```text
T=2^(m_n-1),
U=2^(h_n-1),
F=2^(h_(n-1)+n-3),
Q=2^(n-3+floor((n-2)/2)),
J=2^(n-3).                                                (9.2)
```

For the word `(upper,low-face,high-face)`, the exact atoms are

```text
BBB=J,
BKK=U-J,
KBB=Q-J,
KBK=KKB=F-Q,
KKK=T-U-2F+Q+J,
BBK=BKB=0.                                               (9.3)
```

Both-face symmetry is parity-constant on each gap diagonal; adding upper
symmetry joins the two parity components and makes the diagonal constant.
This gives `Q` and `J`; reflection gives left/right equality and inclusion-
exclusion gives the other atoms.

If `beta_n=U/T`, then

```text
F/T=beta_(n-1),
Q/T=beta_(n-1)^2,
J/T=beta_n beta_(n-1).                                   (9.4)
```

Upper-, low-face-, and high-face-blue are pairwise independent.  They are not
jointly independent: upper blue forces the two face colours equal.  The first
new colour datum across sizes is therefore a pure three-way interaction.

| n | BBB | BKK | KBB | KBK | KKB | KKK |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| 5 | 4 | 4 | 4 | 8 | 8 | 4 |
| 6 | 8 | 24 | 24 | 96 | 96 | 264 |
| 7 | 16 | 240 | 48 | 960 | 960 | 14,160 |

The exact pre-quotient symmetry locates THM-785's black left/right imbalance
in disintegration over unequal node fibres rather than in the raw line tower.

## 10. Mode-B two-end deletion is the defect tower

For `n>=5`, let

```text
p_n=d_Ld_H:X_n->X_(n-2).                                  (10.1)
```

It commutes with complement and reflection, and induces surjective linear maps
on lines, blue subspaces, and defect spaces.  Boundary coordinates give

```text
0 -> F_2^(2n-5) -> L_n -> L_(n-2) -> 0,
0 -> F_2^(n-2)  -> B_n -> B_(n-2) -> 0,
0 -> F_2^(n-3)  -> D_n -> D_(n-2) -> 0.                  (10.2)
```

Surjectivity follows by freely extending an interior word, symmetrically when
restricted to `B_n`.  The kernel dimensions are the differences of the line,
blue, and defect dimensions in (4.4).  Each lower line has `2^(2n-5)` lifts,
and the Mode-B colour channel is triangular:

```text
BB=b_n=2^(n-2)b_(n-2),
BK=0,
KB=(2^(2n-5)-2^(n-2))b_(n-2),
KK=2^(2n-5)k_(n-2).                                      (10.3)
```

A nonzero lower defect forces every lift to remain black: inherited
blackness.  Over a blue lower line, a nonzero defect in the kernel creates
fresh blackness.  Thus black lines have an exact recursive ancestry, not only
a binary colour.

## 11. Loops carry class-sheet holonomy

Let `c_n(t) in I_n` be the ordinary unmerged class.  If `[t]` is a merged loop,
then either

```text
same sheet:       c_n(t)=c_n(a_nt),
converse switch:  c_n(a_nt)=c_n(t)^op != c_n(t).          (11.1)
```

This is endpoint-independent.  A blue loop is necessarily same-sheet; a black
loop may have either holonomy.  Neither loop status nor its sheet status is
hereditary.  Under Mode-B `6->4`, the exact line transitions are

```text
cross/cross 366,  cross/loop 120,
loop/cross   18,  loop/loop    8.                         (11.2)
```

At `n=6` the 24 black loops split into 6 same-sheet and 18 converse-switch
lines; at `n=7` the 114 black loops split `44+70`.

The full finite census is:

| n | tilings | classes | merged | `B[cross,loop]` | `K[cross,loop]` | coloured/plain support | loop holonomy `Bs,Bc,Ks,Kc` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | `1,0` | `0,0` | `1/1` | `0,0,0,0` |
| 4 | 8 | 4 | 3 | `1,1` | `2,0` | `3/3` | `1,0,0,0` |
| 5 | 64 | 12 | 10 | `8,0` | `20,4` | `20/20` | `0,0,4,0` |
| 6 | 1,024 | 56 | 34 | `30,2` | `456,24` | `187/183` | `2,0,6,18` |
| 7 | 32,768 | 456 | 272 | `256,0` | `16,014,114` | `6,126/6,076` | `0,0,44,70` |

The Mode-B line census is:

| map | line degree | blue lifts/blue parent | `BB,BK,KB,KK` | `EE,EL,LE,LL` |
|---:|---:|---:|---:|---:|
| `5->3` | 32 | 8 | `8,0,24,0` | `28,0,4,0` |
| `6->4` | 128 | 16 | `32,0,224,256` | `366,120,18,8` |
| `7->5` | 512 | 32 | `256,0,3840,12288` | `14244,2026,92,22` |

All pullback, torsor, incidence, tensor, defect, parity, colour, loop, and
sheet-transport assertions pass with zero failures.

## 12. Groupoids and the continuation-complete address

Let `P_n` be the finite groupoid of tournaments with a directed Hamiltonian
path.  Its skeleton is `X_n`.  Forgetting the path gives the tournament
groupoid, folding converse gives `M_n`, and quotienting the free complement
action gives `L_n`.  The span (1.4) is a coloured multigraph internal to finite
groupoids.

Terminal path deletion is an honest functor on `P_n`, but is not cartesian
after forgetting the path.  Arbitrary unmarked vertex deletion is instead a
span through `(T,v)` and must retain the `Aut(T)` orbit and stabilizer of `v`.
Deleting an internal path vertex also needs a repaired Hamiltonian path of
`T-v`, because the shortcut can point backward.

Deletion and quotient descent are not invertible and therefore are not arrows
of the fixed-size groupoid.  The graded tower is more precisely a
**groupoid-valued recursion**: each fixed stratum has its reversible symmetry
groupoid, while endpoint deletion, arbitrary deletion spans, and certified
quotients are functors or correspondences between strata.  A well-founded
condensation requires a separately proved complexity decrease.

The information relation proved above is not a total hierarchy.  On the
literal-line domain, write `P <= Q` when carrier `Q` refines carrier `P`.  The
relevant finite lattice contains

```text
node boundary <= colour <= class-sheet packet <= defect packet
       |                                             |
       +-----------> Xi -----------------------------+--> Xi+sheet+defect

independently sorted face transport <= simultaneous coupling <= Xi

Xi+sheet+defect <= exact line ~= (compatible lower lines, coherent phase).
                                                                    (12.1)
```

Several arrows are deliberately absent.  The weighted node kernel and face
matrix aggregate over many lines, so they are linear shadows, not refinements
of one line partition.  A marked tiling is an oriented witness above a line,
not merely the next cell in the lattice.  An owner-labelled metric LRC packet
is a stalk over the combinatorial atlas, not a maximum element on the same
domain.  This typed partial order prevents “more coordinates” from being
confused with “enough information for the next operation.”

The coarsest exact bounded recursive address is naturally Nerode-like: two
marked states agree only if every continuation word in
`{d_L,d_H,a,r}` has the same declared terminal observations.  HYP-6825's
weighted line-WL and primitive face row are strong finite approximations, not
proofs of continuation equivalence for arbitrary `n`.

### 12.1 Finite-exact continuation minimization

This prescription can be computed without guessing an address.  On

```text
L_[3,7]^tot=disjoint union_(3<=k<=7) L_k,                 (12.2)
```

start with equality of a declared current observation `O`.  Repeatedly refine
a cell by the current cells of `r(ell),d_L(ell),d_H(ell)`, using a terminal
symbol when a face is undefined.  The stable partition is the coarsest
refinement of `O` stable under every legal bounded continuation.  A second,
different predicate retains only the unordered multiset of the two face
successor cells.  It still retains `O` and unary reflection, so it is
unordered-successor semantics, not a quotient of current states by
reflection.  In particular, ordered face fields already present in `O` stay
ordered.  Complement is already the identity on `L_k`.

The exact `n=7` minimization is:

| current observation | initial cells | stable labelled faces | stable unordered successors |
|---|---:|---:|---:|
| node boundary | 6,076 | 16,359 | 8,310 |
| + blue/black colour | 6,126 | 16,359 | 8,310 |
| + loop class sheet | 6,128 | 16,359 | 8,310 |
| + exact reflection defect | 8,091 | 16,359 | 8,311 |
| `Xi_n` joint node/face/colour cell | 16,031 | 16,382 | 16,380 |
| `Xi_n` + sheet + defect | 16,270 | 16,382 | 16,380 |
| lower-line pair + coherent phase | 16,384 | 16,384 | 16,384 |
| literal line | 16,384 | 16,384 | 16,384 |

These are partition equalities where the numbers coincide in the nested
rows, not accidental equal cell counts.  Thus colour, sheet, and defect supply
no further distinction after labelled recursive node-boundary refinement.
That common stable partition still has 23 nonsingleton cells: 22 pairs and one
four-line cell, containing 48 lines total (25 excess representatives).
`Xi_n` leaves only two collision pairs.  They are

```text
0x12ca / 0x12cb: lower lines (53,150), epsilon 1/0, defect 0x06a6,
0x146c / 0x146d: lower lines (150,163), epsilon 1/0, defect 0x06a6. (12.3)
```

Each pair consists of the two deck mates from (3.5): the canonical masks
differ by the apex flip, while the coherent common-core phase changes.  Their
reflection defects agree because the apex is fixed by `r_n` and is therefore
annihilated by `1+r_n`.  This makes phase and defect transverse coordinates:
even exact defect ancestry cannot replace the torsor phase.

The four lines form one commuting involution square:

```text
rho_n=(12ca 12cb)(146c 146d),
r_n  =(12ca 146c)(12cb 146d).                            (12.4)
```

All are parallel black cross-lines between `n7-a264`, a non-self-converse node
with two ordinary class-sheet fibres of size 151, and `n7-a270`, a
self-converse node with one fibre of size 57.  The deck involution changes the
marked Hamiltonian-path presentation inside fixed ordinary endpoint classes;
reflection exchanges the two 151-element converse sheets.  Thus the first
failure of recursively refined `Xi_n` is precisely observer-path descent data,
not node, class, colour, sheet, or defect data.

On the rigid 151-presentation endpoint, each deck-pair class isomorphism is
unique.  It fixes the two Hamiltonian-path endpoints and acts as a 5-cycle on
the five interior vertices; the two 5-cycles are conjugate by path reflection:

```text
(0,2,3,5,1,4,6),       (0,2,5,1,3,4,6).                 (12.5)
```

The self-converse endpoint admits three isomorphisms; exactly one fixes both
path ends, and it is the inverse 5-cycle.  No single relabelling maps both
endpoints of the line simultaneously.  The relative presentation holonomy
therefore has order five.  Consequently the first hidden phase is an
endpoint-fixed cyclic relabelling of the interior path, not merely a collision
of canonical codes.

The two common `n=5` cores are the reflected black same-class loops `0x3` and
`0x9` at node 6.  Their compatible one-end loop extensions glue to the four
`n=7` lines.  Conversely, an ordered `Xi_(n+1)` phase collision requires both
of its size-`n` endpoint-face lines to be loops, because the chord-dual deck
move complements each face and reverses its ordered endpoint pair.  The four
residual `n=7` lines are cross-lines, so this branch cannot lift directly to
an `n=8` collision.  Larger phase-blind families may reappear around different
interior cores; collision ancestry is a birth/death process, not a nested
subgraph.

Making only the successor refinement unordered creates two additional `Xi_n`
collision pairs with the face-line order exchanged.  The current `Xi_n`
observation remains ordered; this is not the full reflection quotient of line
states.

The minimization is finite-exact only.  It does not include extension words,
internal-vertex deletion/repair, metric LRC observations, or `n>=8`, and it
does not assert that the displayed partitions stabilize with `n`.

### 12.2 Chord duality and the three-face Cech resolution

Let `A_n` toggle only the apex tile and put

```text
zeta_n=a_n A_n.                                            (12.6)
```

For an apex-zero endpoint, the phase deck mate in (3.5) has apex-zero
representative `zeta_n(t)`.  In the polygon determined by the distinguished
Hamiltonian path, `zeta_n` preserves the `n` boundary edges (the path plus the
apex) and reverses every interior chord.  Moreover

```text
d_L zeta_n=a_(n-1)d_L,       d_H zeta_n=a_(n-1)d_H.        (12.7)
```

Thus the two phase mates always have the same bare endpoint-face lines.  Their
ordered `Xi_n` face-node pairs agree exactly when those lower lines are loops;
the upper-node/class conditions then decide whether a collision survives.
This identifies the missing datum as a chord-duality witness in the marked-
path isomorphism groupoid, modulo endpoint automorphisms and compatible face
restrictions.

Concurrent THM-801 adds the third gap-contraction face `B`.  Its three faces
cover the apex, and their nonempty triple overlap forces the relative endpoint
phases to satisfy the Cech cocycle equation.  Literal three-face line descent
is therefore exact for every `n>=6`: the two-face torsor becomes ordinary
overlap descent data in the finer cover.

On the residual square this resolution is completely visible:

| upper line | high face | gap face | low face | gap incidence |
|---:|---:|---:|---:|---|
| `0x12ca` | `0x35` | `0x115` | `0x96` | black loop `33--33` |
| `0x12cb` | `0x35` | `0x114` | `0x96` | black cross `21--33` |
| `0x146c` | `0x96` | `0x0c3` | `0x0a3` | black loop `33--33` |
| `0x146d` | `0x96` | `0x0c2` | `0x0a3` | black cross `21--33` |

The two endpoint faces are unchanged inside a phase pair; only the gap face
changes.  Here the coherent phase bit is exactly the indicator that the gap
line closes as a loop.  THM-801's node-compressed `Omega_7` already separates
all four lines, without its additional mirror `B2` sidecar.  Hence “phase is
necessary” is cover-relative: it is irreducible for the two-end cover and is
reconstructed as gap-face incidence in the full three-face cover.

## 13. Tournament Analysis and the LRC preservation boundary

The native relation is symmetric, weighted, coloured, and looped.  Forcing the
three sorts into a single binary tournament destroys the theorem's payload.
The computation instead uses information carriers as Tournament Analysis
vertices:

```text
fibre_size, colour_degree, line_support_row, line_weighted_row,
lower_face_support, lower_face_weighted, lower_face_normalized, exact_node.
```

The pairwise observable is separated unordered node pairs.  The switches are
retention and retention per partition cell, with the displayed list as the tie
Hamiltonian path.  At `n=7` both carrier tournaments are transitive, have score
histogram `{0:1,...,7:1}`, zero directed triangles, singleton SCCs, and one
Hamiltonian path; the gauges flip 18 edges.

For the bounded continuation audit the vertices are the eight information
carriers in the rows of the table in Section 12.1, while literal `n=7` lines
are the sample states.  The pairwise observable is the sign of the difference
in separated line-pair utility.  Total separation and separation per cell are
the switches; the table order breaks ties.  Both scalar-ranked tournaments are
transitive with score histogram `{0:1,...,7:1}`, zero directed triangles,
singleton SCCs, and one Hamiltonian path, with 20 edge flips.  The retention
path is `phase,exact,Xi,Xi+defect,node,colour,sheet,defect`; the economy path is
the displayed carrier order.  Transitivity follows from scalar ranking, so
these fingerprints are consistency telemetry rather than a structural
tournament theorem.

The challenged assumption is that vertices must be runners or arcs.  Depending
on the predicate, viable vertices include node classes, Hamiltonian paths,
line instances, defects, gaps, fixed sections, boundaries, wall events,
residues, cover arcs, Fourier modes, matroid circuits, and proof obligations.
Every quotient must name its preserved predicate and destroyed coordinates.

For LRC(14), this theorem does not identify a metagraph node with HYP-6815's
four-coordinate suspension

```text
X_(A,R)={(u,t,c,lambda):u=ct and Phi_(A,R)(u,t)>=lambda}.  (13.1)
```

The three-sorted object is a constructible combinatorial atlas over the
suspension's phase-order strata.  The LRC predicate is nonemptiness of every
integer-slope section at `lambda=1/14`.  A node quotient destroys metric gap
widths, observer phase, runner owners, exact wall side, scale/residue, inverse
winding, endpoint ties, chronology, and carry.  These remain stalk fields
unless fibre-purity, reconstruction, annihilation, or a named residual theorem
proves that a field can be discarded.

The class-sheet bit in (11.1) must not be mistaken for the full LRC holonomy.
The concurrent component-obligation synthesis distinguishes at least

```text
h_class : same/converse tournament sheet of one merged line loop,
h_red   : reduced combinatorial sheet/token return modulo a declared gauge,
Delta_M : metric translation of the component and endpoint phases.          (13.2)
```

A loop may have trivial `h_red` and nonzero `Delta_M`; THM-794 supplies such
repeated packet transport.  THM-796 determines `h_class` only.  Any pullback to
the LRC suspension must therefore attach the other two coordinates rather than
calling a closed finite-state loop metrically closed.

The exact preservation rule is therefore recursive: keep precisely enough
data to make the next intended operation well-defined, recurse on literal
lines and named half-edges before quotienting, refine nodes by their incident
recursive signatures, and attach this finite atlas to the LRC suspension
rather than replacing the suspension by its shadow. ∎
