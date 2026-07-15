---
id: THM-796
title: Three-sorted recursive incidence of tilings, complement lines, and converse-merged tournament nodes
status: PROVED — general two-face tiling pullback, C2 line-phase torsor, colored incidence conservation, and blue/black face recursion; exact node-kernel and non-lumpability census n=3..7
source: codex-2026-07-15-S9
depends_on: [THM-280, THM-345, THM-643, THM-781, THM-793]
related: [THM-477, THM-785, THM-790, HYP-6825, HYP-6870]
verification:
  - 04-computation/three_sorted_metagraph_recursion_codex_S9.py
  - 05-knowledge/results/three_sorted_metagraph_recursion_codex_S9.out
  - 05-knowledge/results/three_sorted_metagraph_recursion_codex_S9.json
  - 03-artifacts/visualizations/tournament-tiling-explorer.html
  - 04-computation/test_tournament_tiling_explorer_line_api_codex_S9.js
---

# THM-796 — the three-sorted recursive metagraph

The merged metagraph has three different sorts.  They must not be collapsed:

```text
tiling endpoints T_n  --free complement quotient-->  lines E_n
       |                                           /        \
       | pi_n                                     / colour   \ endpoint nodes
       v                                         v            v
converse-merged nodes N_n  <---------------- projected coloured multigraph
```

The theorem gives exact forward maps and inverse fibres between these sorts,
then makes their recursion in `n` explicit.  Its main information-boundary
conclusion is that the recursion closes on tilings, becomes a two-sheeted
endpoint-phase torsor on bare lines, and becomes a weighted correspondence—but
not a Markov state—on merged nodes.

## 1. The three sorts and the maps

Put

```text
S_n = {(x,y): 1 <= y < x <= n and x-y >= 2},
M_n = |S_n| = C(n-1,2),
T_n = {0,1}^{S_n}.
```

These are the explorer tilings above the fixed Hamiltonian path
`n -> n-1 -> ... -> 1`.  There are two commuting operations on `T_n`:

```text
kappa_n(t)_(x,y) = 1-t_(x,y),
sigma_n(t)_(n+1-y,n+1-x) = t_(x,y).
```

Here `kappa` flips every non-path tile, while `sigma` is grid reflection.
They are different involutions.  By THM-280, `sigma` induces tournament
converse up to relabelling; `kappa` is the antipodal complement-line operation
inside the fixed-path cube.

Define

```text
E_n = T_n/<kappa_n>,
N_n = tournament isomorphism classes modulo converse,
q_n : T_n -> E_n,
pi_n : T_n -> N_n.
```

A line `e={t,kappa(t)}` is blue when `sigma(t)=t` and black otherwise.  This is
endpoint-independent because `sigma` and `kappa` commute.  For `u in N_n`, set

```text
F_n(u) = pi_n^{-1}(u).
```

THM-781 supplies the exact inverse without scanning the cube:

```text
F_n(u) = union_[T in u] HP(T)/Aut(T),
|F_n(u)| = sum_[T in u] H(T)/|Aut(T)|.
```

Thus a node maps to a set of tilings, never to a fictitious unique tiling.

For a colour `c` and an unordered node pair `{u,v}`, define the exact line
fibre

```text
G_n^c(u,v) = {
  {t,kappa(t)} in E_n :
  {pi_n(t),pi_n(kappa(t))}={u,v} and colour(t)=c
}.
```

These are the inverse fibres of the projected blue/black multigraph.  Parallel
lines and loops remain literal line instances.  The browser explorer now
exposes the executable maps

```text
tilingToMergedNode, mergedNodeToTilings,
tilingToComplementLine, complementLineToTilings,
mergedNodePairToComplementLines.
```

The line descriptor keeps both endpoints even for a loop and records incidence
multiplicity two at its single node.

## 2. Exact two-face pullback of tilings

Let `d_L:T_n->T_(n-1)` delete path endpoint `1`: keep tiles with `y>1` and
shift `(x,y)` to `(x-1,y-1)`.  Let `d_H` delete path endpoint `n`: keep tiles
with `x<n`.  Let

```text
a(t)=t_(n,1)
```

be the apex bit.  Then, for every `n>=4`,

```text
Phi_n(t) = (d_L(t),d_H(t),a(t))
```

is a bijection

```text
T_n  ~=  (T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1},              (2.1)
```

where `(L,H)` belongs to the fibre product precisely when

```text
d_H(L)=d_L(H).                                                (2.2)
```

The common value is the interior tiling on vertices `{2,...,n-1}`.  The
inverse fills the common interior from (2.2), the two disjoint legs from `L`
and `H`, and the only still-missing tile `(n,1)` from the apex bit.  Hence
reconstruction is unique.

The bijection is equivariant in the exact forms

```text
kappa_n(L,H,a) = (kappa_(n-1)L,kappa_(n-1)H,1-a),
sigma_n(L,H,a) = (sigma_(n-1)H,sigma_(n-1)L,a).                (2.3)
```

The compatible pair has

```text
M_(n-2)+2(n-3)=M_n-1
```

free bits.  Formula (2.1) therefore also proves `|T_n|=2^{M_n}`.  Its common
core is exactly THM-793's Mode-B projection `p:T_n->T_(n-2)`; fixing the core
leaves `2(n-3)+1=2n-5` free leg/apex bits, recovering that theorem's
`2^{2n-5}` fibre law.

This is a fixed-path-coordinate pullback.  It is not a pullback of unmarked
tournament classes, and “delete a path endpoint” does not mean “delete a
tournament source or sink.”

## 3. The line tower is a phase torsor

Every complement line has a unique endpoint with apex bit zero.  Using it in
(2.1) gives a canonical bijection

```text
E_n ~= T_(n-1) x_[T_(n-2)] T_(n-1).                           (3.1)
```

This does **not** say that upper lines are a fibre product of lower lines.  A
face commutes with `kappa`, so it gives an endpoint-independent map

```text
bar d_L,bar d_H : E_n -> E_(n-1).
```

For every `n>=5`,

```text
(bar d_L,bar d_H): E_n ->
E_(n-1) x_[E_(n-2)] E_(n-1)                                  (3.2)
```

is uniformly two-to-one.  It is a canonical `C2`-torsor.  In the apex-zero
model (3.1), its deck involution is

```text
rho_n(L,H)=(kappa_(n-1)L,kappa_(n-1)H).                        (3.3)
```

Indeed, equality of the two bare core lines lets one choose coherent core
endpoints in exactly two ways; the choices differ by simultaneous complement.
They reconstruct two distinct upper lines with the same two lower face-lines.
The deck mate has the same upper blue/black colour because `kappa` commutes
with `sigma`.  At `n=4` the empty `T_2` core is degenerate and the corresponding
fibre has size four, which is why (3.2) begins at `n=5`.

This identifies the exact information lost by bare-line recursion: not a new
node scalar, but a simultaneous endpoint phase.  Keeping an oriented endpoint,
or equivalently the apex-zero section, repairs the loss.

The exhaustive audit gives

| n | upper lines | lower-line-pair support | fibre histogram | failures |
|---:|---:|---:|---:|---:|
| 5 | 32 | 16 | `2:16` | 0 |
| 6 | 512 | 256 | `2:256` | 0 |
| 7 | 16,384 | 8,192 | `2:8192` | 0 |

## 4. Projected coloured incidence matrices

Define the directed endpoint kernel

```text
A_n^c(u,v) = #{t in F_n(u): pi_n(kappa(t))=v and colour(t)=c}. (4.1)
```

For every `n`, node pair, and colour:

```text
A_n^c(u,v)=A_n^c(v,u),                                        (4.2)
sum_(v,c) A_n^c(u,v)=|F_n(u)|,                                (4.3)
A_n^c(u,u) is even.                                           (4.4)
```

The involution `t -> kappa(t)` proves (4.2) and pairs the diagonal endpoints
in (4.4).  Consequently the exact projected line multiplicities are

```text
|G_n^c(u,v)| = A_n^c(u,v)       if u != v,
|G_n^c(u,u)| = A_n^c(u,u)/2.                                 (4.5)
```

Equivalently, the node fibre size is the coloured multigraph degree with loops
counted twice:

```text
|F_n(u)| = sum_(v!=u,c)|G_n^c(u,v)| + 2 sum_c|G_n^c(u,u)|.    (4.6)
```

These equations refine the `d=M_n` transport matrix of THM-345 and the
one-level blue/black allocation laws of THM-643.  A simple graph obtained by
replacing every nonempty `G_n^c(u,v)` by one edge destroys exactly the data in
(4.5): parallel multiplicity, two endpoint incidences of a loop, and possibly
colour.

Three loop notions must therefore remain separate:

1. an unmerged class-self line, whose two endpoints have the same ordinary
   tournament class;
2. a merged-node loop, whose endpoints agree only after converse merging;
3. a simple adjacency loop, which records only existence.

For example, at `n=6` there are 8 unmerged class-self lines but 26 merged-node
loops; the latter split as 2 blue and 24 black.

## 5. Recursive node correspondence and its conservation laws

Faces do not descend to functions `N_n->N_(n-1)`.  The honest object is the
weighted correspondence

```text
D_n(u,v) = #{t in F_n(u): pi_(n-1)(d_L(t))=v}.                 (5.1)
```

It obeys

```text
sum_v D_n(u,v)=|F_n(u)|,                                      (5.2)
sum_u D_n(u,v)=2^{n-2}|F_(n-1)(v)|.                           (5.3)
```

For (5.3), a fixed lower tiling has `n-3` free bits on the other exclusive
leg and one free apex bit.  The `d_H` branching matrix equals (5.1): reflection
swaps the two faces by (2.3), while `pi(sigma(t))=pi(t)` after converse merging.

The row `D_n(u,-)` is a recursive node fingerprint.  Key its entries by the
already ordered lower nodes, divide by the gcd of its nonzero entries, and
call the resulting primitive vector `R_n(u)`.  Starting from the transitive-
first order at `n=3`, the lexicographic order of these vectors is an objective
recursive *face order*.  It is a second coordinate, not a replacement for
THM-785's horizontal transitive-to-distributed `C3` flow order.  The output
stores both the primitive vector and its recursive face rank for every node.

Finite exact separation is unexpectedly strong:

| n | nodes | support-row cells | weighted-row cells | primitive-row cells |
|---:|---:|---:|---:|---:|
| 4 | 3 | 2 | 3 | 3 |
| 5 | 10 | 7 | 10 | 10 |
| 6 | 34 | 34 | 34 | 34 |
| 7 | 272 | 264 | 272 | 272 |

Thus at `n=7` support alone leaves eight twin pairs, whereas relative
multiplicities already distinguish all 272 nodes.  Absolute fibre size is not
needed for this finite separation.  This is a verified `n<=7` completeness
statement, not a proof that `R_n` stays injective for every `n`; canonical
converse-orbit code remains the final exact fallback.

This recursive coordinate complements the horizontal flow picture:

```text
THM-785 C3 flow:       where the node lies from transitive to distributed;
primitive D row:       how its whole path-orbit fibre branches into faces;
coloured A rows:       how its tiling endpoints pair into blue/black lines;
canonical orbit code:  final exact tie-breaker if structural rows collide.
```

## 6. Why a complete node fingerprint is still not a recursive state

For `s in F_(n-1)(v)`, define the lift count

```text
L_(u,v)(s) = #{t in F_n(u): d_L(t)=s}.                         (6.1)
```

A node-only Markov recursion would require (6.1) to be constant over every
lower node fibre `F_(n-1)(v)`—the standard strong-lumpability condition.  It
fails immediately and overwhelmingly:

| n | nonuniform nonzero `(u,v)` blocks | all nonzero blocks | max lift range |
|---:|---:|---:|---:|
| 4 | 0 | 5 | 0 |
| 5 | 11 | 19 | 3 |
| 6 | 76 | 112 | 6 |
| 7 | 1,206 | 1,312 | 6 |

For instance, in the `n=7-a026 -> n=6-a01` block, one lower mask has one lift
while another has three.  Therefore an aggregate row can identify the parent
node perfectly while still failing to tell how an individual child tiling
continues.  Recursive reconstruction must retain the lower tiling, its
Hamiltonian-path/automorphism orbit, or an equivalent rooted incident-word
state.  The node is a base address, not the whole stalk.

## 7. Closed blue/black face recursion

Write

```text
f_n=floor((n-1)/2),
r_n=(M_n+f_n)/2,
T=2^{M_n-1},
U=2^{r_n-1},
L=2^{r_(n-1)+n-3},
Q=2^{n-3+floor((n-2)/2)},
J=2^{n-3}.
```

Choose the apex-zero endpoint of each upper line.  Let a three-letter word
record `(upper colour, low-face colour, high-face colour)`, with `B` blue and
`K` black.  For every `n>=4`, the six possible atoms have the exact counts

```text
BBB = J,
BKK = U-J,
KBB = Q-J,
KBK = KKB = L-Q,
KKK = T-U-2L+Q+J,                                            (7.1)
BBK = BKB = 0.
```

Proof of the nontrivial ranks: fix a constant-gap diagonal and enumerate its
tiles by `y=1,...,q`.  Low- and high-face symmetry impose

```text
y ~ q-y,             y ~ q+2-y.
```

Their compositions translate by two.  The equality graph therefore has one
component for odd `q` and two parity components for even `q`.  Summing over
the diagonals and fixing the isolated apex to zero gives
`n-3+floor((n-2)/2)` free bits, hence `Q`.  Upper symmetry adds

```text
y ~ q+1-y,
```

which connects the two parity components.  Every gap diagonal becomes
constant; after fixing the apex there are `n-3` free diagonal values, hence
`J`.  Upper reflection swaps the two faces, so an upper-blue line forces their
colours to agree.  Inclusion-exclusion gives (7.1).

Geometrically, both-face-blue endpoints are parity-constant on every gap
diagonal, while all-three-blue endpoints are constant on every gap diagonal:
they are precisely the apex-normalized Toeplitz, or distance-stationary,
tilings.

There is a sharper probabilistic form.  Let

```text
beta_n = U/T = 2^{r_n-M_n}.
```

Since `r_n+r_(n-1)=M_n`,

```text
beta_n = beta_(n-1) 2^{-floor((n-2)/2)},
L/T = beta_(n-1),
Q/T = beta_(n-1)^2,
J/T = beta_n beta_(n-1).                                    (7.2)
```

Consequently upper-blue, low-face-blue, and high-face-blue are pairwise
independent under the uniform measure on lines.  They are not jointly
independent: upper blue deterministically makes the two face colours equal,
and the triple-blue probability is `beta_n beta_(n-1)`, not
`beta_n beta_(n-1)^2`.  The first new colour datum across sizes is therefore a
pure three-way interaction.

The exact atom counts are:

| n | BBB | BKK | KBB | KBK | KKB | KKK |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| 5 | 4 | 4 | 4 | 8 | 8 | 4 |
| 6 | 8 | 24 | 24 | 96 | 96 | 264 |
| 7 | 16 | 240 | 48 | 960 | 960 | 14,160 |

This exact left/right symmetry belongs to the tiling/line tower before node
projection.  THM-785's observed black left/right drift is therefore a
disintegration effect: it enters when these symmetric line masses are
allocated among unequal node fibres, loops, and transitive-flow positions.

## 8. Finite three-sort census

| n | tilings `T_n` | lines `E_n` | merged nodes `N_n` | coloured node-pair supports | merged loops |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 2 | 1 | 0 |
| 4 | 8 | 4 | 3 | 3 | 1 |
| 5 | 64 | 32 | 10 | 20 | 4 |
| 6 | 1,024 | 512 | 34 | 187 | 26 |
| 7 | 32,768 | 16,384 | 272 | 6,126 | 114 |

Every pullback, reconstruction, core-compatibility, `kappa` naturality,
reflection face swap, apex-zero section, line-torsor, matrix conservation,
diagonal parity, and colour-formula audit has zero failures through `n=7`.

The computation's Tournament Analysis deliberately uses information carriers,
not runners, as vertices.  Its pairwise observable is the number of merged-node
pairs separated by a carrier; the switches are retained pairs and retained
pairs per partition cell; its tie Hamiltonian path is

```text
fibre_size -> colour_degree -> line_support_row -> line_weighted_row
 -> lower_face_support -> lower_face_weighted -> lower_face_normalized
 -> exact_node.
```

At `n=7` the two gauges flip 18 carrier-tournament edges.

## 9. Preservation boundary and the LRC object

This theorem challenges the assumption that the natural vertices must be
runners or arcs.  Here the useful vertices are alternately tiling endpoints,
complement-line instances, merged nodes, lower-face states, and information
carriers.

The quotient ledger is exact:

| operation | preserves | destroys / required repair |
|---|---|---|
| `T_n -> E_n` | antipodal line, colour, two endpoint set | endpoint phase; repair by apex-zero orientation |
| two lower tilings -> two lower lines | face-line pair and core line | simultaneous coherent phase; repair by the `C2` torsor |
| `T_n -> N_n` | tournament converse-orbit class | chosen Hamiltonian path/observer cut and individual lift distribution |
| lines -> weighted coloured node multigraph | all pair multiplicities and loops | endpoint masks/path orbits |
| weighted multigraph -> simple graph | existence of adjacency | multiplicity, loop factor two, and possibly colour |

None of these tournament-only quotients preserves the LRC loneliness predicate,
which is metric and phase-sensitive.  Pulling the address back to the four-
coordinate LRC slope suspension still requires the owner, metric scale,
threshold, endpoint/tie order, inverse-step, and carry/continued-fraction
transport stalk.  THM-778 supplies the analogous lesson on the arithmetic side:
a static node or continued-fraction digit is safe only together with its action
on the labelled fibre.

The comprehensive structural object is therefore not one graph.  It is the
three-sorted incidence diagram together with its endpoint-phase and LRC
transport stalks.
