---
id: THM-4116
title: "Boundary-state gluing and AP odd-shell tree synchronizers"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Ordered boundary extension vectors glue by an
  exact dot product. On a Petersen adjacent two-vertex side with four-edge
  boundary their supports are disjoint, while the analogous K4 side has dot
  product six. Every uncolorable finite simple cubic graph has a parity- and
  Kempe-compressed one-integer adjacent-edge signature. For AP_(2q-1),
  after the complete unequal-v2 reciprocal graph, an added odd-shell graph F
  leaves exactly 2^(c(F)-1) circle-phase sheets. Saturation is equivalent to
  connectivity; its inclusion-minimal repairs are the q^(q-2) labelled
  spanning trees. Deletion robustness is exactly edge connectivity; the
  minimum one-edge-fault-tolerant repairs are cycles when q>=3. This
  synchronizes phase sheets but does not prove LRC safety or LRC(14).
source: codex-snark-apex-boundary-session-20260825
depends_on:
  - THM-4110-sparse-reciprocal-phase-graph-saturation-and-ap13-torsion-tariff
related:
  - THM-261-petersen-root-orthogonality
  - THM-343-H7-impossible
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
  - MISTAKE-507
script: 04-computation/snark_boundary_ap_shell_thm4116.py
output: 05-knowledge/results/snark_boundary_ap_shell_thm4116.out
script_sha256: a52a71c881459a5ffcee55cba183f8c4b8e33bea6b47ed3e2e63bc5cf2ddbe8f
output_sha256: 23a316c682bd89c99f98676e26c35bf7ad349065e6ae9e4ad13a50e28e1bf321
independent_audit_script: 04-computation/snark_boundary_fourpole_thm4116_independent.py
independent_audit_output: 05-knowledge/results/snark_boundary_fourpole_thm4116_independent.out
independent_script_sha256: 33784fcd273cbdeefe3e18569aee48df0acfe0f6f458b684922c6fa18d8bf105
independent_output_sha256: 059a8d93fc40a99d341e76194a68283eef38d0118ab39ebbdef419897d76489b
hash_basis: raw LF bytes
audit: >
  PASS. An independent boundary enumerator checked every simple graph through
  five vertices and every cut (33,866 exact gluing identities), then rebuilt
  the Petersen and K4 supports without the primary recursion. A second,
  committed perfect-matching engine reproduced all 105 edge-side signatures
  in the five-graph cubic atlas. An independent GF(2) incidence-rank path
  reproduced the AP13 histogram, Cayley count, all census totals, hashes, and
  normal/optimized transcripts. The referees found no omitted torsion
  coordinate and confirmed the LRC-safety firewall.
---

# THM-4116 -- boundary-state gluing and AP odd-shell tree synchronizers

**PROVED + VERIFIED-EXACT.**  The useful abstraction of a snark is an empty
boundary-extension pairing, not an arithmetic prime analogy.  Transporting
that exact interface operation to THM-4110's AP torsion sheets exposes the
missing topology: six independent-looking checks do not suffice unless they
connect all seven odd phase bits.

## 1. Ordered boundary states glue by a dot product

Let `G=(V,E)` be a finite loopless graph, let `X` be a subset of `V`, and fix
an order

```text
delta(X)=(e_1,...,e_b)                                      (1)
```

of the edges having one endpoint in `X` and one outside it.  For a word
`sigma in {0,1,2}^b`, let `f_X(sigma)` be the number of assignments of three
labelled colors to

```text
E(G[X]) union delta(X)                                      (2)
```

that restrict to `sigma` on `(1)` and give pairwise distinct colors to all
edges incident with each vertex of `X`.  Define `f_(V-X)` with the same
ordered cut edges.  Then

```text
boxed:
#Col_3(G)=sum_(sigma in {0,1,2}^b) f_X(sigma)f_(V-X)(sigma). (3)
```

Here `#Col_3` counts labelled proper three-edge-colorings.

Indeed, restriction of a coloring gives one extension on each side with the
same boundary word.  Conversely, two side extensions with the same word
agree on every cut edge, and their union is a proper coloring of `G`.
Restriction and union are inverse.  This proves `(3)`.

The order of the boundary edges and the color gauge are load-bearing.  If a
nonempty cubic graph is colorable, all three colors occur at every vertex and
the global `S_3` action is free.  Dividing `(3)` by six therefore counts
global-color-gauge orbits, but `(3)` itself is before that quotient.

### Multipart contraction

The same restriction/union proof is not limited to two pieces.  Let

```text
V=V_1 disjoint_union ... disjoint_union V_k                 (3a)
```

be any vertex partition, let `E_partial` be the set of edges joining
different parts, and order `E_partial`.  For a coloring word
`tau in {0,1,2}^E_partial`, let `f_i(tau|delta(V_i))` count the assignments
to the internal and incident boundary edges that satisfy the edge-coloring
constraint at every vertex of `V_i`.  Then

```text
boxed:
#Col_3(G)=sum_(tau in {0,1,2}^E_partial)
              product_(i=1)^k f_i(tau|delta(V_i)).          (3b)
```

Thus an arbitrary vertex decomposition is an exact finite tensor network.
The tensor support retains extendability and multiplicity; replacing it by a
single colorable/noncolorable bit loses the compatibility obstruction.  The
edge order, repeated incidence of each interpart edge, and common color gauge
remain necessary sidecars.

## 2. Exact Petersen hostile and K4 positive control

Use the Petersen graph on vertices `0,...,9` with edges

```text
(0,1),(0,4),(0,5),(1,2),(1,6),
(2,3),(2,7),(3,4),(3,8),(4,9),
(5,7),(5,8),(6,8),(6,9),(7,9).             (4)
```

For `X={0,1}`, order the cut as

```text
((0,4),(0,5),(1,2),(1,6)).                                (5)
```

The two exact support sets are

```text
A={(a,b,c,d): a!=b, c!=d, {a,b}={c,d}},
B={(a,a,b,b): a,b in {0,1,2}}.                            (6)
```

Every word of `A` has multiplicity one, so `|A|=12` and total left mass is
`12`.  Every word of `B` has multiplicity two, so `|B|=9` and total right
mass is `18`.  The supports in `(6)` are disjoint.  Hence `(3)` gives

```text
#Col_3(Petersen)=0.                                        (7)
```

For the positive control `K4`, take the same side `X={0,1}` and ordered cut

```text
((0,2),(0,3),(1,2),(1,3)).                                (8)
```

Both supports have size and mass `12`, every multiplicity is one, and their
dot product is `6`.  This is one coloring up to the free global `S_3` gauge.

The exact companion computes each side vector by recursive edge assignment
and independently enumerates full graph colorings.  The two paths agree on
both `(7)` and the `K4` value `6`.  Thus this is a literal interface-state
certificate, not a tournament orientation and not a scalar analogy.

### Adjacent-edge parity and Kempe compression

There is a general one-integer snark sidecar hidden in the Petersen control.
Let `G` be a finite simple cubic graph, let `uv` be an edge, and order the
four retained semiedges leaving `{u,v}` as

```text
(uu_1,uu_2,vv_1,vv_2).                                    (8a)
```

The `{u,v}`-side support is always the set `A` in `(6)`, with multiplicity
one: the color of `uv` is the unique third color omitted by each endpoint
pair.

For a coloring of the complementary four-pole, let `k_c` be the number of
semiedges colored `c`, and let `m_c` be the number of internal edges of color
`c`.  Every complementary vertex sees color `c` exactly once, so

```text
|V(G)|-2 = 2m_c+k_c.                                       (8b)
```

A cubic graph has even order; hence every `k_c` is even.  Under global color
relabeling, the 21 parity-admissible four-letter words have four orbits, with
representatives

```text
S=0000,       P=0011,       X=0101,       Y=0110,           (8c)
```

of sizes `3,6,6,6`.  Let `(s_e,p_e,x_e,y_e)` be the per-word extension
multiplicities on these four orbits.  Then

```text
mass(r_e)=3s_e+6(p_e+x_e+y_e),
#Col_3(G)=6(x_e+y_e).                                      (8d)
```

The second identity holds because `A=X disjoint_union Y`.  In particular,
`G` is not three-edge-colorable exactly when `x_e=y_e=0`.

In that case a Kempe involution forces `s_e=p_e`.  Attach a degree-one leaf
to each retained semiedge and take an extension with boundary `0000`.  Its
color-`0/1` subgraph is a union of cycles and two boundary-to-boundary paths.
If a path joined a `u` leaf to a `v` leaf, switching colors `0` and `1` on
that path would produce an `X`- or `Y`-word and hence a coloring of `G`.
Thus the paths pair within the two endpoint blocks.  Switching the unique
`v`-to-`v` path sends `0000` to `0011`.  The same no-cross argument starting
from `0011` supplies the inverse switch.  Therefore

```text
(s_e,p_e,x_e,y_e)=(m_e,m_e,0,0).                           (8e)
```

The union of the three constant-word sectors has `3m_e` extensions.  The
global `S_3` action is free whenever it is nonempty, so `6` divides `3m_e`:
`m_e=2k_e` is even.  Consequently the complement support is empty when
`k_e=0` and is all nine words `{(a,a,b,b)}` otherwise.  The normalized
integer `k_e` is independent of orienting `e` or ordering the two semiedges
inside either endpoint block.

The exact atlas gives

```text
Petersen:    k_e=1 on all 15 edges;
flower J_5:  m_e histogram {4:5, 6:10, 10:10, 12:5};
Blanusa I:   m_e histogram {4:23, 6:4};
Blanusa II:  m_e histogram {4:25, 6:2}.                    (8f)
```

Under the companion's sorted blockwise semiedge convention, the `K4` control
has signature `(0,1,0,1)` on all six edges and six global colorings; an odd
swap inside either block exchanges the `X` and `Y` coordinates.  Thus
uncolorability is essential to `(8e)`.  The flower definition is the
repository's existing `J_5` construction; the two Blanusa definitions are
literal transcriptions of SageMath's official
[small-graph constructors](https://github.com/sagemath/sage/blob/develop/src/sage/graphs/generators/smallgraphs.py#L1291-L1417).
The companion checks every edge of all five graphs.  An independent
perfect-matching enumeration agrees with the recursive edge-coloring path on
all 105 edge sides.

## 3. AP odd-shell component law

Let `q>=2`, put

```text
r=2q-1,            v=(1,2,...,r),
O={1,3,...,2q-1},  T=R/Z.                                  (9)
```

Let `Gamma_neq` be the loopless graph on the speeds `1,...,r` with edge
`ij` exactly when

```text
nu_2(i)!=nu_2(j).                                          (10)
```

For a loopless simple graph `F` on the odd shell `O`, put

```text
Gamma_F=Gamma_neq union F.                                 (11)
```

Use THM-4110's circle-phase kernel and physical orbit

```text
K_Gamma(v)={theta in T^r:
  v_j theta_i-v_i theta_j=0 mod 1 for ij in E(Gamma)},
P_v={t v mod 1:t in T}.                                    (12)
```

Then there is a canonical group isomorphism

```text
K_(Gamma_F)(v)/P_v
  ~= {epsilon in F_2^O:
        epsilon_1=0,
        epsilon_o=epsilon_p for every op in E(F)}
  ~= F_2^(c(F)-1),                                         (13)
```

where `c(F)` is the number of connected components of `F`.  In particular,

```text
boxed: |K_(Gamma_F)(v)/P_v|=2^(c(F)-1).                   (14)
```

### Proof

Take `theta in K_(Gamma_neq)(v)` and set `t=theta_1`.  Every even speed `e`
is adjacent to `1`, so

```text
theta_e=e t.                                               (15)
```

Every odd `o>1` is adjacent to `2`.  Using `theta_2=2t`, its edge equation is

```text
2(theta_o-o t)=0 in T.                                     (16)
```

Therefore there is a unique bit `epsilon_o in F_2` such that

```text
theta_o=o t+epsilon_o/2;             epsilon_1=0.           (17)
```

Conversely, `(15)--(17)` satisfy every unequal-valuation edge.  On an
even--even edge both phases are physical.  On an odd--even edge the residual
is

```text
e epsilon_o/2=0 in T                                      (18)
```

because `e` is even.  Thus subtraction of `t v` identifies
`K_(Gamma_neq)(v)/P_v` with `F_2^(q-1)`, recovering and generalizing the
explicit AP13 half-toggle coordinates of THM-4110.

For an added odd edge `op`, the residual commutator is

```text
(p epsilon_o-o epsilon_p)/2.                               (19)
```

Both speeds are odd, so `(19)` vanishes in `T` exactly when
`epsilon_o=epsilon_p`.  The surviving bits are consequently constant on
each component of `F`.  The component containing speed `1` is fixed to zero;
each other component contributes one free bit.  This proves `(13)--(14)`.
**QED.**

The circle `R/Z` and the full labelled odd-shell incidence are essential.
The theorem concerns added odd--odd edges; same-shell even edges are already
redundant on this particular quotient.

## 4. Minimal synchronizers and the AP13 census

Equation `(14)` gives, for every `q>=2`,

```text
Gamma_F is saturated  iff  F is connected.                 (20)
```

A connected graph on `q` vertices has at least `q-1` edges.  Equality holds
exactly for a tree.  Hence:

```text
minimum added odd-shell constraints = q-1,
inclusion-minimal saturators = labelled trees on O,
number of minimal saturators = q^(q-2).                    (21)
```

The last equality is Cayley's formula.

For AP13, `q=7`.  The complete odd shell has `21` possible edges.  Exhausting
all `binom(21,6)=54,264` six-edge sets gives the exact surviving-sheet
histogram

```text
sheets       1       2       4      8
count    16,807  32,417   5,005     35.                    (22)
```

The first entry equals `7^5`, as `(21)` requires.  An anchor star and an odd
path each leave one sheet.  The six-edge hostile

```text
{1-3,1-5,1-7,1-9,1-11,3-5}                                (23)
```

wastes one edge in a cycle, isolates `13`, and leaves two sheets.  Thus the
number six is necessary but not sufficient: incidence rank, equivalently
connectivity here, is the missing coordinate.

More generally, `m<=5` edges on the seven odd vertices give
`c(F)>=7-m`, hence at least `2^(6-m)` sheets.  The exact minima for
`m=0,...,5` are

```text
64,32,16,8,4,2.                                            (24)
```

### Edge-fault tolerance

The component law also gives a robustness theorem for free.  If
`D subseteq E(F)` is any set of added odd-shell constraints, then

```text
|K_(Gamma_(F-D))(v)/P_v|=2^(c(F-D)-1).                    (24a)
```

Consequently, for `k>=1`, saturation survives deletion of every set of fewer
than `k` added constraints exactly when the edge connectivity of `F` is at
least `k` (with a disconnected graph assigned edge connectivity zero).  In
particular, a tree is a minimum synchronizer but is maximally fragile:
deleting any edge produces exactly two phase sheets.

For `q>=3`, a one-edge-fault-tolerant synchronizer has minimum degree at
least two, so the handshake lemma forces at least `q` added odd-shell edges.
Equality holds exactly when `F` is a labelled `q`-cycle: every degree is then
two and connectedness gives one cycle.  Hence the number of minimum
one-edge-fault-tolerant synchronizers is

```text
(q-1)!/2.                                                  (24b)
```

For AP13 this gives seven constraints and `6!/2=360` labelled cycles.  When
`q=2`, the graph is simple and no one-edge-fault-tolerant synchronizer
exists.  This is a literal transfer of the graph-theoretic “bridge versus
robust gluing” distinction; it is not a transfer of snark noncolorability.

The companion exhausts every odd-shell graph for `q=2,...,6`, checking
`(14)` by a separate bit enumeration (`33,866` graphs and `1,065,508` bit
gates).  It also checks all `27,896` AP13 edge sets of size at most five,
all `54,264` six-edge sets, and all `3,392` unequal-valuation commutator
gates on the `64` explicit AP13 packets.  Normal and optimized runs reproduce
the stored output byte for byte.

## 5. What transfers, and what does not

The exact connection contract is:

```text
source       ordered multipole boundary-color words,
target       AP odd-shell half-toggle bits,
map          matching boundary words / imposing bit equalities,
preserved    exact interface compatibility and saturation,
destroyed    continuous time, odd-primary Smith data, clearance,
             endpoint owners, and internal extension geometry,
sidecar      boundary order and color gauge on the source;
             labelled weighted incidence and physical interval on the target,
test         Petersen/K4 dot products and the complete AP13 rank census
             through six added edges.                                  (25)
```

The native carrier is a relation or reversible action groupoid.  There is no
intrinsic complete antisymmetric relation, so forcing these states into a
tournament would discard the simultaneous compatibility operation.

Most importantly, `(20)` proves synchronization of **circle phases**, not
loneliness.  A physical point can be unsafe (`t=0`) or safe.  Any LRC
consumer must still retain an actual interval, clearance, owner labels,
component ancestry, and the relevant clock.  Therefore `(14)--(24)` sharpen
the physical-entry sidecar of THM-4110 but do not prove LRC(14).
