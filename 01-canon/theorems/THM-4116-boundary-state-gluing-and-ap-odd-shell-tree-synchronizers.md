---
id: THM-4116
title: "Boundary-state gluing and AP odd-shell tree synchronizers"
status: >
  PROVED + VERIFIED-EXACT. Ordered boundary extension vectors glue by an
  exact dot product. On the Petersen two-vertex cut their supports are
  disjoint, while the analogous K4 cut has dot product six. For AP_(2q-1),
  after the complete unequal-v2 reciprocal graph, an added odd-shell graph F
  leaves exactly 2^(c(F)-1) circle-phase sheets. Saturation is equivalent to
  connectivity; its inclusion-minimal repairs are the q^(q-2) labelled
  spanning trees. This synchronizes phase sheets but does not prove LRC
  safety or LRC(14).
source: codex-snark-apex-boundary-session-20260825
depends_on:
  - THM-4110-sparse-reciprocal-phase-graph-saturation-and-ap13-torsion-tariff
related:
  - THM-261-petersen-root-orthogonality
  - THM-343-H7-impossible
  - THM-4105-primitive-reciprocal-phase-descent-and-quantitative-arrival
  - MISTAKE-501
script: 04-computation/snark_boundary_ap_shell_thm4116.py
output: 05-knowledge/results/snark_boundary_ap_shell_thm4116.out
script_sha256: 41f407ddfe435af428f5b4aceeec2267bd1431b8b67555d68c1b141dd5ea050a
output_sha256: 4ebfdcd47e8e8c9d5bf62a0f5aca341ad4d02dffdde91ff4b8f8f9c4b20d0f80
hash_basis: raw LF bytes
---

# THM-4116 -- boundary-state gluing and AP odd-shell tree synchronizers

**PROVED + VERIFIED-EXACT.**  The useful abstraction of a snark is an empty
boundary-extension pairing, not an arithmetic prime analogy.  Transporting
that exact interface operation to THM-4110's AP torsion sheets exposes the
missing topology: six independent-looking checks do not suffice unless they
connect all seven odd phase bits.

## 1. Ordered boundary states glue by a dot product

Let `G=(V,E)` be a finite graph, let `X` be a subset of `V`, and fix an order

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
test         Petersen/K4 dot products and the complete AP13 rank census. (25)
```

The native carrier is a relation or reversible action groupoid.  There is no
intrinsic complete antisymmetric relation, so forcing these states into a
tournament would discard the simultaneous compatibility operation.

Most importantly, `(20)` proves synchronization of **circle phases**, not
loneliness.  A physical point can be unsafe (`t=0`) or safe.  Any LRC
consumer must still retain an actual interval, clearance, owner labels,
component ancestry, and the relevant clock.  Therefore `(14)--(24)` sharpen
the physical-entry sidecar of THM-4110 but do not prove LRC(14).
