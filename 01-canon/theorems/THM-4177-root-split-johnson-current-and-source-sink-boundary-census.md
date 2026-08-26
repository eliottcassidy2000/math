---
id: THM-4177
title: "Root-split Johnson current and source-sink boundary census"
status: >
  PROVED ELEMENTARY ARBITRARY-TENSOR ROOT SPLIT + PROVED CAPACITY-SPECIFIC
  EXPOSED-CURRENT AND ODD-PATH DECOMPOSITIONS + PROVED STRICT UNIVERSAL-
  SOURCE/SINK PADDING DESCENT + FINITE-EXACT COMPLETE q=3..8 PARENT-CLASS
  BY ATTACHMENT-PRESENTATION CENSUS + INDEPENDENTLY AUDITED THROUGH q=6.
  The order-at-least-four no-sink/no-source strict sign law is an OPEN
  CONJECTURE. This theorem does not close the asymmetric order-eleven
  rational bank, exact
  Johnson cosets, or actual response maximizers. Counts are presentations,
  not child orbits or unrooted isomorphism classes.
source: codex-frontier-synthesis-creative-20260826ay
depends_on:
  - THM-002-ocf
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4144-order-eleven-large-homogeneous-module-johnson-centrality
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
  - THM-4168-prime-order-eleven-nontrivial-automorphism-johnson-centrality
  - THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer
  - THM-4172-multideletion-support-tomography-and-same-parity-johnson-holonomy
script: 04-computation/tournament_root_split_source_sink_census_thm4177.cpp
output: 05-knowledge/results/tournament_root_split_source_sink_census_thm4177.out
independent_audit_script: 04-computation/tournament_root_split_source_sink_census_thm4177_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_root_split_source_sink_census_thm4177_independent_audit.out
script_sha256: dcd53fac7254b692fc3ab3c75fe415d8dfeab4cd1dcce459a72ae48e881633dd
output_sha256: 8bb53672d638b576f18103f665f6a29432967a4ac078609a3683c7752126535d
independent_audit_script_sha256: f40ad9ca8648d40aa052a6988d3480eadd11ca031a1e16d6331860a8723c8cb0
independent_audit_output_sha256: 4dad8495c4fdf4a4258bfbfee9c779ba3e0d29956e3e9320573addf073014e45
gentourng_sha256: 89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110
generator_stream_q3_q8_sha256: 0489647b21ed6a8e5ab135f41636685abfa91376224e5d40a9715afe8ed78f86
generator_stream_q3_q6_sha256: 094f9c2f15694ce607004f30549933effeefbe25aa47d579ef2888964c720d8c
hash_basis: Raw file bytes; tracked text artifacts use LF line endings.
primary_audit: >
  PASS. A standalone warning-clean C++17 endpoint/convolution engine checks
  14,522,544 root-capacity identities, 1,823,696 root-split identities, and
  203,206 odd-path capacity coordinates on the complete q=3..8 parent-class
  by attachment-presentation universe. The full q<=8 run passes at Clang O3.
  Complete q<=6 Clang O0, O3, and UBSan streams byte-match.
independent_audit: >
  ACCEPT on its stated q<=6 scope. A standard-library Python referee uses
  Python integers, a separate endpoint implementation, direct disjoint-edge
  enumeration, and direct odd-path recursion. Normal, optimized, and fixed-
  hash-seed complete q=3..6 streams byte-match. It also checks arbitrary-
  tensor root splits and an explicit strong generic-tensor sign hostile.
---

# THM-4177 -- root-split Johnson current and source-sink boundary census

**PROVED ALGEBRA + PROVED CAPACITY IDENTITIES + FINITE-EXACT q<=8 +
INDEPENDENTLY AUDITED q<=6.**

The order-eleven gate of THM-4169 is `G_+>0` and `G_->0`.  This theorem
separates what is formal for every symmetric edge tensor from what uses the
tournament exposure capacity.  It then records a complete small-order
boundary census.  The census strongly suggests an intrinsic source/sink
law from order four onward, but that law remains open in general.  The
order-three directed cycle is a necessary boundary exception: it has no
source or sink and has `D=C=G_+=G_-=0`.

Two type distinctions are load-bearing throughout.

1. If `T=Q+x` is an actual tournament and `c=c(T)`, then
   `c|_Q` is the **corrected restriction of the child tensor**.  It need not
   equal the actual card capacity `c(Q)`.
2. A row `(Q,z)` below is a parent-marked attachment **presentation**.
   Parent automorphisms and forgetting the root can identify rows.  No row
   total is a rooted-orbit or unrooted-child-class count.

## 1. The arbitrary symmetric-tensor root split

Let `T=Q+x`, with old vertex set `V`, and let `z` be any symmetric edge
tensor on the oriented complete graph of `T`.  Write

```text
b_ij=z_ij                         (i,j in V),
a_i=z_ix,
sigma_i=+1                        if i->x,
sigma_i=-1                        if x->i.                 (1)
```

For the old tensor put

```text
d_i=sum_(j!=i)b_ij,
h_i=sum_(i->j in Q)b_ij-sum_(j->i in Q)b_ij,
W=sum_(i<j)b_ij,
A=sum_i a_i,                      S=sum_i sigma_i a_i.     (2)
```

Retain THM-4167's quadratic packet

```text
C(z)=sum_i d_i(z)h_i(z),
D(z)=sum_(e<f, e disjoint f)z_ez_f,
G_eta(z)=D(z)+2 eta C(z),          eta in {+1,-1}.         (3)
```

Thus `G_(+1)=G_+` and `G_(-1)=G_-`.

> **Theorem 1 (root split; arbitrary symmetric tensor).**
>
> ```text
> D(z)=D(b)+sum_i a_i(W-d_i),                            (4)
>
> C(z)=C(b)+sum_i a_i h_i+sum_i sigma_i a_i d_i
>             +sum_i sigma_i a_i^2-A S.                 (5)
> ```
>
> Consequently
>
> ```text
> G_eta(z)=G_eta(b)+sum_i a_i(W-d_i)
>   +2 eta [sum_i a_i h_i+sum_i sigma_i a_i d_i
>                    +sum_i sigma_i a_i^2-A S].          (6)
> ```

### Proof

A disjoint edge pair is either wholly old, or consists of one root edge
`ix` and one old edge avoiding `i`.  This is exactly `(4)`.

At an old vertex the new unsigned and signed degrees are

```text
d_i(z)=d_i+a_i,                   h_i(z)=h_i+sigma_i a_i.
```

At the root they are `d_x=A` and `h_x=-S`.  Expanding
`sum_i d_i(z)h_i(z)+d_xh_x` proves `(5)`, and `(6)` follows.  No
positivity, tournament capacity, primeness, or card identity is used. QED.

There is a useful orientation-resolved form.  Put

```text
P={i:i->x},                       M={i:x->i},
o_i=(d_i+h_i)/2,                  r_i=(d_i-h_i)/2.         (7)
```

Then `(5)` is equivalently

```text
C(z)=C(b)+2[sum_(i in P)a_i o_i-sum_(i in M)a_i r_i
             -sum_(i<j in P)a_i a_j+sum_(i<j in M)a_i a_j].   (8)
```

This fixes the sign convention: old vertices dominating the root lie in
`P`; vertices dominated by the root lie in `M`.

## 2. The capacity-specific root current

Now let `Q` be a tournament and let `T_z=Q+x_z`, with the THM-4169
convention

```text
x->j iff z_j=1.
```

Let `Start_i`, `End_i`, and `Q_ij` be THM-4097's endpoint and ordered
exposed-gap counts of the **actual parent** `Q`.  Thus `Q_ij` counts old
vertex orders with `i` immediately before `j` and every other adjacency
valid.  Put

```text
M={j:z_j=1},                      P=V(Q)-M,
a_i(z)=c_ix(T_z).                                           (9)
```

> **Theorem 2 (affine skew root current).** For every `i`,
>
> ```text
> a_i(z)=Start_i+End_i
>       +sum_(j in M-{i})Q_ij+sum_(j in P-{i})Q_ji.       (10)
> ```
>
> In particular `a_i` is affine and independent of its mutual bit `z_i`:
>
> ```text
> a_i(z)=a_i(0)+sum_(j!=i) beta_ij z_j,
> beta_ij=Q_ij-Q_ji=-beta_ji.                            (11)
> ```

### Proof

Take an exposed child word whose marked gap is `ix`.  If `x` is an
endpoint, deleting it leaves respectively a parent path starting or ending
at `i`, giving `Start_i+End_i`.  Otherwise `x` has one unmarked neighbor
`j`.  A segment

```text
i | x -> j
```

is legal exactly when `j in M`; deleting `x` gives an exposed `i|j` word,
counted by `Q_ij`.  A segment

```text
j -> x | i
```

is legal exactly when `j in P`; deletion gives an exposed `j|i` word,
counted by `Q_ji`.  These cases are disjoint and reversible, proving `(10)`.
Subtracting the `M=empty` value proves `(11)`, whose skewness is literal.
QED.

The uniform corners sharpen to

```text
a_i(0)=o_i(c(Q))+2End_i,          (x is a sink),
a_i(1)=r_i(c(Q))+2Start_i,        (x is a source).        (12)
```

Here `o_i,r_i` use the integer capacity `c=Q+Q^T`.  For example,
`sum_(i->j)Q_ij=H(Q)-End_i`; substitution in `(10)` gives `(12)`.

Skewness supplies exact cut-current corollaries.  For `M={j:z_j=1}`,

```text
sum_i a_i(z)=sum_i a_i(0)+sum_(i in P,j in M)beta_ij,    (13)
sum_(i in M)a_i(z)=sum_(i in M)a_i(0),                  (14)
sum_i a_i(1)=sum_i a_i(0).                              (15)
```

Thus the total root mass changes only by the skew current across the
attachment cut; internal currents cancel.  The antisymmetric `Q` sidecar is
essential.  The symmetric capacity `Q_ij+Q_ji` alone forgets `beta`.

## 3. Odd directed paths are the capacity atoms

For distinct `i,j`, let `P^odd_ij(T)` be the directed vertex-simple paths
whose endpoint set is `{i,j}` and whose number of arcs is odd.  Both directed
endpoint orders are included.  Use `H(empty)=1`.

> **Theorem 3 (odd-path capacity decomposition).**
>
> ```text
> boxed: c_ij(T)=2 sum_(R in P^odd_ij(T)) H(T-V(R)).      (16)
> ```

### Proof

Use the tagged OCF formula of THM-4167.  In one tag, the distinguished odd
cycle through the auxiliary vertex has local segment

```text
j -> x -> i
```

and then a directed path from `i` back to `j` in `T`.  Deleting `x` turns
that cycle bijectively into a directed simple path between `i,j`.  Its
length is `|cycle|-2`, hence odd.  The other tag gives the opposite endpoint
order.

All remaining tagged cycles are an arbitrary vertex-disjoint odd-cycle
collection on `T-V(R)`.  If that collection has `k` cycles, the original
tag has weight `2^(k+1)`.  Summing gives

```text
2 sum_k alpha_k(Omega(T-V(R)))2^k=2H(T-V(R))
```

by THM-002.  The deletion and insertion maps preserve the ambient tag, so
the two endpoint directions are not accidentally identified. QED.

For odd `ell>=1`, define the length layer

```text
c^[ell]_ij=2 sum_(R in P^odd_ij(T), length(R)=ell)
                  H(T-V(R)).                             (17)
```

Then `c=sum_(ell odd)c^[ell]`.  The length-one boundary is especially
simple: the tournament supplies exactly one directed arc on `{i,j}`, so

```text
c^[1]_ij=2H(T-{i,j}).                                   (18)
```

## 4. The exact polarization kernel

For two symmetric tensors `u,v` on the same oriented complete graph, put

```text
B_eta(u,v)=G_eta(u+v)-G_eta(u)-G_eta(v).                 (19)
```

THM-4167 gives

```text
C=2(sum_(common tail)e<f z_ez_f
     -sum_(common head)e<f z_ez_f).                      (20)
```

Hence `B_eta` has the following exact kernel.  Each displayed coefficient
multiplies `u_ev_f+v_eu_f`.

| relation of `e,f` | coefficient in `B_eta` |
|:---|---:|
| same edge | `0` |
| disjoint | `1` |
| common tail | `4 eta` |
| common head | `-4 eta` |
| shared vertex, one entering and one leaving | `0` |

Therefore

```text
G_eta(c)=sum_(ell odd)G_eta(c^[ell])
          +sum_(ell<m, both odd)B_eta(c^[ell],c^[m]).    (21)
```

For `G_+`, the only negative kernel objects are pairs of endpoint arcs with
a common head.  A no-sink hypothesis supplies an outgoing continuation from
that head, but a multiplicity-preserving path-pair injection is not proved.
Section 7 shows why a proof cannot demand every self-layer gate to be
nonnegative.

## 5. A proved strict padding descent

The odd-path atoms do prove one all-order sign reduction.  Let `x triangleright
Q` denote the universal-source extension `x->V(Q)`.  Its root capacities are

```text
a_i=r_i(c(Q))+2Start_i                                  (22)
```

by `(12)`.

> **Lemma 4 (incoming fan injection).** For every `v in V(Q)`,
>
> ```text
> r_v(c(Q)) <= sum_(u->v)a_u.                            (23)
> ```

### Proof

Divide all capacities by two and use `(16)` with the underlying tagged OCF
atoms retained.  Fix an incoming arc `u->v` and an atom of `c_uv/2`, given
by an odd directed path `R` between `u,v` and an OCF collection `Gamma`
off `R`.

- If `R` runs from `u` to `v`, delete its terminal `v` and prepend `x`.
  If `p->v` is the last arc of `R`, the result is an odd root path from
  `x` to `p`, hence an atom of the copy `a_p/2` indexed by `p->v`.
- If `R` runs from `v` to `u`, delete its initial `v` and prepend `x`.
  The result is an odd root path from `x` to `u`, hence an atom of the copy
  `a_u/2` indexed by the original arc `u->v`.

In both cases `v` becomes an unused fixed vertex and `Gamma` is unchanged,
so the OCF weight is preserved.  The images are disjoint: the first old
vertex after `x` dominates `v` in the first case and is dominated by `v` in
the second.  Appending or prepending `v` recovers the original path.  This
is a weight-preserving injection into the disjoint union on the right of
`(23)`. QED.

More precisely, the image consists exactly of target root-path atoms in
which `v` is fixed: `v` lies in neither the root path nor `Gamma`.  Thus
equality in `(23)` holds exactly when every target atom fixes `v`.

Universal-source padding leaves every old--old capacity unchanged: in an
old-edge exposed word the source `x` must be first, and deleting it is a
bijection to the parent exposed word.  Apply `(8)` with `M=V(Q)` and use
`sum_(i<j)a_i a_j=sum_v a_v sum_(u->v)a_u`.  This gives the exact comparison

```text
G_+(x triangleright Q)-G_+(Q)
 =sum_i a_i(W-d_i)
   +4sum_v a_v(sum_(u->v)a_u-r_v).                       (24)
```

Every term in the second sum is nonnegative by `(23)`.  If `|Q|>=3`, then
`a_i>=2H(Q-i)>0`, while `W-d_i` is the positive capacity mass on the old
edges avoiding `i`.  Hence

```text
boxed: G_+(x triangleright Q)>G_+(Q),       |Q|>=3.      (25)
```

There is no equality case in `(25)` for `|Q|>=3`, even if every fan
inequality `(23)` is an equality.  Taking converses gives the exact dual:
universal-sink padding strictly raises `G_-`.

Consequently, a minimum-order counterexample to the proposed implication

```text
|T|>=4 and no sink => G_+>0                                  (26)
```

cannot have a source.  Indeed, delete that universal source; no old vertex
loses an outgoing arc, so the card `Q` still has no sink.  If `|Q|=3`, then
`Q` is the directed cycle and `G_+(Q)=0`, so the strict inequality `(25)`
already gives `G_+(T)>0`.  If `|Q|>=4`, minimality gives `G_+(Q)>0`, and
again `(25)` is a contradiction.  Dually, a minimum-order counterexample to
`|T|>=4 and no source => G_->0` cannot have a sink.  Thus the unresolved
minimal core is genuinely source-and-sink-free.  This is a proved descent,
not a proof of `(26)`.

The capacity hypothesis is indispensable.  On the strong source/sink-free
order-four orientation with label `100111`, assign arbitrary positive edge
weights in lexicographic edge order

```text
(z_01,z_02,z_03,z_12,z_13,z_23)=(1,1,1,1,2,3).
```

Direct calculation gives `(C,D,G_+)=(-4,6,-2)`.  Therefore neither the
root split nor orientation alone implies the desired sign.

## 6. Complete q=3..8 attachment-presentation census

For each `3<=q<=8`, nauty's `gentourng -q q` supplies one representative of
every parent isomorphism class `Q`.  The primary certificate evaluates all
`2^q` patterns `z` relative to that representative and constructs the
**actual child tensor** `c(T_z)` by a fresh endpoint/exposed-gap DP.

The exact construction universe is:

| `q` | parent classes | strong parents | prime parents | `(Q,z)` presentations | strong-child presentations | prime-child presentations |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 1 | 16 | 8 | 0 |
| 4 | 4 | 1 | 0 | 64 | 32 | 19 |
| 5 | 12 | 6 | 3 | 384 | 254 | 126 |
| 6 | 56 | 35 | 15 | 3,584 | 2,771 | 1,599 |
| 7 | 456 | 353 | 197 | 58,368 | 50,668 | 34,160 |
| 8 | 6,880 | 6,008 | 4,008 | 1,761,280 | 1,634,407 | 1,256,103 |

The child columns count presentations, not distinct child classes.  The
same warning applies to all totals below.

### 6.1 Intrinsic child boundary

For `G_+`, the sink strata are:

| `q` | sink parents | child presentations with a sink | `G_+<0` | `G_+=0` | `G_+>0` among sink children | no-sink children with `G_+<=0` |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 6 | 2 | 0 | 4 | 0 |
| 4 | 2 | 20 | 12 | 0 | 8 | 0 |
| 5 | 4 | 76 | 52 | 0 | 24 | 0 |
| 6 | 12 | 440 | 352 | 0 | 88 | 0 |
| 7 | 56 | 4,040 | 3,530 | 0 | 510 | 0 |
| 8 | 456 | 65,248 | 60,216 | 0 | 5,032 | 0 |

The converse census is identical.  Thus, throughout this finite universe,

```text
G_+(T)<=0 => T has a sink,
G_-(T)<=0 => T has a source,                              (27)
```

with no equality rows.  The converses are false: many source/sink children
have positive gates.  For example the transitive-parent uniform sink at
order four has `G_+=16`.

Every strong-child and every prime-child presentation has both gates
strictly positive.  There are zero strong or prime failures in every row of
the first table.

### 6.2 Parent-marked corner equivalences

For every parent in the complete `q<=8` universe, the certificate also
checks the two exact presentation-level equivalences

```text
min_(z!=0)G_+(T_z)<=0                 iff Q has a sink,
min_(z!=0)G_+(T_z)<G_+(T_0)          iff Q has a sink,    (28)
```

and their duals

```text
min_(z!=1)G_-(T_z)<=0                 iff Q has a source,
min_(z!=1)G_-(T_z)<G_-(T_1)          iff Q has a source.  (29)
```

All four disagreement counts are zero at every `q`.  These are finite
parent-existence statements, not an all-order intrinsic classification.

All nonuniform attachments over a strong or prime parent pass both gates.
The exact minimum of `min(G_+,G_-)` on those presentation banks is:

| `q` | strong-parent nonuniform presentations | minimum gate | prime-parent nonuniform presentations | minimum gate | prime-parent/prime-child presentations |
|---:|---:|---:|---:|---:|---:|
| 3 | 6 | 44 | 6 | 44 | 0 |
| 4 | 14 | 640 | 0 | -- | 0 |
| 5 | 180 | 6,476 | 90 | 6,476 | 60 |
| 6 | 2,170 | 54,476 | 930 | 59,508 | 750 |
| 7 | 44,478 | 429,044 | 24,822 | 461,988 | 22,064 |
| 8 | 1,526,032 | 3,320,952 | 1,018,032 | 3,552,196 | 953,904 |

For prime parents the last column also verifies THM-4169's elementary count
`2^q-2-2q` per parent.  It remains a presentation count.

### 6.3 Algebra checks

Across the full primary universe there are

```text
root-capacity formula checks       14,522,544, failures 0,
root-split packet checks            1,823,696, failures 0,
odd-path capacity coordinates         203,206, failures 0.       (30)
```

For an actual child, the root-split check deliberately uses
`b=c(T_z)|_Q`, the corrected restriction.  It never substitutes the actual
card tensor `c(Q)`.  The root-current check separately reconstructs the
actual child root capacities from the actual parent endpoint table via
`(10)`.

## 7. Length-layer census and the first sign obstruction

The primary engine decomposes every actual parent capacity into `(17)` and
checks `(21)`.  On every no-sink parent, it tests all `G_+` self layers and
all distinct-length polarizations; the dual no-source counts are identical.
Each tuple below is `(checks, negative, zero, minimum)`.

| `q` | no-sink parents | self layers | cross layers |
|---:|---:|:---|:---|
| 3 | 1 | `(1,0,1,0)` | `(0,0,0,--)` |
| 4 | 2 | `(4,0,0,12)` | `(2,0,0,20)` |
| 5 | 8 | `(16,0,0,132)` | `(8,0,0,308)` |
| 6 | 44 | `(132,0,0,48)` | `(132,0,0,312)` |
| 7 | 400 | `(1,200,0,0,648)` | `(1,200,0,0,2,376)` |
| 8 | 6,424 | `(25,696,15,3,-3,200)` | `(38,544,0,2,0)` |

The apparent layerwise proof fails first at `q=8`.  All fifteen negative
self gates occur only in the spanning length-seven layer.  Every carrier is
strong, source-free, sink-free, and decomposable; none is prime.  Thirteen
have a module of minimum size two, one minimum size three, and one minimum
size five.  The sharp row is

```text
label=1111011110011111111111111001,
G_+(c^[7])=-3200,
B_+(c^[1],c^[7])=47044,
G_+(c)=1448440.                                         (31)
```

In every one of the fifteen rows, the length-one/length-seven polarization
alone exceeds the negative length-seven self gate.  The frozen output lists
all fifteen labels and the three incident cross terms.  Thus any all-order
proof must permit cross-length compensation.  The observed nonnegativity of
all `q<=8` distinct-length polarizations is **FINITE-EXACT only** and is not
asserted as an all-order conjecture here.

## 8. Open sign law and consequence firewall

The finite data and the padding descent motivate the intrinsic conjecture

```text
OPEN:  |T|>=4 and no sink in T   => G_+(c(T))>0,
OPEN:  |T|>=4 and no source in T => G_-(c(T))>0.          (32)
```

The padding descent reduces a minimum counterexample to the
source-and-sink-free core, with the directed three-cycle handled explicitly
as the zero boundary case.  But neither `(16)`, the kernel table, nor the
finite census supplies the
missing multiplicity-preserving path-pair injection.  Equation `(32)` is not
proved, even for all prime attachments of an arbitrary order-ten parent.

If `(32)` were proved, then every prime order-eleven tournament would pass
THM-4169's strict rational gate `2|C|<D`, since prime tournaments are strong
and have neither a source nor a sink.  That conditional implication does not
promote the present theorem into an order-eleven result.

Even a future proof of `(32)` would establish only rational Johnson
support-floor centrality.  THM-4162's complete `(H,c)` sidecar can compute
each exact Johnson coset, but the layer-dependent map

```text
J_m -> ceil_(a_m+d_m Z)(J_m)
```

is not controlled by `G_+` and `G_-`; distinct anchors and lattices can erase
or reverse a rational gap.  THM-4172 transports homogeneous `C,D` data of
corrected restrictions, not those exact layer anchors or actual card
capacities.  Actual response maximizers are a further, separate object.

Finally, THM-4133 is not a hostile to `(32)`: its order-twelve failure uses
the stricter order-twelve coefficient, while it still has both
`D+2C>0` and `D-2C>0`.

**Later boundary.** THM-4181 computes exact ordinal-sum transfer through two
rooted path-cover parity states and proves its positive remainder on `242,060`
factor-order-at-most-seven presentations. Its all-order remainder remains
**OPEN**, so its strong-core reduction does not change the status of `(32)`.

## 9. Replay

Build and run the full primary certificate with

```bash
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/tournament_root_split_source_sink_census_thm4177.cpp \
  -o /tmp/thm4177-primary

for q in 3 4 5 6 7 8; do gentourng -q "$q"; done \
  | /tmp/thm4177-primary
```

The complete O0 and UBSan control scope is `q<=6`; it byte-matches the O3
stream on that same scope.  The full `q<=8` claim is certified by the O3
run.  This scope distinction is intentional.

Run the independent audit with

```bash
for q in 3 4 5 6; do gentourng -q "$q"; done \
  | python3 -B \
      04-computation/tournament_root_split_source_sink_census_thm4177_independent_audit.py

for q in 3 4 5 6; do gentourng -q "$q"; done \
  | python3 -B -O \
      04-computation/tournament_root_split_source_sink_census_thm4177_independent_audit.py

for q in 3 4 5 6; do gentourng -q "$q"; done \
  | PYTHONHASHSEED=271828 python3 -B \
      04-computation/tournament_root_split_source_sink_census_thm4177_independent_audit.py
```

All three Python streams byte-match the frozen independent output.  The
independent scope is the complete `q<=6` universe, not a second audit of the
`q=7,8` census. **QED for the proved identities, padding descent, and stated
finite universes; `(32)` remains OPEN.**
