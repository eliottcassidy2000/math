---
id: THM-2793
title: "Oriented ramp reference and rooted graceful gap reconstruction"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  One oriented ramp reference upgrades interval Gram tomography
  from an unordered column multiset to the complete ordered
  consecutive-ones matrix: row length and first moment recover both
  endpoints.  For a tree graceful-gap matrix, a vertex fixed at label zero
  removes the remaining cut-complement gauge, and mod-two incidence
  inversion recovers every threshold set, the full graceful labeling, and
  all marked gap-tail decisions.  This is an exact recognition and sidecar
  theorem, not a graceful-label existence theorem.
source: root/oriented-ramp-gap-reconstruction-2026-07-28
depends_on:
  - THM-2787-signed-path-sum-weyl-orbit-and-gap-tail-leaf-insertion
  - THM-2789-interval-gram-tomography-and-graceful-gap-tail-quadratic-detector
related:
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
script: 04-computation/oriented_ramp_rooted_gap_reconstruction_thm2793.py
output: 05-knowledge/results/oriented_ramp_rooted_gap_reconstruction_thm2793.out
script_sha256: bc2c7d650ef70ce95ab50c16f7268d5c4a1184a64662c26f064c491fff0a37bb
output_sha256: 4f29bdb9e1826e5cf749b94380656bba56a3b6eba51104fa817747fdbbae7ca6
hash_basis: LF-normalized bytes
---

# THM-2793 -- oriented ramp reference and rooted graceful gap reconstruction

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2789 proves that the Gram matrix of intervals reconstructs their entire
column multiset but not its order.  The missing coordinate is only
one-dimensional: pair the rows with the oriented position ramp.  Length and
first moment determine an interval's endpoints.

For a tree there is one further binary ambiguity.  Every edge-cut column
determines its vertex side only up to complement.  Fixing the vertex with
label zero removes that gauge simultaneously at every threshold.  The
result is an exact reconstruction and recognition theorem for graceful gap
data.

## 1. One ramp reference restores interval order

Let

```text
X={1,...,N},
C_i=1_(I_i) in {0,1}^X                    (1)
```

for nonempty intervals `I_1,...,I_m`.  Put

```text
G_ij=<C_i,C_j>,
r=(1,2,...,N),
mu_i=<C_i,r>.                              (2)
```

The data `mu_i` are the cross-Gram block against one labelled reference
vector `r`.  Write

```text
ell_i=G_ii=|I_i|.                          (3)
```

If `I_i={L_i,L_i+1,...,R_i}`, then

```text
mu_i=ell_i(L_i+R_i)/2,
R_i-L_i=ell_i-1.                           (4)
```

Therefore

```text
L_i=(2mu_i/ell_i-ell_i+1)/2,
R_i=(2mu_i/ell_i+ell_i-1)/2.               (5)
```

Thus the diagonal of `G` and the single ramp cross-block recover every
ordered row of `C` literally.  Once `(5)` is known, the off-diagonal Gram
entries are redundant for reconstruction; they remain an exact
compatibility check

```text
G_ij=|[L_i,R_i] intersect [L_j,R_j]|.       (6)
```

This upgrades THM-2789's unordered multiset tomography to ordered-matrix
tomography.  It is one labelled reference vector, not one scalar: the
cross-block has one entry per interval row.

## 2. On graceful gaps, the ramp is the edge-sum sidecar

Let `T` be an `m`-edge tree with graceful labeling

```text
f:V(T)->{0,...,m}.                         (7)
```

For the unique edge `e_d={u_d,v_d}` of difference `d`, orient only its
numeric endpoint values:

```text
a_d=min(f(u_d),f(v_d)),
b_d=max(f(u_d),f(v_d))=a_d+d.              (8)
```

Its internal gap row is

```text
I_d={a_d+1,...,b_d} subset {1,...,m}.       (9)
```

Let `mu_d` be the ramp moment of this row.  Then

```text
2mu_d=d(f(u_d)+f(v_d)+1),                  (10)
```

and hence

```text
f(u_d)+f(v_d)=2mu_d/d-1,                   (11)

a_d=(2mu_d/d-d-1)/2,
b_d=(2mu_d/d+d-1)/2.                       (12)
```

So the ramp is exactly the missing edge-sum coordinate beside the prescribed
edge difference.  This links THM-2789's interval Gram to THM-2761's
edge-sum geometry: sum plus absolute difference recovers both endpoint
labels, while the interval row records every threshold between them.

## 3. A root-zero anchor inverts every tree cut

Now retain the tree `T`, the edge-difference bijection

```text
pi:E(T)->{1,...,m},                        (13)
```

the decoded ordered gap matrix `C`, and a vertex `r_0` declared to carry
label zero.

For each threshold `t=1,...,m`, the `t`th column is an edge vector

```text
c_t in F_2^E.                              (14)
```

Reduced incidence over `F_2` is an isomorphism for a tree:

```text
partial_T:F_2^V/<1> -> F_2^E.              (15)
```

Therefore `(14)` has exactly two complementary vertex-side lifts.  Let
`S_t` be the unique lift containing `r_0`:

```text
partial_T 1_(S_t)=c_t,
r_0 in S_t.                                (16)
```

For genuine graceful data,

```text
S_t={v:f(v)<t}.                             (17)
```

In particular

```text
S_1 subset S_2 subset ... subset S_m,
|S_t|=t.                                   (18)
```

Conversely, suppose interval decoding `(5)` is integral and within the
ambient line, `G` passes `(6)`, each row indexed by `d` has length `d`, and
the anchored lifts `(16)` satisfy `(18)`.  Define

```text
f(v)=m-#{t:v in S_t}.                      (19)
```

The flag `(18)` makes `f` a bijection onto `{0,...,m}` with `f(r_0)=0`.
An edge `e_d` crosses exactly the thresholds in row `d`, so

```text
|f(u_d)-f(v_d)|=|I_d|=d.                   (20)
```

Hence `f` is graceful with edge-difference map `pi`.

Equations `(5)--(6)` and `(16)--(20)` give an exact recognition theorem:

> The rooted data `(T,pi,G,mu,r_0)` arise from a graceful labeling with
> `f(r_0)=0` if and only if the decoded rows are valid intervals of lengths
> `1,...,m`, their anchored cut lifts form the ranked flag `(18)`, and the
> Gram compatibility `(6)` holds.

No search over `2^m` cut-side choices remains after fixing `r_0`.

## 4. Every marked gap-tail decision is then explicit

THM-2787 requires a suffix column at a particular gap `t` and the marked
rank equation

```text
|t-iota_t(f(v))|=k.                       (21)
```

THM-2789 recovers how many suffix columns exist but not their positions.
The ramp reconstruction recovers the actual ordered column `c_t`; rooted
cut inversion recovers `f(v)`.  Thus `(21)` and every marked gap-tail
decision are determined exactly by `(T,pi,G,mu,r_0)`.

This identifies the complete controlled-forgetting ladder:

| retained data | exactly recovered | still lost |
|---|---|---|
| `G` and ambient length | incidence-column multiset | column order, cut side |
| `G` plus oriented ramp | ordered gap matrix | complementary vertex side |
| plus root-zero anchor | full graceful labeling and marked tests | existence on a new tree |

The theorem supplies a certificate and reconstruction map.  It does not
construct the required data for an arbitrary tree.

## 5. Both added coordinates are sharp

### 5.1 No ramp

THM-2789's three-point interval families with starts

```text
(0,0,0) and (1,0,0)                     (zero-based starts)          (22)
```

have the same Gram matrix and column multiset but different, non-reflected
column orders.  Their one-based ramp moments are respectively

```text
(1,3,6) and (2,3,6).                       (23)
```

Thus Gram data alone cannot perform ordered reconstruction, while one ramp
separates the hostile exactly.  Pairing with the constant vector would only
repeat the row lengths and would not separate it.

### 5.2 No root anchor

For the one-edge tree, the unique internal gap column is the same for

```text
f=(0,1) and f=(1,0).                       (24)
```

They are the two complementary lifts of the same cut.  Declaring which
vertex carries zero selects exactly one.  Therefore an unrooted gap matrix
cannot recover marked ranks without at least this complement-gauge choice.

## 6. Exact verification

Run

```bash
python 04-computation/oriented_ramp_rooted_gap_reconstruction_thm2793.py
python -O 04-computation/oriented_ramp_rooted_gap_reconstruction_thm2793.py
```

The exact companion uses explicit exceptions and no truth-bearing
assertions.  It checks ramp reconstruction on all `46,233` one-of-each-
length interval families through rank eight.  It then checks all `2,650`
graceful labelings of every nonisomorphic tree through seven vertices:
`15,314` rooted threshold-cut inversions and `140,516` marked gap-tail
decisions.  It also verifies both sharp hostiles.

```text
PROVED HERE (candidate):
  one-ramp endpoint and ordered-row reconstruction;
  graceful ramp/edge-sum identity;
  root-zero mod-two cut inversion;
  exact rooted graceful-data recognition criterion;
  recovery of every marked gap-tail decision;
  Gram-only and anchor-free hostiles;
  exact interval and graceful controls through the stated ranks.

NOT PROVED:
  existence of the certificate for every tree;
  the marked gap-tail extension conjecture;
  the Graceful Tree Conjecture;
  uniqueness without an oriented reference and root gauge;
  an LRC, modular-group, Keller, JC(2), DC(2), knot, or tournament
  consequence.                                                         (25)
```

QED, conditional only on candidate status promotion after independent
hostile audit.
