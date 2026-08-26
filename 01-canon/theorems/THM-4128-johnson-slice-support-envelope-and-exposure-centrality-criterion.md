---
id: THM-4128
title: "Johnson-slice support envelope and exposure centrality criterion"
status: >
  PROVED ELEMENTARY DIRECTED-EXPOSURE QUADRATIC ENVELOPE + SHARP CENTRALITY
  CRITERION + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every tournament ear
  field of order n>=4, the rational THM-4127 support floors on all nonconstant
  cardinality layers form one strictly concave quadratic in the imbalance
  t=n-2m. Its curvature is exactly the disjoint-edge energy and its vertex is
  one raw field-degree correlation divided by that energy. Hence the optimal
  rational layer is the nearest admissible parity-grid point, with sharp
  even/odd centrality criteria. A rational optimizer remains optimal after
  ordinary odd rounding; exact layer-coset rounding and the actual maximizing
  layer require separate sidecars.
source: codex-frontier-synthesis-creative-20260825i
depends_on:
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4127-johnson-slice-hoeffding-variance-and-central-support-dominance
related:
  - THM-4131-strong-tournament-centrality-through-order-eight
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4135-strong-tournament-centrality-complete-order-nine
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - HYP-9029-strong-interval-tiling-law
script: 04-computation/tournament_johnson_slice_support_envelope_thm4128.py
output: 05-knowledge/results/tournament_johnson_slice_support_envelope_thm4128.out
independent_audit_script: 04-computation/tournament_johnson_slice_support_envelope_thm4128_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_johnson_slice_support_envelope_thm4128_independent_audit.out
script_sha256: f4d92806b6df81d9bea5bd0e10123e62c9a98e7493e5b32979e4116e51f0739c
output_sha256: 420808f7876c0fbbfb28875add9968fec08c5a4f7a764f89a12d9852b64692c3
semantic_sha256: d4dcc84bddb0e17848c97f578eaf5554191dea2d6c22841c1c9c14a80420b88f
independent_audit_script_sha256: 618a47772b369194d2d5fbbc60781fd7a35278d83be3f49d678d353bcb1dd437
independent_audit_output_sha256: 4d640604e1f898e45b2773c3acdae88fc81c3d7fe89411876fecaf0816923f51
independent_semantic_sha256: d4dcc84bddb0e17848c97f578eaf5554191dea2d6c22841c1c9c14a80420b88f
hash_basis: raw LF bytes
primary_audit: >
  PASS. A Start/End/exposed-gap subset DP scans every 33,856 labelled
  tournament of orders four through six. It verifies the directed exposure
  carrier, curvature collapse, every layer floor, nearest-grid optimizer,
  centrality criterion, rounding sidecars, and the two first-failure code
  orderings. Literal child DPs agree at codes 2, 20, and 30571. Normal,
  optimized, hash-seeded, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room permutation implementation imports no primary code.
  It enumerates zero-defect Hamilton paths and one-defect orderings, rebuilds
  the directed exposure cut, and independently repeats the complete census
  and named hostiles. Its semantic ledger agrees with the primary ledger;
  normal, optimized, hash-seeded, and frozen outputs byte-match.
---

# THM-4128 -- Johnson-slice support envelope and exposure centrality

**PROVED ELEMENTARY DIRECTED-EXPOSURE QUADRATIC ENVELOPE + SHARP CENTRALITY
CRITERION + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4127 computes the exact variance and support floor on a fixed Johnson
slice. Varying the slice might appear to require `n-1` unrelated
calculations. It does not: the entire family is one strictly concave
quadratic. The missing coordinate is the signed correlation between the
orientation field and the unsigned exposure degrees.

This theorem determines the best rational and ordinary-odd support floor.
It does not claim that the same layer maximizes the actual response, and it
does not absorb THM-4123's layer-specific coset.

## 1. Directed exposure carrier

Let `T` be a tournament on `V={1,...,n}`, `n>=4`, and retain THM-4115 and
THM-4127's notation

```text
F(S)=H+sum_(i in S)h_i+cut_w(S),
sum_i h_i=0,                  w_ij>=0,
d_i=sum_(j!=i)w_ij,           W=sum_(i<j)w_ij.             (1)
```

THM-4114 represents the same response by Start/End/exposed-gap capacities.
Path-flow conservation sharpens that representation. Orient every weighted
edge as the corresponding arc of `T`. If

```text
o_i=sum_(i->j)w_ij,           r_i=sum_(j->i)w_ij,          (2)
```

then

```text
h_i=o_i-r_i,                  d_i=o_i+r_i,                 (3)

F(S)=H+2 sum_(i in S,j notin S,i->j) w_ij.                 (4)
```

To see `(3)`, write `Q(i,j)` for the exposed ordered-gap count. A Hamilton
path contributes one outgoing good gap at every nonterminal occurrence and
one incoming good gap at every noninitial occurrence, so

```text
Start(i)-End(i)
 =sum_(i->j)Q(i,j)-sum_(j->i)Q(j,i).                       (5)
```

Combining `(5)` with `2w_ij=Q(i,j)+Q(j,i)` and the linear coefficient of the
directed cut gives `(3)`; summing an oriented edge's contribution gives
`(4)`. In particular, the raw correlation

```text
C_hd=sum_i h_i d_i
    =sum_i (o_i^2-r_i^2)                                   (6)
```

is the signed exposure-degree tilt.

## 2. One quadratic contains every slice floor

For an interior layer `1<=m<=n-1`, put

```text
q=n-m,                     t=n-2m,
G_n={-(n-2),-(n-4),...,n-4,n-2}.                           (7)
```

Thus `t` runs over the parity grid `G_n`. Retain THM-4127's edge-star and
zero-row-sum residual coordinates `u,z`, and put

```text
D_4=sum_(e<f, e intersect f=empty)w_e w_f,                 (8)

J_0=H+nW/(2(n-1))+||h||^2/(2W)
        +(n-2)||z||^2/(2(n-3)W).                           (9)
```

For odd `n`, `J_0` is a formal central intercept rather than an actual layer.

> **Theorem 1 (quadratic support envelope).** The rational support floor
> `J_m` of THM-4127 is
>
> ```text
> J(t)=J_0+C_hd t/(W(n-2))
>          -D_4 t^2/(W(n-2)(n-3)).                        (10)
> ```
>
> Moreover `D_4>0`, so `(10)` is strictly concave.

Indeed THM-4127 gives

```text
J_m=H+W(n^2-t^2)/(2n(n-1))
       +||h+t u||^2/(2W)
       +((n-2)^2-t^2)||z||^2/(2(n-2)(n-3)W).              (11)
```

The linear coefficient in `(11)` is

```text
<h,u>/W=C_hd/(W(n-2)).                                    (12)
```

Its negative quadratic curvature is `Gamma/(2W)`, where

```text
Gamma=W^2/(n(n-1))+||z||^2/((n-2)(n-3))-||u||^2.          (13)
```

Using the edge ANOVA identities

```text
||z||^2=||w||^2-2W^2/(n(n-1))-(n-2)||u||^2,
(n-2)^2||u||^2=sum_i d_i^2-4W^2/n,                        (14)
```

one obtains the collapse

```text
(n-2)(n-3)Gamma
 =W^2+||w||^2-sum_i d_i^2
 =2D_4.                                                    (15)
```

Substitution proves `(10)`. Finally, choose any Hamilton path
`v_1,...,v_n` of `T`. Its first and third exposed gaps make
`w_(v_1,v_2),w_(v_3,v_4)>=1/2`; they are disjoint, so `D_4>=1/4>0`.
**QED.**

## 3. Exact nearest-layer optimizer

Define the rational exposure vertex

```text
theta=(n-3)C_hd/(2D_4).                                   (16)
```

Completing the square in `(10)` gives

```text
J(t)=J_0+D_4 theta^2/(W(n-2)(n-3))
          -D_4(t-theta)^2/(W(n-2)(n-3)).                  (17)
```

> **Corollary 2 (nearest parity-grid layer).** The maximizing rational
> support floors are exactly the layers whose `t in G_n` minimizes
> `|t-theta|`. There is one maximizing layer except when `theta` lies exactly
> halfway between two consecutive admissible grid points.

Thus all `n-1` rational floors compress to `(W,D_4,C_hd,J_0)` and one clipped
nearest-grid operation. In particular:

```text
n even: the central layer t=0 is optimal
        iff (n-3)|C_hd|<=2D_4;                             (18)

n odd:  at least one central layer t=+-1 is optimal
        iff (n-3)|C_hd|<=4D_4.                             (19)
```

The inequalities include boundary ties. For even order, strict inequality in
`(18)` makes the central layer unique. For odd order, `theta=0` ties the two
central layers; a nonzero `|theta|<2` selects exactly one, while
`|theta|=2` ties that central layer with its next outer neighbor.

Since ordinary odd ceiling is nondecreasing, expanding the rational optimizer
also gives the best ordinary-odd rounded floor across all cardinalities. This
strictly strengthens the **rational** central-only floor whenever `(18)` or
`(19)` fails; its odd ceiling may only tie. It still retains THM-4127's
domination of the full-cube floor.

## 4. Two sidecars that the quadratic does not absorb

THM-4123's exact layer cosets vary with `m`. Their ceilings are not a common
monotone function of `J(t)`, so they can reorder the rational optimizer.
The first labelled example is nonstrong order-six code `140`:

```text
theta=108/77,

t=2: J=2122/21,  exact coset floor=103,
t=4: J=639/7,    exact coset floor=135.                   (20)
```

The `t=4` layer has lattice `102` and actual maximum `135`; it loses
rationally but wins after exact coset rounding. Therefore a bank seeking the
best layer-coset floor must retain the coset sidecar or test additional
layers.

Nor does the optimizer locate an actual maximizing cut. Strong order-six
code `20` has

```text
theta=5/109,
J(0)=10399/109 > J(2)=9671/109,
max_(t=0)F=131 < max_(t=2)F=133.                           (21)
```

The support floor is correctly central while the global response maximum is
not. This preserves THM-4123's hostile boundary.

## 5. Exact finite boundary

The primary audit reconstructs every response from Start/End/exposed-gap
tables for all `33,856` labelled tournaments of orders four through six.
It checks `(3)`, `(10)`, `(15)--(19)`, every rational and rounded floor, and
literal child Hamiltonian counts at codes `2`, `20`, and `30571`.

The exact optimizer histograms are:

| order | labelled / strong | rational optimizer `t` histogram |
|---:|---:|---|
| 4 | `64 / 24` | `(-2,0):8`, `(0):48`, `(0,2):8` |
| 5 | `1,024 / 544` | `(-1):160`, `(-1,1):704`, `(1):160` |
| 6 | `32,768 / 22,320` | `(-2):4,224`, `(0):24,320`, `(2):4,224` |

All labelled order-four and order-five tournaments satisfy the centrality
criterion. At order six, exactly `8,448` fail it, split evenly between
`t=-2` and `t=2`; none is strong. All `22,320` strong labelled order-six
parents are central. This last fact is **FINITE-EXACT**, not an all-order
strong-parent theorem.

The first failure is nonstrong order-six code `2`:

```text
(H,W,D_4,C_hd)=(3,45,375,282),       theta=141/125,
J(2)=629/15 > J(0)=1871/45.                               (22)
```

Ordinary odd rounding ties both at `43`, so `(22)` also separates rational
optimization from strict rounded improvement. Exactly `288` order-six
parents exhibit an exact-coset reordering such as `(20)`; none is strong.

## 6. Information budget and replay

The source is the full tournament-ear cut field. The operation is variation
of the conditioned cardinality, and the preserved rational packet is

```text
(H,W,||h||^2,||z||^2,D_4,C_hd).                           (23)
```

The map to the optimizer is `t=nearest_(G_n)(theta)`. It discards the layer
cosets, the arrangement of values inside a layer, and the location of an
actual maximum. THM-4123's `(a_(T,m),d_(T,m))` is the necessary arithmetic
sidecar for exact lattice rounding; response connectivity remains separate.

Run

```text
python3 -B 04-computation/tournament_johnson_slice_support_envelope_thm4128.py
python3 -B -O 04-computation/tournament_johnson_slice_support_envelope_thm4128.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_johnson_slice_support_envelope_thm4128.py
python3 -B 04-computation/tournament_johnson_slice_support_envelope_thm4128_independent_audit.py
python3 -B -O 04-computation/tournament_johnson_slice_support_envelope_thm4128_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_johnson_slice_support_envelope_thm4128_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte. THM-4131 later
verifies that every rational and exact-coset floor optimizer is central for
every strong tournament of orders four through eight, with a strict central
coset margin; THM-4135/4137 extend this through order ten. THM-4133 **REFUTES**
the all-order extension at order twelve. Order eleven,
actual-maximizer classification, slice intervals, and global `H`-spectrum
completeness remain **OPEN**.
**QED.**
