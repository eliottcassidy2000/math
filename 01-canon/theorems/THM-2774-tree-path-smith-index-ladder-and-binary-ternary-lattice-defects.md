---
id: THM-2774
title: "Tree-path Smith-index ladder and binary/ternary lattice defects"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Every nonsingular full B_m-root frame has elementary two-primary
  Smith form.  A k-edge tree geodesic, however, canonically combines k-1
  D-root differences with its endpoint A-root to give Smith form
  diag(1^(k-1),k), with explicit quotient the signed coordinate sum mod k.
  The quotient is canonically the A_(k-1) weight/root quotient P/Q, and
  dually its torus kernel is the centre mu_k, the k-point diagonal torsion
  fibre.  Off-path coordinate roots extend this to a full ambient frame.
  Thus the first nonstar path P4 carries the first Z/3 defect, and diameter
  at least p supplies a Z/p frame quotient.  These are frame-local
  partial-cube lattice cokernels, not PSL2(Z), graceful existence, Keller,
  or LRC(14).
source: a4-resolvent-next-gate/tree-path-smith-index-2026-07-28
depends_on:
  - THM-2766-quadratic-cubic-pullback-even-sign-kummer-plane-and-weyl-d3-s4
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
script: 04-computation/tree_path_smith_index_ladder_thm2774.py
output: 05-knowledge/results/tree_path_smith_index_ladder_thm2774.out
script_sha256: 6dc6983ccc506da978354a010b7da48e6e43f18dc99c6b3689f3cd014eaea6fe
output_sha256: 9be974965a92b4cb9a5294e53cb86ac878ef888573f928edadacf52980380d10
hash_basis: LF-normalized bytes
---

# THM-2774 -- tree geodesics carry an exact Smith-index ladder

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2770 identifies the graceful-tree arrangement with a type-`B` signed
frame plus one new wall for each long tree geodesic.  The distinction is
arithmetic as well as geometric.  Pure `B` frames have only elementary
binary lattice defects.  A path of length `k` creates a cyclic defect of
exact order `k`.

This gives a literal, limited answer to why two and three recur.  The binary
defect is already present in signed pairs.  The first wall outside the
signed-pair arrangement is a three-edge path and carries `Z/3`.  Longer
paths introduce every later path length as well, so two and three are first,
not exclusive.

## 1. The partial-cube source map

Let `T` be a tree with root `r`, edge set `E`, and root-outward edge
coordinates `y_e`.  Work in the root-zero representative `x_r=0` of vertex
potentials modulo common translation.  Define

```text
iota_r(v)_e=1 if e lies on path(r,v), and 0 otherwise.   (1)
```

Then

```text
x_v=iota_r(v) dot y,                                    (2)
```

which is the inverse incidence formula of THM-2770.  Moreover

```text
support(iota_r(v)-iota_r(u))=path(u,v),
||iota_r(v)-iota_r(u)||_0=dist_T(u,v).                  (3)
```

Thus `iota_r` is the standard isometric partial-cube embedding of a tree,
and every vertex `A`-root has normal

```text
g_uv=iota_r(v)-iota_r(u) in {-1,0,1}^E.                 (4)
```

For a `k`-edge geodesic `P=(e_1,...,e_k)`, reverse individual edge
coordinates if necessary.  This signed-coordinate gauge preserves the
entire `D`-root set and sends the endpoint root `(4)` to

```text
g_P=(1,...,1) in Z^P.                                   (5)
```

Equations `(1)--(5)` are the source map.  The target below is an explicit
finite cokernel of a selected normal frame; it is not the cokernel of the
whole arrangement, which already contains a coordinate basis.

## 2. Pure `B_m` frames are elementary two-primary

A positive `B_m` hyperplane normal is, up to row sign,

```text
e_i,                    e_i-e_j,                    e_i+e_j.  (6)
```

Choose any `m` such normals and let `M` be their `m x m` row matrix.

**Binary-frame theorem.**  If `det M!=0`, then for some `c>=0`

```text
Smith(M)=diag(1^(m-c),2^c),
|det M|=2^c.                                             (7)
```

*Proof.*  Regard a row `e_i` as a half-edge at vertex `i` and a row
`e_i-s e_j` as a signed edge `ij`.  Parallel edges with opposite signs are
allowed, so this is a signed multigraph.  Split it into connected components
and permute rows and columns accordingly.  A nonzero square determinant
forces every component block to have equally many rows and vertices.

Suppose a connected block has `v` vertices, `q` ordinary edges, and `h`
half-edges.  Then `q+h=v`.  Its ordinary edges must connect all `v`
vertices, so `q>=v-1` and hence `h<=1`.

- If `h=1`, the ordinary graph is a tree.  Vertex sign switches make all
  tree rows ordinary differences, and leaf elimination together with the
  half-edge gives Smith form `1^v`.
- If `h=0`, the ordinary graph is unicyclic.  Switch signs along a spanning
  tree.  Leaf elimination reduces the block to its signed cycle, whose last
  determinant is `1-product_cycle(s)`.  A balanced cycle gives determinant
  zero; an unbalanced cycle gives one final invariant factor `2` and all
  other factors `1`.

Taking the direct sum over components proves `(7)`.  In particular, the
smallest binary frame is

```text
C_2=[[1,-1],[1,1]],             Smith(C_2)=diag(1,2).    (8)
```

It is the length-two unbalanced cycle consisting of the two parallel edges
with opposite signs.

Pure signed-pair geometry can create several independent `C2` defects, but
it cannot create odd torsion in a full root frame.

## 3. Every tree path has Smith index equal to its length

On the gauged path coordinates of `(5)`, select the `k-1` lawful `D` roots

```text
e_1-e_2, e_2-e_3, ..., e_(k-1)-e_k                     (9)
```

and the endpoint `A` root `g_P=(1,...,1)`.  Their matrix is

```text
C_k = [ e_1-e_2 ]
      [ e_2-e_3 ]
      [    ...    ]
      [ e_(k-1)-e_k ]
      [ 1 1 ... 1 ].                                    (10)
```

The first `k-1` rows have a unit `(k-1)x(k-1)` minor.  Modulo their span,
all coordinate classes are equal; the last row is `k` times that primitive
class.  Equivalently, the explicit target map

```text
sigma_k:Z^k -> Z/k,        sigma_k(z)=sum_i z_i mod k   (11)
```

kills every row of `(10)`.  Since

```text
det C_k=k,                                               (12)
```

the row lattice is exactly `ker sigma_k`.  Therefore

```text
Smith(C_k)=diag(1^(k-1),k),
Z^k / rowspan_Z(C_k) isomorphic_to Z/k.                 (13)
```

There is a canonical Lie-lattice interpretation of the same quotient.  Put

```text
V={x in R^k:sum_i x_i=0},
Q_A={x in Z^k:sum_i x_i=0},
P_A={x in V:x_i-x_j in Z for all i,j}.                  (13a)
```

Here `Q_A` and `P_A` are the root and weight lattices of `A_(k-1)`.  The
orthogonal projection

```text
pi:Z^k -> V,          pi(z)=z-(sum_i z_i/k)h,
h=(1,...,1),                                             (13b)
```

has kernel `Z h` and image exactly `P_A`.  The difference rows `(9)` map
onto `Q_A`.  Therefore

```text
Z^k/(Q_A+Z h)  isomorphic_to  P_A/Q_A  isomorphic_to Z/k. (13c)
```

The class of `pi(e_1)` generates: all `pi(e_i)` differ by roots, and
`j*pi(e_1)` lies in `Q_A` first at `j=k`.  Thus `(13c)` proves the kernel,
surjectivity, and exact order without a choice of Smith operations.

The dual root-of-unity picture is equally explicit.  On the additive torus
`T=R/Z`, equations `C_k theta=0` first force
`theta_1=...=theta_k=t`, and the last row forces `kt=0`.  Hence

```text
ker(C_k:T^k -> T^k)
 ={(j/k,...,j/k):j in Z/k}.                              (13d)
```

So the Smith defect is literally a `k`-point common-root fibre, not merely
an abstract determinant.  In type-`A` language, `(13c)` is the character
group of the centre `mu_k` of the simply connected `A_(k-1)` group, and
`(13d)` is that centre in diagonal torus coordinates.

This is the promised source-to-target map.  In the original edge gauge, the
rows in `(9)` become `+/-e_i +/-e_(i+1)`, still lawful `D` roots, and the
last row becomes the actual endpoint root `(4)`.

If `T` has `m>k` edges, append the coordinate roots `e_j` for every edge
outside `P`.  After permuting columns this gives a full `m x m` normal frame
with Smith form

```text
diag(1^(m-1),k).                                         (14)
```

Thus the `Z/k` index survives in the full ambient edge lattice.  It remains
frame-local: other arrangement roots are not being quotiented out.

## 4. The first two primes and the complete diameter ladder

The two-edge case `(8)` already lies entirely in `B_m`: its endpoint
two-edge `A` root is also the `D` root `e_1+e_2`.  A three-edge path is the
first time the endpoint root has support outside `B_m`.  For `P4`, it gives

```text
C_3=[[1,-1,0],
     [0,1,-1],
     [1,1,1]],                   Smith(C_3)=diag(1,1,3). (15)
```

Hence the minimal nonstar graceful clutch has, on the same three edge
coordinates,

```text
binary signed-pair defect: Z/2,
ternary long-path defect:  Z/3.                          (16)
```

Dually, `(8)` has the two diagonal torus points `0,1/2`, while `(15)` has
the three diagonal points `0,1/3,2/3`.  Among all `120` triples of the ten
distinct `P4` clutch normals, the complete determinant histogram is

```text
det 0:19,              det 1:73,
det 2:25,              det 3:3.                         (16a)
```

The three index-three frames are the three spanning-tree choices of two
difference roots on the three edge coordinates, together with the long-path
root.  Under THM-2766's half-Hadamard `D3=A3` identification, `h=(1,1,1)` is
one of the four even-sign body diagonals.  `W(D3)=S4` permutes those four
diagonals, and the stabilizer of the marked `h` is the coordinate
`S3=W(A2)`.  Thus the three comparison trees, the `A2` weight/root quotient,
and the ternary fibre are the same marked-rank-three structure.  Marking
`h` selects one `V4`-torsor state; it still does not select modular
generators or an affine monodromy action.

This is the lattice-index shadow of the rank-three `D3=A3` coincidence in
THM-2766.  It is not that theorem's `V4 semidirect S3` action: `(16)` consists
of two selected-frame cokernels, with no multiplication or monodromy action
between them.

Every diameter-`d` path contains a subpath of every length
`1,...,d`.  Consequently `(13)--(14)` give

```text
diameter(T)>=k  implies  T has a full normal frame of index k;
diameter(T)>=p  implies  T has a frame quotient Z/p      (p prime).  (17)
```

Thus primes `2` and `3` are special as the first signed-pair and first
long-geodesic defects.  Trees of larger diameter also exhibit `5,7,...`;
there is no theorem here that only two and three can occur.

## 5. What the index preserves and destroys

The frame construction preserves:

- the chosen tree geodesic and its partial-cube cut coordinates;
- the exact `A` endpoint wall and `D` comparison walls;
- the integral frame index and the explicit cyclic character `(11)`;
- the weight/root quotient `(13c)` and diagonal torsion fibre `(13d)`;
- every prime divisor of the path length.

It destroys the rest of the hyperplane arrangement, polynomial coefficient
signs, a label-grid evaluation, a distinguished modular generator, and any
physical LRC owner/current phase.  In particular:

- `Z/2` and `Z/3` in `(16)` are not the free factors in
  `PSL2(Z)=C2*C3` from THM-2768;
- a frame index does not repair the zero balanced `P4` coefficient of
  THM-2770 and does not prove graceful-label existence;
- a hypothetical transport to an LRC carry lattice would require an actual
  relation/current map, not just the shared prime;
- no quartic field monodromy, Keller exclusion, `JC(2)`, or `DC(2)` follows.

## 6. Exact verification

Run

```bash
python 04-computation/tree_path_smith_index_ladder_thm2774.py
python -O 04-computation/tree_path_smith_index_ladder_thm2774.py
```

The exact companion uses explicit exceptions, no floating point, and no
truth-bearing Python assertions.  It exhausts all `55,041` full `B_m` frames
for `m<=5`; all `26,431` nonsingular frames satisfy `(7)`.  It verifies
`(10)--(13)` for `2<=k<=13`, and on all `5,913` recursive trees through
eight vertices it checks `118,004` geodesic frames, their partial-cube
supports, their full ambient extensions, and the complete shorter-path
ladder inside every diameter.  It also checks the weight/root generator
orders and diagonal torus kernels for `2<=k<=13`, the `D3` body-diagonal
orbit/stabilizer, and the full `P4` minor histogram `(16a)`.  Normal and
optimized runs byte-match the stored transcript.

```text
PROVED HERE (candidate):  elementary two-Smith form for every pure-B frame;
                          partial-cube endpoint-root map;
                          path frame determinant and Smith form;
                          explicit quotient sum mod k;
                          A_(k-1) weight/root quotient realization;
                          exact diagonal k-torsion torus kernel;
                          full ambient extension;
                          P4 first Z/3 defect;
                          all path lengths/prime divisors through diameter.

NOT PROVED:               a global arrangement cokernel;
                          that only primes two and three occur;
                          a modular free-product action;
                          a balanced coefficient or graceful labeling;
                          a Keller/JC(2)/DC(2) result;
                          an LRC(14) relation/current realization.       (18)
```

QED (candidate).
