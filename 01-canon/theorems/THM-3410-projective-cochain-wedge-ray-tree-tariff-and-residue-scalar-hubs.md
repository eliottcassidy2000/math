---
id: THM-3410
title: "Projective cochain wedge, ray-tree tariff, and residue scalar hubs"
status: >
  PROVED analytic integral-wedge and projective-ray tree theorem plus
  PROVED-ELEMENTARY p=2,3,5 scalar-hub families; VERIFIED-EXACT companion;
  INDEPENDENT AUDIT REQUESTED.  Every realized THM-3398 affine cochain is the
  decomposable integral wedge of its owner vector with its mode-centre lift.
  Contracting equal primitive rays gives exact minimum-spanning-tree and
  bottleneck-tree tariffs.  The existing parity, ternary, and five-colour
  half-grid partitions carry explicit nonzero scalar fibres with exact tree
  tariffs.  Fixed-scalar dilation has degree two, but the full hidden scalar
  fibre has sharp degree-three diameter.  No LRC(14) ledger decrement follows.
source: cochain-projective-hub-2026-08-15
audit: self-contained determinant, primitive-ray contraction, strict residue-partition, and scale proofs; 92400 bounded lattice hostiles; 44 direct physical scalar hubs through q=40; q15 half-grid and exceptional-edge regressions; normal/optimized outputs byte-identical; independent immutable-file audit requested
depends_on:
  - THM-3398-general-finite-mode-sheet-cover-cochain
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3408-fixed-zero-additive-order-duality-and-six-core-corridor
  - THM-3409-q15-exceptional-edge-positive-cochain-rigidity-and-leakage-tariff
script: 04-computation/lrc_projective_cochain_wedge_scalar_hub_thm3410.py
output: 05-knowledge/results/lrc_projective_cochain_wedge_scalar_hub_thm3410.out
script_sha256: 5d4ee3abad66c8052a0d78406aa633f29b0bbb110c76438ea2f417774800d7f7
output_sha256: ec8451b180e96287d9dfecd3bfa484d7b4a8e1ac82108694c4656c83968e0179
semantic_sha256: bc6cc2bd14fb566ee7e57adc27291e3dea1085cf542891128cc4feea93d44c51
hash_basis: LF-normalized bytes
---

# THM-3410 -- projective cochain wedge, ray-tree tariff, and residue scalar hubs

**PROVED analytic + PROVED-ELEMENTARY + VERIFIED-EXACT; INDEPENDENT AUDIT
REQUESTED.**

## 1. Inheritance and the change of coordinates

[THM-3398](THM-3398-general-finite-mode-sheet-cover-cochain.md) represents a
physical cyclic-sheet cover by selected consecutive modes and a compatible
integral complete cochain.  [THM-3409](THM-3409-q15-exceptional-edge-positive-cochain-rigidity-and-leakage-tariff.md)
shows that this cochain can be compressed to a faithful spanning tree, but its
tree search is presented there as a finite optimization on one exceptional
edge.  MISTAKE-389 supplies the canonical hostile to forgetting the sidecar:
the same synchronized residue partition may have large nonzero mode-centre
drift.

The corrected near miss is therefore to regard the pair values as fifteen or
more unrelated integers.  They are the two-by-two minors of one integral
two-row matrix.  The least-used sidecar is the **integral mode-centre lift**.

The live concept board is:

| object | representation | invariant / operation | lost coordinate |
|---|---|---|---|
| complete cochain | integral wedge `A wedge u` | primitive projective rays | mode blocks and widths |
| faithful tree | weighted ray-quotient graph | MST / bottleneck threshold | signs and non-tree pairs |
| half-grid packet | integer centering error `E` | `E wedge u=A wedge u` | physical scalar after chart quotient |
| residue partition | `p` colour classes | scalar multiplication `a` | owner-to-colour assignment |
| scale | fixed-ray dilation / full scalar fibre | degrees two / three | which operation changes the scale |

The exact connection contract is:

| field | content |
|---|---|
| source | a physically realized THM-3398 selected-mode packet |
| target | integral columns `(u_i,A_i)`, then a weighted graph of primitive rays |
| map | `A_i=2q u_i x_i`, `P_ij=A_i u_j-u_i A_j`, then contract equal rays |
| preserved | every signed pair value, weighted triangle closure, zero components, both tree tariffs, source-lift gauge, and pure dilation |
| destroyed | selected sheet blocks, interval widths, endpoints, and the cover predicate itself |
| required sidecar | the original modes and widths when one wants a physical cover rather than a tariff |
| cheapest decisive tests | MISTAKE-389's q15 ternary partition and THM-3409's exceptional q15 packet |

The native carrier is an integral projective line arrangement with weights.
There is no need to force it into a tournament: zero determinants are genuine
parallel-ray components, while absolute values alone discard the cyclic ray
order.

## 2. Every realized complete cochain is an integral wedge

Fix `q>=2`.  Let a realized packet have positive distinct owners
`u_0,...,u_(r-1)` and containing lifted mode centres

```text
x_i=(n_i+h_i/(2q))/u_i,             n_i in Z.          (1)
```

Put

```text
A_i=2q u_i x_i=2q n_i+h_i in Z.                        (2)
```

Then THM-3398's pair cochain is exactly

```text
P_ij=2q u_i u_j(x_i-x_j)
    =A_i u_j-u_i A_j.                                  (3)
```

Thus `P=A wedge u` is a decomposable integral two-form.  In particular its
weighted triangle closure

```text
u_k P_ij-u_j P_ik+u_i P_jk=0                           (4)
```

is the Laplace relation among the two-by-two minors of the matrix with columns
`(u_i,A_i)`.  A spanning tree is faithful because its determinant differences
recover all column slopes `A_i/u_i` up to one common additive constant; `(4)`
then recovers every non-tree pair.

Changing the real source lift by an integer `b` sends

```text
A_i -> A_i+2qb u_i.                                    (5)
```

This is one common integral unimodular shear of the columns.  It preserves
their contents, primitive determinants, the complete cochain, and every tariff.
Hence the construction is independent of the representative of physical time
modulo one.

## 3. Primitive rays and exact contraction of the tree problem

For each column define

```text
c_i=gcd(u_i,A_i)>0,
rho_i=(s_i,b_i)=(u_i/c_i,A_i/c_i),
gcd(s_i,b_i)=1,             s_i>0.                    (6)
```

Use the signed wedge convention

```text
[rho_i,rho_j]=b_i s_j-s_i b_j.                         (7)
```

Equations `(3)` and `(6)` give the exact factorization

```text
P_ij=c_i c_j [rho_i,rho_j].                            (8)
```

Two owners lie on the same rational ray exactly when `P_ij=0`.  Partition the
owners into equal-ray classes `C`.  Put

```text
mu_C=min_(i in C)c_i.                                  (9)
```

For distinct ray classes `C,D`, the primitive determinant
`Delta_CD=|[rho_C,rho_D]|` is a positive integer, and `(8)` implies

```text
min_(i in C,j in D)|P_ij|=mu_C mu_D Delta_CD.          (10)
```

Give the complete class graph this edge weight.  Then, for this packet,

```text
tau_1 = minimum-spanning-tree weight of (10),
tau_infinity = least threshold at which the graph (10) is connected. (11)
```

Here the tariffs minimize over spanning trees on all original owners.  Given
any such tree, contract each equal-ray class, discard loops, and delete cycle
edges.  What remains contains a class tree whose edges are among the paid
cross edges, and every one costs at least `(10)`.  Conversely choose one pair
attaining `(10)` for every edge of a class tree and add zero trees inside the
classes.  This lifts the class tree without extra cost.  The same argument
works for the largest edge, so it proves both identities.

If there are `K>=2` classes and their minima in increasing order are
`mu_(1)<=...<=mu_(K)`, integrality of every `Delta_CD` gives the universal
content tariffs

```text
tau_1 >= mu_(1) sum_(j=2)^K mu_(j),
tau_infinity >= mu_(1)mu_(K).                          (12)
```

Indeed root any class tree at a class of minimum content.  The edge entering
class `j` costs at least `mu_(1)mu_(j)`.  The edge entering a maximum-content
class proves the second bound.  When `K=1`, the cochain and both tariffs are
zero.  The first inequality need not be sharp: the exceptional q15 control in
Section 7 has lower bound `8` and exact tariff `10`.  Its bottleneck bound `3`
is sharp.

## 4. Two exact scales

Under THM-3398's pure degree-`lambda` dilation,

```text
(q,u_i,x_i) -> (lambda q,lambda u_i,x_i/lambda).       (13)
```

Equation `(2)` gives

```text
A_i -> lambda A_i,
c_i -> lambda c_i,
rho_i -> rho_i,
P_ij -> lambda^2 P_ij.                                (14)
```

Therefore both tariffs scale by `lambda^2`, and

```text
tau_1/q^2,             tau_infinity/q^2                (15)
```

are exact pure-dilation invariants of the projective packet.

There is a second integral frame at a synchronized half-grid physical time
`c`, meaning `H_i=2q u_i c in Z` for every active owner.  Define the centering
error

```text
E_i=A_i-H_i=2q u_i(x_i-c) in Z.                        (16)
```

Since `H_i/u_i=2qc` is common,

```text
P_ij=E_i u_j-u_i E_j.                                 (17)
```

If the selected mode has THM-3398 width `w_i`, physical membership in its open
interval gives the exact bound

```text
|E_i|<w_i/7.                                           (18)
```

Thus MISTAKE-389's missing coordinate is an integral bounded error vector, not
an endpoint ambiguity.  Vanishing cochain means that the error columns lie on
one rational ray; half-grid integrality alone does not force that ray to be
zero or common across owners.

## 5. Scalar hubs over the residue partitions

The even quadratic ladder extends to the existing ternary and five-colour
half-grid minima.  Let `p in {2,3,5}`, put `q=pd`, and use the following owner
words:

```text
V_2(d)=(1,2d-1),                                      d>=4;

V_3(3)=(1,5,7),
V_3(d)=
 (1,2d-1,2d-2), d=0 mod 3,
 (1,2d,  2d-1), d=1 mod 3,                            d>=5 odd,
 (1,2d-2,2d  ), d=2 mod 3;

V_5(d)=(1) union ({2d-2,...,2d+2}\5Z),
                                      d>=5 odd, 3 not|d. (19)
```

The last line is written in increasing order after the initial `1`.  In every
case `V_p(d)\{1}` contains one representative of every nonzero residue modulo
`p`.  Define

```text
Delta_p(d)=max(1,max_(v in V_p(d),v!=1)|2d-v|).        (20)
```

Thus `Delta_2=1`, `Delta_5=2`, and `Delta_3` is one or two.  For every integer
`a>=1` satisfying

```text
gcd(a,p)=1,                    7a Delta_p(d)<q,         (21)
```

take

```text
c_a=a/(2dq),                   U_(p,d)=d V_p(d).        (22)
```

Then the owners in `(22)` partition `Z/qZ` into its `p` residue classes.  The
owner `v=1` covers class zero; owner `v!=1` covers

```text
ell == -a v^(-1)  (mod p).                             (23)
```

Their unique containing single-atom mode centres and integral lifts are

```text
x_1=0,                         A_1=0,
x_v=a/(pdv),                   A_v=2ad       (v!=1).   (24)
```

Consequently the entire cochain is

```text
P_(1,v)=-2ad^2,
P_(v,w)=2ad^2(w-v)             (v,w!=1),               (25)
```

and its exact tree tariffs are

```text
tau_1=2a(p-1)d^2=2a(p-1)q^2/p^2,
tau_infinity=2ad^2=2aq^2/p^2.                          (26)
```

### Proof of the partition and tariffs

Every `v` in `(19)` is a unit modulo `p`, so owner `dv` is constant on sheet
classes modulo `p`.  For the root, the centre error in `(16)` is `-a`.  For a
nonroot owner, the proposed centre in `(24)` has error

```text
E_v=a(2d-v).                                           (27)
```

Condition `(21)` and `(18)` put the physical time strictly inside each
single-atom interval.  The dangerous class is `(23)`, obtained from
`v ell == -a mod p`.  Since the danger arc has length `1/7` and the `p` phase
classes have spacing `1/p>=1/5`, no second class can fire.  Multiplication by
the unit `a` permutes the nonzero classes, proving the exact partition.

Equation `(25)` is direct substitution of `(24)` into `(3)`.  Every root edge
has absolute weight `2ad^2`; every nonroot edge has this weight times the
positive integer `|w-v|`.  Any tree on `p` vertices therefore costs at least
`p-1` copies of `2ad^2` and has maximum edge at least `2ad^2`.  The root star
attains both bounds, proving `(26)`.

For `p=2`, `(22)` is exactly

```text
q=2d, c_a=a/q^2, U=(q/2,q(q-1)/2),                    (28)
```

so `(25)` recovers the even ladder `P=-aq^2/2`.  At `p=3,d=5,a=1`, it gives
MISTAKE-389's q15 packet `(5,40,50)` with

```text
A=(0,10,10), P=(-50,-50,100),
(tau_1,tau_infinity)=(100,50).                         (29)
```

The same mechanism now supplies genuine scalar fibres over the ternary and
five-colour partitions, not only parity.

## 6. The sharp family-diameter exponent is three

For fixed `a`, equation `(26)` has the same degree-two homogeneity as the
dilation law `(14)`.  These scalar-hub packets are not asserted to be literal
pure dilates: their owner words vary with `d`.  The chart quotient nevertheless
forgets `a`, and `(21)` permits `a` of order `q` along unbounded subfamilies.
Hence the same unlabeled residue partition supports tariffs of order `q^3`.
More precisely,

```text
tau_infinity/q^3 < 2/(7 Delta_p(d)p^2),
tau_1/q^3 < 2(p-1)/(7 Delta_p(d)p^2).                  (30)
```

These are uniform on the displayed scalar fibres.  The exponent is sharp:
in the parity family take `q=14m` and `a=2m-1`.  Then `(21)` holds and

```text
tau_infinity/q^3=tau_1/q^3=a/(2q) -> 1/14.            (31)
```

Thus `q^-2` is the correct normalization for **one fixed projective scalar
under pure dilation**, but it is not compact over the whole scalar fibre.
Uniform family diameter has order `q^3`.  This separates two operations that
the bare block quotient makes look identical.

## 7. Exact and hostile controls

The standard-library companion performs two independent kinds of replay.

1. On `92,400` integral matrices with `2<=r<=5`, distinct positive owners
   through seven, and lifts in `{-2,-1,0,1,2}`, it checks `(3)--(12)` exactly.
   There are `118` one-ray cases and `92,282` multiray cases.  The clean
   contraction hostile

   ```text
   u=(2,4,3,6), A=(1,2,-1,-2)                         (32)
   ```

   has four owners but two rays; both exact tariffs are `5` after the two
   zero components are contracted.
2. Starting from the original strict danger predicate rather than predicted
   masks, it reconstructs every admissible scalar hub through `q=40`: exactly
   `44` instances, split `30/11/3` for `p=2/3/5`.  It independently finds the
   containing centres, checks the residue partition, forms `(3)`, and
   enumerates all labelled trees.  At `q=39,p=3,a=3`, failure of the unit
   condition leaves only `13` of `39` sheets.  At `q=35,p=5,a=3`, failure of
   the strict radius leaves only `21` of `35` sheets.

For THM-3409's exceptional edge, the canonical positive representative has

```text
u=(1,7,8,9,11,13),
A=(0,-2,-2,-3,-3,-4),
c=(1,1,2,3,1,1),                                      (33)
```

and primitive rays

```text
(1,0),(7,-2),(4,-1),(3,-1),(11,-3),(13,-4).           (34)
```

The ray MST gives `(tau_1,tau_infinity)=(10,3)`, reproducing enumeration of
all `6^4=1296` trees.  The content bounds `(12)` are `(8,3)`.  Dilations by
`1,2,3,7` preserve the normalized tariffs

```text
(tau_1/q^2,tau_infinity/q^2)=(2/45,1/75).             (35)
```

This simultaneously supplies a positive control for the projective compiler
and a hostile to claiming that the first content bound is always an equality.

Run

```text
python3 04-computation/lrc_projective_cochain_wedge_scalar_hub_thm3410.py
python3 -O 04-computation/lrc_projective_cochain_wedge_scalar_hub_thm3410.py
```

The two outputs must match the stored artifact byte for byte.  The companion
contains no floating literal, random choice, external solver, network call, or
optimization-sensitive assertion.

## 8. Information boundary and nonclaims

The projective graph is a faithful compiler for the affine cochain and its
tree tariffs.  It is not a cover classifier: it forgets blocks, widths,
strict endpoints, and the common Helly intersection.  Conversely, the bare
block partition, owner quotient orders, and additive-order densities forget
the scalar `a` and cannot control even the `q^-2` normalized tree tariff.
THM-3408's fixed-zero fractional dual therefore neither sees nor contradicts
the positive-cochain scalar hubs.

The displayed owners are unrestricted positive transverse owners, often far
outside `{1,...,14}`.  No map into a surviving refined LRC row, wall weight,
owner chronology, or safe-mass ledger is proved.  The theorem supplies a new
exact sidecar and a sharp obstruction to quotient-only transport; it makes no
LRC(14) ledger decrement and does not prove LRC(14).

**QED.**
