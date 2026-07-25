---
id: THM-2330
title: "Partition-lattice Gordian lift spectrum and exact merge cocycle"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED, WITH CITED KNOT INPUTS.
  For a labelled knot packet, a partition of its labels, and a
  target knot, the minimum marked product distance over all target
  factorizations is attained. Its excess over ordinary Gordian distance is
  exactly the obstruction to a compartment-preserving geodesic lift. These
  lift costs decrease under partition coarsening; their cover drops form an
  exact nonnegative coboundary, compose by min-plus convolution on disjoint
  packets, and admit a stable homogenization. At the unknot, the construction
  is exactly the partition transform of the higher subset defect from
  THM-2248. Brittenham--Hermiller's T(2,7)-mirror example gives a finite and
  stable sole-merge drop of at least one, and their cited large-gap family
  makes finite sole-cover drops unbounded; no exact value of the connected
  sum's unknotting number is claimed. Exact conical word-metric controls show
  both that the root defect complex does not determine target-conditioned
  lifting and that every fixed-arity root-defect truncation, hence every
  graph or tournament built only from original pair defects, can miss the
  final partition-refinement drop. The abstract controls are not asserted
  to be knot examples.
source: codex-2026-07-25-partition-gordian-lift
depends_on:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2308-mirror-double-nakanishi-floor-and-sharp-stable-mixture-profile
  - THM-2317-alexander-fibre-fan-and-prime-marked-gordian-quotient
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
script: 04-computation/partition_gordian_lift_spectrum_thm2330.py
output: 05-knowledge/results/partition_gordian_lift_spectrum_thm2330.out
script_sha256: cbf30bed327fe095ad019a2531d7fa10d7e06a5ba9413515ff04f780e50e4717
output_sha256: 2df58cc7462dba503236253de1729d853067b12b02546a19be43ace07303e9b9
hash_basis: working-tree bytes (LF)
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
  - "Horst Schubert, Die eindeutige Zerlegbarkeit eines Knotens in Primknoten, Sitzungsberichte der Heidelberger Akademie der Wissenschaften 1949/3, 57--104."
---

# THM-2330 -- the missing object is a partition-indexed lift spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED, WITH CITED KNOT
INPUTS.**

THM-2317 identifies connected-sum nonadditivity as strict contraction under
the map which forgets labelled factor compartments. Its root defect records
only the two extreme states:

```text
every factor separate  --->  every factor forgotten.
```

There is a canonical object between those extremes. It asks, for every
partition of the labels and every target knot, how cheaply a path can reach
the target while respecting exactly those remaining compartments. The
result is a monotone function on the partition lattice. Its differences are
the precise costs of forgetting successive factor walls.

## 1. The target-indexed lift cost

Let `I` be a nonempty finite label set and let

```text
x=(K_i)_(i in I)
```

be a labelled packet of oriented knots. For `A subset I`, write

```text
K_A=#_(i in A) K_i,
K_emptyset=U.                                        (1)
```

Let `Pi(I)` be the lattice of set partitions of `I`, ordered by refinement:

```text
pi <= rho
```

means that every block of `pi` is contained in a block of `rho`. Denote the
discrete and indiscrete partitions by `0hat` and `1hat`.

For `pi in Pi(I)`, form the Cartesian Gordian product

```text
Gamma_pi=product_(B in pi) Gamma_G
```

with product distance and connected-sum realization

```text
d_pi(y,z)=sum_(B in pi)d_G(y_B,z_B),

Sigma_pi(y)=#_(B in pi)y_B.                          (2)
```

The prescribed starting point is

```text
x^pi=(K_B)_(B in pi).                               (3)
```

For a target knot `J`, define its factorization fibre and the
**partition-indexed Gordian lift cost**

```text
F_pi(J)={y in Gamma_pi:Sigma_pi(y)=J},

Lambda_x(pi;J)
 =min_(y in F_pi(J)) d_pi(x^pi,y).                  (4)
```

The minimum exists. The fibre is nonempty: put `J` in one compartment and
`U` in all the others. All distances are nonnegative integers, so
well-ordering already gives attainment. More concretely, Schubert prime
decomposition makes `F_pi(J)` finite: a factorization is obtained by
distributing the finitely many prime-summand occurrences of `J` among the
labelled blocks of `pi`.

Define the **geodesic-lift obstruction**

```text
Omega_x(pi;J)
 =Lambda_x(pi;J)-d_G(K_I,J).                        (5)
```

> **Geodesic-lift theorem.** For every `pi` and `J`,
>
> ```text
> Omega_x(pi;J)>=0.                                 (6)
> ```
>
> Equality holds if and only if some Gordian geodesic from `K_I` to `J`
> lifts, edge for edge, to a path in `Gamma_pi` beginning at `x^pi` and
> ending at some factorization in `F_pi(J)`.

### Proof

THM-2317 proves that `Sigma_pi` maps every marked product edge to a genuine
Gordian edge. Hence it preserves the length of every chosen marked path and
is nonexpansive on endpoint distances:

```text
d_G(K_I,J)
 <=d_pi(x^pi,y)                    for y in F_pi(J). (7)
```

Minimizing gives (6).

If equality holds, choose a minimizing endpoint `y` and a product geodesic
from `x^pi` to `y`. Its image has length

```text
Lambda_x(pi;J)=d_G(K_I,J),
```

so the image is an unmarked geodesic. Conversely, the length of any
edge-for-edge lift of an unmarked geodesic is `d_G(K_I,J)`, forcing the
reverse inequality in (6). QED.

Only existence of one lift is asserted. Other geodesics may fail to lift.
The datum depends on the labelled packet and partition, not only on the
isotopy class of `K_I`.

## 2. Coarsening and the exact merge cocycle

If `pi<=rho`, merge all `pi`-coordinates inside each block of `rho`:

```text
m_(pi,rho):Gamma_pi -> Gamma_rho,

(m_(pi,rho)(y))_C=#_(B in pi, B subset C)y_B.       (8)
```

This graph map satisfies

```text
m_(pi,rho)(x^pi)=x^rho,
Sigma_rho after m_(pi,rho)=Sigma_pi.                (9)
```

> **Partition monotonicity.** If `pi<=rho`, then
>
> ```text
> Lambda_x(rho;J)<=Lambda_x(pi;J),
> Omega_x(rho;J)<=Omega_x(pi;J).                    (10)
> ```

Indeed, (8) sends every endpoint in `F_pi(J)` to one in `F_rho(J)`.
Connected-sum nonexpansivity inside each new block gives

```text
d_rho(x^rho,m_(pi,rho)(y))<=d_pi(x^pi,y).           (11)
```

Minimize over `y`.

For comparable partitions define the **merge drop**

```text
c_x^J(pi,rho)
 =Lambda_x(pi;J)-Lambda_x(rho;J)       (pi<=rho).   (12)
```

Then

```text
c_x^J(pi,rho)>=0,

c_x^J(pi,tau)
 =c_x^J(pi,rho)+c_x^J(rho,tau)
                    whenever pi<=rho<=tau.          (13)
```

Thus `c` is the exact nonnegative coboundary of the potential `Lambda` on
the partition poset. In particular, the sum of cover drops along any chain
depends only on its endpoints, and every partition-lattice diamond obeys
the corresponding square identity.

At the indiscrete partition there is one coordinate and one possible target
tuple, so

```text
Lambda_x(1hat;J)=d_G(K_I,J),
Omega_x(1hat;J)=0.                                  (14)
```

Consequently

```text
Omega_x(pi;J)=c_x^J(pi,1hat).                       (15)
```

The obstruction is exactly the total distance lost while all remaining
compartment walls are forgotten.

## 3. The root slice is the partition transform of higher defect

At `J=U`, Schubert cancellation says

```text
F_pi(U)={(U)_(B in pi)}.
```

Therefore

```text
Lambda_x(pi;U)=sum_(B in pi)u(K_B),                 (16)

Omega_x(pi;U)
 =sum_(B in pi)u(K_B)-u(K_I).                       (17)
```

Let THM-2248's labelled subset defect be

```text
Delta(A)=sum_(i in A)u(K_i)-u(K_A).                 (18)
```

Substitution in (17) gives the exact partition transform

```text
Omega_x(pi;U)
 =Delta(I)-sum_(B in pi)Delta(B).                   (19)
```

Hence the root slice contains no information beyond the full weighted
subset-defect function. The new information in (4) is the target variable
and its path-lifting meaning.

For the discrete partition,

```text
Omega_x(0hat;U)
 =sum_i u(K_i)-u(K_I)
 =Delta(I).                                         (20)
```

If `rho` covers `pi` by merging exactly two blocks `A,C in pi`, then

```text
c_x^U(pi,rho)
 =u(K_A)+u(K_C)-u(K_(A union C)).                   (21)
```

Thus each root cover drop is THM-2176's symmetric pair defect, but evaluated
on the two *current composite blocks*. A tournament on the original labels
does not retain those recursively formed vertices.

## 4. Exact min-plus composition

Let `I,H` be disjoint finite label sets, with packets `x` on `I` and `z` on
`H`, and partitions `pi in Pi(I)`, `eta in Pi(H)`. Their disjoint-union
partition has no block meeting both sides.

> **Min-plus convolution law.**
>
> ```text
> Lambda_(x disjoint_union z)(pi disjoint_union eta;J)
>  =min_(A#C=J)
>      [Lambda_x(pi;A)+Lambda_z(eta;C)].             (22)
> ```

The minimum ranges over ordered connected-sum factorizations of `J`.

### Proof

Any endpoint factorization of `J` on the left groups into

```text
A=#_(B in pi)y_B,
C=#_(D in eta)y_D,
A#C=J.                                              (23)
```

Its cost is the sum of its `I`- and `H`-coordinate costs, hence is at least
the corresponding right side of (22). This proves `>=`.

Conversely, choose a minimizing factorization `A#C=J` and minimizing
endpoints for the two terms. Their disjoint union is an admissible endpoint
on the left with exactly the summed cost, proving `<=`. All minima are
attained by the integer-valued argument following (4). QED.

In min-plus notation, if functions on knots are convolved by

```text
(f tensor g)(J)=min_(A#C=J)[f(A)+g(C)],
```

then (22) says

```text
Lambda_(x disjoint_union z)^(pi disjoint_union eta)
 =Lambda_x^pi tensor Lambda_z^eta.                  (24)
```

This is an exact compositional law for the scalar lift spectrum. It does not
recover the actual minimizing endpoint factorization or path; those remain
the witness sidecar.

## 5. Stable lift spectrum

Write `nK=#^n K`. For positive integers `n`, put

```text
a_n(pi;J)=Lambda_(nx)(pi;nJ),                       (25)
```

where `nx` replaces every `K_i` by `nK_i`. If `y^(n)` and `y^(m)` are
minimizing endpoint factorizations, their coordinatewise connected sums are
an admissible factorization for `a_(n+m)`. Joint nonexpansivity gives

```text
a_(n+m)(pi;J)<=a_n(pi;J)+a_m(pi;J).                 (26)
```

Fekete's lemma therefore defines the **stable lift cost**

```text
Lambda_x^hash(pi;J)
 =lim_(n->infinity)a_n(pi;J)/n
 =inf_(n>=1)a_n(pi;J)/n.                            (27)
```

Likewise

```text
d_G^hash(A,J)
 =lim_(n->infinity)d_G(nA,nJ)/n                    (28)
```

exists. Define

```text
Omega_x^hash(pi;J)
 =Lambda_x^hash(pi;J)-d_G^hash(K_I,J).              (29)
```

Dividing the finite inequalities by `n` and taking limits proves

```text
Omega_x^hash(pi;J)>=0,

pi<=rho
 implies
 Lambda_x^hash(rho;J)<=Lambda_x^hash(pi;J),         (30)

Omega_x^hash(pi;J)
 =lim_(n->infinity)Omega_(nx)(pi;nJ)/n.             (31)
```

Stable merge drops are again potential differences and hence obey (13).
At the two boundary slices,

```text
Lambda_x^hash(1hat;J)=d_G^hash(K_I,J),              (32)

Lambda_x^hash(pi;U)=sum_(B in pi)u_hash(K_B),       (33)

Omega_x^hash(pi;U)
 =sum_(B in pi)u_hash(K_B)-u_hash(K_I).             (34)
```

Equation (33) again uses the singleton root fibre. Stable zero means only
that the lift obstruction is sublinear; it need not produce a finite
geodesic lift.

## 6. The Brittenham--Hermiller pair

Let

```text
K=T(2,7),             Kbar=mirror(K),
x=(K,Kbar),           J=U.                          (35)
```

The cited one-body values are

```text
u(K)=u(Kbar)=3,
```

and Brittenham--Hermiller prove

```text
u(K#Kbar)<=5.                                       (36)
```

The two-label partition lattice has a single cover. Equations (14) and (16)
give

```text
Lambda_x(0hat;U)=6,
Lambda_x(1hat;U)=u(K#Kbar)<=5,                     (37)

c_x^U(0hat,1hat)
 =Omega_x(0hat;U)
 =6-u(K#Kbar)>=1.                                  (38)
```

Thus their shortcut is exactly the positive cost of forgetting the sole
factor wall. Equation (38) does **not** assert that the drop is exactly one:
the exact value of `u(K#Kbar)` is not imported.

The signature lower bound is sharp on every self-power of `K` and `Kbar`,
so

```text
u_hash(K)=u_hash(Kbar)=3.
```

Repeating the cited path in (36) gives `u_hash(K#Kbar)<=5`. Therefore

```text
Omega_x^hash(0hat;U)
 =6-u_hash(K#Kbar)>=1.                              (39)
```

The bypass persists linearly after homogenization, but (39) still supplies
no exact stable diagonal.

There is also no universal upper bound on a single partition-cover drop.
Brittenham--Hermiller's cited Corollary 1.6 gives, for every `N>=1`,
infinitely many knots `L` satisfying

```text
u(L#mirror(L))<=2u(L)-N.                            (39a)
```

For the two-label packet `(L,mirror(L))`, the unique cover therefore has

```text
c_(L,mirror(L))^U(0hat,1hat)
 =2u(L)-u(L#mirror(L))>=N.                          (39b)
```

This is an unbounded **finite** merge-drop statement. No stable lower bound
for these additional families is asserted.

## 7. Two exact hostile controls

The following controls are abstract conical metric monoids, not knots.
They delimit what can follow from the metric-monoid laws alone.

### 7.1 Root defect does not determine a target lift

Let `M=N^3`, let `e_1,e_2,e_3` be the standard basis, and compare two
translation-invariant word metrics inherited from `Z^3`.

The first uses only `+/-e_i`, each of cost one. The second adds

```text
v=e_1+e_2-e_3
```

as a generator of cost one. For `x=(e_1,e_2)`, both metrics have

```text
ell(e_1)=ell(e_2)=1,
ell(e_1+e_2)=2.                                    (40)
```

Thus their entire root subset-defect data on this packet agree and vanish.
At target `J=e_3`, however, the ordinary metric has

```text
d(e_1+e_2,e_3)=3,
Lambda_x(0hat;e_3)=3,
Omega_x(0hat;e_3)=0,                               (41)
```

while the metric with `v` has

```text
d(e_1+e_2,e_3)=ell(v)=1,
Lambda_x(0hat;e_3)=3,
Omega_x(0hat;e_3)=2.                               (42)
```

Indeed the only two ordered factorizations of `e_3` into two elements of
`N^3` are `(e_3,0)` and `(0,e_3)`, and in the second metric

```text
ell(e_1-e_3)=ell(e_2-e_3)=2.                       (43)
```

The new generator gives a global shortcut which cannot be assigned to
either nonnegative compartment. This is the smallest target-conditioned
reason that (19) cannot replace the full spectrum.

### 7.2 Every fixed-arity root-defect truncation can miss the terminal merge

Fix `r>=3`, put `M=N^r`, and let

```text
g=e_1+...+e_r.
```

Compare the ordinary `l_1` metric with the restriction of the weighted word
metric on `Z^r` having generators

```text
+/-e_i of cost 1,          +/-g of cost r-1.        (44)
```

For a subset `S` of size `s`, direct minimization over the net number `k`
of `g`-letters gives

```text
ell(e_S)
 =min_(k in Z)
   [|k|(r-1)+s|1-k|+(r-s)|k|].                     (45)
```

Hence

```text
ell(e_S)=s                  for S proper,
ell(g)=r-1.                                        (46)
```

The ordinary `l_1` metric instead has `ell(g)=r`, while agreeing on every
proper subset. Because `N^r` is conical, every factorization of zero is
coordinatewise zero. For the packet `(e_1,...,e_r)`, (16) therefore gives

```text
ordinary l_1:
  Lambda(pi;0)=r and Omega(pi;0)=0 for every pi;

metric (44):
  Lambda(pi;0)=r and Omega(pi;0)=1 for pi!=1hat,
  Lambda(1hat;0)=r-1 and Omega(1hat;0)=0.           (47)
```

Every cover drop is zero until the last two blocks are merged, when the
drop is one. In particular all original pair defects are zero in both
metrics, yet the terminal merge differs. More generally the two metrics
agree on every interaction of arity less than `r`.

Thus no tournament or weighted graph built only from the original pair
defects can recover the partition-refinement spectrum in the axiomatic
category of conical commutative nonexpansive metric monoids. This is not a
claim that an actual higher-order knot packet exists; it is an obstruction
to proving pairwise completeness from the formal laws.

## 8. Source, map, loss, and repair

| source object | map / target | preserved | information destroyed | repair sidecar |
|---|---|---|---|---|
| partitioned Cartesian Gordian product `Gamma_pi`, based at `x^pi` | connected-sum realization `Sigma_pi` into the ordinary Gordian graph | every chosen local crossing path, its length, total endpoint, relabelling, disjoint-union composition | which compartment owns each crossing and which factorization of `J` is reached | retain `pi`, the endpoint tuple, and the marked path |
| distance from `x^pi` to `F_pi(J)` | scalar `Lambda_x(pi;J)` | optimal compartment-respecting cost | minimizing endpoint and path multiplicity | lift-witness fibre |
| partition spectrum | root slice `J=U` | full subset defect through (19) | every non-root target and its lifting obstruction | target variable `J` |
| full partition lattice | original-label pair graph or tournament | binary root drops | composite-block vertices and the coarsening stage where mixing first pays | merge cocycle on all comparable partitions |
| finite spectrum | stable spectrum | linear asymptotic rate | bounded and sublinear lift obstruction | finite `Omega_n` profile |

The faithful object remains the marked path category and its endpoint
factorization. `Lambda`, `Omega`, and the merge cocycle are exact scalar
compressions adapted to the particular question "when does a shortest path
stop respecting these factor walls?"

## 9. Scope and failure boundary

1. The construction refines a labelled factorization, not knot isotopy of
   the total connected sum alone.
2. At the unknot it is exactly the transform (19) of THM-2248, not a second
   independent root invariant.
3. A positive `Omega` proves nonliftability of every geodesic through the
   chosen compartments; it does not identify the crossing changes of an
   unmarked shortcut.
4. Stable zero means sublinear loss, not a finite lift.
5. The word-metric controls prove an axiomatic no-go, not the existence of a
   higher-order knot symbiont.
6. The Brittenham--Hermiller application uses only `u(K#Kbar)<=5`; neither
   its exact finite value nor its exact stable value is claimed.
7. Symmetric merge drops have no intrinsic orientation. A tournament becomes
   lawful only after a separate asymmetric observable is supplied.

## 10. Exact companion

The companion enumerates set partitions and target factorizations in the
two conical word-metric controls. It independently checks:

```text
root-slice agreement but target-lift disagreement;
agreement through every proper arity and the final r-body drop;
partition monotonicity and all cover drops for r=3,...,6;
the partition diamond/coboundary identities;
the disjoint-packet min-plus convolution on a finite target bank.
```

Reproduce with

```bash
python3 04-computation/partition_gordian_lift_spectrum_thm2330.py
python3 -O 04-computation/partition_gordian_lift_spectrum_thm2330.py
```

Both transcripts must equal the stored output byte-for-byte after LF
normalization. The companion verifies the exact finite controls, not any
unknown knot value. The ordinary run owns the runtime assertions; the
optimized run is transcript-parity evidence because Python removes
`assert` statements under `-O`.
