---
source: codex-2026-07-25-planar-chirp-quotient-resolution
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Restricting the
  eta-average in planar chirp
  tomography to a character subgroup Xi does not recover an individual
  Gram edge: it recovers exactly the sum over one derivative-coset of the
  annihilator Xi-perp. For F_(p^2) and a one-dimensional F_p chirp line,
  this is a p-point affine-line sum. All p+1 projective chirp directions
  reconstruct the edge field by one affine-plane incidence identity;
  omitting any direction has a nonzero zero-mean hostile, and two
  transverse directions admit a genuine-signal checkerboard hostile.
  At p=13, one Bockstein line therefore has resolution 13. THM-2337 makes
  all 14 directions algebraically available as derived functionals on each
  chosen planar graph, so every off-diagonal edge has an exact
  coefficient formula; this is not an independently supplied physical
  pair-probe. The exact remaining boundary is aligned singleton
  location/energy. Pair-polarized graph chirps sharpen this to one
  common-target atom. Same-frequency graph polarization cannot locate that
  atom; cross-frequency coherence can, and is the precise missing anchor.
  The averaged graph-channel difference energy is positive exactly when a
  nonzero target survives, and its sum over graphs dominates THM-2337's
  word-support energy by the sharp factor 1/169. No LRC profile is excluded.
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2355-quadratic-grouped-current-repair
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
script: 04-computation/lrc14_planar_chirp_quotient_resolution_probe.py
output: 05-knowledge/results/lrc14_planar_chirp_quotient_resolution_probe.out
script_sha256: 99efb2277febbf938c2af4081b50b0e093192f97b83827ca58cfac812f6a74c8
output_sha256: a82d26dac963f9895d245e9042e36bb7d7d9a4c4881127f6b5eb9738068d5ac5
hash_basis: working-tree bytes (LF)
---

# Planar chirps have an exact quotient-resolution ladder

THM-2356 identifies a full planar target--jet transform that recovers every
off-diagonal Gram entry. PROVED THM-2337 supplies its coefficient-level LRC
input as a joint target/first-jet table; after a field identification, the
table splits into the resulting planar graphs. The chirp intensities are then
derived quadratic functionals of those coefficients, not independently
measured scalar-cover observables. The correct refinement is exactly which
aggregate each projective subbank returns.

The answer is a derivative-coset sum. In the `F_13^2` specialization this
is finite affine Radon data, and its resolution loss can be computed
without estimates.

## 1. Subgroup inversion

Let `G,A` be finite abelian groups of common order `q`, and let
`phi:G->A` be planar:

```text
D_h phi(u)=phi(u+h)-phi(u)
```

is bijective for every `h!=0`. For a signal `z:G->C`, write

```text
F_z(eta,chi)=sum_x z(x) eta(phi(x)) chi(x),

E_z(eta,chi)=|F_z(eta,chi)|^2.                     (1)
```

Fix a subgroup `Xi<=Ahat`, and let

```text
Xi_perp={a in A: eta(a)=1 for every eta in Xi}.    (2)
```

For `h!=0` define the restricted inverse

```text
T_Xi(h,y)
 =1/(q |Xi|) sum_(eta in Xi) sum_(chi in Ghat)
    E_z(eta,chi)
    conj(chi(h)) conj(eta(D_h phi(y))).             (3)
```

Then

```text
T_Xi(h,y)
 =sum_{u:
       D_h phi(u)-D_h phi(y) in Xi_perp}
    z(u+h) conj(z(u)).                              (4)
```

This is exact. Expanding `E_z` and summing in `chi` first forces
`x-u=h`. The subgroup orthogonality relation

```text
sum_(eta in Xi) eta(a)
 =|Xi| if a in Xi_perp, and 0 otherwise            (5)
```

then gives (4). Since `D_h phi` is bijective, every aggregate in (4)
contains exactly

```text
|Xi_perp|=q/|Xi|                                   (6)
```

Gram edges.

The endpoints of the ladder are familiar:

```text
Xi={1}:          the full lag-h autocorrelation;

Xi=Ahat:         the single edge z(y+h)conj(z(y)).  (7)
```

No positivity or asymptotics enter this statement.

## 2. One F_13 jet line has resolution thirteen

Let `p` be an odd prime, take `K=F_(p^2)`, and fix a nontrivial
additive character `psi:F_p->C^*`. Define

```text
Psi_K(t)=psi(Tr_(K/F_p)(t))
```

and put

```text
phi(x)=x^2/2,

eta_a(t)=Psi_K(a t),
```

and choose a one-dimensional `F_p` subspace

```text
L=F_p a_0 <= K.
```

The corresponding chirp subgroup has `p` characters. Since

```text
D_h phi(u)-D_h phi(y)=h(u-y),                      (8)
```

formula (4) becomes

```text
T_L(h,y)
 =sum_(u: Tr(a_0 h(u-y))=0)
    z(u+h)conj(z(u)).                               (9)
```

The constraint is one affine `F_p` line containing `p` points. Thus for
`p=13` a single projective Bockstein/chirp subbank has exact resolution
`13`, not atomic resolution.

This conclusion is invariant under the chosen additive bases. Changing
bases changes which affine-line direction is called vertical; it does not
change the fibre size.

## 3. The full projective bank and its sharp boundary

For fixed `h`, transport the edge field through the planar derivative:

```text
g_h(v)
 =z(u+h)conj(z(u)),       u=(D_h phi)^(-1)(v).      (10)
```

The one-dimensional chirp subgroups are the `p+1` projective directions
in the dual affine plane. Their restricted inversions are precisely the
sums of `g_h` over all affine lines of the corresponding direction.

Let `R_L g(line_L(v))` denote the line sum through `v`. In an affine plane,
each other point lies with `v` on exactly one of the `p+1` lines through
`v`. Hence

```text
sum_(L in P^1(F_p)) R_L g(line_L(v))
 =p g(v)+sum_w g(w),                                (11)
```

and therefore

```text
g(v)
 =1/p [
    sum_(L in P^1(F_p)) R_L g(line_L(v))
    -sum_w g(w)
   ].                                               (12)
```

The total `sum_w g(w)` is itself obtained by summing the `p` parallel
line sums in any one direction. So the `p+1` direction bank reconstructs
every edge exactly.

This is sharp for arbitrary edge fields. If one dual direction `L_0` is
omitted, choose a nonzero zero-mean function `c` on `F_p` and put

```text
g(v)=c(l_0(v)).
```

Every retained affine line meets each `l_0`-level once, so every retained
line sum is zero; the omitted direction sees `g`.

Even two transverse directions fail on a genuine signal. To state the
hostile in its physical start coordinate, write

```text
gtilde_h(u)=z(u+h)conj(z(u)).
```

For the quadratic derivative, multiplication by `h` pulls projective
line directions bijectively between the `v`-plane in (10) and the
`u`-plane. Choose two chirp directions whose pulled-back normals are the
row and column normals. At `p=13`, put the lag

```text
h=(3,3)
```

and prescribe the four nonzero edge values

```text
gtilde(0,0)=+1,   gtilde(0,1)=-1,

gtilde(1,0)=-1,   gtilde(1,1)=+1.                  (13)
```

Set `z(u)=1` at those four starts and
`z(u+h)=gtilde(u)` at the four disjoint endpoints, with `z=0` elsewhere.
Then `gtilde_h` is exactly the displayed checkerboard and its row and
column sums vanish. Transport by `D_h phi` gives the same hostile for the
two selected chirp directions. Thus two transverse jet axes do not become
sufficient merely because the edge field came from a scalar signal.

## 4. LRC connection contract

THM-2337 gives the joint table `A(q,z)` on `G x B`. After chosen
identifications with `F_169`, define its `c`-th planar-graph signal by

```text
Z_c(q)=A(q,q^2/2+c).                                (14)
```

The joint table is canonical, while these graphs depend on the chosen
additive identifications and resulting field multiplication. For any such
choice, the `169` graphs partition the table. On graph `c`, a linear target/jet
character restricts as

```text
Psi_K(bq+az)=Psi_K(ac) Psi_K(bq+a q^2/2),
```

and the constant `Psi_K(ac)` cancels from intensity. Thus every projective
subbank in this theorem, including the full fourteen-direction bank, is an
exact derived functional of THM-2337's coefficients. This is algebraic
availability, not a claim that an independent physical pair-probe
observable was already present in the scalar cover.

```text
source:
  THM-2337's joint target/jet coefficient array, grouped along a chosen
  planar graph and transformed by a target-jet character subgroup;

target:
  off-diagonal grouped-current Gram edges;

map:
  character inversion in the target variable and in the selected jet
  subgroup;

preserved:
  lag h, derivative-coset address, and the exact sum of all edge currents
  in that coset;

destroyed:
  the individual source point inside a coset of Xi_perp;

exact invoice:
  q/|Xi| edges per returned aggregate;

F_169 boundary:
  one F_13 jet line returns 13-edge line sums; all 14 projective
  directions return atoms by affine-plane inversion;

needed sidecar:
  a phase-coherent cross-frequency target anchor, an aligned singleton
  energy/location theorem, or a proof that some nonzero graph signal has
  support at least three, so its diagonal follows from the recovered
  off-diagonal Gram entries;

cheapest decisive hostile:
  test every proposed two-axis service against the genuine-signal
  checkerboard (13), and every claimed full-bank landing against a
  singleton graph signal at q=0.                    (15)
```

Because the graphs partition a nonzero joint Abel array, at least one
`Z_c` is nonzero. The full projective bank reconstructs all its
off-diagonal edges. What it cannot do uniformly is locate a singleton:
every one-point signal has the same graphwise chirp-intensity table,
independent of its target. Hence the remaining zero-fibre debt is
cross-frequency coherence/aligned singleton energy, not construction of
the coefficient-level planar grouping.
No scalar profile is excluded.

## 5. Cross-graph polarization leaves one common-point hostile

There is a further exact compression of the singleton boundary. For graph
`c`, write

```text
F_c(a,b)=sum_q Z_c(q) Psi_K(bq+a q^2/2).            (16)
```

This is the graph-offset-calibrated amplitude. The raw restriction of the
joint character carries the known factor `Psi_K(ac)`; remove it before
polarization, or multiply the raw cross product for graphs `c,d` by
`Psi_K(-a(c-d))`.

Suppose derived pair-polarized graph intensities at one common chirp
setting `(a,b)` are retained. The four real
measurements

```text
|F_c+F_d|^2, |F_c-F_d|^2,

|F_c+iF_d|^2, |F_c-iF_d|^2
```

recover the complex cross product

```text
C_(c,d)(a,b)=F_c(a,b)conj(F_d(a,b)).                (17)
```

The same planar inversion as THM-2356, now bilinear rather than
sesquilinear in one signal, gives for every `h!=0`

```text
Z_c(y+h)conj(Z_d(y))
 =1/q^2 sum_(a,b) C_(c,d)(a,b)
    Psi_K(-bh-a(hy+h^2/2)).                         (18)
```

Consequently exactly one of the following holds:

```text
(A) some recovered cross edge with h!=0 is nonzero;

(B) every nonzero graph signal Z_c is supported at one common target q_0.
                                                               (19)
```

Indeed, the self-products `c=d` show that each graph in case (B) has
support at most one. If two nonzero graphs were supported at distinct
targets, their cross product would give a nonzero edge with the difference
of those targets as `h`.

Case (A) already lands a nonzero target, because two distinct target
positions cannot both be zero. Thus graph-pair polarization collapses the
entire residual to a **common target atom**.

If target landing still fails, that common point is `q_0=0`, and the
joint table has the exact vertical-tensor form

```text
A(q,z)=1_(q=0) B(z),       B(z)=A(0,z).             (20)
```

Conversely every such tensor has
`Z_c(q)=B(c)1_(q=0)` on every chosen graph. Taking THM-2333's formal
full-term target hostile `delta_0(q)` and any nonzero/full-support jet
profile `B` shows the boundary is sharp at the group-algebra level.
For `B=1`, all `169` graph signals survive and all are zero-target
singletons. This is not asserted to come from canonical interval weights.

That last atom is a genuine stopping boundary for same-frequency
graph-polarized intensities.
If

```text
Z_c(q)=alpha_c 1_(q=q_0)
```

for every graph, then every same-frequency graph-linear combination has

```text
|sum_c lambda_c F_c(a,b)|^2
 =|sum_c lambda_c alpha_c|^2,                       (21)
```

independent of `q_0` (the common chirp phase factors out). Including the
known graph-offset phases `Psi_K(ac)` changes the coefficients
`lambda_c alpha_c` but still does not expose `q_0` at that common `(a,b)`.

This stopping rule must not be extended across frequencies. The single
plus-interference family already distinguishes `q_0=0` from `q_0!=0`:

```text
|F_c(0,0)+F_c(0,b)|^2
 =|alpha_c|^2 |1+Psi_K(bq_0)|^2.                   (22)
```

It is even in `q_0`, so it does not by itself distinguish `q_0` from
`-q_0`. Full quadrature polarization across frequencies recovers

```text
D_b
 =F_c(0,0)conj(F_c(0,b))
 =|alpha_c|^2 Psi_K(-bq_0),                         (23)
```

from the four intensities with `F_c(0,0) +/- F_c(0,b)` and
`F_c(0,0) +/- iF_c(0,b)`. Nondegeneracy of the trace pairing makes
`b -> D_b/D_0` determine `q_0` uniquely. Thus the sharp next service is
precisely one phase-coherent cross-frequency target anchor, or
equivalently one aligned singleton energy—not another same-frequency
graph polarization.

There is also an exact averaged detector which does not try to name the
target. Put `F_c(b)=F_c(0,b)` and

```text
D_c^diff
 =1/169 sum_b |F_c(b)-F_c(0)|^2
 =sum_(q!=0)|Z_c(q)|^2
   +|sum_(q!=0)Z_c(q)|^2.                          (24)
```

The second equality is Parseval plus character orthogonality. Consequently
`D_c^diff>0` if and only if graph `c` has a nonzero-target coefficient.
This pins the remaining LRC service more economically than full tomography:
retain a lawful phase ratio between the graph channels `b` and `0`, or a
pair-twist measurement of their difference. Separate channel intensities
do not determine (24). At the THM-2337 coefficient level it is an exact
derived quadratic functional; no claim is made that scalar-cover dynamics
already expose it as a physical observable.

## 6. The graph detector dominates the older word-support energy

The new graph-channel debt is not independent of THM-2337's word-support
quadratic form. Let

```text
A_0(q)=sum_z A(q,z),        D_graph=sum_c D_c^diff.
```

Because `c=z-q^2/2` is a bijection in `z` for each fixed `q`, expanding
(24) gives

```text
D_graph
 =sum_(q!=0,z)|A(q,z)|^2
   +sum_c |sum_(q!=0)A(q,q^2/2+c)|^2.              (25)
```

The second term is nonnegative. Cauchy in the `169` jet coordinates
therefore yields

```text
D_graph
 >=1/169 sum_(q!=0)|A_0(q)|^2.                     (26)
```

For each of THM-2337's three word-support masks `P_sigma`, the mask is
zero at `q=0` and takes only the values zero and one. Its exact landing
energy

```text
E_sigma=sum_q P_sigma(q)|A_0(q)|^2
```

hence satisfies the quantitative bridge

```text
D_graph >= E_sigma/169.                            (27)
```

The constant is sharp for arbitrary joint arrays. Fix two distinct
nonzero targets `q_1,q_2` and a nonzero scalar `alpha`, and set, for every
graph address `c`,

```text
A(q_1,q_1^2/2+c)= alpha,
A(q_2,q_2^2/2+c)=-alpha,
```

with all other entries zero. Every graph sum in the second term of (25)
cancels, while equality holds in both Cauchy inequalities for a mask
containing `q_1,q_2`.

Thus any positive lower bound for the older word-support energy
automatically pays for the graph-channel detector. Conversely,
`D_graph` may remain positive when the uncoloured marginal `A_0(q)`
cancels, so the jet-resolved graph route is strictly more informative.
This identifies one shared target rather than two competing ones:
prove positivity of a canonical target-support quadratic form, using the
jet resolution only where the uncoloured marginal loses phase.

## 7. Exact verification

Run

```text
python3 04-computation/lrc14_planar_chirp_quotient_resolution_probe.py
python3 -O 04-computation/lrc14_planar_chirp_quotient_resolution_probe.py
```

The companion checks all `28,392` ordered distinct point pairs in
`F_13^2`, the unique affine direction through each pair, exact
reconstruction from all fourteen direction families, a nonzero hostile
after omitting one family, and the genuine-signal two-direction
checkerboard. Normal, optimized, and stored transcripts agree after LF
normalization.

The separate independent THM-2356 audit companion additionally checks
`80` exact finite-field instances of (27), including the sharp
two-target equality hostile.
