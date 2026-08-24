# Hopf/S6 manuscript intake and audit ledger -- 2026-08-24

**Status gate.**  The unsigned, undated direct PDF
[*The (3,4,infinity) modular family of 2-tori, completed at its three special
points, is a complex structure on S6*](https://alpo.ge/s6.pdf) is treated here
as a **MANUSCRIPT CLAIM / UNDER AUDIT**.  It was downloaded on 2026-08-24 as a
108-page PDF whose author/title metadata fields are blank.  Its creation and
modification dates are the non-informative `1980-01-01` epoch; the body gives
no author, date, or version.  Its SHA-256 is

```text
283bba102dd1d5dc346af81b28145bdaaea6654398d5032e76e97bafb9a858f2.
```

The title is transcribed from page 1.  The host is Levent Alpöge's personal
domain, but the homepage neither lists nor links this file.  Attribution to
Alpöge is therefore a plausible inference, not explicit bibliographic data in
the artifact.  No peer-review or publication status is inferred.  This file
separates four levels which must not be collapsed:

```text
PROVED IN REPO
FINITE-EXACT FROM DISPLAYED DATA
READ / NO IMMEDIATE CONTRADICTION
MANUSCRIPT CLAIM / OPEN REFEREE OBLIGATION.                  (1)
```

In particular, this repository does **not** currently record the Hopf problem
as solved.  The historical baseline is Agricola--Bazzoni--Goertsches--
Konstantis--Rollenske,
[*On the history of the Hopf problem*](https://arxiv.org/abs/1708.01068),
which states the classical question: whether `S6` admits an integrable complex
structure.  The new manuscript is potentially frontier-changing, so every
global consequence below stays conditional until specialist verification.

## 1. Primary-source bundle

1. **S6 manuscript:** direct PDF above, all 108 pages rendered and inspected;
   displayed integer data independently transcribed into the exact audit below.
2. **Published 2020 comparison theorem:** Campana--Demailly--Peternell,
   [*The algebraic dimension of compact complex threefolds with vanishing
   second Betti numbers*](https://arxiv.org/abs/1904.11179), published in
   *Compositio Mathematica* 156 (2020), 679--696,
   [DOI](https://doi.org/10.1112/S0010437X19007802).
3. **1996 predecessor:** Campana--Demailly--Peternell,
   [*The algebraic dimension of compact complex threefolds with vanishing
   second Betti number*](https://arxiv.org/abs/math/9607215).
4. **Exact repo companion:**
   [integer monodromy and Smith audit](../../04-computation/hopf_s6_triangle_monodromy_snf_audit_20260824.py)
   with [frozen output](../results/hopf_s6_triangle_monodromy_snf_audit_20260824.out).
5. **Exact gauge/filling companion:**
   [centralizer and elliptic-action audit](../../04-computation/hopf_s6_gauge_centralizer_filling_audit_20260824.py)
   with [frozen output](../results/hopf_s6_gauge_centralizer_filling_audit_20260824.out).
6. **Exact completed-overlap companion:**
   [`b=2` centralizer/cusp Cech audit](../../04-computation/hopf_s6_b2_completed_centralizer_cech_audit_20260824.py)
   with [frozen output](../results/hopf_s6_b2_completed_centralizer_cech_audit_20260824.out).
7. **Unconditional local theorem:**
   [THM-3955, node cotangent normalization kernel and conductor
   torsion](../../01-canon/theorems/THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion.md),
   extended at triple crossings by
   [THM-3957](../../01-canon/theorems/THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel.md).
8. **Published symmetry comparison:** Huckleberry--Kebekus--Peternell,
   [*Group Actions on `S6` and complex structures on `P3`*](https://arxiv.org/abs/math/9812076),
   proves that a complex `S6` is not almost homogeneous and has
   `h^0(X,T_X)<=2`.  Its proof uses the 1996 Campana--Demailly--Peternell
   constant-meromorphic-function conclusion, so it is a load-bearing
   comparison gate rather than an independent resolution of the manuscript's
   claimed algebraic-dimension-one conflict.

The S6 manuscript is a source under examination, not a proved dependency of
THM-3955, THM-3957, or of any established repo theorem.

## 2. The manuscript's claimed construction

The construction starts from the orbifold triangle group

```text
Delta(3,4,infinity)=<g1,g2 | g1^3=g2^4=1>
```

and a rank-four integral representation with matrices `T1,T2`.  Put
`T0=(T1*T2)^(-1)=I+N`, with `N^2=0` and `rank N=2`.  On the dual lattice the
monodromies are

```text
A1=(T1^(-1))^t,   A2=(T2^(-1))^t,   M0=(T0^(-1))^t.     (2)
```

The claimed geometric spine is:

```text
equivariant period functions tau, mu, beta
  -> complex two-torus family over the punctured orbifold curve
  -> logarithmic-transform fillings at orbifold orders 3 and 4
  -> infinite-fan toric filling at the unipotent cusp
  -> translations by local sections during the three gluings
  -> a compact complex threefold X fibred over P1.       (3)
```

At the cusp the claimed central fibre `W` is the quotient of the prescribed
`A2` triangulation: its normalization is the degree-six del Pezzo surface,
opposite hexagon sides are identified, the double locus has three rational
curves, and there are two triple points.  Thus `W` is nonnormal.  At the two
elliptic points the reduced fibres are claimed bielliptic surfaces with
multiplicities 3 and 4.

The local-section translations are compressed to three integers
`(ell0,ell1,ell2)`.  The manuscript claims

```text
p=12 ell0-4 ell1-3 ell2,
pi1(X)=Z/|p|,
(ell0,ell1,ell2)=(0,1,-1) gives p=-1.                   (4)
```

For this choice it then claims integral homology equal to `S6`, invokes the
simply-connected homology-sphere recognition step, and concludes that `X` is
diffeomorphic to `S6`.  It also claims algebraic dimension one.

## 3. Dependency spine and nonfinite gates

The manuscript supplies two homology computations, but they share local
geometric inputs and therefore are not fully independent certificates.

### Route A: collapse and Mayer--Vietoris

```text
toric cusp quotient and collapse model
  -> homology of W and its collar
  -> special-fibre inclusion maps
  -> Mayer--Vietoris maps
  -> integral homology of X.                            (5)
```

### Route B: nearby cycles and Leray

```text
integral nearby-cycle local systems
  -> local invariant and coinvariant lattices
  -> special-point contributions and differentials
  -> Leray groups and extensions
  -> integral homology of X.                            (6)
```

The first genuinely nonfinite gates are:

1. proper discontinuity and compactness for the infinite toric fan quotient;
2. correctness of the cusp collapse model after the self-identifications;
3. integral nearby cycles for the non-simple-normal-crossing, nonnormal fibre;
4. equality of the geometric attaching maps with the displayed integral
   matrices, including orientations and meridians;
5. the final smooth-manifold and recognition hypotheses.

Exact Smith computations verify consequences of the displayed matrices.  They
cannot verify that those matrices are the correct maps for the constructed
analytic space.

## 4. Exact finite audit

The companion script uses exact `ZZ` arithmetic and no floating point.  Its
verified universe and conclusions are:

| object | status | exact conclusion |
|---|---|---|
| displayed `T1,T2,T0` | **FINITE-EXACT** | determinants one; orders 3 and 4; `T0=I+N`, `N^2=0`, `rank N=2` |
| dual monodromy | **FINITE-EXACT** | matrices (2) and `A1*A2*M0=I` |
| cusp invariants | **FINITE-EXACT** | invariant ranks on exterior degrees `0..4` are `(1,2,4,2,1)` and the images are saturated |
| lattice coinvariants | **FINITE-EXACT** | `SNF([A1-I | A2-I])=(1,1,1,0)`, hence one free displayed coinvariant |
| degree-two coinvariants | **FINITE-EXACT** | Smith diagonal `(1,1,1,1,1,0)`; primitive covector `(0,0,1,6,0,0)` |
| cusp pushout matrix | **FINITE-EXACT** | rank 4, Smith diagonal `(1,1,1,1,0,0)` |
| oriented conductor quotient | **FINITE-EXACT GIVEN THE STATED QUOTIENT** | `H_0..H_4(W)` ranks `(1,2,4,2,1)`, all free; `chi(W)=2` |
| one reversed branch | **FINITE-EXACT HOSTILE** | rank rises to 5 and Smith diagonal becomes `(1,1,1,1,2,0)` |
| special-surface map `alpha2` | **FINITE-EXACT** | rank 3, Smith diagonal `(1,1,1,0)`, primitive left annihilator `(4,2,3,2)` |
| displayed twists | **FINITE-EXACT** | `A1*v1=v1`, `A2*v2=v2`, `(ell1,ell2)=(1,-1)` |
| chosen clutch presentation | **FINITE-EXACT** | Smith diagonal `(1,1)` for `p=-1` |
| hostile clutch control | **FINITE-EXACT** | Smith diagonal `(1,7)` for `(0,1,1)` and `p=-7` |
| twist gauge quotient | **FINITE-EXACT** | `Z^3/ker(p)=Z` via signed `p`; admissible image is the units modulo 12 |
| simultaneous centralizer | **FINITE-EXACT** | `C_GL4(Z)(T1,T2)={+/- (I+bE14)}` |
| marked elliptic fillings | **FINITE-EXACT GIVEN THE DISPLAYED AFFINE ACTIONS** | order 3 permits every integral centralizer shift; order 4 permits exactly even shifts |

For coprime integers `a,b`, the same two-generator presentation is

```text
R(a,b;ell)= [[a,-ell1],[b,ell2-b ell0]],
p=ab ell0-b ell1-a ell2,
det R=-p.                                                (7)
```

The script checks 54,320 admissible tuples with `2<=a,b<=15`,
`|ell0|<=3`, `|ell1|,|ell2|<=6`.  Every checked full-rank presentation is
cyclic of order `|p|`.  The finite experiment is only a hostile control; the
general statement follows directly when
`gcd(ell1,a)=gcd(ell2,b)=1`: the entries of `R` have gcd one and
`gcd(p,ab)=1`, so its first Smith divisor is one and the second is `|p|`.

The native response increments are

```text
Delta_ell0 p=ab,   Delta_ell1 p=-b,   Delta_ell2 p=-a. (8)
```

For coprime `(a,b)`, the exact kernel is

```text
ker(phi)=Z(1,a,0) direct_sum Z(1,0,b).
```

Indeed, reducing `ab d0=b d1+a d2` modulo `a` and `b` gives
`d1=a r`, `d2=b s`, and then `d0=r+s`.  Thus `p` is noninjective on raw
labelled triples.  This does **not** show that it loses an analytic gluing
parameter.  Under the manuscript's Seifert coordinates
`e0=ell0`, `beta1=-ell1`, `beta2=-ell2`, the two kernel generators are exactly
the standard presentation changes

```text
(e0,beta1,beta2) -> (e0+1,beta1-a,beta2),
(e0,beta1,beta2) -> (e0+1,beta1,beta2-b).
```

Consequently `p` is complete for the displayed Seifert data modulo these
standard gauge moves.  Whether the corresponding holomorphic fillings are
also gauge-equivalent is **OPEN**.  The analytic gluing lattice must first be
quotiented by global section/fibre-translation gauge; no kernel representative
may be promoted to a moduli sidecar before that calculation.

For `(a,b)=(3,4)` this quotient is exact without a finite scan.  If

```text
k1=(1,3,0),        k2=(1,0,4),        q=(0,-1,1),
```

then `det[k1 k2 q]=1`, `p(k1)=p(k2)=0`, and `p(q)=1`.  Hence
`Z^3/<k1,k2>~=Z` via signed `p`.  Under the admissibility conditions
`3 does not divide ell1` and `ell2` odd, `p` runs exactly through the units
modulo `12`, with decoder

```text
(ell1 mod 3,ell2 mod 4) : (1,1)->5, (1,3)->11,
                          (2,1)->1, (2,3)->7.           (8c)
```

This is completeness for oriented labelled Seifert presentation data.  The
displayed fundamental group and homology see only `|p|`; no theorem here says
that `|p|` classifies six-manifolds or analytic threefolds.

Assuming only the manuscript's stated orientation-preserving identification
of opposite boundary curves, the conductor pushout itself can be audited
without degeneration theory.  The normalization `dP6`, its six-curve hexagon,
and the quotient double locus have Euler characteristics `6,6,2`, so
`chi(W)=2`.  The degree-two pushout map has rank four and saturated Smith
factors.  Its long exact sequence gives

```text
H_k(W;Z)=(Z,Z^2,Z^4,Z^2,Z),       k=0,...,4,            (8a)
pi1(W)=Z^2.                                                (8b)
```

The fundamental-group statement uses the actual hexagon attaching word, a
commutator in the free group on two double-locus loops, not merely its zero
abelianization.  Reversing the relative orientation of one paired branch
changes the Smith response by a `Z/2` factor while leaving cell counts and
Euler characteristic unchanged.  Thus branch degrees are a load-bearing
sidecar; this computation still does not verify that the analytic fan quotient
realizes the stated oriented topological quotient.

## 5. Complex-geometric audit ledger

The following items survived a first full reading with **NO IMMEDIATE
CONTRADICTION**, which is not proof:

1. the explicit period-function construction and lattice action in Sections
   2--3;
2. the claimed free/proper action on a sufficiently small fan neighborhood;
3. the finite cyclic actions used for the order-3 and order-4 logarithmic
   transforms;
4. the collar gluing description and the dependence of the displayed
   topological presentation on the three twist coordinates;
5. the local ambient differential whose fibre-differential image becomes
   conductor-supported torsion in Section 10.

The priority **OPEN REFEREE OBLIGATIONS** are:

1. audit the half-weight modular square root, its automorphy multiplier, and
   the claimed `O(-1)` torsor class;
2. rebuild the infinite-fan analytification and proper-discontinuity estimates
   independently rather than only following the manuscript;
3. verify that the claimed map `X -> P1` is a global holomorphic algebraic
   reduction, then audit the multiple-fibre invariant module and the no-carry
   step used to deduce algebraic dimension one;
4. derive cusp specialization, cup products, and integral nearby cycles from
   a model which retains all branches of the nonnormal fibre;
5. check relative duality and base change in the final cohomological argument;
6. trace the generic-fibre and filling torus translations through the three
   clutches, recording exactly which symmetries fail to globalize; and
7. run a specialist topology audit of the collapse and recognition maps.

Until these gates close, statements about compactness, smoothness, homology,
or a complex structure on `S6` remain **MANUSCRIPT CLAIM / UNDER AUDIT**.

## 6. The CDP20 conflict: what is and is not established

Campana--Demailly--Peternell 2020 propagates a twisted-differential vanishing
through singular fibres.  In the published argument, Proposition 2.4 reduces
the global vanishing problem to all fibres, and the conormal discussion treats
vanishing of torsion-free fibre differentials as sufficient for the ambient
restriction.  The S6 manuscript says an ambient cotangent section can remain
nonzero while its image in fibre differentials is torsion supported on the
conductor and hence dies after the torsion-free quotient.

[THM-3955](../../01-canon/theorems/THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion.md)
proves the exact local algebra independently.  For

```text
B=k[x,y]/(xy),
E=B dx direct_sum B dy,
T=torsion Omega_(B/k),
K=ker(E -> Omega_(B/k)/T),
```

one has

```text
T=k*(x dy)=k*(-y dx),
Omega_(B/k)/T ~= Omega_(Btilde/k),
0 -> B*d(xy) -> K -> T -> 0,
K/B*d(xy) ~= k.                                        (9)
```

This is characteristic-free.  It establishes a genuine local implication
failure: conormal vanishing plus torsion-free vanishing need not imply ambient
vanishing.  It does **not** establish that the class globalizes on the claimed
threefold, that every CDP20 hypothesis is met, or that the published theorem's
conclusion is false.  Those are separate global checks.

The typing is important: the map
`Omega_B/T -> Omega_Btilde` is an isomorphism.  The extra conductor kernel is
in the preceding ambient map `E -> Omega_B/T`, not between torsion-free and
normalized differentials.

The manuscript's fibre also has two claimed triple points, so the node model
is not the end of the local audit.  [THM-3957](../../01-canon/theorems/THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel.md)
proves for `B=k[x,y,z]/(xyz)` that

```text
T ~= Btilde/B,                  Ann(T)=(xy,xz,yz),
0 -> T -> Omega_B -> Omega_Btilde
  -> k[x]dx direct_sum k[y]dy direct_sum k[z]dz -> 0.   (9a)
```

Here torsion has generic rank one on each double axis and fibre dimension two
at the triple point.  Unlike the node, `Omega_B/T -> Omega_Btilde` is not
surjective; its cokernel is a full differential module on each conductor axis.
This stronger local correction is characteristic-free and still makes no
globalization claim.

The manuscript also challenges a first-homology count inherited from the 1996
paper by inserting nontrivial monodromy coinvariants.  The displayed lattice
coinvariant calculation is exact, but its global use remains coupled to the
geometric attaching-map audit.

## 7. Repo connections and strict firewalls

### Conductor three-site ledger -- productive

The exact connection to
[MISTAKE-466](../../01-canon/MISTAKES.md) and
[THM-3944](../../01-canon/theorems/THM-3944-repeated-factor-double-torus-one-place-square-conductor-collapse.md)
is the operation

```text
original nonnormal regular locus
  -> normalization minus conductor preimage
  -> full normalization.                                (10)
```

Each arrow can erase data.  The needed sidecar is the conductor ideal, branch
labels, restriction/extension maps, and the actual consumer.  THM-3944 carries
a `mu_3` Kummer class; THM-3955 carries point-supported node torsion; THM-3957
adds axis-supported torsion and a normalization cokernel.  Their payloads are
not interchangeable.

### LRC affine local systems -- precise program, not a proof

The manuscript's lawful reusable operation is "monodromy plus translations,
then coinvariants."  Applied to the existing LRC carrier, the target program
is to place the sheet transitions of
[THM-773](../../01-canon/theorems/THM-773-prime-seven-sheet-monodromy-and-tournament-fibre.md)
over the labelled cell complex of
[THM-2047](../../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md),
quotient only diagonal translation, and retain owner, metric phase, and legal
next-event data as the stalk.  The tournament is a coarse base only.

This is not a triangle-group or modular-cusp transfer; MISTAKE-233 forbids
that analogy.  No LRC(14) consequence is claimed.

### Planar Jacobian boundary incidence -- methodological

[THM-3950](../../01-canon/theorems/THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow.md),
[THM-3951](../../01-canon/theorems/THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry.md)
and the incoming audited
[THM-3952](../../01-canon/theorems/THM-3952-minimal-mobius-internal-split-carriers-are-four-critical-colors-and-nonentry.md),
[THM-3953](../../01-canon/theorems/THM-3953-rationally-split-hidden-cubic-ramification-triangle-nonentry.md), and
[THM-3954](../../01-canon/theorems/THM-3954-extra-common-debt-creates-A-singularities-and-nonunibranch-residual-boundary.md)
show that normalization must retain a boundary-incidence multigraph before an
affine-plane conclusion is drawn.  The later results add four critical colors,
a forbidden boundary triangle, and a normal `A`-singularity whose residual
prime has two normalization addresses.  These obstructions are specific to an
`A2` open surface.  The S6 cusp fibre instead has a cyclic/self-identified
conductor incidence pattern.  The shared rule is "normalization plus incidence
sidecar"; there is no direct Jacobian theorem transfer.
[THM-3956](../../01-canon/theorems/THM-3956-split-hidden-cubic-integrality-and-repeated-root-trichotomy.md)
is now **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED** and closes
the split hidden-cubic lane by integrality and a repeated-root trichotomy.  It
reinforces the boundary/predicate discipline but supplies no Hopf theorem.
[THM-3958](../../01-canon/theorems/THM-3958-one-hidden-root-principal-different-and-pure-power-boundary.md)
is also **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**: its
one-hidden-root lane closes through a principal-different/forbidden-unit gate.
No claim from the live THM-3959 or THM-3960 namespaces is used here; consult
their frontmatter rather than caching their rapidly changing status.

### Explicitly rejected vocabulary bridges

- order `3,4` monodromy has no map to Mahler's base `3/2` orbit;
- the determinant `p` has no additive `A+B=C`, radical, or height, hence no
  ABC consequence;
- tori, logarithmic transforms, and monodromy do not instantiate IUT's typed
  copies or disputed height comparison;
- Hopf lattice data has no intrinsic binary observable, so forcing a
  tournament would be cosmetic;
- unstable homotopy groups of spheres do not by themselves test integrability
  of an almost-complex structure.

## 8. Reproduction and promotion rule

From the repository root:

```text
python 04-computation/hopf_s6_triangle_monodromy_snf_audit_20260824.py
python -O 04-computation/hopf_s6_triangle_monodromy_snf_audit_20260824.py
```

The output should match
[the frozen result](../results/hopf_s6_triangle_monodromy_snf_audit_20260824.out).
The LF-byte SHA-256 values at promotion are

```text
script  588cc30430474b088b27b3a124a2d4257786a1a38ed656b831839e5753cf203d
output  a7e8e94071a58825fde77dcaa4298614c8499e18d1f5aeabd884976efae443bd
```

Before any promotion of the manuscript's main conclusion, rerun the exact
audit, close the six nonfinite gates in Section 5, and obtain independent
specialist verification of the analytic construction and topology.

## 9. Exact extensions found by the 2026-08-24 audit

### The irreducible one-Euler-fibre grammar isolates torus dimension two

[THM-3991](../../01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md)
proves independently of the manuscript that a rank-`n`, lattice-periodic,
unimodular toric central fibre divided by an index-`d` translation sublattice
has

```text
chi(W)=d*n!,                 number of components=d.       (15)
```

Only top-simplex orbits contribute Euler characteristic, and their number is
forced by normalized volume.  If this is the only nonzero-Euler fibre of a
compact torus fibration with homology `S^(2n+2)`, then `d*n!=2`; the solutions
are `(n,d)=(1,2),(2,1)`, and irreducibility `d=1` forces `(n,d)=(2,1)`.
Thus no same-form example with `n>=3` satisfies the sphere Euler condition.
Extra singular
fibres, fixed-point quotients, singular cones, or a different degeneration
grammar are the exact escape routes.

### The correct discrete moduli problem is a gauge quotient

Equation `(8)` must be interpreted after quotienting the twist lattice by
global fibre translations and changes of local section.  The immediate
calculation is a Cech `H^1` problem for the sheaf of torus sections, not the
selection of a preferred representative in `ker(phi)`.  Until that quotient
is computed, “same `p` gives distinct complex threefolds” is unsupported.

### Narrow the representation search through a rational-elliptic spine

Remark 3.12 of the manuscript identifies the `tau`-spine with a rational
elliptic surface whose three Kodaira fibres are

```text
IV*, III, I1,                 8+3+1=12.                   (16)
```

This suggests a sharper search than enumerating arbitrary coprime `(a,b)`
pairs and rank-four matrices.  First classify three-fibre rational elliptic
surfaces with the desired elliptic orders and cusp width; then derive the
rank-four extension and Abel--Jacobi coordinate from Gauss--Manin data.  For
orders `3,4` with a multiplicative cusp, the elementary Euler filter also
permits `IV+III+I5`; the primitive width-one cusp singles out the displayed
`IV*+III+I1` spine.

### The integral centralizer gives a marked local filling filter

The frozen exact companion computes the simultaneous integral centralizer

```text
C_GL4(Z)(T1,T2)={+/- (I+b E14): b in Z}.                 (17)
```

The element `I+b E14` shifts the manuscript's `beta` coordinate, and hence
the punctured-family constant `c0`, by an integer.  On the dual lattice put
`Q_b=I+bE41`.  The displayed period identity is

```text
Pi_(c0+b)=Pi_(c0) Q_b,          Q_b A_j=A_j Q_b.         (18)
```

Thus `c0 mod Z` is the exact unipotent-centralizer orbit for the marked
punctured group family.  It is **not** the completed-filling answer.

The elliptic affine actions supply a new exact local filter.  With primitive
invariant vectors `epsilon_j,delta` and invariant covectors `psi_j`, the
displayed lattices satisfy

```text
Q_b(ell_j epsilon_j+c_j delta)-(ell_j epsilon_j+c_j delta)
  =b ell_j delta,
psi_1(delta)=3=m_1,             psi_2(delta)=2, m_2=4.  (19)
```

Translation conjugacy at order three is therefore automatic for every
integer `b`; at order four it holds exactly when `b` is even because `ell2`
is odd.  The negative centralizer sign cannot extend: the invariant covector
`gamma` sees `-2ell_j`, which is not divisible by `m_j`.  Consequently every
base- and fibre-preserving isomorphism arising from this marked centralizer
must satisfy

```text
c0'-c0 in 2Z.
```

For a general order-four translation coefficient `c2`, the first filled
marked address is

```text
kappa=[c0-c2] in C/(2Z);                               (20)
```

the manuscript's displayed `c2=0` gives `c0 mod 2Z`.  A one-`delta` change in
the order-four vector is an exact hostile: it preserves `ell2` and `p` but
changes the invariant-covector residue by `2 mod 4`.

### The even generator has zero completed Cech class

**Status: VERIFIED-EXACT FROM DISPLAYED DATA + PROVED CONDITIONAL
IMPLICATION.**  Assume the manuscript's analytic pieces and transition maps
exist as displayed; assume their compatibility with the punctured torus
exponential sequence `0 -> Lambda -> E -> T -> 0`; and assume that
zero-winding differences of the three filled local isomorphisms are exactly
the Cech cocycles of the complete logarithm sheaf identified in Proposition
9.23 as the completed vertical-translation bundle

```text
0 -> O -> V -> O(-1) -> 0.                              (21)
```

Under those hypotheses the marked centralizer shift `b=2` extends across all
three fillings.  More generally every `b=2k` extends, and the
centralizer-generated completed orbit subgroup is exactly `2Z` within the
prescribed marked affine class.  This is a conditional classification of that
marked centralizer action, not of the threefolds under arbitrary
biholomorphisms.

First, the cusp contributes no hidden analytic modulus.  Replacing `c0` by
`c0+b` for an integer `b` gives

```text
Pi_(c0+b)=Pi_(c0) Q_b,
C_(c0+b)(t)=C_(c0)(t)+b E21.                            (22)
```

The integral automorphism `Q_b` fixes `Lambda_tor` pointwise and acts as the
identity on `Lambda/Lambda_tor`.  For
`lambda_bar=(a,d) in Z^2`, the second change in `(22)` sends it to the integral
vector `(0,ba)`, so its componentwise exponential is one.  Consequently every
`c_lambda(t)`, every `Psi_lambda`, the cusp quotient `N0`, the exponential
overlap map `E0`, and `G0` are literally unchanged for integral `b`.

For `b=2k`, rewrite both elliptic fillings through this common punctured
family.  With basis `(gamma_hat,u_hat,w_hat,delta_hat)`, the two affine
numerators change by

```text
d1= 2k delta_hat,                d2=-2k delta_hat.
```

Translations by

```text
r1=k(1/12,-1/2,0,0),            r2=k(1/6,0,0,0)
```

give exact conjugacies because

```text
d1/3=(I-A1)r1+k u_hat,
d2/4=(I-A2)r2-k w_hat.                              (23)
```

Their translations on the punctured overlaps have lattice cocycles

```text
kappa1=-k u_hat,                 kappa2=k w_hat.
```

Since `A1 A2 A0=I`, the triangle relation

```text
kappa1+A1 kappa2+A1 A2 kappa0=0
```

forces `kappa0=k w_hat`.  The manuscript's toric shear supplies the missing
cusp map explicitly: `B0(-k u_bar)=-k w_hat`, so, because
`s o g0=s-1`, `Phi_(-k u_bar)` has boundary lift `-k s w_hat` and cocycle
`k w_hat`.  It commutes with the `Psi_lambda` action and descends across
`N0`.  The three `kappa_j` define a class of the period local system on
`B^circ=P1 minus {p0,p1,p2}`.  The universal-cover vector bundle `E` of the
torus family is coherent on this Stein curve, so Cartan B gives
`H^1(B^circ,E)=0`; the exponential sequence therefore realizes that period
class by a global translation of the punctured family.  Subtracting it from
the three displayed local maps leaves zero winding at every overlap.  Thus
the full **discrete** overlap obstruction is zero.

After these winding classes are removed, the overlap translations admit
holomorphic logarithms in `V`.  Equation `(21)` and

```text
H^1(P1,O)=H^1(P1,O(-1))=0
```

give `H^1(P1,V)=0`.  The additive Cech cocycle is therefore a coboundary;
exponentiating the commuting vertical fields patches the three local maps to
a base- and fibre-preserving isomorphism `X_c0 -> X_(c0+2k)`.  Conversely the
order-four invariant covector gives residue `2b ell2 mod 4`, which is zero for
even `b` and is `2 mod 4` for odd `b` because `ell2` is odd.  Hence no odd
marked centralizer shift in this prescribed affine class extends.  This
completes the promised generator
test and makes `(20)` the exact orbit address for the marked completed
centralizer action.

The conclusion above remains conditional on the displayed analytic
construction and `(21)`.  It verifies no compactness, topology, homology,
diffeomorphism, or complex structure on `S6`.

The exact companion/frozen-output hashes are

```text
script  c093b38c979bfae763a8b76531315c9e44b4d542a163e8c0af5fafd927dd0f10
output  d594fe56f10b815c5d6799bb7f1c444dfb114d90a1c16e779b683e79a1c733c1
```

The completed-overlap companion/frozen-output hashes are

```text
script  4e8c731cfb1d1eeeb535eac45782dad574ebbe92b22507a30061cdabf095c7b8
output  4542332b4d21e3ee9b2a57dbfa14d557d7179d7352c154fa151b8ef67cfd8c5a
```
