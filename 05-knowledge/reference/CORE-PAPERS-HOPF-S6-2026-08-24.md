# Hopf/S6 manuscript intake and audit ledger -- 2026-08-24

**Status gate.**  The direct PDF
[*The (3,4,infinity) modular family of 2-tori, completed at its three special
points, is a complex structure on S6*](https://alpo.ge/s6.pdf) is treated here
as a **PREPRINT CLAIM / UNDER AUDIT**.  It was downloaded on 2026-08-24 as a
108-page PDF with SHA-256

```text
283bba102dd1d5dc346af81b28145bdaaea6654398d5032e76e97bafb9a858f2.
```

The title is transcribed from page 1.  No peer-review or publication status is
inferred.  This file separates four levels which must not be collapsed:

```text
PROVED IN REPO
FINITE-EXACT FROM DISPLAYED DATA
READ / NO IMMEDIATE CONTRADICTION
PREPRINT CLAIM / OPEN REFEREE OBLIGATION.                    (1)
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
5. **Unconditional local theorem:**
   [THM-3955, node cotangent normalization kernel and conductor
   torsion](../../01-canon/theorems/THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion.md),
   extended at triple crossings by
   [THM-3957](../../01-canon/theorems/THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel.md).

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

For `(a,b)=(3,4)`, the move `(1,3,0)` lies in the kernel of the evaluator.
Consequently `|p|` is a lawful ordinal for one selected response, but is not
a lossless address for the gluing triple or the analytic threefold.

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
5. the local conductor-supported differential mechanism in Section 10.

The priority **OPEN REFEREE OBLIGATIONS** are:

1. audit the half-weight modular square root, its automorphy multiplier, and
   the claimed `O(-1)` torsor class;
2. rebuild the infinite-fan analytification and proper-discontinuity estimates
   independently rather than only following the manuscript;
3. verify the multiple-fibre invariant module and the no-carry step used in
   the algebraic-dimension computation;
4. derive cusp specialization, cup products, and integral nearby cycles from
   a model which retains all branches of the nonnormal fibre;
5. check relative duality and base change in the final cohomological argument;
6. run a specialist topology audit of the collapse and recognition maps.

Until these gates close, statements about compactness, smoothness, homology,
or a complex structure on `S6` remain **PREPRINT CLAIM / UNDER AUDIT**.

## 6. The CDP20 conflict: what is and is not established

Campana--Demailly--Peternell 2020 propagates a twisted-differential vanishing
through singular fibres.  In the published argument, Proposition 2.4 reduces
the global vanishing problem to all fibres, and the conormal discussion treats
vanishing of torsion-free fibre differentials as sufficient for the ambient
restriction.  The S6 manuscript says this passage loses a class supported on
the conductor of its nonnormal cusp fibre.

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
[THM-3960](../../01-canon/theorems/THM-3960-natural-one-parameter-cubic-normal-monogenic-closure.md)
is now independently audited and closes the full natural globally monogenic
one-parameter cubic family. Independently audited
[THM-3959](../../01-canon/theorems/THM-3959-centered-degree-five-root-map-color-closure.md)
closes the centered degree-five color rows, while
[THM-3961](../../01-canon/theorems/THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt.md)
classifies normality in the irreducible arbitrary-`q` monogenic grammar.
Proved THM-3962 closes all coefficient-constant cylinders, including their
conductor debts, and proved THM-3963 closes the moving scalar `P^2` debt.
Proved THM-3964 closes the displayed scalar-coefficient graph-root family;
proved THM-3965 rules out a one-place discriminant in its exact constant-
`(g,h)` unit-ideal family; THM-3966/3968 supply the class/Euler and canonical-
different boundary invoices. Proved THM-3967 closes every irreducible
`deg_P(q)<=2` row, and proved THM-3969 closes the collision-free `Xi in k^*`
affine-`P` graph packet. THM-3970 proves the exact osculation/gcd reframe,
not a normalization closure. THM-3971 proves the all-`m` determinantal
exact-volume/no-Darboux near miss. THM-3972 gives the nonconstant squarefree-
`Xi` blowup normalization and closes `a=t`, constant-`(c,r)`; general
canonical-compatible rows remain open. THM-3973 gives an exact-volume
completion passport and a globally nonmonogenic finite cubic whose
ramification curve meets the affine-plane open, preventing Keller. THM-3974
closes the two-by-three cell and its
transpose at every height, plus the wider height-two seven-piece floor.
THM-3975 adds the two-color/plinth ledger, a finite-free tower of exact rank
`2ell+1` for every `ell>=floor(n/2)`, and the all-height rational no-mate
theorem for every nonconstant `f in k(p)`.
THM-3976 computes the rational-compression pseudoplane intersection and only
its internal degree/support floors. THM-3977 closes the lowest simultaneous
cusp/arm seam by critical points and generic-fibre residues, including both
tested formal corrections. THM-3978 gives the exact linear-seam response
ideals: incompatible invariant constants block global algebraization even
though each completed color has a mate. THM-3979 proves all-order formal cusp
lifting; THM-3981 then proves its canonical quadrature genus-two and
transcendental, including logarithmic nonentry on every exceptional scalar
slice. Alternative gauges, finite Keller entry, and unrestricted Darboux
entry remain open.
None of these supplies a Hopf dependency here.

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
```

The output should match
[the frozen result](../results/hopf_s6_triangle_monodromy_snf_audit_20260824.out).
The companion deliberately rejects `python -O`: its gates are Python
assertions, so optimized execution would erase the audit (MISTAKE-476).
The LF-byte SHA-256 values at promotion are

```text
script  3748afedf8af898b05cd9514970c453637d58ab9a69629a76d94ef14775df6b8
output  a7e8e94071a58825fde77dcaa4298614c8499e18d1f5aeabd884976efae443bd
```

Before any promotion of the manuscript's main conclusion, rerun the exact
audit, close the six nonfinite gates in Section 5, and obtain independent
specialist verification of the analytic construction and topology.
