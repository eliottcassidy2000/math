---
id: THM-2695
title: "Secondary Kummer Bockstein, Picard divisibility spectrum, and prime-alignment boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent immutable promotion review
  pending.  For an irreducible equivariant squareclass plane on a normal
  complex variety, divisor parity and the mu4-to-mu2 coefficient Bockstein
  give an exact three-level all-or-nothing spectrum: ramified, quasi-etale
  but mu4-nonliftable, or mu4-liftable.  More generally the secondary
  obstruction image is Pic(U)[ell]/ell Pic(U)[ell^2].  THM-2681, the toric
  d^2=abc normalization, and a C3-action on G_m^2 realize the three levels.
  The THM-2657/2688 class 7 is the same central-extension lifting mechanism
  after a 13-primary pushout, but no nonzero 2-primary/13-primary coefficient
  map exists.  The Q8 class x^2+xy+y^2 is a further binary-A4 lift
  obstruction, not a Keller obstruction without a separate spin/Q8 sidecar.
source: root-2026-07-28-secondary-kummer-bockstein
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2686-coprime-galois-invariant-kummer-vanishing-and-a4-standard-plane-gate
related:
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2688-affine-facet-normalization-vertical-sphere-diagonal-simplex-and-lens-quotient
script: 04-computation/jacobian_secondary_kummer_bockstein_thm2695.py
output: 05-knowledge/results/jacobian_secondary_kummer_bockstein_thm2695.out
script_sha256: 7f0f3381036be3123ba0c50cc626ecbe09db24817ecc9dd04398ebdd4c02ebce
output_sha256: ad1c5a1f59a42d0baaf2eaf8f30aeda79afa215837155b318b51d5d9a28cc3b6
hash_basis: LF-normalized bytes
---

# THM-2695 -- the secondary Kummer obstruction is Picard divisibility

**PROVED CANDIDATE + VERIFIED-EXACT; independent immutable promotion review
pending.**

THM-2685 gives the codimension-one parity test for extending a generic
Kummer packet quasi-etale.  A surviving packet has a second, different
question: does its `mu_2` torsor lift through `mu_4`?  The answer is measured
exactly by divisibility of its two-torsion Picard class.  The two tests form a
three-level spectrum, not one blended residue.

## 1. Two consecutive equivariant obstruction maps

Let `X` be a normal integral finite-type complex variety, put

```text
U=X_reg,                  K=C(X),
```

and let a finite group `Q` act by complex automorphisms of `X`.  Let `W` be
a `Q`-stable irreducible plane in `K^*/K^{*2}`.  The first map is the
codimension-one Kummer residue

```text
partial:W -> direct_sum_E F_2,
partial_E([f])=v_E(f) mod 2,                              (1)
```

where `Q` permutes the prime divisors `E`.  If `partial|_W=0`, purity over
the regular scheme `U` identifies the generic plane with a plane in
`H^1_et(U,mu_2)`.  The second map is the connecting homomorphism

```text
beta_4:H^1_et(U,mu_2) -> H^2_et(U,mu_2)                  (2)
```

belonging to the exact coefficient sequence

```text
1 -> mu_2 -> mu_4 -> mu_2 -> 1.                          (3)
```

The complex action fixes `mu_2` and `mu_4`, so both `(1)` and `(2)` are
`Q`-equivariant.  Their kernels on the irreducible plane are therefore zero
or the whole plane.  Exactly one of the following three levels occurs:

```text
level 0: ker(partial|_W)=0;
level 1: partial|_W=0 and beta_4|_W is injective;
level 2: partial|_W=0 and beta_4|_W=0.                    (4)
```

There is no intermediate one-dimensional kernel.  This remains true for the
standard `F_2[S3]` plane even though characteristic two divides `|S3|`:
irreducibility, not semisimplicity, is the input.

## 2. Exact Picard-divisibility formula

Write

```text
lambda_2:H^1_et(U,mu_2) -> Pic(U)[2]=Cl(X)[2]            (5)
```

for the Kummer line-bundle class.  Then

```text
beta_4(alpha)=0
iff
lambda_2(alpha) belongs to 2 Pic(U)[4].                  (6)
```

More precisely,

```text
H^1_et(U,mu_2) / image(H^1_et(U,mu_4))
  ~= Pic(U)[2] / 2 Pic(U)[4]
  ~= image(beta_4).                                      (7)
```

The same proof works for every prime `ell` invertible on `U`:

```text
image(beta_(ell^2)) ~= Pic(U)[ell] / ell Pic(U)[ell^2],  (8)
```

where `beta_(ell^2)` belongs to
`1->mu_ell->mu_(ell^2)->mu_ell->1`.

### Proof

Compare the two etale Kummer rows:

```text
0 -> O(U)^*/O(U)^{*4} -> H^1(U,mu_4) -> Pic(U)[4] -> 0
              |                 |               | times 2
0 -> O(U)^*/O(U)^{*2} -> H^1(U,mu_2) -> Pic(U)[2] -> 0.  (9)
```

Normality gives `Pic(U)=Cl(X)`.  The left vertical arrow is surjective: a
unit modulo squares is the image of the same unit modulo fourth powers.  The
long exact sequence of `(3)` says that `beta_4(alpha)=0` exactly when `alpha`
comes from `H^1(U,mu_4)`.

Necessity in `(6)` follows from the right square.  Conversely, suppose
`lambda_2(alpha)=2M` with `M in Pic(U)[4]`.  Lift `M` to a `mu_4` torsor.
Its image differs from `alpha` by a unit squareclass, and that difference
lifts through the surjective left arrow.  This proves `(6)`.  Taking
cokernels in `(9)` and using exactness at `(2)` gives `(7)`.  Replacing two
by `ell` proves `(8)`.

At the generic field the secondary class always vanishes: the map

```text
K^*/K^{*4} -> K^*/K^{*2}
```

is surjective by the same representative.  Thus `beta_4` is not another
valuation of `f`.  It detects failure to choose the fourth-root torsor
globally etale, equivalently the Picard divisibility quotient `(7)`.

Every global-unit squareclass is at level 2.  In the class-group branch of
THM-2685, the standard plane is at level 2 precisely when its projected plane
in `Cl(X)[2]` lies in `2Cl(X)[4]`; otherwise `(2)` embeds another standard
plane in `H^2_et(U,mu_2)`.

## 3. All three levels occur

The same irreducible `C3` standard plane occurs at every level of `(4)`.

### Level 0: ordered-root completion

On THM-2681's ordered-root completion, THM-2685 computes the restored divisor
rows

```text
(1,0),             (0,1),             (1,1).             (10)
```

They have rank two.  The generic unit-chart plane therefore does not extend
quasi-etale.

### Level 1: the sharp toric A4 carrier

For

```text
X=Spec C[a,b,c,d]/(d^2-abc),                             (11)
```

one has

```text
O(X)^*=C^*,                 Cl(X)=(Z/2)^2,               (12)
```

with `C3` cycling the three nonzero classes.  The standard `V4` torsor is
quasi-etale, but

```text
2Cl(X)[4]=0.                                                (13)
```

Thus `beta_4` is injective on its standard plane.  This is a concrete
`A4`-carrier calculation: the sharp quasi-etale hostile from THM-2655 already
carries a nonzero secondary Kummer Bockstein.

### Level 2: a unit plane

On `X=(G_m)^2`, define

```text
sigma(x,y)=(y,(xy)^(-1)).                                (14)
```

This automorphism has order three and acts irreducibly on the unit
squareclasses `[x],[y]`.  Each class lifts through `(3)` by the corresponding
fourth-root torsor, so `beta_4` vanishes on the whole standard plane.

A class-group level-2 control is given by a `C3`-stable `(Z/4)^2`: its
two-torsion standard plane is exactly twice its four-torsion.  The theorem
does not require a geometric realization of this additional control.

## 4. Exact relation to the odometer lens Bockstein

For a principal `G`-cover classified by `a in H^1(B,G)` and a central
extension

```text
0 -> A -> G_tilde -> G -> 0,                             (15)
```

the connecting class `delta(a) in H^2(B,A)` is exactly the obstruction to
lifting the cover to `G_tilde`.  The two applications are

```text
Kummer:   G=C2,   G_tilde=C4,        A=C2,       delta=beta_4;
odometer: G=C13,  G_tilde=C_(13^6),  A=C_(13^5), delta=7. (16)
```

All coefficient actions in `(16)` are trivial.  For the THM-2657 quotient
map `k|->2k mod13`, the lift of `1` is `7`; thirteen turns give kernel
coordinate `7`.  For `C4->C2`, the lift is `1`; two turns give the nonzero
kernel element.

Push the odometer kernel `C_(13^5)` onto `C13`.  The middle group becomes
`C169`, so the coarse nonvanishing in THM-2688 is literally the standard
one-step 13-Bockstein rescaled by

```text
2^(-1)=7 mod13.                                           (17)
```

The remaining four 13-adic kernel digits record physical odometer depth,
not another coarse cohomology class.

There is no nonzero homomorphism between any of the 2-primary and 13-primary
coefficient groups in `(16)`.  Hence no naturality square transports the
number `7` to `(1)` or `(2)`.  The common preserved predicate is only
**failure of a specified coefficient lift**.  The lens construction forgets
divisors, units, class groups, and the quartic standard plane; the Kummer
construction forgets carry digits, physical translations, and the lens
cover.

THM-2686 gives the sharp action boundary.  For `C3` acting on a 2-primary
module, `H^1(C3,-)=H^2(C3,-)=0`, so the action linearizes and
`V4 semidirect C3=A4` is split.  The odometer deck and kernel are both
13-primary, so its extension can remain nonsplit.  The lens clutch therefore
cannot become an `A4` group-extension obstruction.

The other live repo uses of *Bockstein* remain separately typed.  THM-2571
uses the Smith cokernel of a Cayley difference, THM-2580 uses its Hasse-moment
filtration, and THM-2337 uses a divided relation-response sidecar.  They share
an integral-lift failure predicate, but there is no coefficient-sequence map
from those objects to `(2)` or `(16)`.

## 5. Quaternionic refinement and its stopping boundary

Let `x,y` be the two characters of the standard `V4` plane.  The unique
nonzero `C3`-invariant class in `H^2(V4,F_2)` is

```text
q=x^2+xy+y^2=beta_4(x)+xy+beta_4(y).                     (18)
```

Over `C`, the degree-one coefficient Bockstein is cup-square.  The mixed cup
term `xy` is extra data: the individual secondary Bocksteins do not determine
the quaternionic obstruction.  The class `(18)` restricts nontrivially to
all three nonzero cyclic subgroups, so it classifies

```text
1 -> C2 -> Q8 -> V4 -> 1.                                (19)
```

Adjoining the order-three automorphism gives the binary tetrahedral extension
`2.A4`.  Pulling `(18)` to a resolvent regular locus is precisely the
obstruction to lifting the `V4` torsor to `Q8`.

This is not a Keller obstruction without a new geometric sidecar.  Neither a
quartic cover nor the Keller equation currently requires a `Q8`, `2.A4`, or
spin lift.  On `(11)`, the free locus upstairs is affine three-space minus the
three fixed axes and is simply connected.  Its regular quotient has
fundamental group `V4`, so a `Q8` lift would require an impossible section
`V4->Q8`.  Thus nonzero secondary and quaternionic Bocksteins are compatible
with the complete known abstract `A4` carrier.

For a new resolvent normalization, the cheapest exact audit is:

1. compute the divisor-parity matrix on the two Kummer generators;
2. compute `2:Cl(X)[4]->Cl(X)[2]` together with the unit squareclasses;
3. only if an independent binary-lift sidecar exists, compute `(18)`.

Without such a sidecar, step 3 classifies a lift failure but excludes no new
`A4` or `S4` candidate.  No quartic Keller exclusion, degree bound, `JC(2)`,
general Jacobian, Dixmier, or LRC conclusion follows.

## 6. Exact companion and review status

Run

```bash
python 04-computation/jacobian_secondary_kummer_bockstein_thm2695.py
python -O 04-computation/jacobian_secondary_kummer_bockstein_thm2695.py
```

Both LF-normalized transcripts byte-match

```text
05-knowledge/results/jacobian_secondary_kummer_bockstein_thm2695.out.
```

The companion uses optimization-stable explicit guards.  It checks the full
cyclic kernel and generator-lift residues, all wrap cocycles, the `C169`/`C4`
pushouts, cross-prime Hom boundary, the `C3`-fixed part of
`H^2(V4,F_2)`, all three `Q8` line restrictions, the THM-2681 parity rank,
the exponent-two and four-divisible Picard controls, and the standard-plane
orbit.

The abstract proof, three controls, and Q8/no-transfer boundaries have passed
an independent hostile audit.  That audit also required replacing a
tautological cyclic check by the present kernel/image and complete wrap-class
checks.  The canonical optimization-stable artifact and frontmatter remain a
**PROVED CANDIDATE** until the main-agent immutable review and status
promotion.  The declared hashes use LF-normalized bytes.

QED (candidate).
