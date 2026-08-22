---
id: THM-3175
title: "First-Witt Hensel wedge and infinitesimal Pluecker covariance"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A lift-independent first-Witt determinant obstructs a common lift of two
  residual polynomial branches and transforms by the determinant of every
  polynomial two-row frame.  Applied to the fixed factorial residuals, it
  separates all ten eligible resultant-exception primes through offset six,
  without neighboring-prime transport.
audit: >
  Two independent derivations proved lift independence, the first-order
  Hensel criterion, polynomial-frame covariance, and the determinant-wall
  boundary.  An independent immutable audit also rederived the slope-zero
  integrality argument and the sharp zero-derivative hostile.  Independent
  O(p) recurrence-and-automatic-differentiation runs reproduced every
  displayed exceptional-prime jet and wedge.  Fresh normal and optimized
  replays byte-match the stored transcript and declared LF hashes.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
related:
  - THM-2580-hasse-bockstein-carry-tower-and-salem-local-unit-boundary
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
  - THM-3153-four-step-prime-resonance-second-euclidean-newton-separation
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3170-five-step-prime-resonance-euclidean-newton-holotopy
script: 04-computation/factorial_first_witt_hensel_pluecker_thm3175.py
output: 05-knowledge/results/factorial_first_witt_hensel_pluecker_thm3175.out
script_sha256: 50d08b2514811a256b5f4a9695e0c9f45ed9ccef6057434c4d0883857eb8b88e
output_sha256: 910a636c6f1f57b1bb97a8d48d8d01e5fa0d788c3d8436b40b8bc03e045b141b
hash_basis: LF-normalized bytes
---

# THM-3175 -- first-Witt Hensel wedge and infinitesimal Pluecker covariance

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem isolates a reusable first-order replacement for a failed
residual resultant.  It is stated over a DVR because the mechanism is not
special to factorial moments.  The name *first-Witt* refers only to the first
two residue layers; no general Witt-vector theorem is assumed.

## 1. The first-Witt wedge

Let `O` be a Henselian DVR with uniformizer `pi`, fraction field `K`, and
residue field `k`.  Let `F,G in O[v]`, let `a in k` be a common root of their
reductions, and choose any lift `atilde in O`.  Define

```text
f0=F(atilde)/pi mod pi,             g0=G(atilde)/pi mod pi,
f1=Fbar'(a),                        g1=Gbar'(a),

omega_a(F,G)=f0 g1-g0 f1 in k.                              (1)
```

The divisions are lawful because `a` is a common residual root.

Then:

1. `omega_a(F,G)` is independent of `atilde`.
2. For `X=(F,G)^t`, `M(v) in Mat_2(O[v])`, and `Xnew=M(v)X`,

   ```text
   omega_a(Xnew)=det(Mbar(a)) omega_a(X).                    (2)
   ```

3. If `(f1,g1)!=(0,0)`, the pair has a common lift modulo `pi^2`
   reducing to `a` if and only if `omega_a(F,G)=0`.
4. If both derivatives are units, the two individual first Hensel digits are

   ```text
   z_F=-f0/f1,                       z_G=-g0/g1,
   omega_a(F,G)=f1 g1(z_G-z_F).                         (3)
   ```

5. If the residual gcd is exactly the simple factor `v-a` and `omega_a` is
   nonzero, then `F,G` have no common factor belonging to their slope-zero
   Newton segments over `K`.

Thus `(1)` is the missing infinitesimal sidecar at a height-zero collision.

## 2. Proof of the invariant and lift criterion

Replace `atilde` by `atilde+pi h`.  Taylor expansion modulo `pi^2` sends

```text
(f0,g0) -> (f0,g0)+h(f1,g1).                                (4)
```

Taking the determinant with `(f1,g1)` proves lift independence.

A prospective common lift has the form `atilde+pi z` modulo `pi^2`.  Its two
equations are

```text
f0+z f1=0,                         g0+z g1=0.                (5)
```

If the derivative vector is nonzero, this two-by-one linear system is
solvable exactly when its determinant `(1)` vanishes.  When both derivatives
are units, solving each equation separately gives `(3)`.  The nonzero-
derivative hypothesis is sharp: if `f1=g1=0`, then `omega=0` identically,
while `(5)` is solvable only when `f0=g0=0`.

For the common-factor assertion, pass if necessary to a finite extension of
the completed DVR.  A monic common factor belonging to the slope-zero
segment remains monic and degree-preserving on reduction.  Since its
reduction divides the simple gcd `v-a`, it is linear and supplies a common
first lift of `a`.  Equation `(5)` would force `omega=0`, a contradiction.

## 3. Polynomial-frame covariance

Write

```text
M(v)=[[m11,m12],[m21,m22]].                                 (6)
```

At the common residual root, the value-over-`pi` column transforms by
`Mbar(a)`.  Differentiating `MX` introduces the extra column `M'X`, but
`X(a)=0 mod pi`, so this column vanishes in the residue field.  The derivative
column therefore transforms by the same `Mbar(a)`.  Taking the determinant
of the two transformed columns proves `(2)`.

If `det(Mbar(a))` is a unit, wedge nonvanishing is an invariant of the entire
Euclidean frame.  At a determinant wall only the forward identity survives;
the converse is false.  Multiplying one row by `pi` kills the transformed
wedge even when the original wedge is nonzero.  This is the exact boundary
between a flat two-row holotopy and a lossy projected chart.

## 4. Factorial-moment consumers

For the factorial functional `L(t^j)=j!`, put

```text
d=p+s,
A=A_(p+s-2),                       B=A_(p+s-1),
A_n(v)=L((d-t+v t^2)^n).                                  (7)
```

THM-3148 gives the complete height-zero residual pair

```text
Abar=s F_(s-2,s),                 Bbar=s F_(s-1,s),         (8)
```

and the fixed resultant

```text
delta_s=Res(F_(s-2,s),F_(s-1,s)).                          (9)
```

For `s=3,4,5,6`, THM-3148's untruncated-degree condition is
`p>2(s-1)`.  The eligible prime divisors of `(9)` therefore give exactly ten
residual exceptions (in particular, the factor `p=7` of `delta_6` is below
the eligibility boundary).  Every residual gcd is simple linear.  Values
and derivatives modulo `p^2` are computed in `O(p)` ring operations from the
division-free THM-3124 recurrence

```text
M_(n+1)=d^n(d-n-1)
       +2(n+1)(2n+1)v M_n
       +n(n+1)(1-4dv)M_(n-1),                              (10)
```

and its formal derivative, beginning with `M_0=1` and `M_1=d-1+2v`.

The following rows display

```text
(s,p,a; A(a)/p,B(a)/p,A'(a),B'(a); omega,z_A,z_B) mod p.   (11)
```

```text
(3,29,28; 19,24,6,1; 20,21,5),

(4,7,5; 3,5,3,3; 1,6,3),
(4,4547,3243; 2552,3588,4280,1741; 3739,3671,3132),

(5,11,2; 10,6,4,4; 5,3,4),
(5,20747,5645; 6993,17070,8474,18464; 7075,13673,8568),
(5,249721,222768;
 129614,136893,167574,84554; 39549,155348,15040),

(6,139,113; 16,34,32,28; 55,69,108),
(6,3767,1641; 1566,2007,2749,2460; 131,2466,1042),
(6,12041,9004; 6042,5530,1522,7790; 10951,1515,463),
(6,807241,22489;
 341729,683313,364516,409175; 439007,414823,557512).       (12)
```

Every displayed derivative is a unit and every displayed wedge is nonzero.
Equations `(3),(5)` therefore show that the two simple residual branches have
different first Hensel digits in all ten cases.  Hence every eligible
resultant exception through offset six has no common slope-zero lift.

For offsets four and five this is an intrinsic second proof of the large
zero-face exclusions in THM-3153 and THM-3170.  Their neighboring-prime
transports remain valid independent controls, but are no longer necessary
for the zero face.  The offset-six rows are the zero-face input to THM-3176;
this theorem alone says nothing about its positive Newton slopes.

## 5. Exact evidence

Run

```text
python 04-computation/factorial_first_witt_hensel_pluecker_thm3175.py
python -O 04-computation/factorial_first_witt_hensel_pluecker_thm3175.py
```

and compare byte-for-byte with the declared stored output.  The companion
uses integer and modular arithmetic only.  It independently recomputes
`delta_3,...,delta_6` by Bareiss resultants and their exact factorizations;
constructs every fixed residual gcd; checks all ten roots, first jets,
derivative scalings, Hensel digits, and nonzero wedges; and supplies explicit
lift-independence, polynomial-frame covariance, and determinant-wall
controls.  The largest recurrence has `807246` steps and constant state.  No
floating point, random sampling, or fitted formula is used.

## 6. Connection contract and scope

The source object is the two-column first jet

```text
[(F/pi,G/pi);(F',G')] mod pi.                              (13)
```

Its target is the exterior coordinate `omega`.  A unit polynomial frame
preserves nonvanishing and destroys the individual root label, phase, and
higher Hensel digits.  At determinant walls even nonvanishing can be lost.
When `omega=0` and the derivative vector is nonzero, only a common lift modulo
`pi^2` is proved; an actual common root requires the higher-Witt tower.

THM-2580 and THM-2585 are related Hasse--Bockstein analogues: their
triangular chart laws also preserve a first carry obstruction.  They are not
dependencies.  In particular, `(1)` does not manufacture LRC's missing
semantic endpoint, positive physical observable, or owner alignment.

For Gaussian moments, the lemma offers a lawful way to refine a collided
finite-field face by one arithmetic jet.  It does not prove that an arbitrary
radial support supplies such a prime, a simple residual face, or a useful
unit frame.  No arbitrary-offset induction, arbitrary support, `NC(2)`,
`GMC(2)`, or `LRC(14)` conclusion is asserted here.

**QED.**
