---
id: THM-3178
title: "Resultant tangent cone, conormal norm, and squarefree first-Witt factorization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Along any degree-r common residual factor, the resultant vanishes to order
  at least r and its first possible coefficient is the norm of one intrinsic
  conormal Pluecker coordinate.  This formula does not require the common
  factor to be squarefree.  On the squarefree locus it factors into the
  branchwise first-Witt Hensel wedges, and polynomial row frames act by their
  determinant in the collision algebra.
audit: >
  Two independent algebraic derivations obtained the tangent sign and norm
  factor.  Exact split, irreducible, and repeated-root controls reproduce the
  formula through common-factor degree three; separate controls verify lift
  independence, branchwise first-Witt recovery, polynomial-frame norm
  scaling, a determinant wall, and the non-full-gcd boundary.  The canonical
  pure-integer companion replays normally and under optimization against its
  stored transcript and declared LF hashes.  The immutable audit rederived
  every sign and hypothesis and repaired the cubic consumer's characteristic-
  three degree-loss boundary before promotion.
source: root/multiscale-newton-flag/2026-08-02
depends_on: []
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2946-full-macaulay-maximal-minor-gcd-and-chart-free-resultant
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3175-first-witt-hensel-wedge-and-infinitesimal-pluecker-covariance
script: 04-computation/resultant_tangent_cone_first_witt_norm_thm3178.py
output: 05-knowledge/results/resultant_tangent_cone_first_witt_norm_thm3178.out
script_sha256: a2b283a0174b9e327dff811ff0da219f16eac27ac8e3272df470da9c8e16a6c8
output_sha256: 81c44bf184a5ee3e781328db1a6d337e840936babf305b3f5f339fcb90364d8e
hash_basis: LF-normalized bytes
---

# THM-3178 -- resultant tangent cone, conormal norm, and squarefree first-Witt factorization

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The first-Witt wedge of THM-3175 is the degree-one case of a larger
determinant-line identity.  The correct object is not a list of chosen roots:
it is one exterior coordinate in the conormal algebra of the resultant
divisor.  This formulation continues across repeated-root seams, while the
individual Hensel branches do not.

## 1. Conormal tangent-cone theorem

Let `O` be a Henselian DVR with uniformizer `pi` and residue field `k`.  Let

```text
F,G in O[v],                     deg F=m, deg G=n,          (1)
```

and assume their leading coefficients are units, so reduction preserves both
degrees.  Suppose that for a monic `H in k[v]` of degree `r>=1`,

```text
Fbar=H f,                        Gbar=H g,
gcd(f,g)=1.                                                   (2)
```

Thus `H` is the full monic residual gcd.  No squarefreeness of `H` and no
coprimality between `H` and `fg` is assumed.  Put

```text
A=k[v]/(H).                                                   (3)
```

Choose any monic lift `Htilde in O[v]` and let `I=(pi,Htilde)`.  The regular
conormal module `I/I^2` is free over `A` on `[pi],[Htilde]`.  There are unique
classes `phi_F,phi_G in A` such that

```text
[F]=phi_F[pi]+f[Htilde],       [G]=phi_G[pi]+g[Htilde]
                                                   in I/I^2.
```

Define the conormal Pluecker coordinate

```text
theta=phi_F g-phi_G f in A.                                  (4)
```

Then `theta` is independent of `Htilde`, and

```text
v_pi(Res(F,G)) >= r,                                         (5)

pi^(-r) Res(F,G) mod pi
 =(-1)^[r(m-r+1)] Res(f,g) Norm_(A/k)(theta).                (6)
```

Here `(6)` means the reduction of the lawful quotient in `(5)`, and the
resultant convention is

```text
Res(P,Q)=lc(P)^deg(Q) product_(P(alpha)=0) Q(alpha).          (7)
```

Since `Res(f,g)` is nonzero and multiplication by an element of a finite
`k`-algebra has nonzero determinant exactly when that element is a unit,

```text
v_pi(Res(F,G))=r                 iff theta is a unit in A.    (8)
```

If `theta` is not a unit, `(6)` vanishes and the resultant has valuation at
least `r+1` (or is zero).  Formula `(6)` itself remains a polynomial identity
even when `gcd(f,g)!=1`; in that case its right side is zero and `(8)` is not
available.  This is the sharp full-gcd boundary.

## 2. Proof of the tangent formula

First, changing the monic lift by

```text
Htilde -> Htilde+pi u
```

sends

```text
phi_F -> phi_F-u f,              phi_G -> phi_G-u g.
```

The two added terms cancel in `(4)`, proving that `theta` is intrinsic.

Exact monic division by `Htilde` gives

```text
F=Htilde ftilde+pi Phi_F,          G=Htilde gtilde+pi Phi_G,
```

where `ftilde,gtilde` reduce to `f,g` and `Phi_F,Phi_G` reduce to
`phi_F,phi_G`.  Introduce an independent parameter `eps` in those two exact
identities.  At `eps=0` the pair has the exact common factor `Htilde`, so its
resultant is divisible by `eps^r`; at `eps=pi` it is the original resultant.
After division by `pi^r` and reduction modulo `pi`, only the coefficient of
`eps^r` survives.  It is therefore enough to calculate the residue
deformation

```text
F_eps=Hf+eps phi_F,               G_eps=Hg+eps phi_G.         (9)
```

This proves, rather than assumes, that the tangent coefficient depends only
on the conormal classes.

Over the universal characteristic-zero coefficient ring, it is enough to
prove the coefficient identity on the Zariski-dense locus where `H`, `f`, and
`g` split into simple, pairwise-disjoint roots.  Both sides are integral
polynomials in all coefficients, so the universal identity extends across the
discriminant and incidence walls and then specializes to every residue field.

Let `a` be a root of `H`.  The nearby root of `F_eps` is

```text
alpha_a=a-eps phi_F(a)/(H'(a)f(a))+O(eps^2),
```

and hence

```text
G_eps(alpha_a)=-eps theta(a)/f(a)+O(eps^2).                 (10)
```

The other roots of `F_eps`, those reducing to roots of `f`, contribute
`Res(f,Hg)` at order zero.  Multiplying `(10)` over the `r` roots of `H`
gives

```text
[eps^r]Res(F_eps,G_eps)
 =(-1)^r Norm_(A/k)(theta) Res(f,Hg)/Res(H,f).               (11)
```

Now

```text
Res(f,Hg)=Res(f,H)Res(f,g),
Res(f,H)=(-1)^[r(m-r)]Res(H,f).                              (12)
```

Equations `(11)--(12)` are exactly `(6)`.  In particular all coefficients
below `eps^r` vanish, proving `(5)`.  The polynomial-continuation argument is
why `(6)` survives when `H` has repeated roots even though `(10)` no longer
has a branchwise meaning.

Equivalently, the Sylvester map at `eps=0` has an `r`-dimensional kernel and
cokernel.  Its first transverse map is multiplication by `theta` on `A`, so
the determinant tangent cone is `Norm(theta)`; `(12)` is the complementary
determinant and its orientation sign.  This is the determinant-line form of
the same proof.

## 3. Squarefree factorization into first-Witt branches

Assume now that `H` is squarefree.  Then `A` is etale and `H'` is a unit in
`A`.  Define

```text
Omega=H' theta in A.                                         (13)
```

Over a residue-field splitting extension, let `a` run through the roots of
`H`.  Lift `a` to a root of `Htilde` in the unramified Hensel extension.  At
that lift,

```text
F/pi=phi_F(a),                   G/pi=phi_G(a),
Fbar'=H'(a)f(a),                 Gbar'=H'(a)g(a).
```

Therefore

```text
Omega(a)
 =(F/pi)Gbar'-(G/pi)Fbar'
 =omega_a(F,G),                                             (14)
```

the lift-independent first-Witt wedge of THM-3175.  Since

```text
Norm_(A/k)(H')=(-1)^[r(r-1)/2] Disc(H),                     (15)
```

formula `(6)` becomes

```text
pi^(-r) Res(F,G) mod pi
 =(-1)^[r+(m-r)r+r(r-1)/2]
   Res(f,g) Disc(H)^(-1) Norm_(A/k)(Omega).                 (16)
```

Thus on the squarefree locus:

```text
v_pi(Res(F,G))=r
 iff every branch wedge omega_a is nonzero
 iff no residual branch has a common first Hensel lift.      (17)
```

The derivative vector at each branch is nonzero because `H'(a)` is nonzero
and `f,g` cannot both vanish at `a`.  Hence the final equivalence follows
from the two-equation first-order Hensel criterion, including branches where
one of the two derivatives vanishes.

Squarefreeness is load-bearing only for `(13)--(17)`, not for `(6)`.  The
exact repeated-cubic control has `theta` a unit and resultant valuation
exactly three, while `Disc(H)=Norm(Omega)=0`.  It is a sharp witness that a
repeated collision cannot be decomposed into ordinary first-Witt root labels.

## 4. Polynomial-frame covariance and determinant walls

Let `X=(F,G)^t` and `M(v) in Mat_2(O[v])`.  Reducing `M` in the collision
algebra gives `M_A in Mat_2(A)`.  Multiplication in `I/I^2` sends both
conormal columns through the same matrix, so

```text
theta(MX)=det(M_A) theta(X),                                (18)

Norm(theta(MX))=Norm(det(M_A)) Norm(theta(X)).              (19)
```

If `det(M_A)` is a unit, exact transverse order is preserved.  Formula `(6)`
may be applied to the transformed pair whenever its leading degrees and
full-gcd hypotheses also remain valid.  If `det(M_A)` is not a unit, the
converse fails: a frame can erase one collision component and kill the norm.
The companion realizes this on the split quadratic algebra with
`det(M_A)=v-1`.  This is the multi-branch version of THM-3175's determinant
wall.

The source, target, and loss are now explicit:

```text
two-row conormal frame  --exterior determinant--> theta in A
                         --finite norm-----------> resultant tangent.
```

The first arrow retains every residual component but forgets the row gauge;
the second forgets which component failed.  A zero norm needs a component
selector or the next conormal/Witt layer.  This is precisely the distinction
between a flat Pluecker holotopy and its lossy scalar projection.

## 5. Exact evidence

Run

```text
python 04-computation/resultant_tangent_cone_first_witt_norm_thm3178.py
python -O 04-computation/resultant_tangent_cone_first_witt_norm_thm3178.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
integer arithmetic and exact rational forward-difference interpolation only.
It computes each deformation resultant at enough integer values to recover
the entire polynomial in `eps`; no fitted degree below the universal
`m+n` bound is assumed.

The exact bank contains common-factor degrees `1,2,3`, split and irreducible
squarefree algebras, and a repeated cubic.  It checks every lower resultant
coefficient, the sign in `(6)`, lift independence of `(4)`, the discriminant
conversion `(15)--(16)`, three branchwise wedges, polynomial-frame covariance
and norm scaling, a one-component determinant wall, and a non-full-gcd
hostile whose first nonzero order rises from one to two.  A separate cubic
`1+2` control recovers THM-2598's exact residual-quadratic discriminant gate.
There is no floating point, random sampling, symbolic black box, or
assertion-sensitive test.

## 6. Cubic `1+2` discriminant consumer

THM-2598's homogeneous quartic boundary supplies the integral cubic

```text
S_A(U)=U^3+2pU^2+(p^2-4Ar)U-Aq^2.                          (20)
```

Take `A` as the DVR uniformizer, assume the residue characteristic is outside
`{2,3}`, and assume `p` is a unit.  (In characteristic three, `S_A'` loses
degree, so the preserved-degree hypothesis of Section 1 fails.)  At `A=0`,

```text
S_0=U(U+p)^2,
gcd(S_0,S_0')=H=U+p,
f=U(U+p),                         g=3U+p.                   (21)
```

For the pair `(S_A,S_A')`, the conormal coordinate at `U=-p` is

```text
theta=omega_(-p)(S_A,S_A')=-2p(4pr-q^2),
Res(f,g)=-2p^2.                                             (22)
```

Since a monic cubic has `Disc(S)=-Res(S,S')`, formula `(6)` gives

```text
Disc(S_A)/A mod A
 =-2p^2 omega_(-p)(S_A,S_A')
 =4p^3(4pr-q^2).                                           (23)
```

Thus the residual-quadratic gate `4pr-q^2` is exactly the first-Witt
transversality coordinate of the double resolvent root, rather than an
isolated discriminant coefficient.  It distinguishes a visible simple
transposition from a higher collision in the resolvent divisor.  This is a
reinterpretation of THM-2598's already-proved coefficient identity; it does
not remove that theorem's `V4`/index or Keller-transfer boundary.

## 7. Factorial consumer and scope

For `r=1`, `H=v-a` and `H'=1`, so `(16)` reduces to

```text
pi^(-1)Res(F,G) mod pi=(-1)^m Res(f,g) omega_a(F,G).         (24)
```

Thus THM-3175's wedge is literally the transverse coordinate of the simple
resultant divisor.  For a future fixed-offset face with several simple common
roots, `(16)` replaces a root-by-root product calculation with one norm in
the etale collision algebra.  Formula `(18)` then transports that certificate
through every unit Euclidean frame without choosing or matching roots.

The unit-leading-coefficient hypothesis is essential to this exact resultant
chart.  In the factorial resonance, the full high-degree polynomials lose
degree modulo the prime, so `(6)` cannot be applied directly to their full
resultant.  One must first extract a monic slope-zero Weierstrass factor (as
in THM-3175); no exact valuation for the unnormalized full factorial resultant
is claimed.

The norm records whether all collision components are transverse, but not
which component carries a physical phase or semantic label.  It supplies no
prime, no simple/squarefree face, and no arbitrary-support observer.  No
arbitrary-offset induction, `NC(2)`, `GMC(2)`, or `LRC(14)` conclusion is
asserted here.

**QED.**
