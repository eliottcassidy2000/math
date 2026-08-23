---
id: THM-3823
title: "Nonlinear cubic etale surface has compact Euler characteristic eight"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over C, the smooth
  integral affine etale surface U=Spec S[A/D,omega/D] of THM-3811 has compactly
  supported Euler characteristic 8.  Its branch normalization is P1 minus
  three points and has Euler characteristic -1; four normalization points are
  identified at the origin and the other two singular branch values are
  unibranch, giving branch Euler characteristic -4.  Fibrewise deletion of
  ramification gives chi_c(U)=3*5-7=8.  Hence no finite etale morphism
  A2_C -> U exists.  This Euler mechanism alone leaves a nonfinite,
  nonproper degree-d atlas with signed Euler sheet debt 8d-1 and at least two
  divisorial boundary primes.  THM-3841 later excludes every dominant plane
  morphism to this surface, while THM-3845 gives an independent Keller-degree
  contradiction.  No general Jacobian claim follows.
source: root / nonlinear-cubic Euler-atlas obstruction lane, 2026-08-23
audit: >
  THREE INDEPENDENT HOSTILE AUDITS PASS (jc_sparse_direct_search, root and
  thm3791-hostile-audit, 2026-08-23).  The first audits reconstructed the
  normalization, three punctures, four origin preimages, two unibranch triple
  values, sheet stratification and signed Euler debt.  The strongest audit
  independently rederived the
  branch incidence, irreducibility, rational inverse, all normalization poles,
  four vertex branches and tangent slopes, complete derivative pullbacks, two
  unibranch triple values, companion collision locus, constructible fibre
  strata, and universal finite-degree Euler obstruction were independently
  rederived.  Good-reduction point counts and two bad-reduction hostiles were
  retained as controls.  The audit caught and repaired a non-load-bearing
  bounded-degree script gate; the canonical scripts now use universal integer
  divisibility.  Primary and independent assertion-free companions have 34
  and 58 active gates; the earlier 16-gate packet remains a hostile control.
  Normal and optimized streams match frozen outputs after LF normalization.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3578-zariski-main-boundary-rank-and-sheet-debt
related:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3845-nonlinear-cubic-keller-atlas-total-degree-contradiction
script: 04-computation/jc2_nonlinear_cubic_etale_surface_euler_thm3823.py
output: 05-knowledge/results/jc2_nonlinear_cubic_etale_surface_euler_thm3823.out
script_sha256: 11922872c37d824cd22da0875bca48044afe7267f92ad2c67be3fd5d4088c815
output_sha256: 1007441d08053e669ca299ca6022c6b2462682ef3fa54e1ebe08e3850fb64b30
semantic_sha256: e2fcbcad4c643b3a0e360789ef5ad169e0f5d43eefd5eb4546426b6df1ac5101
independent_script: 04-computation/jc2_nonlinear_cubic_etale_surface_euler_independent_audit_thm3823.py
independent_output: 05-knowledge/results/jc2_nonlinear_cubic_etale_surface_euler_independent_audit_thm3823.out
independent_script_sha256: 20a5587387ed15a9acef37755d732729415a64e9e55cb45757abe2b9ad60e2ad
independent_output_sha256: 6da6f3e471778376d0b4d9f69c374bf847acaad5512d6eada588df08fdcfbf89
independent_semantic_sha256: 14808b30725ffc7fcc363302dbe76e91e8ba073ed332d8d62802b17f30344b40
compact_script: 04-computation/jc2_nonlinear_cubic_compact_euler_thm3823.py
compact_output: 05-knowledge/results/jc2_nonlinear_cubic_compact_euler_thm3823.out
compact_script_sha256: 713e7a925e4ae42852dc22c5db766510ff0cca74ce6487e224bc883c1659c5e5
compact_output_sha256: ff8f7c8ac1d35b46a4a4f237e3f36c1024fcb9f26a5fa0a782dff16b555e204a
compact_semantic_sha256: ecf52db47db9533b4bf1eec5c23bea5a0e2a91e33a21e6c941863e0cba937cfe
hash_basis: raw LF bytes
---

# THM-3823 -- compact Euler characteristic excludes finite plane atlases

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over the
complex numbers.  Let

```text
Xbar=Spec S -> A2_(A,C),       E=ramification divisor,
U=Xbar minus E=Spec S[A/D,omega/D]                       (1)
```

be the normal finite cubic completion and smooth integral affine etale
surface of THM-3811.  Then

```text
chi_c(U(C))=8.                                           (2)
```

Consequently no finite etale morphism `A2_C -> U` exists.

## 1. The branch normalization

The discriminant branch `Gamma` of the cubic completion is parametrized by

```text
R=(q-3)(q+1)(q+2),             D=3q^2+7,
A(q)=-2q^2R/D^2,               C(q)=-qR/D.              (3)
```

For the explicit discriminant `Delta(A,C)`, the double-root incidence
equations

```text
A(q^3+7q+3)-C(C+q^2)=0,
A(3q^2+7)-2Cq=0                                      (4)
```

satisfy

```text
Res_q((4))=-A Delta.                                    (5)
```

Moreover `Res_q(D,q^3+7q+3)` is nonzero,
`Delta(0,C)=-4C^5`, and `Delta` is one reduced irreducible factor.  Thus (3)
covers the whole reduced irreducible branch, not merely one component.  On a
dense open, the double root

```text
T_double=A^2(27A-9C^2+7C)/(2(C^2-21A^2))               (6)
```

satisfies `q=T_double/A`, giving the rational inverse.

The two roots of `D` are distinct.  At each, `A(q),C(q)` have pole orders two
and one, respectively; their numerators are coprime to `D`.  At infinity the
pole orders are one and two.  These are all poles, so the affine normalization
is

```text
Gamma_tilde=P1_C minus {two roots of D, infinity},
chi_c(Gamma_tilde)=2-3=-1.                              (7)
```

## 2. Complete singular-preimage census

The four distinct parameters

```text
q=-2,-1,0,3                                             (8)
```

map to the origin.  Their tangent slopes `C/A` are

```text
-19/4, -5, infinity, 17/3,                              (9)
```

so they are four genuine normalization branches.

Pulling the two branch derivatives back along (3) gives

```text
Delta_A = 4q^3R^3(3q^2-7)^3(q^3+7q+3)/D^6,
Delta_C =-4q^4R^3(3q^2-7)^3(q^3+21q+12)/D^7.           (10)
```

The two terminal cubics in (10) are coprime, and the common numerator divisor
is exactly

```text
q^3R^3(3q^2-7)^3.                                      (11)
```

It is coprime to `D`.  Hence (8) and the two roots of `3q^2-7` are every
affine singular preimage.  The latter map to two distinct nonzero triple-root
values.  At either one, `q=C/(3A)`, so each has exactly one normalization
preimage and makes no Euler correction.  Only the four-to-one origin changes
(7), by `4-1`.  Therefore

```text
chi_c(Gamma)=-1-(4-1)=-4.                              (12)
```

## 3. Fibrewise Euler integration

On the normalization, the cubic factors as

```text
f(T)=(T-Aq)^2(T-(C-2Aq)).                               (13)
```

The gap between the simple companion and the ramified root is

```text
(C-2Aq)-Aq=qR(3q^2-7)/D^2,                             (14)
```

and the derivative at the companion is the square of (14).  Thus the
companion collides with the ramified point exactly at the origin and the two
triple-root values.  Let `S_0` be this three-point set.

The finite cubic completion has three points over `A2 minus Gamma`.  After
deleting `E`, one simple companion remains over every point of
`Gamma minus S_0`; all geometric points over `S_0` lie on `E`, so those fibres
are empty.  The constructible strata and their contributions are therefore

```text
base stratum             chi_c(base)   fibre size   contribution
A2 minus Gamma                 5            3             15
Gamma minus S_0               -7            1             -7
S_0                            3            0              0. (15)
```

Additivity and finite-cover multiplicativity on each stratum give

```text
chi_c(U)=15-7=8,                                         (16)
```

proving (2).

## 4. Finite-atlas obstruction and exact boundary

The surface `U` is integral and hence connected.  A nonempty finite etale
morphism

```text
phi:A2_C -> U                                             (17)
```

has open-and-closed image, so it is surjective of a constant positive integer
degree `d`.  On complex analytic spaces it is a finite covering, whence

```text
1=chi_c(A2_C)=d chi_c(U)=8d,                             (18)
```

impossible for every `d>=1`.

Finiteness is load-bearing.  A dominant etale finite-type map is quasi-finite,
but without properness sheets may escape at the boundary.  The nonfinite
etale open immersion

```text
D(x)=G_m x A1 -> A2                                      (19)
```

has compact Euler characteristics `0` and `1`, so the finite-cover equation
does not extend to the nonproper lane.  Therefore this Euler argument forces
any dominant etale polynomial map `A2_C -> U` into the nonfinite, nonproper
lane.  It does not itself obstruct that lane; THM-3841 and THM-3845 close it
by independent later mechanisms.

## 5. Every surviving atlas has two further boundary debts

Let `psi:A2_C -> U` be a dominant etale morphism of generic degree `d`.  The
constructible defect function

```text
delta_psi(u)=d-#psi^-1(u)                                      (20)
```

obeys THM-3578's exact compact-Euler integral.  Since `chi_c(U)=8`,

```text
integral_U delta_psi d chi_c=d chi_c(U)-chi_c(A2)=8d-1.        (21)
```

This is a signed Euler-weighted identity; it does not say every stratum has
positive Euler contribution.  Independently, THM-3811 gives

```text
Cl(U)=Pic(U)=Z^2.                                               (22)
```

In the finite normalization supplied by Zariski Main, THM-3578 therefore
forces at least two irreducible divisorial boundary primes, and the sum of
the generic divisorial sheet debts is at least two.  A hypothetical atlas
must simultaneously satisfy this class-group invoice, the Euler debt (21),
and THM-3822's intrinsic `SL2`/punctured-arm equations.  These are necessary
conditions for the surviving nonproper lane, not a contradiction.

## 6. Verification

The primary companion verifies (3)--(18) in 34 active gates.  The independent
companion separately checks branch coverage and inverse, tangent branches,
the complete singular divisor, two triple images, companion fibres, five good
finite-field point-count controls, two bad-reduction hostiles, and the
nonproper boundary in 58 active gates.  Both scripts are assertion-free and
their normal and optimized streams equal the frozen outputs.  The additional
16-gate `jc2_nonlinear_cubic_compact_euler_thm3823.py` companion separately
replays the normalization and Euler ledger.  **QED.**
