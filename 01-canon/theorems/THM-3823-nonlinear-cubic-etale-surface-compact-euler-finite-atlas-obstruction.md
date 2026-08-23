---
id: THM-3823
title: "Nonlinear cubic etale surface has compact Euler characteristic eight"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over C, the
  THM-3811 affine etale surface U has compactly
  supported Euler characteristic 8.  Hence no finite etale morphism
  A2_C -> U exists.  Every dominant etale plane atlas is necessarily
  nonfinite; if its generic degree is d, its exact signed Euler sheet debt is
  8d-1, and its finite envelope has at least two divisorial boundary primes.
  This Euler argument alone does not obstruct a nonproper dominant atlas;
  THM-3845 subsequently excludes every polynomial plane atlas of this
  surface.  No general planar Jacobian claim follows.
source: root / nonlinear-cubic Euler-atlas obstruction lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (jc_sparse_direct_search and root,
  2026-08-23).  The audits independently reconstructed the normalization as
  P1 minus two
  denominator roots and infinity, checked pole orders and absence of
  cancellations, and used uniqueness of the double/triple root to verify
  that only q=-2,-1,0,3 are identified.  They checked the two singleton
  nonzero triple fibres and the exact three-sheet/simple-companion
  stratification, obtaining chi_c(U)=8.  The finite-cover contradiction and
  both THM-3578 consequences were rederived with their proper/nonproper and
  signed-Euler qualifications intact.  A direct singular-locus elimination
  separately recovered exactly the origin and two triple-root targets and
  confirmed that signed Euler debt is not the unweighted divisorial debt.
  No repair was found.  The 16-gate exact companion checks the
  rational branch normalization, all three normalization punctures, the
  four distinct origin preimages, the two affine nonorigin triple-root
  parameters, and the complete compact-Euler ledger.  The proof separately
  checks the projection strata against THM-3811 and the boundary-rank and
  Euler-debt consequences against THM-3578.  Normal and optimized replay
  byte-match the frozen output, both raw hashes match, and the script has no
  inactive Python assert.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3578-zariski-main-boundary-rank-and-sheet-debt
related:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
  - THM-3600-danielewski-arm-plane-atlas-singular-shear-and-no-filling
  - THM-3845-nonlinear-cubic-keller-atlas-total-degree-contradiction
script: 04-computation/jc2_nonlinear_cubic_compact_euler_thm3823.py
output: 05-knowledge/results/jc2_nonlinear_cubic_compact_euler_thm3823.out
script_sha256: 713e7a925e4ae42852dc22c5db766510ff0cca74ce6487e224bc883c1659c5e5
output_sha256: ff8f7c8ac1d35b46a4a4f237e3f36c1024fcb9f26a5fa0a782dff16b555e204a
semantic_sha256: ecf52db47db9533b4bf1eec5c23bea5a0e2a91e33a21e6c941863e0cba937cfe
hash_basis: raw LF bytes
---

# THM-3823 -- the nonlinear cubic surface has Euler mass eight

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over `C`.
Let

```text
pi:U -> A2_(A,C)                                                   (1)
```

be the smooth affine etale surface and degree-three quasi-finite projection
constructed in THM-3811.  Then

```text
chi_c(U)=8.                                                        (2)
```

In particular, there is no finite etale morphism `A2 -> U`.  More generally,
every dominant etale morphism `psi:A2 -> U` is nonfinite.  If its generic
degree is `d`, then

```text
d chi_c(U)-chi_c(A2)=8d-1                                         (3)
```

is its exact signed Euler sheet debt in the sense of THM-3578, and every
finite envelope of `psi` has at least two divisorial boundary primes.

## 1. Euler characteristic of the branch

THM-3811 parametrizes the irreducible discriminant curve `Gamma` by

```text
R(q)=(q-3)(q+1)(q+2),
A(q)=-2q^2R(q)/(3q^2+7)^2,
C(q)=-qR(q)/(3q^2+7).                                             (4)
```

This is the normalization away from the explicitly audited singular fibre.
Its affine normalization is

```text
Gamma_tilde = P1 minus ({3q^2+7=0} union {infinity}).              (5)
```

Indeed the quadratic `3q^2+7` has two distinct roots.  Neither numerator in
`(4)` vanishes there, while `A` and `C` have poles of orders one and two at
infinity.  Thus

```text
chi_c(Gamma_tilde)=2-3=-1.                                        (6)
```

Exactly four distinct normalization values

```text
q=-2,-1,0,3                                                       (7)
```

map to the square-zero origin.  Away from that target, the double root of
the cubic is unique; at a triple root, the root is again unique.  Hence
there are no other identifications in the normalization map.  Identifying
the four points `(7)` to one lowers compactly supported Euler characteristic
by three.  The two cuspidal triple-root values are bijective normalization
fibres and do not alter Euler characteristic.  Therefore

```text
chi_c(Gamma)=-1-4+1=-4.                                           (8)
```

The exact companion verifies `(4)` directly against the discriminant of
THM-3811, checks all puncture cancellations and degrees, and checks that the
four origin parameters are squarefree.  It also verifies that
`3q^2-7=0` gives two distinct affine parameters disjoint from `(7)` and
from the normalization punctures; these are exactly the two nonzero
triple-root values of THM-3811.

## 2. The degree-three projection stratifies as 3 plus 1

Over `A2 minus Gamma`, the finite cubic completion is etale of degree three,
and removing the ramification divisor removes nothing.  Thus

```text
chi_c(A2 minus Gamma)=1-(-4)=5,
chi_c(pi^-1(A2 minus Gamma))=3*5=15.                              (9)
```

Over a nonexceptional point of `Gamma`, the double root lies on the deleted
ramification divisor and the unique simple companion remains in `U`.  The
companion therefore gives one sheet over

```text
Gamma^o=Gamma minus {origin, two triple-root values}.              (10)
```

The fibre of `U` is empty at each of the three deleted target values: the
origin has only the square-zero vertex, and each triple-root fibre has only
the ramified root.  Consequently

```text
chi_c(Gamma^o)=-4-3=-7,
chi_c(U)=15+(-7)=8,                                               (11)
```

which proves `(2)`.

This calculation also explains why the answer is larger than the generic
degree.  The branch itself has negative Euler characteristic: replacing the
three generic sheets by one companion sheet over `Gamma^o` adds rather than
subtracts Euler mass.

## 3. Finite atlases are impossible

If a finite etale morphism `f:A2 -> U` had degree `d>=1`, then over the
complex topology it would be a finite covering.  Compactly supported Euler
characteristic is multiplicative for finite coverings, so `(2)` would give

```text
1=chi_c(A2)=d chi_c(U)=8d,                                       (12)
```

an impossibility.  A proper etale plane atlas is also impossible, since a
proper quasi-finite morphism is finite.

This is not a contradiction to the open source-plane problem.  A
noninvertible Keller map is expected to be nonproper, and a dominant etale
morphism of affine surfaces is only quasi-finite in general.

## 4. Every surviving atlas has a large boundary invoice

Let `psi:A2 -> U` be any dominant etale morphism of generic degree `d`.
THM-3578 applies because `U` is a normal smooth affine surface.  Its
constructible fibre-defect function

```text
delta_psi(u)=d-#psi^-1(u)                                        (13)
```

satisfies the exact Euler integral

```text
integral_U delta_psi d chi_c
  =d chi_c(U)-chi_c(A2)=8d-1.                                   (14)
```

This is a signed Euler-weighted identity, not a claim that every stratum
contributes positively.  Separately, THM-3811 gives

```text
Cl(U)=Pic(U)=Z^2.                                                 (15)
```

The Zariski-main boundary-rank inequality of THM-3578 therefore forces at
least two irreducible divisorial boundary primes in every finite envelope of
`psi`, and the sum of the generic divisorial sheet debts is at least two.

Together with THM-3822, a hypothetical atlas must now pay three different
invoices:

```text
geometric arm:  every component of h=0 carries a nonconstant unit;
divisor class:  at least two finite-envelope boundary primes;
Euler mass:     exact signed sheet debt 8d-1.                     (16)
```

The useful constructive question is whether one nonproperness passport can
pay all three invoices simultaneously.  Standard finite covers and standard
hyperbolic arms cannot.

## 5. Exact replay and scope

Run

```bash
python3 04-computation/jc2_nonlinear_cubic_compact_euler_thm3823.py
python3 -O 04-computation/jc2_nonlinear_cubic_compact_euler_thm3823.py
```

Both executions must byte-match

```text
05-knowledge/results/jc2_nonlinear_cubic_compact_euler_thm3823.out
```

The companion reports `CHECKS=16` and `RESULT=PASS`.  It contains no inactive
Python `assert`.  This Euler mechanism rules out finite or proper plane
atlases only; THM-3845 later excludes the nonfinite case by an independent
degree argument.  **QED.**
