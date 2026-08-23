---
id: THM-3824
title: "Rule 30 fixed-division tariff and physical phase separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every period
  n, old image depth r, and fixed division exponent e, the largest sublattice
  on which A->2^(-e)T_n A is integral and remains equal modulo T_n^r Z^n is
  explicit.  Its quotient separates coordinatewise integrality from one
  odd-core carry, with all e=0, r=0, n=1, pre-core, and post-core boundaries.
  On the true nonnegative Rule 30 phase ray, the exact two-slot free defect
  is strictly increasing at every fixed scale, so no same-scale exact Smith
  collision realizes a nonzero tariff class.  The fixed branch is not a
  congruence for the nonlinear exact-gap stratum; every Rule 30 prize remains
  open.
source: root + rule30_next_operation + rule30_tariff_audit / incoming-signal extension session, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (rule30_tariff_audit, 2026-08-23).  A second
  companion constructs the actual image lattices by column HNF, the induced
  integral action, and K/J by unimodular Smith kernels.  It checks 1,134 full
  tariff filtrations and 213,052 exhaustive residue definitions, and finds
  the n=2 exact-gap hostile that repairs the scope.  The primary, physical,
  and independent companions contain no Python assert gate; normal and
  optimized raw LF streams match all frozen transcripts.
depends_on:
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
  - THM-3804-rule30-all-period-amplitude-lattice-smith-law
related:
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
script: 04-computation/rule30_fixed_division_tariff_thm3824.py
output: 05-knowledge/results/rule30_fixed_division_tariff_thm3824.out
script_sha256: 16e402d570b05aaca6402872edc1a84daa8866d609df4d37ba9433d1e597879d
output_sha256: 9223c5026576dff9cd20c9445ed1f9e5a3b4af52d6e803cc54a9bffd793482b9
semantic_sha256: aed0155a46398446366eb19751960b3276ea144777a483d50fa68cd2589b815e
physical_script: 04-computation/rule30_physical_sibling_defect_thm3824.py
physical_output: 05-knowledge/results/rule30_physical_sibling_defect_thm3824.out
physical_script_sha256: 15dd9d342497230b2393128c8f8f7cc13a8f2441b8e2e4d19c1f0c825b59437c
physical_output_sha256: d6d35abdb005763a2f7f444f3269b0b029fddcd1a8385ae1844ce46454117537
physical_semantic_sha256: bb16b41f50e1d12b847dfbf9696bd7e061626163a339d68a98d30cd4162a18d5
independent_script: 04-computation/rule30_fixed_division_tariff_independent_audit_thm3824.py
independent_output: 05-knowledge/results/rule30_fixed_division_tariff_independent_audit_thm3824.out
independent_script_sha256: de1a2013c942ab5b0a5bc4dc5f294d0adc6a057b3e604f4461655d56a6ebec6b
independent_output_sha256: 9bff73508609556be5f41ea32b232c23587935410ab4c637d9d1c2d9b0bcb95b
independent_semantic_sha256: 5e28eecc3d2178fa35d9993e209b1a12ae76d7b00355b6f58e6d61c0ed655827
hash_basis: raw LF bytes
---

# THM-3824 -- exact fixed-division tariff and the physical phase boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The theorem
below applies the next native amplitude operation to THM-3804's Smith
quotient.  It proves no center-column periodicity, balance, density, or
complexity result.  Every Rule 30 prize remains **OPEN**.

## 1. The fixed linear branch

For `n>=1`, on coordinates indexed by `Z/nZ`, retain

```text
(T_n A)_j=A_(2j)+A_(2j+1).                            (1)
```

Fix `r,e>=0` and put

```text
L=T_n^r Z^n,
D_e={A in Z^n:T_nA in 2^e Z^n},
N_e(A)=2^(-e)T_nA             on D_e.                 (2)
```

The largest translation sublattice `H<=L` for which `H<=D_e` and
`N_e(H)<=L` is

```text
K=L intersect T_n^(-1)(2^eL).                         (3)
```

Indeed, the two requirements on `H` are together exactly
`T_nH<=2^eL`, hence `H<=K`; conversely `K` has both properties.  Thus
`L/K` is the exact extra tariff needed to make equality modulo `L` a
congruence for the **fixed** linear branch `N_e` on its divisibility domain.

Put

```text
J=L intersect D_e.                                    (4)
```

Then

```text
K <= J <= L,
L/J = failure of coordinatewise integral division,
J/K = post-integrality Smith carry,
L/K = combined fixed-division tariff.                 (5)
```

## 2. All-period group law

Write

```text
n=2^a m,          m odd,
b=min(a,r),       d=n/2^b,
s=max(0,r-a).                                          (6)
```

When `e=0`, one has `K=J=L`, so all three quotients in `(5)` are trivial.
Assume `e>=1`; omit every `Z/1` summand.

Before the odd core, `r<a`, so `d=n/2^r` is even, and

```text
L/K ~= (Z/2^e)^(d/2),
L/J = L/K,                 J/K=0.                     (7)
```

At and after the odd core, `r>=a`, so `d=m` is odd, and

```text
L/K ~= (Z/2^e)^(d-1) direct-sum Z/2^(e-1),            (8)

L/J ~= (Z/2^e)^(d-1)
       direct-sum Z/2^max(e-s-1,0),                   (9)

J/K ~= Z/2^min(s,e-1).                                (10)
```

These formulas include every boundary.  If `r=0<n` with `n` even, `(7)`
has `d=n`; if `n` is odd, `(8)--(10)` apply for every `r>=0`.  At `n=1`,
the total tariff is `Z/2^(e-1)`.  At the odd-core boundary `s=0`, the carry
in `(10)` is trivial.

### Proof

THM-3804 gives the exact factorization

```text
T_n^r=E_(n,d) T_d^s Q_(n,b),                          (11)
```

where `Q_(n,b)` is onto and `E_(n,d)` is injective.  Therefore

```text
L=E T_d^s Z^d.                                        (12)
```

The parametrization `phi_s(z)=E T_d^s z` is injective.  This is immediate
when `s=0`; if `s>0`, then `d` is odd and `det T_d=+/-2`, so `T_d^s` has no
kernel.  Since `T_nE=ET_d`, the condition `phi_s(z) in K` is

```text
E T_d^(s+1)z=2^e E T_d^s w        for some w in Z^d. (13)
```

Cancel the injective `E`, then the common injective map `T_d^s` when
`s>0`.  No integral inverse or determinant division is used.  Equation
`(13)` becomes

```text
T_dz in 2^e Z^d,                                     (14)
```

and hence

```text
L/K ~= image(T_d mod 2^e).                            (15)
```

For even `d`, THM-3804 gives

```text
SNF(T_d)=diag(1^(d/2),0^(d/2));                       (16)
```

for odd `d`, it gives

```text
SNF(T_d)=diag(1^(d-1),2).                             (17)
```

Unimodular Smith changes remain invertible modulo `2^e`, so `(15)--(17)`
prove `(7)--(8)`.

For `J`, condition `(4)` pulls back to

```text
T_d^(s+1)z in 2^e Z^d.                                (18)
```

On the odd core,

```text
SNF(T_d^(s+1))=diag(1^(d-1),2^(s+1)),                 (19)
```

which proves `(9)`.  Finally `T_d^k 1=2^k 1`.  The constant-vector subgroup
killed modulo `2^e` has the full kernel size forced by `(19)`, so it is the
entire cyclic kernel.  Quotienting the `k=s+1` kernel by the `k=1` kernel
proves `(10)` and rules out a hidden nonsplit extension.

## 3. Exact-gap hostile: the branch label is external

The theorem is deliberately about divisibility by `2^e`, not preservation
of the nonlinear stratum where `e` is the exact common valuation.  At

```text
n=2, r=0, e=1,
x=(1,1),          delta=(2,0) in K,          x'=x+delta=(3,1),
```

one has

```text
T_2x=(2,2)      with exact gap 1,
T_2x'=(4,4)     with exact gap 2.                      (20)
```

Thus translation by `K` can change the adaptive branch.  An application to
a physical exact-gap slice must separately retain the gap owner and the
all-odd normalized type.

The smallest nontrivial total carry is the complementary control

```text
n=2, r=3, e=2,       delta=(4,4) in L,
N_e(delta)=(2,2) notin L,       L/K ~= Z/2.            (21)
```

## 4. The physical exact free defect separates phase

Retain THM-3512's packed Rule 30 orbit

```text
R_0=1,
R_(t+1)=R_t xor ((2R_t) or (4R_t)),
q=2^m,
U_m(t)=(R_(t+q)-R_t)/2^(v_m),                         (22)
```

where `v_m` is phase-independent and every `U_m(t)` is odd.  For a true
sibling pair define its exact Smith free coordinate

```text
Delta_m(t)=U_m(t+q)-U_m(t).                            (23)
```

Then, for every `m,t>=0`,

```text
boxed: Delta_m(t+1)>Delta_m(t)>0.                     (24)
```

Consequently no two distinct nonnegative phases at one fixed scale have the
same exact `Q_(2,r)` Smith shadow for any old quotient depth `r>=1`; their
integer free coordinates already differ.  Hence no same-scale physical
sibling pair realizes a nonzero class of the ambient tariff by an exact
Smith collision.

### Proof of `(24)`

Induction in `(22)` shows that `R_k` has top bits at positions `2k` and
`2k-1`, no higher bit, and a permanent low bit.  Thus, for `k>=1`,

```text
(3/2)4^k < R_k < 2*4^k.                               (25)
```

Put `c=4^q` and clear the fixed positive denominator:

```text
M_t=2^(v_m)Delta_m(t)
   =R_(t+2q)-2R_(t+q)+R_t.                            (26)
```

The two sides of `(25)` give

```text
M_t     < (2c^2-3c+2)4^t,
M_(t+1) > (6c^2-16c+6)4^t.                           (27)
```

Their difference is positive because

```text
4c^2-13c+4>0                     for c>=4.            (28)
```

The same bounds give

```text
M_t>(3c^2/2-4c)4^t>0.                                (29)
```

Dividing by `2^(v_m)` proves `(24)`.

The separation in `(24)` is an unbounded encoding of absolute chronology,
not a finite-state observer.  Cross-scale comparisons, finite reductions,
off-ray sections, larger physical blocks, the projective owner, and center
chronology remain outside the theorem.

## 5. Exact verification

The primary companion checks `40` effective Smith forms, `48` odd-core
power Smith forms, `1,296` power-image/kernel gates, and `3,960` all-period
tariff instances.  The physical companion checks `35,842` strict
inequalities through `m=8,t=4096` as a control of `(24)`.

The independent companion instead constructs the actual lattices.  It
checks `162` image HNFs, `1,134` full `K/J` filtrations, `213,052`
exhaustive residue definitions, all boundary types, and hostile `(20)`.

Run

```bash
python3 -B 04-computation/rule30_fixed_division_tariff_thm3824.py
python3 -B -O 04-computation/rule30_fixed_division_tariff_thm3824.py
python3 -B 04-computation/rule30_physical_sibling_defect_thm3824.py
python3 -B -O 04-computation/rule30_physical_sibling_defect_thm3824.py
python3 -B 04-computation/rule30_fixed_division_tariff_independent_audit_thm3824.py
python3 -B -O 04-computation/rule30_fixed_division_tariff_independent_audit_thm3824.py
```

Each normal and optimized raw LF stream equals its named frozen output.
