---
id: THM-3514
title: "U_full owner-boundary K4 x F13 endpoint factorization and Walsh spectrum"
status: >
  PROVED STRUCTURAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the fixed
  THM-3479 U_full refined endpoint bank, removing only the H guard exposes
  39 common guard atoms.  The speed-13 owner-in factor annihilates all 13
  middle atoms and leaves 26 boundary atoms.  Restoring the translated guard
  reconstructs the complete frozen 13^3 bank exactly.  Pair expansion has
  117 chamber/drift types, of which exactly 52=(2*2*13) survive; every
  surviving type is nonzero separately for q_H, q_q5, and their difference.
  The resulting 4x13 bridge table has rank four, all four F_2^2 Walsh rows are
  nonzero at every drift, and every Walsh row has all thirteen drift-Fourier
  modes nonzero.  This is a Cartesian endpoint factorization and an undirected
  K4 carrier, not a tournament, THM-2471 ancestry support map, physical
  current, absolute H1 class, row exclusion, or LRC(14) theorem.
source: codex/ufull-owner-boundary-independent-audit/2026-08-16
audit: >
  Independent central tensor contraction of materialized endpoint atom
  tables.  The audit imports the promoted THM-3479 endpoint primitives but no
  candidate guard/bucket companion, derives the owner and guard boundary law
  directly from PATTERN_E, uses a separately written two-pointer atomizer,
  and compares 78 direct guarded endpoint controls with guard restoration.
depends_on:
  - THM-3479-literal-half-twist-relation-current-two-transplant-certificate
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3505-r5-projective-slope-pencil-and-anchored-haar-basis
script: 04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
output: 05-knowledge/results/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.out
script_sha256: f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc
output_sha256: 7684cdb6bb1641780977d9e2def3753d802bf026e126dbc5351bd8f8ddebd906
semantic_sha256: d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6
hash_basis: LF-normalized UTF-8 bytes; exact finite-field semantic ledger
---

# THM-3514 -- the fixed `U_full` endpoint bridge has full `K4 x F_13` spectrum

**PROVED STRUCTURAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Fixed endpoint bank and owner boundary

Use THM-3479's `U_full` tuple.  In current-coordinate order its speed word is

```text
w=(1,183,27,131,53,313,13,2197,742586),               (1)
```

with guard coordinate `H` at speed one and owner coordinate `c1` at speed
thirteen.  The refined endpoint shifts are

```text
ell(alpha,beta,tau)
 =(tau,-alpha,-beta,0,0,0,0,alpha,beta)
                                      in F_13^9.       (2)
```

Let `E_(alpha,beta,tau)` be the Boolean endpoint set built from
`PATTERN_E`.  Let `E^circ_(alpha,beta)` be the set obtained by removing
**only** its `H:guard_safe` factor and retaining its owner-in and all seven
out factors.

Put

```text
s=91t=7a+r,             a in F_13,       0<=r<7.       (3)
```

The owner factor in `PATTERN_E` is

```text
distance(13t,Z)<1/14.                                  (4)
```

Since `13t=a+r/7` modulo one, its literal strict support in each sheet is

```text
[0,1/2) union (13/2,7).                                (5)
```

The THM-3479 interval engine represents the same indicator almost everywhere
by

```text
[0,1/2) union [13/2,7).                                (6)
```

Equations `(5)` and `(6)` differ at exactly thirteen singleton left
endpoints, one `r=13/2` point per sheet.  The engine integrates indicator
functions by endpoint differences, so this half-open choice changes no
endpoint integral.  Its owner-only primitive is exactly the cyclic family

```text
[14a-1,14a+1)     in units T_DEN/182,     a in F_13.  (7)
```

This direct source-level comparison fixes the boundary convention rather
than inferring it from the computed output.

## 2. The 39 common guard atoms and exact restoration

Partition the residual coordinate into the half-open chambers

```text
C_L=[0,1),        C_M=[1,6),        C_R=[6,7).          (8)
```

Thus `r=1` is assigned to `C_M` and `r=6` to `C_R` in the interval-engine
representative.  The `39=13*3` atoms are

```text
I_(a,C)={t:s=7a+r and r in C}.                         (9)
```

The `H` guard has unsafe window

```text
[-13-7tau,13-7tau) modulo 91.                         (10)
```

It is constant almost everywhere on each atom.  Write its safe value as
`g_C(a+tau)`.  Direct reduction of `(10)` gives the danger-sheet arcs

```text
D_L={12,0,1},
D_M={11,12,0,1},
D_R={11,12,0}.                                        (11)
```

By `(5)`--`(7)`, every middle atom of `E^circ_(alpha,beta)` is empty.  Only

```text
F_13 x {C_L,C_R}                                      (12)
```

remains, giving exactly 26 active owner-boundary atoms before any Fourier
contraction.

Let `AX` and `BY` denote the two linear endpoint functionals used in
THM-3479.  Define

```text
A_(alpha,beta)(a,C)
  =AX(E^circ_(alpha,beta) intersect I_(a,C)),

B_(alpha,beta)(a,C)
  =BY(E^circ_(alpha,beta) intersect I_(a,C)).          (13)
```

Linearity and `(10)` give the exact restoration formulas

```text
AX(E_(alpha,beta,tau))
 =sum_(a,C) g_C(a+tau) A_(alpha,beta)(a,C),

BY(E_(alpha,beta,tau))
 =sum_(a,C) g_C(a+tau) B_(alpha,beta)(a,C).            (14)
```

The independent audit checks `(14)` against direct guarded endpoint
evaluation for six `(alpha,beta)` pairs and all thirteen `tau`, hence 78
direct controls.  Materializing all 169 unguarded atom tables and restoring
the guard reconstructs all `13^3=2197` frozen values with digest

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682. (15)
```

The normalized inverse values in the certified split field are

```text
A^-(q_H)  =320618948602619577408,
A^-(q_q5) =503604956476841920373,
A^-(q_H)-A^-(q_q5)
           =389266878372286537904 !=0.                (16)
```

These agree exactly with THM-3479's frozen refined bank.

## 3. Relative drift leaves exactly 52 supported types

For a left atom `(a,C)` and right atom `(b,D)`, put

```text
d=b-a in F_13.                                        (17)
```

If `zeta` is the fixed primitive thirteenth root, define the guard-pair
kernel

```text
K_(C,D,d)(k)
 =sum_(tau in F_13)
   g_C(tau) g_D(d+tau) zeta^(-k tau).                 (18)
```

Translation gives

```text
sum_tau g_C(a+tau)g_D(b+tau)zeta^(-k tau)
 =zeta^(ka) K_(C,D,d)(k).                             (19)
```

All `39^2=1521` pair kernels are nonzero at each of the two frequencies used
here: `k=1` for `q_H` and `k=0` for `q_q5`.  Expand the product in `(14)`,
retain the common-sheet phase in `(19)`, include the `alpha,beta` inverse
characters, and group by `(C,D,d)`.  This produces

```text
3*3*13=117                                             (20)
```

endpoint drift types for each of `q_H`, `q_q5`, and their difference.

Exactly the `5*13=65` types with `C=C_M` or `D=C_M` vanish.  Their abstract
guard kernels are nonzero; they vanish because the actual owner-supported
endpoint atom is empty.  Every remaining type in

```text
{C_L,C_R}^2 x F_13                                   (21)
```

is nonzero separately for `q_H`, `q_q5`, and the bridge.  Therefore the
support count is sharply

```text
52=2*2*13.                                            (22)
```

Summing the 52 bridge buckets recovers the final value in `(16)`.  Restricting
instead to equal sheets, equal chambers, or equal guard atoms gives three
different values, none equal to the full bridge.  Equality of these geometric
addresses is therefore not a hidden ancestry relation.

## 4. Full `K4 x F_13` bridge spectrum

Encode `C_L=0` and `C_R=1`.  For every drift, the four chamber-pair states

```text
(0,0),(0,1),(1,0),(1,1) in F_2^2                     (23)
```

are the vertices of a canonical undirected `K4` carrier.  Let `D_ij(d)` be
the bridge bucket in state `(i,j)` and drift `d`.  The `4 x 13` matrix

```text
D=(D_ij(d))                                           (24)
```

has row rank four in the certified split field.  Its four Walsh rows are

```text
W_0=D_00+D_01+D_10+D_11,
W_L=D_00+D_01-D_10-D_11,
W_R=D_00-D_01+D_10-D_11,
W_X=D_00-D_01-D_10+D_11.                              (25)
```

Every entry of every row in `(25)` is nonzero.  Moreover, for every
`W in {W_0,W_L,W_R,W_X}` and every `k in F_13`,

```text
sum_(d in F_13) W(d) zeta^(-kd) !=0.                  (26)
```

Thus the bridge has complete support in

```text
widehat(F_2^2) x widehat(F_13).                       (27)
```

Nonzero reduction modulo the certified split prime proves nonvanishing of
the corresponding cyclotomic quantities.  Rank four modulo that prime proves
that a characteristic-zero `4 x 4` minor is nonzero.

The four objects in `(23)` are vertices of a complete undirected carrier.
No intrinsic binary observable or orientation has been supplied, so they do
not form a tournament.  The Walsh basis records the constant, two coordinate,
and mixed/checkerboard channels; it does not manufacture directed arcs.

## 5. The root-label debt reduces to one common torsor gauge

There is an elementary equivariant alignment, but it is weaker than physical
root realization.  The guard pullback uses the regular translation

```text
a -> a+tau.                                             (28)
```

On THM-2471's first-collision chart,

```text
w_u=(y+u)/13,
w_u+tau/13=w_(u+tau),                                  (29)
```

so the root labels carry the same regular `+tau` action.  Every equivariant
bijection between two regular `F_13` torsors is

```text
u=a+c,                    c in F_13.                    (30)
```

Indeed, equivariance gives `f(a)=f(0)+a`; hence there are exactly thirteen
choices and no canonical one.  Changing `c` multiplies a primitive
frequency-`k` coefficient by `zeta^(kc)`, so it cannot change nonvanishing.
If one common gauge is used on both endpoint legs, then

```text
u_R-u_L=a_R-a_L=d.                                    (31)
```

Independent gauges instead change the drift by `c_R-c_L`.  Thus a genuine
common-ancestry construction would have to carry one common gauge.

Equations `(28)`--`(31)` align thirteen abstract labels only.  The chamber
bit remains, no guard atom is identified with a physical root node, and no
base, owner, word, source, or horizon support is supplied.

## 6. Cartesian pairing is not THM-2471 ancestry

The construction takes the Cartesian product of the independently
marginalized endpoint tables in `(13)`.  It does not construct a support
predicate on THM-2471's Boolean stalk variables

```text
(y,u,v,a,b,e')                                        (32)
```

and does not retain a common base, collision root, owner branch, word,
source sheet, or horizon.  The relative sheet drift `(17)` is an exact guard
covariance coordinate, not a map into the THM-2471 ancestry fibre product.

Consequently `(14)`--`(31)` prove endpoint spectral closure and abstract
label alignment only.  They do
not prove any of

```text
common ancestry or a root-character bispectrum realization;
a lawful Boolean or physical current;
a nonzero absolute H^1 class;
a grouped exact-address coefficient C(a;X,m);
an all-unit projector B(q);
exclusion of a scalar row;
LRC(14).                                               (33)
```

In particular, the endpoint edge differences remain compatible with the
`B^1`/zero-absolute-`H^1` boundary already recorded in THM-3479.

## 7. Exact audit package

Run from the repository root:

```text
python -B 04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
python -B -O 04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py
```

Both executions reproduce

```text
05-knowledge/results/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.out
```

exactly.  The semantic ledger is

```text
d52c9f0a56c14a83e1e6b175c7b725314c99f09d44509bc8582847a5857f7da6.
```

The audit imports no candidate guard-sheet or drift-bucket script.  It uses
THM-3479's pinned endpoint primitives as the established parent mechanism,
then independently materializes the atom tables and performs the guard
restoration, pair-kernel construction, central drift contraction, rank test,
Walsh transform, and drift Fourier transform.
