---
id: THM-4288
title: "A23 partial-normalization relative differential and etale base-change obstruction"
status: >
  PROVED FORMAL-LOCAL + FINITE-EXACT AUDIT PASS, RELATIVE TO THM-4272/4280/4284
  only for the A23 and degree-shell interpretation. Every strict intermediate
  A_m subset A_ell, 1<=ell<m, has relative-differential length
  min(m-ell,2ell) and dual-number closed fibre, so it cannot be a full
  Cartesian base change of an etale Keller map. The final triple-point
  blowdown has the exact one-dimensional defect c_C'-c_R'-c_E'. This isolates
  common-tangent descent but does not construct the missing Cartesian square,
  prove ambient regularity, or prove JC(2).
source: root/continuation-frontiers-20260829b/2026-08-29
depends_on:
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4280-integral-three-channel-fat-contact-observer-and-sharp-five-jet-bound
  - THM-4284-a23-conductor-defect-and-degree-shell-first-character-nondescent
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction
primary_script: 04-computation/jc23_a23_partial_normalization_differential_thm4288.py
primary_output: 05-knowledge/results/jc23_a23_partial_normalization_differential_thm4288.out
primary_script_sha256: a7ef08760646f3e757112373ef456d9a5e0fc8413fe21a90860f4c44e7720bbc
primary_output_sha256: 723e9b457bae9bc53db0859fa06486622cc7659fae23ebbe59bf04aecf899a8a
hash_basis: raw LF bytes
audit: >
  PASS. An adversarial algebra audit independently checked the A-algebra
  presentation, distinguished the A- and B-annihilators, and restricted the
  ramification claim to 1<=ell<m. A standard-library finite-jet checker over
  F_1000000007 reconstructs the truncated pair rings, differential quotient,
  annihilator intersection, dual-number fibre, and unique triple-point
  functional. Ordinary and optimized outputs are byte-identical. The Keller
  statement is explicitly only a Cartesian-base-change obstruction.
---

# THM-4288 -- A23 partial-normalization relative differential and etale base-change obstruction

**PROVED FORMAL-LOCAL + FINITE-EXACT AUDIT PASS. THE KELLER CARTESIAN SQUARE,
RAW DESCENT, AND JC(2) REMAIN OPEN.**

## 1. Strict partial normalizations

Let `k` be a field of characteristic different from two, put `R=k[[b]]`, and
for `j>=1` define

```text
A_j={(f_0,f_1) in R direct-sum R : f_1-f_0 in b^j R}
   ~=k[[b,q_j]]/(q_j(q_j-b^j)).                         (1)
```

Fix

```text
m>=2,                    1<=ell<m,                    n=m-ell,
A=A_m,                   B=A_ell.                       (2)
```

Inside `R direct-sum R`, write

```text
t=(0,b^ell),             q=(0,b^m)=b^n t.              (3)
```

Then

```text
B=A[t]
 ~=A[T]/(T^2-b^ell T, b^n T-q),                        (4)
B/A~=b^ell R/b^m R,
length_k(B/A)=m-ell.                                    (5)
```

The relative Kähler module is

```text
Omega_(B/A)~=B/(b^n,2t-b^ell) dt.                       (6)
```

Put

```text
d=min(m-ell,2ell).                                      (7)
```

There are exact identifications

```text
B/(b^n,2t-b^ell)~=R/(b^d),
length_k Omega_(B/A)=d,                                 (8)

Ann_B Omega_(B/A)=Fitt_B^0 Omega_(B/A)
                  =(b^n,2t-b^ell),                     (9)

Omega_(B/A)~=A/(b^d,q),
Ann_A Omega_(B/A)=Fitt_A^0 Omega_(B/A)=(b^d,q).         (10)
```

The distinction between `(9)` and `(10)` is essential: `2t-b^ell` need not
belong to `A`.

## 2. Proof of the differential formula

Every element of `B` has branch difference in `b^ell R`, so it can be written
as a diagonal element plus an `A`-multiple of `t`. This proves `B=A[t]`.
The displayed quotient in `(4)` eliminates `q` as `b^nT` and reduces to

```text
R[T]/(T^2-b^ell T)=A_ell,                              (11)
```

which proves that `(4)` has no hidden kernel.

For an `A`-derivation `D:B->M`, set `v=D(t)`. Differentiating the two
relations in `(4)` gives

```text
b^n v=0,                   (2t-b^ell)v=0.              (12)
```

Conversely these equations define an `A`-derivation: the ambiguity in writing
an element as `a+ct` has coefficient in `(b^n,q)`, and `qv=t b^n v=0`.
This proves `(6)`.

Modulo `2t-b^ell`, one has `t=b^ell/2`; substituting in
`t(t-b^ell)=0` gives `b^(2ell)=0`. Together with `b^n=0` this proves `(7)--(9)`.
The map `A->R/(b^d)` sends `b` to `b` and `q` to zero. Its kernel is
`(b^d,q)`, proving `(10)`.

The finite map `Spec B->Spec A` is subintegral, but not unramified. At its
closed point,

```text
Omega_(B/A) tensor_B k ~= k dt,                         (13)
```

and, more visibly,

```text
B tensor_A k=B/(b,q)~=k[t]/(t^2).                       (14)
```

Thus every strict intermediate partial normalization in `(2)` has a
dual-number contact fibre.

The scope `1<=ell<m` is sharp. At `ell=m` the map is the identity. The full
normalization corresponds informally to `ell=0`; its relative Kähler module
vanishes, although the map is still non-etale because it is not flat. Hence
the unqualified claim that every proper normalization is ramified would be
false.

## 3. The exact etale base-change firewall

For any base change `A->C` and any point of `B tensor_A C` above the contact,
base change of Kähler differentials and `(13)` leave a one-dimensional
relative cotangent fibre. Therefore `A->B` cannot be, even completed-locally
at the contact, the **full pullback in a Cartesian square** of an etale map.

A planar Keller map is etale, so this gives a precise conditional firewall:

```text
if a Keller boundary construction identified A->B with its full
Cartesian base change, then ell=m.                       (15)
```

No such Cartesian identification is known. A graph image, saturation,
normalization, selected component, or map into a pullback is not itself a
base change of the Keller morphism. Thus `(15)` does not exclude a Keller pair;
it identifies a theorem that a future compactification argument would have
to prove before invoking etaleness.

## 4. A23 specialization and the differential boundary

THM-4272 gives the raw contact `A_12`, and THM-4284 gives actual saturated
graph types with `ell in {1,2,4}`. Equations `(5)` and `(8)` become

| `ell` | graph type | conductor loss `12-ell` | `length Omega_(A_ell/A_12)` |
|---:|---:|---:|---:|
| 1 | `A_1` | 11 | 2 |
| 2 | `A_3` | 10 | 4 |
| 4 | `A_7` | 8 | 8 |

For `ell=1`, `dt` pulls to `(0,db)` on the normalization, so modulo the
diagonal branch cotangent it is literally the anti-diagonal direction. For
`ell>=2`, it pulls to zero on normalized branch cotangents and is instead a
singular/conormal direction.

The table supplies another firewall: relative-differential length does not
measure the full conductor loss. At `ell=1`, two differential lengths see an
eleven-dimensional function-gluing defect. A cotangent obstruction can prove
non-unramifiedness without reconstructing all raw descent jets.

## 5. Exact terminal triple-point defect

Start from the two branches `q=0` and `q=b^m`. After `m-1` common-tangent
blowups, put

```text
z=q/b^(m-1).                                             (16)
```

Immediately before the final blowup that makes the reduced total transform
nodal, the local reduced divisor is the ordinary plane triple point

```text
D=Spec k[[b,z]]/(b z(z-b)),                             (17)
E: b=0,                 R: z=0,                 C: z=b.
```

Let

```text
V={(e(z),r(b),c(b)):e(0)=r(0)=c(0)}                    (18)
```

be the value-gluing ring in the normalization. Define

```text
lambda(e,r,c)=c'(0)-r'(0)-e'(0).                        (19)
```

Then restriction from the plane has the exact sequence

```text
0 -> O_D -> V --lambda--> k -> 0.                       (20)
```

In particular, the seminormal/value-gluing quotient retains exactly one
first-order class beyond the plane triple point.

To prove `(20)`, let the common value be `a` and set

```text
F_0(b,z)=e(z)+r(b)-a.                                   (21)
```

It has the required restrictions on `E` and `R`. Its error on `C` is

```text
h(b)=c(b)-e(b)-r(b)+a,
h(0)=0,                       h'(0)=lambda(e,r,c).       (22)
```

If `lambda=0`, write `h=b^2H` and use

```text
F=F_0+bzH(b).                                           (23)
```

This has all three required restrictions. Conversely the chain rule for an
ambient linear form `alpha b+beta z` gives branch derivatives

```text
(e'(0),r'(0),c'(0))=(beta,alpha,alpha+beta),            (24)
```

so every plane restriction kills `lambda`. The hostile triple `(0,0,b)` has
equal values and `lambda=1`, proving surjectivity.

## 6. Last-blowup descent and the actual-shell hostile

Blow up the origin in `(17)`. The three strict branches meet the new central
exceptional `P^1` at three distinct nodes. Maps from this reduced nodal star
to a separated target glue by equality of node values; their first derivatives
need not obey `(24)`.

Translate a common target value on `E_0` to the origin and use its formal
logarithm. If the maps on `E`, `R`, and the exceptional tree are constant and
the map on `C` has

```text
L_C(b)=c_1 b+c_2 b^2+...,                               (25)
```

then descent through the last blowdown is equivalent to

```text
lambda=c_1=0
iff
L_C in b^2 k[[b]].                                      (26)
```

Thus the last triple-point blowdown detects precisely common tangent, or
graph type at least `A_3`. The raw filtration

```text
A_1 superset A_2 superset ... superset A_m,
A_j/A_(j+1)~=b^jR/b^(j+1)R                              (27)
```

has `m-1` one-dimensional conductor layers, but `(26)` identifies only its
first character; it does not prove full order-`m` descent.

THM-4280 proves that every actual degree-34 or degree-42 `C_0->E_0` class has
`c_1!=0`. Put such a map on the positive-genus component and the common
constant on every rational and exceptional component of the reduced resolved
curve. This is a genuine resolved-curve morphism with matching node values,
but `(26)` proves that it does not descend across the last blowdown.

This hostile shows that the following data are insufficient, even together:

- agreement at every resolved node;
- constancy on the rational exceptional tree;
- properness of the elliptic target;
- membership in the actual `E_0` Hom lattice; and
- the correct degree-34/42 shell.

It is a local resolved-data hostile, not an assertion that this morphism is a
Keller response.

## 7. A smaller positive bridge

There are two precise sufficient hypotheses for `(26)`.

1. If the resolved map extends to an `E_0`-valued morphism on a surface
   neighbourhood of the final exceptional `P^1`, it descends through that
   blowup. Every `P^1->E_0` is constant; after shrinking into an affine target
   neighbourhood, `sigma_* O_(Bl_p S)=O_S` descends the coordinate functions.
2. For common-tangent descent alone, it suffices that the three pullbacks of
   the invariant differential be restrictions of one ambient-regular Kähler
   one-form. Its tangent samples obey `(24)`; constants on `E` and `R` force
   the `C` sample to vanish.

A dualizing or logarithmic differential is not thereby enough: neither has
been identified with an ambient Kähler form here. THM-3955/3957 show why
normalization and torsion-free passage can erase conductor-supported
cotangent data. The open Keller bridge is therefore not generic properness;
it is an ambient-neighbourhood extension, or at least promotion of the Keller
residue from dualizing/logarithmic to ambient-regular Kähler type on this
terminal chart.

## 8. Cross-frontier mechanism and audit

THM-4286 and this theorem instantiate one operation-level mechanism without
transferring mathematics between LRC and the Jacobian problem:

| lane | coarse quotient | next operation | forgotten coordinate | exact sidecar |
|---|---|---|---|---|
| LRC | old inactive signature | replace deleted masks | responder activity/margin | debt responses; on index 396 one exact margin completes the Boolean profile |
| A23 | resolved node values | blow down final exceptional | `c_C'-c_R'-c_E'` | ambient first-order circuit |

In both cases the quotient is a valid observation but not a congruence for the
next native operation. The connection is methodological; there is no functor
or theorem transferring an LRC deck to a plane curve.

The exact checker

```text
04-computation/jc23_a23_partial_normalization_differential_thm4288.py
```

constructs truncated pair rings over `F_1000000007`, checks the differential
quotient, the `A`-annihilator intersection, the dual-number fibre, and the
rank-one triple functional. It freezes the A23 cases `ell=1,2,4`, the hostile
near-identity case `ell=11`, and the basic `m=2,ell=1` control. Its output is

```text
05-knowledge/results/jc23_a23_partial_normalization_differential_thm4288.out
```

Reproduce with

```bash
python3 -B 04-computation/jc23_a23_partial_normalization_differential_thm4288.py
python3 -B -O 04-computation/jc23_a23_partial_normalization_differential_thm4288.py
```

The formal proof, not the finite-field audit, establishes the theorem. Nothing
here constructs the Cartesian square in `(15)`, proves an ambient Keller
extension, gives raw `A_23` descent, enters the exact-`M=12` seam, or proves
`JC(2)` or `DC(2)`.
