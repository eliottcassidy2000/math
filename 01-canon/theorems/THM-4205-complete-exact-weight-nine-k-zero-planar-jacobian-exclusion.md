---
id: THM-4205
title: "Complete exact-weight-nine K-zero planar-Jacobian exclusion"
status: >
  PROVED RELATIVE TO THM-4147/4176/4192/4209 + VERIFIED-EXACT +
  INDEPENDENTLY SOURCE-AUDITED. On the live b=d=0 reduced (2,3) seam, the
  entire exact-weight-nine K=0 coefficient face is empty of nonautomorphic
  planar Keller pairs. The new Y-only row eta=0,zeta!=0 has source residual
  length 20 and affine critical length 24 for arbitrary Phi,Theta, without
  any projected-discriminant hypothesis. Its complete packet
  (11,8,2,2,2,1), prime pure-cubic owner, and strict finite/full deficits
  close the last coefficient row. THM-4192, THM-4209, and THM-4176 supply the
  exhaustive P-only, off-anti mixed, and anti-diagonal rows. Entry, other
  reduced cells, M>=10, JC(2), and DC(2) remain OPEN.
source: codex-jc-kzero-w9-complete-owner-audit-20260826
depends_on:
  - THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion
  - THM-4176-complete-repeated-top-wall-planar-jacobian-exclusion
  - THM-4192-complete-p-only-k-zero-planar-jacobian-exclusion
  - THM-4209-generic-mixed-b-k-zero-exact-weight-nine-planar-jacobian-exclusion
related:
  - THM-4155-generic-y-only-delta-zero-weight-nine-planar-jacobian-exclusion
script: 04-computation/jc23_k0_y_only_complete_exclusion_thm4205.py
output: 05-knowledge/results/jc23_k0_y_only_complete_exclusion_thm4205.out
independent_audit_script: 04-computation/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.py
independent_audit_output: 05-knowledge/results/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.out
script_sha256: 4bc1b18c3ba84a140fa8c8682dbc5cc11baf047a33b222bdb8b0e1ce123e6529
output_sha256: 32117d44d900c97aa6066188d0deb924c93041c8a0bac3762a9a65594fa90682
independent_audit_script_sha256: 672d2b89b2609c77b0838b3b04f08289f65c8879ef84412d7f8adf0e26d10d3d
independent_audit_output_sha256: 310ffe235a1397bc8b335e967e12c291e8cf02d9143eaae6695f9feab8b5a263
hash_basis: raw LF bytes
primary_audit: >
  PASS. The exact source certificate specializes the complete Y-only source
  before elimination, proves the polynomial (A,B) gradient identities and
  Hessian bridge, excludes p=0, t=0, and source-infinity loss, obtains
  p^6 R_20 with parameter-uniform endpoint units, restores both universal
  fibres, and reconstructs every Newton face, the packet, cubic owner, and
  response inequality.
independent_audit: >
  ACCEPT. A separate source-coordinate referee uses the alternative pair
  (A,C_0), proves its localized Hessian bridge, obtains p^8 R_20 with the same
  endpoint units, checks a disjoint normalized projection, independently
  reconstructs the polygon and packet, and retains the cubic residue field
  as one prime owner. Normal, optimized, and fixed-hash streams byte-match.
---

# THM-4205 -- complete exact-weight-nine K-zero planar-Jacobian exclusion

**PROVED RELATIVE TO THM-4147/4176/4192/4209 + VERIFIED-EXACT +
INDEPENDENTLY SOURCE-AUDITED; JC(2) AND DC(2) REMAIN OPEN.**

## 1. Theorem, exact face, and inheritance

Work over `C` in the live exact-`M=9`, `b=d=0` reduced `(2,3)` seam. Put

```text
P=T+X^2T^2,                         Y=XTP,
K=2848/45-(7/6)Delta.
```

On `K=0`, necessarily

```text
Delta=5696/105!=0.                                      (1)
```

The complete normalized source is

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +Phi P^2Y+(5696/105)P^4+Theta P Y^2
  +eta P^3Y+zeta Y^3.                                  (2)
```

There is no ellipsis in `(2)`. Its only weight-nine terms are the last two,
so exact weight nine is equivalent to

```text
(eta,zeta)!=(0,0).                                      (3)
```

> **Theorem.** For every `Phi,Theta,eta,zeta in C` satisfying `(3)`, the
> source `(2)` is not the normalized source of a nonautomorphic planar Keller
> pair in the inherited reduced seam.

No critical-resultant discriminant, projected-root squarefreeness, or
separation at `T=-1/6` is assumed.

The inheritance pass is:

- closest proved mechanism: THM-4147's complete packet and prime-carrier
  finite/full response obstruction;
- canonical hostile: two distinct Morse source points may have one projected
  coordinate, so a projected resultant may repeat without source
  nonreducedness;
- corrected near miss: MISTAKE-519, where THM-4209 called an anti-diagonal
  specialization open despite THM-4176's older universal quantifier;
- least-used sidecar: the unspecialized source critical ideal together with
  the residue field of the pure cubic boundary point.

The connection contract is

```text
source:       full critical ideal on the rational (s,p) chart
target:       affine Morse fixed-sheet count and monodromy response
map:          source ideal -> saturated resultant -> restored coordinate fibres
preserved:    intersection length, Hessian reducedness, cubic residue field
destroyed:    s-coordinate and ordering of geometric carrier conjugates
sidecar:      source-infinity rows, p/t exclusions, universal fibres, owner field
decisive test: endpoint units plus the two Hessian bridges.               (4)
```

## 2. Exhaustive coefficient partition

Condition `(3)` has the following pairwise-disjoint exhaustive partition:

| row | coefficient condition | proved owner |
|---|---|---|
| `P` | `zeta=0, eta!=0` | THM-4192 |
| `Y` | `eta=0, zeta!=0` | Sections 3--6 below |
| `G` | `eta*zeta*(eta+zeta)!=0` | THM-4209 |
| `A` | `eta*zeta!=0, eta+zeta=0` | THM-4176 |

For row `A`, `zeta=-eta` and `(1)` gives `eta*Delta!=0`, exactly the complete
repeated-top scope of THM-4176. In particular its repeated-tangent
`Theta=-Delta` rows are already included. Thus only row `Y` requires new
work.

## 3. Complete Y-only source and two critical ideals

Set `eta=0`, `zeta!=0`, and use

```text
s=XT,                  p=T+s^2,                  t=p-s^2,

G=-s^2/(2t)+H,
H=-3p+(8/3)p^2-(1376/135)p^3+Phi*s*p^3
  +(5696/105)p^4+Theta*s^2*p^3+zeta*s^3*p^3.          (5)
```

Define

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2 H_p,
B=(C_0+sA)/t^2.                                       (6)
```

Direct cancellation makes `A,B` polynomials and gives

```text
t^2 G_s=pA,                    2t^2 G_p=t^2B-sA.       (7)
```

On the source critical ideal,

```text
p det D(A,B)=2t^2 det Hess(G)                 mod (A,B). (8)
```

The independent pair `(A,C_0)` instead satisfies

```text
t^2 G_s=pA,                    2t^2 G_p=C_0,
p det D(A,C_0)=2t^4 det Hess(G)              mod (A,C_0) (9)
```

in the localization `C[s,p,1/(pt)]`. The audit verifies `(9)` by expanding
the exact product-rule correction terms, not by silently cancelling `p` or
`t`.

For `(A,B)` there is no hidden point on either excluded divisor:

```text
p=0: A=-s, B=-6;
t=0: A=-s, and the only candidate s=p=0 has B=-6.      (10)
```

At source infinity the two independent pairs have leading rows

```text
(A,B):   (3zeta p^2,9zeta p^2),
(A,C_0): (3zeta p^2,6zeta p^2).                       (11)
```

Hence no finite `p!=0` point is lost at `s=infinity`.

## 4. Exact affine critical length without projection simplicity

Specializing the source first and eliminating `s` gives the two exact
identities

```text
Res_s(A,B)=p^6 R_20(p),
Res_s(A,C_0)=p^8 R_20'(p),                             (12)
```

with identical endpoint units

```text
R_20(0)=R_20'(0)=1,259,712 zeta^3,

[p^20]R_20=[p^20]R_20'
  =(40,870,620,168,192/1,225) zeta^7.                 (13)
```

They are nonzero for every `Phi,Theta` as soon as `zeta!=0`. Equations
`(10)--(13)` therefore give residual source intersection length `20`, counted
with multiplicity, uniformly on the whole Y-only face.

For a hypothetical Keller realization, the inherited target-Hessian
congruence and either `(8)` or `(9)` make every source critical point Morse.
Thus the source critical scheme is reduced even when several points share a
`p`- or `T`-coordinate. Projection discriminants play no role.

The rational chart collapses exactly the two inherited coordinate fibres.
The normalized audit restores

```text
T=0:    X^2=-6, G=0,   det Hess(G)=+6;
T=-1/6: X^2= 6, G=1/2, det Hess(G)=-6.                 (14)
```

There are two Morse points on each fibre, so the exact affine critical length
is

```text
L_Y=20+2+2=24.                                         (15)
```

## 5. Complete Y packet and the owner sidecar

For `Q=q^-1`, clearing the denominator of `G=q` gives

```text
F_Q(s,p)=(s^2-p)(1-QH)-Q s^2/2.                       (16)
```

Its Newton polygon and Pick data are

```text
Newt(F_Q)=conv{(0,1),(2,0),(5,3),(3,4),(0,5)},
(2Area,boundary,g_Pick)=(30,10,11).                   (17)
```

The four nonvertical primitive faces are

```text
s^2(1-Q/2)-p,
s^2(1-Q/2-Q*zeta*s^3p^3),
Q p^3s^3zeta(p-s^2),
Q p^4((5696/105)p+zeta*s^3).                          (18)
```

In particular `Phi,Theta` do not occur on a boundary face. The vertical edge
consists of affine `s=0` points and is not a source-puncture packet. The
complete nonvertical packet is

```text
e_Y=(11,8,2,2,2,1),
n_Y=sum e_Y=26,                 delta_Y=sum(e-1)=20.   (19)
```

As in THM-4147, the finite critical scheme and THM-3827 imply geometric
connectedness. The Newton genus is at most `11`, while the Keller residue
identity and defect `20` force genus at least `11`. Equality makes `(19)`
complete: there is no hidden affine ramification or unrecorded genus loss.

At `K=0` the unique nonrational boundary equation is

```text
zeta W^3=q-1/2.                                        (20)
```

It is separable over `C(q)` and irreducible because `q-1/2` has valuation one,
which cannot be the valuation of a cube. Thus `(20)` is one prime cubic
closed point, not three independently assignable geometric points.

THM-4147's finite-separable-carrier transport keeps this owner quantifier. If
the boundary image is finite, its residue field is intermediate in a prime
degree-three extension. THM-4120 rules out the rational finite alternative,
so all three conjugates form one degree-three horizontal carrier. If the
response is at the origin, the whole closed point responds there. No arbitrary
subset of conjugates is admitted.

## 6. Strict response contradictions

For a finite carrier response, removing its three index-two punctures gives

```text
n=20,                         beta=3.                  (21)
```

If at least one handle is nonidentity, the sharp merger capacity is

```text
2n-L_Y-1+beta=40-24-1+3=18<n-1=19.                   (22)
```

If both handles are identities, their total carrier index is only
`beta=3<19`. Both finite alternatives fail. For a full response, `n=26`, and
the commutator-overlap bound gives

```text
ind([X,Y])<=2(n-L_Y)=2(26-24)=4<delta_Y=20.            (23)
```

Thus row `Y` contains no nonautomorphic planar Keller pair. Together with the
exhaustive table in Section 2, this proves the theorem.

## 7. Boundary, nonconsequences, and replay

The hostile specializations are correctly typed:

```text
zeta=0, eta!=0:     P-only, already THM-4192;
eta=zeta=0:         exits exact weight nine;
eta*zeta!=0,
  eta+zeta=0:       complete anti-diagonal, already THM-4176.          (24)
```

The theorem closes the full exact-weight-nine `K=0` **coefficient face inside
the inherited reduced seam**. It does not prove entry into that seam, exclude
other reduced cells, treat `K!=0` outside the named earlier chambers, treat
`M>=10`, or prove `JC(2)` or `DC(2)`.

Replay the exact artifacts from the repository root:

```bash
python3 -B 04-computation/jc23_k0_y_only_complete_exclusion_thm4205.py \
  | diff -u 05-knowledge/results/jc23_k0_y_only_complete_exclusion_thm4205.out -
python3 -O -B 04-computation/jc23_k0_y_only_complete_exclusion_thm4205.py \
  | diff -u 05-knowledge/results/jc23_k0_y_only_complete_exclusion_thm4205.out -
PYTHONHASHSEED=4205 python3 -B \
  04-computation/jc23_k0_y_only_complete_exclusion_thm4205.py \
  | diff -u 05-knowledge/results/jc23_k0_y_only_complete_exclusion_thm4205.out -

python3 -B \
  04-computation/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.py \
  | diff -u \
  05-knowledge/results/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.out -
python3 -O -B \
  04-computation/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.py \
  | diff -u \
  05-knowledge/results/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.out -
PYTHONHASHSEED=271828 python3 -B \
  04-computation/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.py \
  | diff -u \
  05-knowledge/results/jc23_k0_y_only_complete_exclusion_independent_audit_thm4205.out -
```
