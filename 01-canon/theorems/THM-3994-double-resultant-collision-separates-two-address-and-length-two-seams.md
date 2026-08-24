---
id: THM-3994
title: "Double resultant collisions separate two-address and length-two seams"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the first-height
  constant-(c,r) affine-P graph
  of THM-3972, the two discriminant seams have the same scalar symptom but
  different completed base geometry. At cr=3 the double resultant root is
  supported at two distinct transverse basepoints in one fibre. At 4cr=3 it
  is supported at one curvilinear length-two base scheme with completed local
  ideal (L,V^2). The latter Rees graph has an A1 chart V^2=LZ and local class
  group Z/2. Resultant multiplicity therefore records total intersection
  length but not its distribution among addresses or the singularity of the
  completed graph.
source: root + frontier_transfer_scout / completed-zero-locus seam, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (thm3994_hostile_audit, 2026-08-24). Every
  resultant, discriminant, factorization, Jacobian, scheme-length, Rees-chart,
  normality, local-class-group, scope, hash, and dependency claim was rebuilt.
  The displayed inverse from (D,V) to (t,v) was added so the Z/2 claim concerns
  the algebraic local ring itself, not merely its completion. Normal and
  optimized companion executions are byte-identical.
depends_on:
  - THM-3969-affine-p-graph-debt-relative-p1-normalization
  - THM-3972-simple-collision-affine-p-graph-blowup-normalization
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3990-componentwise-harmonic-obstruction-and-repair-quotient
script: 04-computation/jc2_affine_p_double_collision_thm3994.py
output: 05-knowledge/results/jc2_affine_p_double_collision_thm3994.out
script_sha256: d77ef7469f6c705b9fece263dd9f6a07892d2960625d0ddb4b94c18a8c0cbed0
output_sha256: 240b2cdc0c6b7c5e209fbf9c17c2e7fb8d4e40fdd54f1c51a97e8b84f194d213
hash_basis: raw LF bytes
---

# THM-3994 -- a double resultant does not determine the collision geometry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. In the first-height
constant-`(c,r)` row of THM-3972, assume `r!=0` and put

```text
D=1-tv^3,                    N=3r+3v+cv^2.              (1)
```

These are the pole and numerator rows on the finite chart of
`S=A1_t times P1_v`. Their resultant is

```text
Xi=Res_v(D,N)=c^3+27r^3t^2+27(1-cr)t,                  (2)
disc_t(Xi)=-27(cr-3)^2(4cr-3).                         (3)
```

The squarefree theorem excluded the two seams in `(3)`. They now separate
exactly as follows.

## 1. The base ideal is the completed zero-locus sidecar

In THM-3972 the rational graph map is represented by `[P_h:D_h]`, with

```text
P_h=r^2D_h+VN_h.                                       (4)
```

Both seams below have `c!=0`, so `[V:W]=[1:0]` is not a basepoint: there
`N_h=cV^2!=0`. Every finite common zero used below also has `v!=0`.
Consequently `(4)` gives the equality of local base ideals

```text
(P_h,D_h)=(N_h,D_h).                                   (5)
```

Thus the scheme-theoretic common zero of `D,N`, rather than the scalar
resultant alone, determines the local graph closure through its Rees algebra.

## 2. The seam cr=3 consists of two reduced addresses

Set `c=3/r`. Then

```text
Xi=27(r^3t-1)^2/r^3.                                   (6)
```

At its unique collision parameter `t_0=1/r^3`, put

```text
Q=v^2+rv+r^2.
```

Direct factorization gives

```text
N=(3/r)Q,
D|_(t=t_0)=-(v-r)Q/r^3.                                (7)
```

The discriminant of `Q` is `-3r^2`, so `Q` has two distinct roots. Neither
root is zero or `r`. At either root the Jacobian of `(D,N)` in `(t,v)` is

```text
J_b=D_t N_v=-v^3(3+6v/r),                              (8)
```

and it cannot vanish: `2v+r=0` is incompatible with `Q=0`. Equivalently,
the norm of `(8)` modulo `Q` is `27r^6!=0`. The base scheme over `t_0`
therefore consists of two distinct reduced transverse points. Locally the
graph is the ordinary blowup at each point and is smooth there.

Hence the order-two zero in `(6)` is the sum of two intersection lengths
`1+1`, carried by two different `v`-addresses in the same fibre.

## 3. The seam 4cr=3 is one curvilinear length-two centre

Set `c=3/(4r)`. Then

```text
Xi=27(8r^3t+1)^2/(64r^3).                              (9)
```

Its common zero is supported at

```text
t_0=-1/(8r^3),                 v_0=-2r.                (10)
```

Write `T=t-t_0` and `V=v-v_0`. Exact expansion gives

```text
N=3V^2/(4r),
D=8r^3T+(3/(2r))V + terms of total order at least 2.   (11)
```

Since

```text
det partial(D,V)/partial(T,V)=8r^3!=0,                 (12)
```

`L=D,V` are actual algebraic local coordinates. Indeed `V=v+2r`, so at the
point `V-2r=v` is a unit, and the inverse change is

```text
v=V-2r,                    t=(1-L)/(V-2r)^3.            (13)
```

Equations `(5)` and `(11)` therefore identify the local base ideal and its
completed quotient as

```text
I=(L,V^2),
k[[L,V]]/I = k[[V]]/(V^2).                             (14)
```

Thus the double root in `(9)` is one curvilinear scheme of length two, not two
reduced addresses.

The local graph closure is `Proj Rees(I)`. On the chart where the homogeneous
coordinate of `L` is nonzero it has equation

```text
V^2=LZ.                                                (15)
```

This chart has one ordinary `A1` singularity at `L=V=Z=0`; its other Rees
chart is smooth. The explicit inverse `(13)` identifies the algebraic graph
local ring itself with

```text
k[L,V,Z]_(L,V,Z)/(V^2-LZ).
```

This ring is normal, and the usual quadratic-cone or toric divisor computation
gives

```text
Cl(O_(graph,0)) = Z/2.                                 (16)
```

Equivalently, `(15)` is the invariant ring
`k[s^2,st,t^2]=k[s,t]^(mu_2)`, and its minimal resolution has one exceptional
curve with intersection matrix `[-2]`, whose Smith cokernel is `Z/2`.

## 4. Exact survivor and scope

Both seams have an order-two scalar resultant and total base intersection
length two. They differ in every sidecar relevant to normalization:

```text
cr=3:   two supports, lengths 1+1, smooth graph charts;
4cr=3:  one support, length 2, one normal A1 graph singularity.            (17)
```

Therefore neither resultant order nor the number of collision fibres is a
normalization-place count. A global analysis of either seam must retain the
labelled exceptional, pole, and ramification incidences before taking a class
group or Smith quotient.

The local `Z/2` in `(16)` is not a global class-group computation and gives
no `3`-torsion obstruction by itself. This theorem does not resolve the
global seam, construct a Keller map, or prove `JC(2)`. It closes only the
scheme-theoretic local geometry hidden by the double resultant. **QED
candidate.**

## Reproduction

```bash
python3 04-computation/jc2_affine_p_double_collision_thm3994.py
python3 -O 04-computation/jc2_affine_p_double_collision_thm3994.py
sha256sum 04-computation/jc2_affine_p_double_collision_thm3994.py \
  05-knowledge/results/jc2_affine_p_double_collision_thm3994.out
python3 agents/check_docs.py
```
