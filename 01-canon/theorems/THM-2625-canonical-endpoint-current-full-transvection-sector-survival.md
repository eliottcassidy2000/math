---
id: THM-2625
title: "Canonical endpoint current: full transvection-sector survival"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the explicit
  canonical typed non-cover control of THM-2309/THM-2349, with owner c1,
  word {a}, clock two, and marked triangle (X,m,Y)=(13,1,742599), the
  THM-2334 marked current in THM-2620's left-deep endpoint allocation is
  nonzero in all 28,561 ordered endpoint cells.  Every one of its 2,185
  admissible (q,Delta) determinant-sector aggregates is nonzero, including
  all 2,016 nondegenerate parabolic sectors; the only empty formal sectors
  are the structurally impossible q=0, Delta!=0 cells.  Both endpoint DFTs
  and all target marginals have full 169-element support.  This is an exact
  dual-field cyclotomic certificate on one complex unrestricted mod-thirteen
  current, not positivity, chronology, a covering-row theorem, all-91-unit
  survival, seven-clock gluing, a row exclusion, or LRC(14).
source: codex-2026-07-28-canonical-endpoint-current
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
related:
  - THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice
script: 04-computation/lrc14_canonical_endpoint_current_thm2625.py
output: 05-knowledge/results/lrc14_canonical_endpoint_current_thm2625.out
script_sha256: eb89eb4753f95b00ba902e1ff7c691326da62ae52501817180b81c13ed6bc62f
output_sha256: 69b809444259b31e0c80ce0e4189505de455be1bd4f220b0500c00dced0f491f
hash_basis: working-tree bytes (LF)
---

# THM-2625 -- the canonical current occupies every endpoint transvection sector

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2620 proves that an allocated endpoint pair `(L,R)` contains a determinant
coordinate erased by the target difference `q=L-R`, and that every
`q!=0, Delta!=0` fibre is a thirteen-edge parabolic graph.  That theorem is
structural: it does not say that a THM-2334 current survives in any such
fibre.  This theorem supplies the missing existence statement on the
canonical exact control.  Unexpectedly, the current does not merely hit one
fibre: it fills every endpoint cell and every admissible determinant sector.

The conclusion remains coefficient-side and same-event.  It is a precise
bridge from the relation current to the parabolic grammar, not the physical
chronological intertwiner still required by LRC(14).

## 1. Typed control and marked triangle

Use the canonical THM-2309 typed row

```text
W=(1,14,27,40,53,66,13,2197,742586)
  =(H,q1,q2,q3,q4,q5,c1,c2,c3),                         (1)

nu_13(c1,c2,c3)=(1,3,5).
```

It is a typed **non-cover** control; no scalar covering hypothesis is made.
Take owner `c1`, target coordinates `a=c2`, `b=c3`, the first positive
THM-2349 word `sigma={a}`, and delayed clock

```text
k=2,                         R=13^2=169.                 (2)
```

The exact Boolean sets have measures

```text
mu(E_1)=1882176/28589561,
mu(Q_(1,{a}))=143103830843/5727632650740,
mu(E_1 intersect T^(-2)Q_(1,{a}))=21376087/17907461390>0. (3)
```

Use the marked frequency triangle

```text
(X,m,Y)=(13,1,742599),             Y=X+m c3.             (4)
```

Both `X` and `Y` have 13-adic valuation one and `gcd(m,91)=1`.

## 2. Split the THM-2334 coefficient before taking target difference

Let

```text
Ghat=L_13^perp/<W>,                       |Ghat|=169.     (5)
```

The canonical quotient representatives are

```text
v1=(0,12,0,0,0,0,0,1,0),
v2=(0,0,12,0,0,0,0,0,1),
ell(alpha,beta)=alpha v1+beta v2,                       (6)

tau(ell)=(ell_a,ell_b)=(alpha,beta) in F_13^2.
```

Thus `tau` is a bijection.  In the cyclotomic ring containing all rational
endpoint phases define

```text
P_ell=zeta_13^(m ell_b) (2 pi i X) Ehat_Q^ell(X),

Q_ell=conj((2 pi i Y) Ehat^ell(Y)).                     (7)
```

These are cyclotomic integers: multiplication by the displayed frequencies
turns each interval Fourier coefficient into a finite signed endpoint sum.
They descend **separately** through the representative gauge.  Under
`ell -> ell+sW`, the left factor in (7) gains
`zeta_13^(s(X+m c3))=zeta_13^(sY)=1`, while the conjugate right factor gains
`zeta_13^(-sY)=1`, since `13|Y`.  This is stronger than descent of their
product alone.

Define unnormalized inverse endpoint transforms

```text
Lstar(l)=sum_ell P_ell zeta_13^(-<tau(ell),l>),

Rstar(r)=sum_ell Q_ell zeta_13^(+<tau(ell),r>).          (8)
```

The plus sign on `Rstar` is load-bearing: `Q_ell` is the negative-character
transform produced by conjugation.  In THM-2620's left-deep allocation,

```text
L=pi(u+D beta+m e_c3),                  R=pi(v),

q=L-R.                                                   (9)
```

Writing

```text
C0=delta_hat(m)/(4 pi^2 X Y) != 0,
```

the complete allocated joint current is

```text
J(l,r)=C0 Lstar(l)Rstar(r)/169^2.                       (10)
```

Equation (10) is exact after fixing this allocation gauge; it is not an
intrinsic physical endpoint pair.

## 3. Determinant-sector refinement and Radon compatibility

For `q=l-r` and `Delta=det(l,r)=det(q,r)`, put

```text
Sstar(q,Delta)
 =sum_(r: det(q,r)=Delta) Lstar(r+q)Rstar(r).            (11)
```

Also define the unnormalized THM-2334 target aggregate

```text
Anum(q)=sum_ell P_ell Q_ell
                 zeta_13^(-<tau(ell),q>).               (12)
```

Fourier orthogonality gives, identically and before computation,

```text
sum_Delta Sstar(q,Delta)
 =sum_r Lstar(r+q)Rstar(r)
 =169 Anum(q).                                          (13)
```

After multiplying by `C0/169^2`, equation (13) is exactly the THM-2334
target marginal

```text
A(q)=C0 Anum(q)/169.                                    (14)
```

Thus this refinement commutes with the lawful present/bare-endpoint Radon
pushforward; no new marginal has been substituted for the canonical one.

## 4. Exact full-support theorem

For the control (1)--(4),

```text
support(Lstar)=169,                    support(Rstar)=169,

support(J)=13^4=28,561,                support(Anum)=169. (15)
```

Among the `169*13=2,197` formal `(q,Delta)` labels, exactly twelve have empty
index sets:

```text
q=(0,0),                         Delta=1,...,12.          (16)
```

Every other sector sum in (11) is nonzero:

```text
support(Sstar)=2,185.                                      (17)
```

In particular,

```text
Sstar(q,Delta)!=0
   for every q!=0 and every Delta in F_13,               (18)
```

and all

```text
168*12=2,016                                             (19)
```

nondegenerate `q!=0,Delta!=0` parabolic fibres of THM-2620 carry nonzero
current.  Each such fibre has thirteen endpoint edges; all `26,208` of these
individual coefficients are nonzero because (15) gives full joint support.
The degenerate admissible sector `(q,Delta)=(0,0)` also survives.

The separate computation of (12) proves every `A(q)` nonzero.  This does not
follow merely from (17), since nonzero determinant sectors could cancel when
summed over `Delta`; both statements are checked independently and (13) is a
control tying them together.

## 5. Exact characteristic-zero certificate

All unnormalized quantities in (7)--(12) lie in

```text
Z[zeta_N],                   N=50,334,435,734,703,120.   (20)
```

The companion constructs two exact specializations

```text
p1=352341050142921841,          zeta_N -> h1=435817657216,

p2=956354278959359281,          zeta_N -> h2=153943385426666320. (21)
```

Explicit Lucas certificates prove both `p_i` prime from the complete
factorization of `p_i-1`; exact-order tests prove `h_i` has order `N`.
Therefore each map in (21) is a ring homomorphism
`Z[zeta_N] -> F_(p_i)`.  A nonzero image rigorously proves that the original
cyclotomic integer is nonzero; no floating-point inference is used.

Both fields give the support counts (15)--(19).  With coordinates ordered
lexicographically, the SHA-256 digests of `(Lstar,Rstar,Sstar)` are

```text
p1:
33d77c87c2970bf99bd1f0f139e74b8eef70677ba28844391f6d681ebf2e6c69
755b046e39d71fbe83295dd7e43727eaa30a5220f359e2d389268d124e6a8f27
d7e20ba35d1b25c1761b03c297f2449d74b29be0defbc7496b10eff8b9458b41

p2:
97471c735336996fab845da55ebbca7a19786297ebe5dd9f7f5026a0a816c9ee
07bc100112b3cf83cfddec6a6354df5caf7d5a12edd7f9d10c921ff08153271b
bb09105d839e5abdfadf7b4e1959bcc8e5a4af680a8dbe21ec485ab7bb7ec8f2. (22)
```

The companion also checks every inverse DFT, separate `P/Q` representative
descent on hostile gauges, all 169 identities (13), the structural sector
histogram, full joint support, and the exact digest serialization.  Normal
and optimized Python runs byte-match the stored output.  An independent
referee rebuilt the DFT and sector loop, reran both fields, recovered every
count and identity, and audited the allocation and conjugation signs.

## 6. Exact boundary of the result

This theorem proves survival of one actual marked-current coefficient bank in
the complete abstract endpoint plane.  It proves none of the following:

```text
nonnegative or sign-coherent endpoint weights;
all-91-unit projected survival;
uniformity over typed rows or any covering row;
cross-clock ancestry or endpoint gluing;
an endpoint-to-physical-root intertwiner;
a seven-clock projective closure;
a row exclusion or LRC(14).
```

The current is complex, unrestricted mod thirteen, same-event, and tied to
the allocation (9).  The next LRC obligation is therefore no longer “find a
nonzero parabolic endpoint sector” on this control.  It is to preserve one of
the now-proved surviving sectors through positivity, chronology, semantic
root identification, and adjacent-clock ancestry without taking the Radon
quotient that erases its endpoint origin.

QED.
