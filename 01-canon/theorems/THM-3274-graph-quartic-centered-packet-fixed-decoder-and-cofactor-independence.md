---
id: THM-3274
title: "Graph-quartic centered-packet fixed decoder and cofactor independence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the actual fixed-plus-escaping-C3 graph-quartic
  anatomy of THM-3201, rescaling the THM-3271 centered-norm packet by the
  base-field scalar pi^(3m) makes its reduction equal to
  (W-A^3)(W-7A^3/27)^3 when 3 does not divide m, and to
  W(W+2A^3/27)^3 when 3 divides m.  Thus in both Newton lanes the fixed
  packet component is the unique simple residual factor and is recoverable,
  without a sheet mark, from the scalar packet characteristic polynomial.
  The fixed spectral denominator has respective leading constants
  (20A^3/27)^3 and (2A^3/27)^3.  The THM-3273 covariant has normalized
  residue zero in the first lane and 27A^12/8 in the second: C sees the
  moving packet collision, not fixed/moving separation.  Conversely, an
  exact norm-one cofactor gauge preserves the quartic, packet, collision
  covariant, fixed projector, total cofactor norm, and cofactor valuations
  while changing the pointed physical-Jacobian ratio from 1 to lambda^4.
  Hence the centered packet now recovers the fixed scalar but cannot recover
  the true chain-rule cofactor or Keller equality.  This is a local decoder
  and an information-independence theorem, not a C3 exclusion or JC(2)
  theorem.
source: jc-graph-packet-cofactor-2026-08-02
audit: >
  An independent audit rederived both depressed leading quartics, packet
  factors, normalized projector denominators, and covariant residues from
  THM-3201's actual root asymptotics.  It also checked the tame cyclic
  hostile's derivative stars, norm-one gauge, unchanged invariant list, and
  lambda-four pointed ratio.  Fresh normal and optimized replays byte-match.
depends_on:
  - THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate
  - THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector
  - THM-3273-quartic-centered-norm-packet-collision-covariant-and-discriminant-factor
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
related:
  - THM-3272-tame-pure-c3-centered-norm-packet-integral-fixed-projector
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
script: 04-computation/jc_graph_quartic_packet_fixed_decoder_cofactor_hostile_thm3274.py
output: 05-knowledge/results/jc_graph_quartic_packet_fixed_decoder_cofactor_hostile_thm3274.out
script_sha256: c265202808047c5126a99bf4465ff2ad5aa7135e5b55b5f902986449348dc164
output_sha256: 774659e0d1c77f6bf15df295526895168024ff61680f3c00a04418c312f08993
hash_basis: LF-normalized bytes
---

# THM-3274 -- the graph packet finds the fixed scalar, but not its cofactor

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and exact scope

Retain the actual fixed-plus-escaping-`C3` graph-quartic setup of
[THM-3201](THM-3201-c3-local-resolvent-splitting-and-matching-newton-gate.md),
Section 6.  Thus the completed splitting field has a Kummer uniformizer
`pi`, the fixed source sheet is finite, and the three moving source sheets
have

```text
X_i=A zeta^(-im) pi^(-m)+O(pi^(-m+1)),       m>=1,      (1)
```

where `A` is a nonzero residue-field element.  The base field contains
`pi^3`; in THM-3201's notation,

```text
s=pi^(3m)=(t/tau)^m is in K*.                            (2)
```

Let

```text
f(X)=X^4+pX^2+qX+r                                      (3)
```

be the depressed graph quartic, let `xi=X mod f`, and use the centered-norm
packet of [THM-3271](THM-3271-universal-quartic-centered-norm-packet-and-local-singleton-projector.md):

```text
N(X)=-(20X^3+18pX+27q)/27,
nu=N(xi),
P_N(Z)=Norm_(K[X]/(f)/K)(Z-nu).                         (4)
```

The normalized scalar characteristic polynomial

```text
P_tilde(W)=s^4 P_N(s^(-1)W)                             (5)
```

lies in `K[W]`.  Its roots are the four integral normalized packet values
`sN(z_j)`.  The theorem identifies its complete leading reduction and asks
two separate questions:

1. does the unmarked scalar polynomial recover the fixed packet component?;
2. does that recovery also identify the physical Jacobian cofactor attached
   to the fixed component?

The answer is respectively **yes** and **no**.

## 2. Off resonance: the covariant vanishes but the fixed gap survives

Suppose `3` does not divide `m`.  The moving leading trace in `(1)` is zero.
After depression and multiplication of the roots by `pi^m`, the four leading
roots are

```text
0, A, A zeta, A zeta^2.                                (6)
```

Their monic polynomial and centered packet map are

```text
f_0(X)=X(X^3-A^3)=X^4-A^3X,
N_0(X)=A^3-(20/27)X^3.                                 (7)
```

The first root in `(6)` is the fixed sheet.  Hence its packet value is
`A^3`, whereas all three moving roots have packet value `7A^3/27`.  It
follows that

```text
bar(P_tilde)(W)
  =(W-A^3)(W-7A^3/27)^3,                               (8)

fixed-minus-moving gap =20A^3/27.                      (9)
```

In particular the fixed factor is simple even though the complete packet
reduction is not squarefree.

This lane also separates two meanings of "packet collision."  Substitution
of `(7)` in the collision covariant of
[THM-3273](THM-3273-quartic-centered-norm-packet-collision-covariant-and-discriminant-factor.md)
gives

```text
bar(s^4 C(p,q,r))=C(0,-A^3,0)=0.                       (10)
```

Equivalently, `w(C)>-12m`.  The zero in `(10)` records the three equal
**moving** packet initials.  It does not record a fixed/moving collision,
because `(9)` is nonzero.  Thus `C` being a unit is sufficient for total
packet separability but is not necessary for decoding the fixed component.

## 3. On resonance: the tame singleton model reappears

Suppose `3|m`.  The three moving initials in `(1)` coincide.  Depression
turns the four scaled leading roots into

```text
-3A/4, A/4, A/4, A/4.                                 (11)
```

Their polynomial is

```text
f_1(X)=(X+3A/4)(X-A/4)^3
      =X^4-(3A^2/8)X^2+(A^3/8)X-3A^4/256.             (12)
```

The fixed value of `N_1` is zero and each moving value is `-2A^3/27`.
Consequently

```text
bar(P_tilde)(W)=W(W+2A^3/27)^3,                        (13)

fixed-minus-moving gap =2A^3/27,                       (14)

bar(s^4 C(p,q,r))=27A^12/8.                            (15)
```

This is exactly the separated singleton/triple geometry behind THM-3272,
now derived at the pole scale of the actual THM-3201 graph quartic.

## 4. The unmarked fixed decoder in both lanes

Let `n_0=N(z_f)` be the packet value of the finite fixed sheet, and let

```text
H(Z)=product_(i=0)^2 (Z-N(z_i)),
delta_0=H(n_0).                                        (16)
```

Equations `(8)` and `(13)` say that `bar(P_tilde)` has a unique simple
factor and one coprime cubic factor.  Hensel factorization therefore lifts
that linear factor uniquely.  The lifted root is `s n_0`, so `n_0` is
recoverable from `P_N` without retaining a sheet mark.  The corresponding
quartic-algebra projector is the THM-3271 spectral idempotent

```text
e_f=H(nu)/H(n_0).                                      (17)
```

Moreover the normalized denominator has the exact leading residues

```text
bar(s^3 delta_0)=
  (20A^3/27)^3 =8000A^9/19683,   if 3 does not divide m;
  ( 2A^3/27)^3 =   8A^9/19683,   if 3 divides m.        (18)
```

Both are nonzero.  Thus `(17)` is integral in the packet order after the
base-field scaling `(2)`.  This repairs a genuine information question left
open by the earlier scalarization ledger: for the actual graph asymptotics,
the fixed **packet scalar** is not the missing coordinate.

The repair is deliberately narrow.  The scalar packet need not recover the
fixed source root itself, and `(17)` does not attach an independently given
branchwise cofactor to the decoded idempotent.

## 5. Exact norm-one cofactor hostile

The latter omission is intrinsic.  Let `k` be a characteristic-zero field
containing `zeta_3`, put

```text
K=k((t)),             L=K[y]/(y^3-t),
g(X)=X^3-t,           f(X)=(X-1)g(X).                  (19)
```

The cubic is irreducible by its `t`-valuation and `L/K` is a tame cyclic
`C3` extension.  The fixed/cubic derivative stars from
[THM-3064](THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md)
are

```text
D_0=g(1)=1-t,
D_C=(y-1)g'(y)=3y^2(y-1).                              (20)
```

Begin with the Keller-compatible abstract cofactor packet

```text
c_0=D_0^(-1),                   c=D_C^(-1),
J_0=c_0D_0=1,                  J_C=cD_C=1.             (21)
```

For any base unit `lambda in k*`, make the equivariant twist

```text
c_0' =lambda^(-3)c_0,           c'=lambda c.            (22)
```

It preserves the total four-sheet cofactor norm because

```text
c_0' Norm_(L/K)(c')
 =lambda^(-3)lambda^3 c_0 Norm_(L/K)(c)
 =c_0 Norm_(L/K)(c).                                   (23)
```

It also preserves every cofactor valuation and the product of the four
physical Jacobian values.  The quartic `f`, its centered packet, `P_N`, `C`,
the fixed projector, and all derivative stars are entirely unchanged.  For
reference, depression of `(19)` gives the packet-collision covariant

```text
C=-27(125t-1)(8000t^2-475t+8)/64,
bar(C)=27/8,                                             (24)
```

so this hostile is safely off the packet-collision boundary.

Nevertheless `(22)` changes the physical packet to

```text
J_0'=lambda^(-3),               J_C'=lambda,
rho'=J_C'/J_0'=lambda^4.                               (25)
```

The pointed Keller decoder is therefore

```text
Delta_J'=Norm_(L/K)(rho'-1)=(lambda^4-1)^3.             (26)
```

At `lambda=2`, all invariants listed above agree with `(21)`, while

```text
rho'=16,                       Delta_J'=15^3=3375 !=0.  (27)
```

Thus even the augmented data

```text
(f, N, P_N, C, fixed spectral projector,
 derivative stars, total cofactor norm, cofactor valuations)              (28)
```

do not determine sheetwise Keller equality.

## 6. Exact frontier consequence

The graph-packet ledger is now sharper:

```text
unmarked scalar centered packet
  -> fixed packet scalar and its idempotent:       recovered in both lanes;

fixed packet scalar + symmetric cofactor product
  -> true physical fixed/cubic Jacobian ratio:     not recovered.          (29)
```

The least remaining sidecar is the true chain-rule cofactor in the decoded
`1+3` decomposition, equivalently THM-3064's pointed element `rho` or shifted
norm `Delta_J`, together with whatever affine-owner incidence proves that
this cofactor comes from the same polynomial graph.  Another symmetric root
covariant cannot supply that information: `(22)` fixes all such root data.

The hostile does **not** claim that both abstract cofactor packets are
realized by polynomial maps.  Its role is information-theoretic: any future
JC argument using only `(28)` has forgotten a coordinate on which Keller
equality genuinely depends.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc_graph_quartic_packet_fixed_decoder_cofactor_hostile_thm3274.py
python3 -O 04-computation/jc_graph_quartic_packet_fixed_decoder_cofactor_hostile_thm3274.py
```

Both modes must reproduce
`05-knowledge/results/jc_graph_quartic_packet_fixed_decoder_cofactor_hostile_thm3274.out`
byte for byte.  The companion independently checks both normalized packet
factorizations, both fixed/moving gaps and projector denominators, the two
covariant residues, the norm-one gauge, the `lambda^4` ratio, and the exact
`lambda=2` shifted norm.
