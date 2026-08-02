---
id: THM-3083
title: "Exceptional binary-point/ternary-direction S4 tomography clutch"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING IMMUTABLE INDEPENDENT
  HOSTILE AUDIT.  After choosing the exceptional identification
  AGL_2(F_2)=S4=PGL_2(F_3), the mean-zero four-point binary module is an
  S4-equivariant similarity of the sum-zero part of the four rank-one even
  ternary line channels.  Binary Walsh characters become signed 2+2
  contrasts of ternary directions, not individual lines.  Integrally the
  ternary channel lattice has quotient F3^3, while radial-plus-standard
  splitting adds Z/4; the composite Smith form is 1,3,3,12.  This is a
  representation/lattice clutch, not a physical quartic, tree, Keller,
  owner, or LRC intertwiner.
source: root-exceptional-two-three-tomography-2026-08-01
depends_on:
  - THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law
related:
  - THM-2387-degree-eighteen-h4-elliptic-three-isogeny-atlas
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2971-discriminant-cover-edge-orientation-sextic-algebra-intertwiner
  - THM-2996-prime-modular-affine-defect-trichotomy-and-spherical-quartic-uniqueness
  - THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch
  - THM-3067-tetrahedral-modular-two-three-flag-quotient-and-origin-loss
  - THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan
script: 04-computation/exceptional_binary_point_ternary_direction_s4_tomography_clutch_thm3083.py
output: 05-knowledge/results/exceptional_binary_point_ternary_direction_s4_tomography_clutch_thm3083.out
script_sha256: dcc5de253bdee3c9a9e21a8ea2b6441e86c667a2131d7b66700c02714708ea35
output_sha256: bf969955d11175a27bb448f7286c1622ac9f9372aee57dc58781d6493321a7ef
hash_basis: LF-normalized bytes
---

# THM-3083 -- exceptional binary-point/ternary-direction S4 tomography clutch

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING IMMUTABLE INDEPENDENT
HOSTILE AUDIT.**

## 1. Inheritance and exact claim

[THM-3076, finite affine-plane line
tomography](THM-3076-finite-affine-plane-line-quotient-tomography-and-p2-three-view-law.md)
proves that the `p+1` line averages on `F_p^2` split into a constant channel
and `p+1` centered direction channels.  At `p=2` the latter are three
one-dimensional Walsh channels.  At `p=3` they have dimension two, but parity
`x -> -x` cuts each of them to dimension one.  Simultaneously,

```text
AGL_2(F_2) ~= S4 ~= PGL_2(F_3),                            (1)
```

and both groups have a faithful four-point action.  This section proves the
precise point--direction clutch created by these two exceptional facts.

Let

```text
B=F_2^2,                         D=P^1(F_3).                (2)
```

After choosing an isomorphism in `(1)` and a compatible equivariant bijection

```text
phi:D -> B,                                                   (3)
```

there is an `S4`-equivariant similarity

```text
Phi:k[B]_0 -> H_0,                                           (4)
```

where `k[B]_0` is the mean-zero binary permutation module and `H_0` is the
sum-zero part of the four even ternary line channels.  Its Gram multiplier is
`18`.  A binary Walsh character maps to a signed `2+2` contrast of the four
ternary directions.  It does **not** map to one ternary line.

The integral clutch is not saturated.  Its exact successive defects are

```text
Z R direct_sum H_(0,Z)  subset  H_Z  subset  E_(Z,even)^0,

H_Z/(Z R direct_sum H_(0,Z))       ~= Z/4,
E_(Z,even)^0/H_Z                   ~= (F_3)^3,              (5)

Smith[(Z R direct_sum H_(0,Z)) -> E_(Z,even)^0]
                                      = diag(1,3,3,12).     (6)
```

Thus the rational module match exposes distinct binary and ternary integral
defects of total index `108`.  It is not a literal identification of the two
trees or a physical quartic decoder.

## 2. The exceptional four-point identification

Write

```text
B={00,10,01,11},
D={infinity,0,1,2}.                                         (7)
```

Use the following generators on `D`:

```text
S(z)=-1/z,             R(z)=-1/(z+1),             J(z)=-z. (8)
```

On `B`, take

```text
S_B(x)=x+(1,1),
R_B(x)= [[0,1],[1,1]] x,
J_B(x)= [[0,1],[1,0]] x+(1,1).                             (9)
```

Then `S_B,R_B` generate `A4=V4 semidirect C3`, `J_B` supplies an odd
reflection, and the three generators generate all of `S4`.  The explicit
bijection

```text
phi(1)=00,       phi(infinity)=10,
phi(0)=01,       phi(2)=11                                  (10)
```

intertwines `(8)` and `(9)`.  The companion enumerates the resulting
24-element group on both sides.  Once these generator pairings are fixed,
`(10)` is the unique equivariant bijection.

The choice matters.  Without an exceptional-group isomorphism there is no
absolute map from binary points to ternary directions.  A polarity is needed
only if one further identifies Fourier annihilator lines in the dual ternary
plane with primal lines; it is not needed for the physical line averages
below.

## 3. Even ternary line tomography

Let `V=F_3^2`, and work over a field `k` of characteristic different from
`2,3`.  For `L in D`, let `P_L` average translations along `L`, let `P_V`
average all of `V`, and put

```text
E_L=P_L-P_V.                                                (11)
```

Let

```text
W^+={f in k^V:f(-x)=f(x)}.                                 (12)
```

The full image of `E_L` consists of mean-zero functions on the quotient
`V/L`.  Parity identifies its two nonzero quotient classes, so

```text
im(E_L|W^+)=k h_L,             h_L=3 1_L-1.                (13)
```

Direct counting gives

```text
<h_L,h_M>=18 delta_(L,M),
sum_(L in D) h_L=9 delta_0-1=(8,-1,-1,-1,-1,-1,-1,-1,-1).
                                                                  (14)
```

Consequently

```text
W^+=k 1 direct_sum direct_sum_(L in D) k h_L.              (15)
```

The scalar ambiguity in a lift from `PGL_2(F_3)` to `GL_2(F_3)` is `+/-I`.
It acts trivially on `W^+`, so `PGL_2(F_3)=S4` genuinely acts on `(15)` and
permutes the four `h_L`.

Define

```text
H_0={sum_L a_L h_L:sum_L a_L=0},
Phi(f)=sum_(L in D) f(phi(L)) h_L.                          (16)
```

For `f,g in k[B]_0`, `(14)` yields

```text
<Phi(f),Phi(g)>=18<f,g>.                                   (17)
```

Equivariance follows from `(10)` and the permutation law for the `h_L`.
Both sides of `(16)` are the standard irreducible `[31]` module, so after
the choice `(3)` this intertwiner is unique up to scalar.

## 4. The three binary channels become matching contrasts

For nonzero `c in B`, put

```text
chi_c(x)=(-1)^(c dot x).                                   (18)
```

Its positive and negative fibres each have two points.  Transporting them by
`phi^(-1)` gives one of the three bipartitions of the four directions.  Thus

```text
Phi(chi_c)=H_c
 =sum_(chi_c(phi(L))=+1)h_L-sum_(chi_c(phi(L))=-1)h_L.     (19)
```

These are the three signed perfect-matching contrasts, not four individual
line channels.  This is the exact way the binary and ternary grammars share
one `S4` object.

There is also an exact quartic reading, but no new quartic invariant.  Label
the roots of a depressed quartic by the affine torsor `B`, and let `r:B->k`
be the root-value function.  Then `sum r=0`, and with normalized Walsh
coefficient

```text
rhat(c)=(1/4)sum_x chi_c(x)r(x),                            (20)
```

the matching-pair sum `s_c` satisfies

```text
rhat(c)=s_c/2,
Phi(r)=sum_(c!=0) (s_c/2)H_c,
u_c=s_c^2=4 rhat(c)^2.                                     (21)
```

The last quantity is exactly the cubic-resolvent root already owned by
[THM-2598](THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary.md).
So `(21)` linearly realizes its square-root channels, while squaring loses
their `V4` phase.  The pair

```text
T^4-7T^2+6T,                 T^4-7T^2-6T                  (22)
```

has matching-square multiset `{1,4,9}` in both cases and opposite square-root
orientation.  This is the sharp phase hostile.

## 5. Exact integral defects

Let `E_(Z,even)^0` be the even integer functions on `F_3^2` with total sum
zero.  Record their values `a_L` on the four nonzero antipodal pairs.  The
value at zero is forced to be

```text
f(0)=-2 sum_L a_L,                                         (23)
```

so these pair values identify `E_(Z,even)^0` with `Z^4`.

Let

```text
H_Z=direct_sum_L Z h_L.                                    (24)
```

If `c=(c_L)` is its coefficient vector, its four pair values are

```text
a=(3I_4-J_4)c.                                             (25)
```

The determinantal divisors of `3I_4-J_4` are

```text
1,1,3,9,27,                                                (26)
```

and hence

```text
Smith(3I_4-J_4)=diag(1,3,3,3),
E_(Z,even)^0/H_Z ~= (F_3)^3.                               (27)
```

Put

```text
R=sum_L h_L,
H_(0,Z)={sum_L c_Lh_L:sum_Lc_L=0}.                         (28)
```

In the coefficient lattice, the constant line plus the `A3` augmentation
lattice has index four.  In the pair-value coordinates `(23)`, columns for
`R,h_I-h_2,h_0-h_2,h_1-h_2` have determinantal divisors

```text
1,1,3,9,108                                                (29)
```

and Smith invariants `(1,3,3,12)`.  This proves `(5)--(6)`.  Equivalently,
on the saturated standard target `A3`, the map `(16)` sends the binary
integral augmentation lattice to `3A3`, with cokernel `(F_3)^3`.

This ternary standard defect is different from THM-3045's single trivial
`F3` defect in the six-edge lattice.

## 6. Why this coincidence stops at two and three

At `p=2`, all points are parity-fixed, there are three directions, and every
centered line channel has rank one.  For odd `p`, the even part of a centered
line channel has rank

```text
(p-1)/2.                                                   (30)
```

It is rank one exactly at `p=3`, when there are four projective directions.
For `p>3`, both the rank and the number `p+1` of directions are too large for
the four-point `S4` clutch.  Thus `(1)--(4)` isolate a genuine exceptional
two/three phenomenon rather than the first case of an all-prime tree
identification.

## 7. Stopping boundaries

The theorem deliberately retains the following losses.

1. The radial vector `R=(8,-1^8)` is nonzero.  Four ternary line channels are
   not three binary directions; only their sum-zero standard part matches.
2. Binary characters map to signed `2+2` combinations, not to single lines.
   A line-by-line tomography intertwiner has the wrong type.
3. The map is linear, not multiplicative.  Its natural extension sends a
   binary point basis vector `delta_(phi(L))` to `h_L`; but
   `delta_(phi(L))^2=delta_(phi(L))` while `h_L^2!=h_L`.
4. THM-2387 has the four `F3` projective lines and their Weil switching class,
   but no Galois-equivariant matching to quartic sheets or Cardano roots.  If
   such a common torsor were supplied, equality of sign characters would be a
   conditional tautology, not a new discriminant identity.
5. Nothing here selects an affine origin, a sign of `s_c`, a Feuerbach point,
   Farey flank, physical word/current, owner, Keller graph order, or LRC
   carrier.  THM-2971's primitive discriminant-cover coordinate is not
   reconstructed.

## 8. Exact evidence

Run

```bash
python 04-computation/exceptional_binary_point_ternary_direction_s4_tomography_clutch_thm3083.py
python -O 04-computation/exceptional_binary_point_ternary_direction_s4_tomography_clutch_thm3083.py
```

Both executions byte-match the stored transcript.  The companion uses only
explicit `require` gates.  It checks all 24 paired group elements and the
unique generator-fixed bijection; the complete ternary projector algebra on
the even basis; every group/channel action; Gram similarity and the three
Walsh matchings; both Smith computations; the radial hostile; and both
quartic phase controls.

```text
PROVED HERE:       chosen exceptional S4 point--direction identification;
                   even p=3 line-channel decomposition and Gram;
                   S4 similarity of the two standard modules;
                   binary Walsh to signed matching-contrast formula;
                   exact index-27, index-4, and Smith-108 lattice defects;
                   quartic square-root reinterpretation and phase hostile.

NOT PROVED:        canonical identification without a chosen S4 gauge;
                   linewise or multiplicative tomography intertwiner;
                   new quartic/resolvent/discriminant information;
                   physical Cardano, Weil, owner, current, or tree carrier;
                   Keller, JC(2), DC(2), GMC, or LRC(14).                  (31)
```

QED.
