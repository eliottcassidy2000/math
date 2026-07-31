# Literal absent/present allocation on the retained LRC sheet

Status: **scratch exact result; no tracked canon change**.

## Verdict

There is a literal carrier-absent amplitude on the selected labelled Boolean
ancestry sheet.  The factor `13` in its central harmonic is derivable, but it
must not be its definition.

The resulting same-sheet Boolean square is nonzero in both certified fields
and has a nonzero coefficient mixed face.  It does **not**, however, induce
the endpoint steps `(0,1)` and `(12,0)`.  Its intrinsic endpoint steps are
both zero.  Consequently its determinant mixed face is zero, while the
determinant-one parallelogram in the pending witness is a separate
transported-sheet/carried-sample control.

## 1. The object that must stay fixed

Let

```text
I=[142004992589460,142005019034340)
J=I+431933040
```

be the selected source sheet and its target copy.  Let `A_a` and `B_b` be
the source and target one-sided physical carriers at primal carrier twists
`a,b in F_13`.  The twist unit is

```text
u=T/13=22910530602960.
```

Rebuilding the selected cell `(s,t,e)=(0,4,1)` gives `240` source support
pieces and `241` target support pieces.  Exact half-open intersection with
the full one-sided carriers gives

```text
I intersect A_0 = I,       I intersect A_a = empty for a!=0,
J intersect B_0 = J,       J intersect B_b = empty for b!=0.       (1)
```

This is stronger than the corresponding statement for the translated unit
atom: `(1)` uses the full physical one-sided carriers.

The pending carried bank instead evaluates the translated sheets `I_a` and
`J_b` themselves.  That replacement no longer keeps the selected ancestry
copy fixed and therefore cannot define a same-sheet allocation toggle.

## 2. Direct omitted-carrier construction

For an endpoint address representative `ell`, let

```text
p_ell = Phi_L(1_I times E_ell),
q_ell = Phi_R(1_J times E_ell),                              (2)
```

where `Phi_L,Phi_R` are the exact `x_sweep` and endpoint-sum functionals and
`E_ell` is the inherited endpoint mask.  Equation `(2)` is evaluated
directly from the interval machinery with no carrier factor and no carrier
twist variable.

To compare it with a carried bank, extend the absent value constantly over
the omitted primal twist:

```text
p_abs(a)=p_ell,                 q_abs(b)=q_ell.              (3)
```

For the literal present state, retain the sheet selector and multiply by the
carrier indicator.  Equation `(1)` gives

```text
p_pre(a)=1_(a=0) p_ell,         q_pre(b)=1_(b=0) q_ell.      (4)
```

Under the repo's unnormalized thirteen-point Fourier convention,

```text
hat(p_abs)(k)=13 p_ell 1_(k=0),   hat(p_pre)(k)=p_ell,
hat(q_abs)(l)=13 q_ell 1_(l=0),   hat(q_pre)(l)=q_ell.       (5)
```

Thus at the central harmonic `k=l=0`, after the separate endpoint-address
DFTs, the laws hold pointwise on the full endpoint planes:

```text
P_0(L)=13 p(L),       P_1(L)=p(L),
Q_0(R)=13 q(R),       Q_1(R)=q(R).                          (6)
```

This proves why `13*dual[address,0]` appears: selector idempotence identifies
the **raw** direct absent value with the twist-zero present value, and the
factor `13` then comes from the constant extension `(3)`.  Neither equality
is a definition of carrier absence.

## 3. The true four corners

At any fixed endpoint origin, put `w=p(L)q(R)`.  Equations `(6)` give the
literal allocation square

```text
(v_00,v_10,v_01,v_11)=w(169,13,13,1).                      (7)
```

Its Segre--Hadamard coordinates are

```text
(D_0,D_1,D_2,D_3)=w(196,168,168,144).                      (8)
```

In particular

```text
D_3=(P_0-P_1)(Q_0-Q_1)=144w,                               (9)
D_0D_3=D_1D_2.
```

The calculation is not special to `13`.  For any finite twist group of order
`m`, whenever a retained sheet meets the present carrier only at the identity
twist, the same argument gives the rank-one central profile

```text
w(m^2,m,m,1)
```

and Hadamard profile

```text
w((m+1)^2,m^2-1,m^2-1,(m-1)^2).                         (9b)
```

Thus the Boolean mixed coordinate is universally `(m-1)^2w`; the LRC case
is `m=13`.

The exact fixed-origin values are:

```text
field 352341050142921841
  (P0,P1,Q0,Q1)
    =(161934023005863791,175075409527953449,
      300649489018742460, 50230041473974177)
  corners
    =(346389504894336138,351883238969953710,
      351883238969953710,216790045382338969)
  (D0,D1,D2,D3)
    =(209922877787817004,129599459511997169,
      129599459511997169,211754122479689528)

field 956354278959359281
  (P0,P1,Q0,Q1)
    =(155560747474790541, 11966211344214657,
      879469786637801433,361914377113479889)
  corners
    =(17884554354397974,295638590014756546,
      295638590014756546,96307143767239679)
  (D0,D1,D2,D3)
    =(705468878151150745,877931689546517576,
      877931689546517576,479268797051483842).
```

All displayed factors, corners, and Hadamard coordinates are nonzero.

### Primal common atom versus central Fourier square

The nonzero value in `(9)` must not be moved before the carrier transform.
At a fixed endpoint value `w=pq`, the four **primal** allocation arrays are

```text
B(a,b)=w,
P(a,b)=1_(a=0)w,
Q(a,b)=1_(b=0)w,
H(a,b)=1_(a=b=0)w.                                      (9a)
```

The sole primal point at which all four are co-supported is `(a,b)=(0,0)`.
There the four-vector is `(w,w,w,w)` and its pointwise Boolean mixed face is
zero.  The allocation face

```text
Omega=B-P-Q+H
```

equals `w` on the `12*12` points with `a,b!=0` and zero otherwise.  Summing
this face at central carrier harmonic gives `144w`, exactly `(9)`.

Therefore `(7)--(9)` are a genuine **central transformed allocation K4**,
not a nonzero pre-Fourier common-atom face.  They avoid the mixed-character
one-corner hostile because all four central transforms are nonzero, but they
still do not satisfy THM-2772's strongest requirement that a nonzero mixed
face already live on one fourfold co-supported primal atom.

## 4. THM-2763 representative gauge and the marked-sheet sidecar

The delta at twist zero in `(4)` is not an invariant unmarked coordinate.
THM-2763's representative relation is

```text
(ell,a,b) ~ (ell+W,a+1,b-1).                              (10)
```

The literal constructor does satisfy the required **marked** covariance.
Changing representative transports the endpoint masks and the selected sheet
together:

```text
source: (ell,I,delta_0) -> (ell+W,T I,delta_1),
target: (ell,J,delta_0) -> (ell+W,T J,delta_12).           (11)
```

Direct exact evaluation gives equality coefficientwise at all `169`
endpoint addresses, on both endpoint sides, in both certified fields.  The
carrier-mask support vectors undergo the same cyclic reindexing:

```text
delta_0 -> delta_1                 on the source,
delta_0 -> delta_12                on the target.          (12)
```

The central carrier harmonic is the orbit sum over the twist coordinate.
The changes `h -> h+1` and `h -> h-1` merely reindex that sum, so the `k=0`
and `l=0` coefficients in `(6)` descend.  What descends is the coefficient
together with the sidecar

```text
(marked ancestry sheet, twist origin, endpoint representative). (13)
```

The phrase “fixed sheet” therefore means a fixed point of this marked torsor,
not the same physical interval in every representative chart.  If the
sidecar is forgotten and `I,J` are held physically fixed while
`ell -> ell+W`, coefficient covariance fails already at address `(0,0)`;
exactly `27/169` rows fail on each endpoint side in either field.  The first
hostiles are

```text
field 352341050142921841:
  source 254455016269350867 -> 0,
  target 231164267889491750 -> 0;

field 956354278959359281:
  source 318932490657369324 -> 0,
  target 630230755085920022 -> 0.                         (13a)
```

Thus
the four-state allocation is globally typed on the marked-sheet pullback but
is only a distinguished-`e_0` local chart after forgetting that sidecar.

## 5. Sharp endpoint-translation no-go

The literal toggle changes the carrier-allocation state while retaining the
same endpoint origin.  Equation `(6)` gives the complete covariance:

```text
P_1(L)=13^(-1)P_0(L),       Q_1(R)=13^(-1)Q_0(R).          (14)
```

Thus its intrinsic shift, character, and scalar triples are

```text
source: ((0,0),(0,0),13^(-1)),
target: ((0,0),(0,0),13^(-1)).                            (15)
```

The smallest failed identity in the proposed symplectic bridge occurs before
any determinant:

```text
L_1-L_0=(0,0), not (0,1);
R_1-R_0=(0,0), not (12,0).                                (16)
```

Accordingly all four literal flag corners have the same address

```text
(L,R,q,Delta)=((0,0),(0,0),(0,0),0),
```

and the determinant mixed face is

```text
det((0,0),(0,0))=0, not 1.                                (17)
```

The nonzero coefficient face `(9)` therefore does not map to the geometric
determinant face by the proposed allocation-to-translation identification.

As a hostile comparison, replacing `I` and `J` by all their translated
carrier sheets gives the pending carried banks.  Exhaustion of all endpoint
shifts, characters, and scalars finds no affine covariance from the literal
bare bank to either replacement bank in either certified field.  Hence the
steps `(0,1)` and `(12,0)` cannot be recovered by retyping that replacement
as a toggle.

There is also a cheaper structural witness.  On each endpoint side and in
each field, the inverse endpoint-address support has size

```text
literal absent/present: 81/81,
transported-sheet replacement: 121.                         (18)
```

Endpoint translation and character multiplication only translate/modulate
the inverse-address bank, so they preserve its support cardinality.  The
inequality `81!=121` rules out every translation/character/scalar covariance
before any value comparison.

## 6. Consequence and boundary

This proves a useful positive/negative pair:

```text
positive:
  one literal retained sheet has a nonzero central transformed Boolean
  allocation coefficient face in two exact-order fields;

negative:
  its sole fourfold co-supported primal point has zero mixed face, and
  that allocation acts in the carrier-twist coordinate, not as a
  determinant-one endpoint-origin translation.
```

It supplies neither the root-deck functional `mu -> c_j` nor a seven-clock
Cech correction.  It does not imply row exclusion or LRC(14).

## 7. A filtered/Bockstein lead

The universal profile has a potentially useful, but currently only
conditional, interpretation:

```text
(v00,v10,v01,v11)=w(13^2,13,13,1),
D3/v11=(13-1)^2=144 == 1 mod 13.                          (19)
```

If `w` is viewed as a unit in the relevant localization, the four corners
carry the integral `13`-adic filtration profile

```text
(2,1,1,0).                                                (20)
```

Reduction modulo `13` sends the projective corner vector to
`(0,0,0,1)`.  This is exactly THM-2772's virtual one-corner hostile, so simply
using `(19)` as the root-deck map is not lawful: it forgets three allocation
vertices.

Here, however, that one-corner reduction has a canonical integral lift,
marked-sheet provenance, and the filtration `(20)`.  A plausible next object
is therefore the two-dimensional Bockstein/associated-graded (or Rees)
sidecar of this allocation square.  Such a sidecar would have to distinguish
the lifted profile `(20)` from an arbitrary post-Fourier one-corner vector
before the normalized unit `144 mod 13=1` could be used to pay the Cech
correction.  No such root-deck map is proved here.

Reproduction:

```bash
python3 .scratch/lrc_absent_present_constructor/audit.py
python3 -O .scratch/lrc_absent_present_constructor/audit.py
```
