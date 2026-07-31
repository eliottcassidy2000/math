# Independent hostile audit: atom-dependent common sheet

**Witness source under audit**

```text
.scratch/lrc_atom_dependent_four_corner/atom_dependent_four_corner.py
stable corrected LF SHA-256:
8f064ff8922a787211cb66ad7f7a9cdfe5a87ded84100d7dbd39e264e387482f
```

**Verdict:** promote a sharply split result, not the originally suggested
unified determinant-current bridge.

1. **Positive:** one literal THM-2584 inverse-branch sheet survives as a
   source/target common physical atom, including all simultaneous factors,
   semantic carry data, and all four transformed central
   bare/source/target/both marginals at one fixed endpoint origin.  Both the
   full orbit coefficient and the literal zero-twist mask pass the carrier
   representative gauge.  The transformed central `D3` is nonzero in both
   certified fields and has the universal ratio described in Section 4.
2. **Negative common-atom gate:** before twist summation, the sole address
   supporting all four vertices is `(a,b)=(0,0)`, where their values are
   equal and `D3=0`.  The nonzero transformed `D3` is supported on the
   `12*12` bare-only complement, so it is a virtual marginal face rather
   than THM-2772's pre-Fourier common-atom face.
3. **Positive but separate:** four carried-array endpoint samples on the
   orbit obtained by replacing the fixed sheet by each of its thirteen
   translates have nonzero coefficient `D3`, while the chosen geometric
   determinant-label mixed face is `det(s,t)=1`.
4. **Negative bridge:** on the literal fixed sheet the allocation toggle has
   zero endpoint step and scalar `1/13`, so its determinant mixed face is
   zero.  The transported-sheet orbit bank has no shift--character--scalar
   covariance with the bare bank anywhere on either full 169-point endpoint
   fibre in either field.  Thus neither construction induces the proposed
   determinant-one endpoint translations.

No seven-clock/root-deck map, Cech correction, row exclusion, or LRC(14)
conclusion follows.

## 1. Multiplicity was successfully unfolded

The old interval profile retained only weight

```text
27,581,135,604.
```

The companion and an independent constructor replay unfold it as

```text
U inverse-branch labels: 966,606
V inverse-branch labels:  28,534
product:                  27,581,135,604.
```

At the source midpoint and its `+SHIFT` target, the complete ordered label
sets coincide.  The independent packed-label digests are

```text
U: d7b599a50d8647ac7957149d1e84588b3602fdb897dc923432102ce0e3a51cf5
V: 7072f55b7c54d7c32f2b8ac87ed2aa59f091b3dfdd3d1f193778e9a75c93d00f
```

The concrete path

```text
(a,b,eprime)=(59162,26,56658)
```

is active throughout the full selected source interval and its translated
target.  Independent inverse-label event enumeration places both intervals
inside one chamber with positive left/right slacks:

```text
U_Q: (1114492399020, 2185428144900)
U_E: (3723575735880,17550030017520)
V_E: (1442440654200,19186496489160).
```

Changing `b` to `25` or `27`, `a` to `59161`, or `eprime` to `56657`
deactivates the corresponding nearby path.  Hence the path check is not a
restatement of positive aggregate weight.

This repairs the main multiplicity objection.  The physical pullback is

```text
(x,a,b,eprime) -> (x+SHIFT,a,b,eprime),
```

not an anonymous choice of one scalar unit from a weighted profile.

## 2. Physical and semantic atom typing

The selected cell is

```text
(s,t,clock)=(0,4,1),
I=[142004992589460,142005019034340),
weight(I)=27581135604.
```

The unit sheet is `1_I` with the inverse-branch path retained as a sidecar.
It is contained simultaneously in every native and adjacent source/target
factor in order

```text
(E3,clock,q1,q2,c2,c3).
```

Thus its label is intrinsically `common_after_all_factors`, not an ordered
first-failure artefact.  Its unit semantic values are

```text
source carry 12: (root-zero part 0, root-one part 103478815440)
target carry 6:  (root-zero part 0, root-one part 103478815440).
```

The all-567 control still finds `193` nonempty source wings and `193`
nonempty target wings, with every opposite-wing intersection empty.  The
new sheet lies in `A intersect B`; it does not evade the old no-go by
reusing `A\B` and `B\A`.

## 3. Fixed-sheet versus transported-sheet present states

The outer selector `1_xi` and carrier factor must be kept conceptually
separate.  On the chosen sheet,

```text
1_xi C_source,0 = 1_xi,
1_xi C_target,0 = 1_xi.
```

Therefore direct omitted-carrier endpoint evaluation

```text
integral 1_xi * endpoint factors
```

equals the present zero-twist raw evaluation

```text
integral 1_xi C_0 * endpoint factors.
```

The absent carrier has no twist variable, so its raw value is extended
constantly over the thirteen formal twists.  The unnormalized central
Fourier coefficient is consequently

```text
13 * dual[address,0].
```

This is a derived atomwise idempotence identity.  It would be circular to
*define* physical absence as the zero-twist present coefficient without
first proving `1_xi C_0=1_xi`.

The same outer selector must remain in the present states.  If `C_a` is the
carrier translated by `aT/13`, then the selected interval is so short that

```text
1_xi C_0=1_xi,
1_xi C_a=0 for every a!=0.
```

The pre-correction companion instead evaluated `C_a` on its own translated
support and summed the resulting thirteen endpoint functionals.  That
operation replaces `xi` by `T_a xi`; it is a transported-sheet orbit bank,
not the present allocation of the fixed atom.  Retaining the ancestry label
as a name while moving its physical interval does not repair this type
change.

The literal mask descent is explicit.  In the distinguished section the
present twist mask is `delta_(a=0)`.  Under THM-2763's representative gauge
it transports covariantly:

```text
(xi,a=0,ell) -> (T xi,a=1,ell+W)
```

on the source side; the target sends `delta_0` to `delta_12`.  The companion
checks these mask permutations and all `169` source plus `169` target
`h=0` coefficient rows in each field.  The selected one-term `k=0` sum is
therefore representative invariant when the mask origin is retained as a
sidecar.

## 4. Central transformed marginal square and its raw hostile

At

```text
r=0, k=l=0, L=R=(0,0), q=(0,0), Delta=0,
```

write `p` and `q` for the source and target endpoint functionals of the
fixed sheet at twist zero.  The literal fixed-sheet factors are

```text
(P0,P1,Q0,Q1)=(13p,p,13q,q).
```

Consequently, over every field in which thirteen is a unit,

```text
(v00,v10,v01,v11)=pq*(169,13,13,1),
(D0,D1,D2,D3)=pq*(196,168,168,144).
```

These are the **separately summed central marginals**, not four raw values
on one twist address.  Before the two twist sums, with `w` denoting the
common endpoint factor, the exact arrays are

```text
B(a,b)=w,
P(a,b)=delta_(a=0) w,
Q(a,b)=delta_(b=0) w,
H(a,b)=delta_(a=0)delta_(b=0) w.
```

Their sole fourfold co-supported address is `(0,0)`, where

```text
(B,P,Q,H)=(w,w,w,w),
B-P-Q+H=0.
```

At each of the `12*12=144` addresses with `a,b!=0`,

```text
(B,P,Q,H)=(w,0,0,0),
B-P-Q+H=w.
```

Therefore the central transformed face `144w` is exactly the sum of 144
bare-only virtual faces.  The central rank-one square is numerically
nonzero, but it does not pass THM-2772's “common atom before Fourier
allocation” gate.

For `p=352341050142921841`:

```text
(P0,P1,Q0,Q1)=
(161934023005863791,
 175075409527953449,
 300649489018742460,
 50230041473974177)

(v00,v10,v01,v11)=
(346389504894336138,
 351883238969953710,
 351883238969953710,
 216790045382338969)

(D0,D1,D2,D3)=
(209922877787817004,
 129599459511997169,
 129599459511997169,
 211754122479689528).
```

For `p=956354278959359281`:

```text
(P0,P1,Q0,Q1)=
(155560747474790541,
 11966211344214657,
 879469786637801433,
 361914377113479889)

(v00,v10,v01,v11)=
(17884554354397974,
 295638590014756546,
 295638590014756546,
 96307143767239679)

(D0,D1,D2,D3)=
(705468878151150745,
 877931689546517576,
 877931689546517576,
 479268797051483842).
```

Every displayed factor, product, and Hadamard coordinate is nonzero.  In
both fields

```text
D3=(P0-P1)(Q0-Q1) != 0,
D0 D3-D1 D2=0.
```

Here `D3`, not the zero Pluecker residual, is the nonzero **central marginal**
mixed coordinate.  Its mechanism is Boolean orthogonality on the fixed
selector.  The strongest literal common-atom result is instead the raw
no-go above: fourfold co-support forces equality and zero mixed face.

## 5. Carried endpoint parallelogram is a different square

The four values in the original incoming note are

```text
(Pcarried(L),Pcarried(L+s),Qcarried(R),Qcarried(R+t)),
```

where each carried bank already sums endpoint functionals over the thirteen
translated copies `T_a xi`.  They are neither
`(bare,present,bare,present)` nor four states on literal fixed `xi`.  Their
coefficient `D3` values are

```text
170114988031260853,
757304766188814060
```

in the two fields.  Independently, with

```text
L=R=(0,0), s=(0,1), t=(12,0),
```

the endpoint labels are

```text
q:     (0,0),(0,1),(1,0),(1,1)
Delta: 0,0,0,1,
```

so

```text
delta_Delta=det(s,t)=1.
```

These are different types:

```text
carried-sample D3 = coefficient Mobius face,
delta_Delta       = determinant-label Mobius face.
```

Pluecker zero proves only rank-one decomposability and cannot identify
either face with the other.

## 6. Exhaustive allocation-to-translation no-go

For the **transported-sheet orbit bank**, for every one of the `169`
endpoint shifts, every one of the `169` characters, and one globally
consistent scalar, the companion tests

```text
Ppresent(x)=c chi_u(x) Pbare(x+s)
Qpresent(x)=d chi_v(x) Qbare(x+t)
```

on all `169` endpoint points.  Both covariance sets are empty for both
certified fields.  The search includes the proposed steps and both signs
because it exhausts all shifts.

Therefore neither the original carried-sample quartet nor the separately
formed translated bare/orbit-present candidate is a physical realization
of the fixed-sheet allocation toggle.

For the **literal fixed sheet**, the covariance is instead immediate:

```text
Ppresent_fixed(x)=(1/13)Pbare_fixed(x),
Qpresent_fixed(x)=(1/13)Qbare_fixed(x).
```

Its induced endpoint translations are

```text
s=t=(0,0),
det(s,t)=0.
```

This explains the bridge failure more strongly than an empty search: the
fixed-atom **central marginal** is nonzero, but its raw common-support face
is zero and its intrinsic endpoint action is scalar with zero geometric
movement.

## 7. Minimal missing map

The first missing object is not yet the endpoint translation map.  It is a
non-idempotent **local carrier action on the same atom**.  Indeed, for any
fixed atom of weight `w` and Boolean source/target carrier values
`C_S,C_T in {0,1}`, the raw allocation vector is

```text
(B,P,Q,H)=w*(1,C_S,C_T,C_S C_T),
```

and therefore

```text
B-P-Q+H=w(1-C_S)(1-C_T).
```

Thus a pure indicator allocation has zero mixed face wherever both carriers
are present; its mixed face is supported precisely where both are absent.
This selected witness is not an accident but the sharp Boolean mechanism.
THM-2772's desired common-present face must use amplitudes that change
nontrivially on the same atom--for example a local character, phase, signed
weight, cocycle, or another non-idempotent sidecar--rather than carrier
membership alone.

Only after such a local action exists is the next bridge a map

```text
tau_S,tau_T on the retained path-labelled sheet
  -> endpoint operations on (L,R)
```

whose endpoint images commute, have independent normalized steps, satisfy
`P1(L)=chi_S P0(L+s)` and `Q1(R)=chi_T Q0(R+t)` fibrewise, and transport the
path, semantic labels, endpoint order, and representative gauge together.

The current literal indicator supplies only the zero-step scalar map, while
the transported orbit has no affine shift--character--scalar map at all.  A
future successful construction must therefore add local amplitude data,
use an atom whose carrier translates overlap with a nontrivial weight, or
retain a larger endpoint/translation groupoid in which moving the sheet is
itself part of the typed atom.  Fitting two nonzero samples after choosing
`s,t` is not such a map.

## 8. Replay status

The independent ancestry/arithmetic companion is

```text
.scratch/lrc_atom_dependent_four_corner_audit/audit.py
```

and passes ordinary/optimized byte comparison.  The stable corrected full
witness source keeps the outer atom selector in every present twist, audits
its gauge mask, and passes both an ordinary replay and a separately launched
optimized replay.  Their stdout transcripts are byte-identical:

```text
lines:   86
bytes:   10862
SHA-256: a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f
ordinary:
  .scratch/lrc_atom_dependent_four_corner/atom_dependent_four_corner.normal.out
optimized:
  .scratch/lrc_atom_dependent_four_corner_audit/witness.optimized.out
```

Both end in `ALL EXACT CHECKS PASSED`.  Replays of pre-correction hash
`5e1235171c95453cd5a728a3ec1793b1d9c7e6b7de8e1f085fe23b2e349542a9`
remain evidence only for transported-orbit side controls, not for the
fixed-sheet allocation conclusion.
