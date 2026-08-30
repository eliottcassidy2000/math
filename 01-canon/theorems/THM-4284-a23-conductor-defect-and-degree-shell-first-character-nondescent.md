---
id: THM-4284
title: "A23 conductor defect and degree-shell first-character nondescent"
status: >
  PROVED FORMAL-LOCAL RELATIVE TO THM-4272/4280 + FINITE-EXACT AUDIT PASS.
  The saturated graph of two branches of A_(2m-1), whose formal-log
  difference has order s, is exactly A_(2ell-1) for ell=min(m,s). Full raw
  descent is ell=m, whereas first-character descent is only ell>=2. At the
  raw A23 contact, actual Hom ramification 1,2,4 gives graph types A1,A3,A7;
  every degree-34/42 class gives A1 and loses eleven conductor lengths. Thus
  common-tangent preservation would already suffice for the intended
  contradiction, but it is not proved for the Keller graph. JC(2) remains
  open.
source: codex-continue-frontiers-20260829
depends_on:
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4280-integral-three-channel-fat-contact-observer-and-sharp-five-jet-bound
related:
  - THM-4279-four-channel-formal-log-hasse-observer-for-e0-hom-at-fat-contact
primary_script: 04-computation/jc23_a23_conductor_graph_spectrum_thm4284.py
primary_output: 05-knowledge/results/jc23_a23_conductor_graph_spectrum_thm4284.out
primary_script_sha256: 2029c54c17b2c6b9a57982e13644767c3105faab47dd3f13906ddea65373a723
primary_output_sha256: 9ead08564cd036220996f2cc70e20c50530daebbb485a19ac7031bd1ca696f1e
hash_basis: raw LF bytes
audit: >
  PASS. Exact Groebner elimination independently recovers the four saturated
  affine-modification presentations at orders 1,2,4,12. A truncated
  normalization-pair audit freezes conductor colengths, and hostile controls
  verify the unsaturated z-line, anti-diagonal cotangent cokernel, and
  unchanged W-saturation of the transverse contact divisor.
---

# THM-4284 -- A23 conductor defect and degree-shell first-character nondescent

**PROVED FORMAL-LOCAL RELATIVE TO THM-4272/4280 + FINITE-EXACT AUDIT PASS.
THE KELLER APPLICATION AND `JC(2)` REMAIN OPEN.**

## 1. The graph-image theorem

Let `k` be a characteristic-zero field, put `R=k[[b]]`, and define

```text
A_m=k[[b,q]]/(q(q-b^m)),                       m>=2.    (1)
```

The normalization embedding is

```text
A_m -> R direct-sum R,
b   |-> (b,b),                       q |-> (0,b^m).    (2)
```

It identifies

```text
A_m={ (f_0,f_1):f_1-f_0 in b^m R },
(R direct-sum R)/A_m ~= R/(b^m)                              (3)
```

where the quotient map is branch difference. Thus the full nonreduced
gluing obstruction is the branch-difference class modulo `b^m`.

Take two formal branch maps to the same point of a smooth one-dimensional
target. After choosing a target coordinate, and after applying the formal
logarithm when the target is `E_0`, write their translated series as

```text
L_0,L_1 in bR,                  Delta=L_1-L_0.          (4)
```

Let

```text
s=ord_b(Delta),                 ell=min(m,s),            (5)
```

with `s=infinity` if `Delta=0`. The completed coordinate ring of the
saturated scheme-theoretic image of the two branch graphs is

```text
B_Delta=A_m[(L_0,L_1)]^ = A_ell.                       (6)
```

Consequently,

```text
full raw descent       iff Delta in b^m R iff ell=m;
first channel vanishes iff Delta in b^2 R iff ell>=2.   (7)
```

For `ell<m`, this partial normalization is the exact affine modification

```text
A_ell=A_m[q/b^(m-ell)]                                  (8)
```

inside the total quotient ring. Its defect and analytic type are

```text
length_k(A_ell/A_m)=m-ell,
A_ell has two-branch type A_(2ell-1).                   (9)
```

Thus preserving the first formal-log character is a length-two/common-
tangent condition. It is strictly weaker than preserving the complete
length-`m` conductor overlap.

## 2. Proof of the graph-image theorem

Every diagonal pair `(f,f)` belongs to `A_m`. Subtracting `(L_0,L_0)` from
the new graph coordinate shows that the completed ring in `(6)` is generated
over `A_m` by `(0,Delta)`. Its possible branch differences form the ideal

```text
(b^m,Delta)=(b^ell).                                    (10)
```

Indeed, when `s<m`, write `Delta=b^s u` with `u` a unit. Multiplication by
the diagonal unit `(u^(-1),u^(-1))` gives `(0,b^s)`, and diagonal
multiplication then gives every pair whose difference lies in `b^sR`. When
`s>=m`, the new element already lies in `A_m`. Equations `(3)` and `(10)`
prove `(6)--(7)`.

Under `(2)`,

```text
q/b^(m-ell) |-> (0,b^ell),                              (11)
```

which proves `(8)`. The quotient in `(9)` is the ideal quotient
`b^ellR/b^mR`. Finally, translating the equation
`q(q-b^ell)=0` by `q-b^ell/2` gives, up to a unit and a formal coordinate
rescaling, `y^2-b^(2ell)=0`; this is type `A_(2ell-1)`.

The use of the formal logarithm does not change the graph image: in
characteristic zero it is an invertible formal target coordinate with inverse
the formal exponential.

## 3. The raw A23 spectrum and the shorter sufficient bridge

THM-4272 identifies the raw contact with `(1)` at `m=12`, hence with
`A_23`. Pair a constant map on its rational branch with a nonconstant actual
global map `C_0->E_0` on the other branch. THM-4280 sharpens the possible
formal-log orders to

```text
s in {1,2,4}.                                            (12)
```

Equations `(5)--(9)` give the exact graph-singularity dictionary

```text
map ramification       1        2        4
saturated graph        A_1      A_3      A_7
conductor loss         11       10       8.             (13)
```

In particular, THM-4280 proves that every degree-`34` or degree-`42` class
has nonzero first character `c_1`. Its branch difference has `s=1`, so its
saturated graph is `A_1` and the projection to the raw `A_23` source loses
exactly eleven conductor lengths.

This reverses the useful proof obligation. To contradict such a degree class,
a future Keller argument need not prove full equality on `12Q`. It would
suffice to establish that the **saturated** Keller graph retains the common
tangent, equivalently that its branch difference is in `b^2R`. That would
force `c_1=0`, contradicting THM-4280. No such common-tangent preservation is
proved here.

## 4. Cotangent form of the missing sidecar

At the closed point, normalization induces

```text
m_(A_m)/m_(A_m)^2 ->
 (m_R/m_R^2) direct-sum (m_R/m_R^2),
db |-> (db,db),                   dq |-> (0,0).          (14)
```

The image is the diagonal cotangent line and the cokernel is the
anti-diagonal line. The two branch pullbacks of the target cotangent have the
same image in degree one exactly when the coefficient of `b` in `Delta`
vanishes. Hence the missing first-channel sidecar is precisely

```text
gr_b^1(R/(b^m)),                                        (15)
```

the anti-diagonal branch cotangent/conductor class.

Pointwise agreement after resolution retains only the residue-field value;
it need not retain `(15)`. This is the exact quotient loss behind the open
resolved-to-raw step.

## 5. Minimal saturation hostile

Set `Delta=b` and introduce its graph coordinate `z`. The two graph branches
have ideals

```text
I_0=(q,z),                   I_1=(q-b^m,z-b).           (16)
```

Their saturated union is

```text
I_0 intersect I_1=(q-b^(m-1)z, z(z-b)).                 (17)
```

By contrast, the cleared equations alone give

```text
K=(q-b^(m-1)z, q(q-b^m)).                               (18)
```

Substitution of the first generator into the second gives only

```text
b^(2m-2) z(z-b)=0.                                      (19)
```

Thus `K` has a parasitic full `z`-line over `b=q=0`, whereas

```text
K:b^infinity=(q-b^(m-1)z,z(z-b)).                       (20)
```

For `m=12`, the new regular coordinate is exactly `z=q/b^11`; the correct
saturation is `A_1`, not the raw `A_23` graph. The quotient becomes regular
on the eleventh common-tangent blowup chart. This describes common-tangent
depth only; an embedded-resolution convention may perform another blowup to
make the whole total transform simple normal crossings.

Composing the branch coordinate `b` with the elliptic formal exponential
gives an `E_0`-valued version of the same hostile. Equations `(18)--(20)` are
the precise MISTAKE-455 mechanism: cleared graph equations neither remove the
parasitic exceptional component nor prove that the projection is the raw
source after correct saturation.

## 6. A transverse saturation criterion

There is one positive route toward the missing bridge. In

```text
S=k[[b,W]],                         g=b^m-W,             (21)
```

one has

```text
(g):W^infinity=(g).                                       (22)
```

Indeed, `S/(g)~=k[[b]]`, and `W` maps to the nonzerodivisor `b^m`. Therefore,
if a branch-log difference `Delta in S` is genuinely regular and belongs to
`(g)` after inverting `W`, it already belongs to `(g)` before inversion.
Specialization then gives

```text
Delta(b,0) in b^m k[[b]].                                (23)
```

The unresolved hypothesis is regularity in the unsaturated raw two-variable
ring; a resolved or cleared presentation does not provide it.

The tangent/conormal split is also exact. If

```text
Delta=alpha b+beta W+terms of total order at least two,  (24)
```

then restriction to `W=b^m` has first `b`-coefficient `alpha`, independent
of `beta`. At the origin `dg=-dW`, whereas the contact tangent is generated
by `db`. Thus the first `b`-channel is tangential and a bare `W`-jet is
conormal. This proves formal independence only; it does not identify any
physical transverse Jacobian from THM-4265 with `beta`.

## 7. Exact audit and scope

The optimization-safe script

```text
04-computation/jc23_a23_conductor_graph_spectrum_thm4284.py
```

uses exact Groebner elimination of the localization variable to recover
`(20)` and its order-`2`, order-`4`, and order-`12` analogues. A separate
truncated normalization-pair calculation freezes the colengths, while hostile
controls retain the unsaturated `z`-line, cotangent cokernel, and transverse
`W`-saturation. Its stored output is

```text
05-knowledge/results/jc23_a23_conductor_graph_spectrum_thm4284.out
```

Reproduce with

```bash
python3 -B 04-computation/jc23_a23_conductor_graph_spectrum_thm4284.py
python3 -B -O 04-computation/jc23_a23_conductor_graph_spectrum_thm4284.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_a23_conductor_graph_spectrum_thm4284.py
```

The proved connection ledger is

```text
source:              two normalized branches of the raw A_(2m-1) contact;
target:              their saturated graph image in a formal target chart;
map:                 branch difference Delta modulo b^m;
preserved predicate: conductor-overlap length and first target character;
destroyed data:      anti-diagonal jets discarded by resolved node values;
restoring sidecar:   gr_b^1(R/(b^m)), or the full difference modulo b^m;
cheapest hostile:    Delta=b and z=q/b^(m-1).            (25)
```

Nothing here proves that the Keller response has a saturated raw graph,
retains the common tangent, or is regular in the ring `(21)`. It supplies no
wall exclusion, exact-`M=12` entry, `JC(2)`, or `DC(2)` conclusion.
