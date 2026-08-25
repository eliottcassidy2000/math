---
id: THM-4039
title: "Exceptional-quartic actual J3/J4 lift and scalar-gate rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The frozen
  THM-3688 actual-target certificate has exact third-stage debt
  Lambda(D_3)=0 over the common irreducible quartic field.  THM-3737
  therefore supplies actual F_4,G_4 with J_3=0, and THM-3677's universal
  fourth-order identity then supplies the J_4 scalar gate and actual
  F_5,G_5.  Thus all four exceptional embeddings admit one stagewise actual
  lift through J_4.  The induced next-cokernel map from every admissible
  L_0-kernel is zero, and the coupled F_4/F_5 choice response at J_5 is also
  zero, with the earlier frozen data through F_3,G_3 fixed.  This theorem
  does not evaluate J_5; downstream THM-4043 reaches J_6, and THM-4046
  reaches J_7 and proves the sharp J_8 obstruction.  No conclusion outside
  that exceptional family, and no JC(2) conclusion, is asserted here.
source: jc2-double-zero-rebuild-20260824 / retained-jet continuation, 2026-08-24
audit: >
  PASS -- an independent full-polynomial reconstruction recovered the
  canonical normal sidecars, degrees 367,364,746 for delta(F_2),delta(G_2),
  D_3, the direct i+j=4 identity, and the same zero quartic-field scalar.
  A separate retained-jet implementation propagated value, x, q, and mixed
  x-q data through all canonical target expressions.  Ten good split fibres
  across p=137,163,179,193 supplied hostile controls.  Independent algebra
  audits checked the adjacent and coupled two-step rigidity arguments and
  their representative/sidecar scope.  Normal and optimized production
  executions byte-match the frozen transcript.
depends_on:
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3687-russell-cylinder-exceptional-quartic-actual-j0-lift
  - THM-3688-russell-cylinder-exceptional-quartic-actual-j1-j2-lift
  - THM-3703-russell-cylinder-exceptional-quartic-sagbi-module
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
related:
  - THM-3651-russell-cylinder-degree-seven-double-zero-sixth-order-closure
  - THM-4034-exceptional-quartic-global-conductor-degree-178
  - THM-4043-exceptional-quartic-shifted-stable-identities-and-j6-lift
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_j3_lift_rigidity_thm4039.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_j3_lift_rigidity_thm4039.out
script_sha256: abd5bef4b306710cbc0e5ff76fc7a7b7bb580eb296dd5ad3312e39ea2df36993
output_sha256: 6f20e7d1efa07a71d530596be59066c6220be22af74be144a40368568b70fa8d
hash_basis: raw LF bytes
---

# THM-4039 -- the exceptional actual lift reaches `J_4`

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The scalar
gate left open by THM-3688 vanishes exactly, not numerically.  This produces
actual target coefficients through `J_4`.  At the same time, the apparent
kernel freedom at one stage is proved unable to tune the next scalar gate;
even the coupled two-step `F_4,F_5` choice cannot tune `J_5`.

All rings below have characteristic zero.

## 1. Frozen actual certificate and stable coefficients

Retain the irreducible field, restriction algebra, and actual leading pair
of THM-3683, THM-3687, and THM-3703:

```text
K=Q[alpha]/(72783360 alpha^4-77822208 alpha^3-28419741 alpha^2
                                      +7849770 alpha-1276420),

S=K[B,C,E] subset K[x],
L_0(P,Q)=C'Q-E'P,
C'G_1-E'F_1=1.                                           (1)
```

THM-3688 freezes actual target representatives through order three and
proves the full polynomial identities

```text
J_0=1,                         J_1=J_2=0.                (2)
```

For an actual representative `H`, write

```text
gamma_k(H)=gamma(partial_q^k H).
```

The pulled source coefficients are

```text
a_m=sum_(0<=k<=m/2) gamma_k(F_(m-2k))/k!,
b_m=sum_(0<=k<=m/2) gamma_k(G_(m-2k))/k!,                (3)
```

and

```text
J_n=sum_(i+j=n+1)(j a_i'b_j-i a_i b_j').                (4)
```

In particular, put

```text
a_2=F_2+gamma_1(F_0),       b_2=G_2+gamma_1(G_0),
a_3=F_3+gamma_1(F_1),       b_3=G_3+gamma_1(G_1),

a_4^*=gamma_1(F_2)+gamma_2(F_0)/2,
b_4^*=gamma_1(G_2)+gamma_2(G_0)/2.                       (5)
```

Here `F_i,G_i` in formulas denote their frozen restrictions when no
confusion is possible.  Directly expanding `(4)` gives

```text
J_3=D_3+4L_0(F_4,G_4),                                  (6)

D_3=4L_0(a_4^*,b_4^*)
    +(3F_1'b_3-F_1b_3')
    +2(a_2'b_2-a_2b_2')
    +(a_3'G_1-3a_3G_1').                                (7)
```

The normal sidecars in `(5)` are attached to the frozen canonical target
expressions.  Replacing them by arbitrary restriction lifts would change
the question.

## 2. Exact cancellation of the `J_3` scalar

THM-3737 proves the cutoff-free image identity

```text
L_0(S^2)=ker Lambda,
Lambda(H)=5H(-1)/18-H(0)+13H(1)/18.                     (8)
```

The production companion propagates the four exact local coordinates

```text
gamma(H), d_x gamma(H), gamma(partial_q H),
d_x gamma(partial_q H)                                  (9)
```

at `x=-1,0,1` through the entire THM-3703 canonical target grammar.  It
checks every canonical address and restriction against the frozen THM-3688
certificate, reconstructs the required jets of `F_2,G_2,F_3,G_3`, and
checks `(7)` independently against the direct sum `(4)` at all three points.

Each of the four displayed summands in `(7)` has nonzero `Lambda` value.
For the canonical comma-separated four-coordinate field serialization their
exact lengths and hashes are

```text
summand  characters  sha256
04          46444    bcbe6bd74afcdf59f2003eb10325b021e6a073a785fdeb72524566c22f3176b6
13          92834    a2decb387cde9f3d294702c20d9f76d79909c8490a4782ec69442e7f816130b2
22          46462    f41828e4b1f257e1e42688221e0f2e001787e13244082d9c5119071164cb390c
31          92840    c52bbb16e38bb4962fe0ffeaeb15107c120d409c85714d8f0b1faaf994b5d03d
                                                                    (10)
```

Nevertheless they cancel exactly in `K`:

```text
boxed: Lambda(D_3)=(0,0,0,0).                           (11)
```

This is not a one-embedding numerical cancellation.  An independent
full-polynomial path reconstructs `gamma_1(F_2),gamma_1(G_2)` with degrees
`367,364`, obtains `deg D_3=746`, checks the complete polynomial identity
`(7)`, and recovers `(11)`.  Ten good-reduction fibres over four split primes
give additional hostile controls but are not used to infer exact vanishing.

## 3. Actual continuation through `J_3` and `J_4`

By `(8)` and `(11)`, `-D_3/4` lies in `L_0(S^2)`.  Hence there are
restrictions `F_4,G_4 in S`, and therefore actual target representatives,
for which

```text
J_3=0.                                                   (12)
```

No bounded solve or untyped normalization polynomial is being substituted
here: the equality `(8)` is precisely the actual-target image theorem.

The exceptional polynomial is `Q_r` on THM-3677's zero-fourth-debt
parabola.  That theorem proves for every arbitrary target two-form the
universal identity

```text
J_0=1 and J_2=0       ==>       Lambda(J_4)=0.           (13)
```

It uses neither `J_1` nor `J_3`.  After choosing `(12)`, first take the
provisional actual coefficients `F_5=G_5=0` and apply `(13)` to that
truncated pair.  Its `J_4` debt therefore lies in `ker Lambda`.  The general
new term is `5L_0(F_5,G_5)`, so a second application of `(8)` supplies actual
`F_5,G_5` that cancel the complete debt and give

```text
boxed: J_0=1,              J_1=J_2=J_3=J_4=0.           (14)
```

Because `(14)` is proved over the irreducible quartic field, it holds after
each of its four complex embeddings.

## 4. The adjacent kernel has zero cokernel response

The positive lift does not arise from a freely tunable scalar at every
stage.  First identify the exact kernel of `L_0` on actual restrictions:

```text
ker(L_0|S^2)
 ={(C'U,E'U):U in K[x], C'U in S, E'U in S}.            (15)
```

Indeed, for a kernel pair `(P,Q)`, the unique potential is

```text
U=G_1P-F_1Q,                  P=C'U, Q=E'U.              (16)
```

At the retained points every element of `S` has one common value, while

```text
C'=(3,3,3),                 E'=(-9,4,9).                (17)
```

If both members of `(15)` lie in `S`, the first row in `(17)` makes the
three values of `U` equal; the second row then forces that common value to
be zero:

```text
U(-1)=U(0)=U(1)=0.                                     (18)
```

For every `H in S`, the chain rule at the common target point gives

```text
(H'(-1),H'(0),H'(1)) in
span_K{(3,3,3),(-9,4,9)}=ker Lambda.                    (19)
```

Since `P'=C''U+C'U'`, equations `(17)--(19)` imply

```text
Lambda(U')=0.                                           (20)
```

Now let `m>=2`.  With all earlier actual representatives fixed, changing
`(F_m,G_m)` by `(P,Q)` preserves `J_(m-1)` and changes the next debt by

```text
T_m(P,Q)=mF_1'Q-F_1Q'+P'G_1-mPG_1'
        =U'+A_mU,                                       (21)

A_m=C''G_1-mC'G_1'+mF_1'E'-F_1E''.                     (22)
```

Equation `(18)` kills `A_mU` at the retained points and `(20)` kills the
remaining term.  Therefore

```text
boxed: Lambda(T_m(P,Q))=0             for every m>=2.   (23)
```

Thus all actual solutions of one stable equation, with the earlier data
fixed, give the same immediately following scalar gate.  Changing the
actual representative inside the same restriction fibre is also invisible
to that immediate gate because only `gamma_0(F_m),gamma_0(G_m)` occur there.
The normal sidecar first occurs one equation later.

## 5. Even the coupled `F_4,F_5` choice cannot tune `J_5`

For the frozen data through `F_3,G_3`, the next scalar has a stronger
choice-invariance.  Change an `F_4,G_4` solution by the kernel pair `(15)`
and put

```text
v=(U'(-1),U'(0),U'(1)).                                 (24)
```

Both derivative rows `P'=3v` and `Q'=E'v` belong to the tangent plane
`ker Lambda`.  The two independent equations

```text
Lambda(v)=0,                       Lambda(E'v)=0         (25)
```

have intersection exactly the constant line, so

```text
v=lambda(1,1,1).                                        (26)
```

Let `delta` denote one normal `q` derivative of an actual representative.
The retained vertical compiler rows are

```text
n_B=delta(B)=(0,0,0),
n_C=delta(C)=(2,0,-2),        n_E=delta(E)=(-2,4,-2).   (27)
```

The ordinary retained triple makes the vertical value of a target function
depend only on its retained derivative row: both `B'` and `n_B` vanish,
while the `C',E'` rows form the tangent plane.  Consequently `(26)` gives
the retained-row identities

```text
delta(P)=lambda n_C,             delta(Q)=lambda n_E,   (28)
```

independently of the chosen actual representatives.  Moreover

```text
C'n_E-E'n_C=12(1,1,1).                                 (29)
```

The direct change in the `J_5` debt caused by `(P,Q)` is

```text
4a_2'Q-2a_2Q'+2P'b_2-4Pb_2'
       +6(C'delta(Q)-delta(P)E').                       (30)
```

At the retained points `P=Q=0`.  The frozen identity `J_1=0`, together
with `F_1=0,G_1=1/3` there, gives

```text
3b_2-a_2E'=-F_1'/6.                                    (31)
```

Equations `(26),(29)--(31)` show that both parts of `(30)` have zero
`Lambda`: the first is `-lambda F_1'/3`, and the second is a constant row.

The change `(21)` at `m=4` must also be compensated by changing
`F_5,G_5` to retain `J_4=0`.  Since `U=0` and `U'=lambda(1,1,1)` at the
retained points, its retained value is the constant row
`lambda(1,1,1)`.  The compensating equation is therefore

```text
L_0(Delta F_5,Delta G_5)=h,       h_ret=-(lambda/5)(1,1,1).   (32)
```

Write a compensation using the Bezout pair:

```text
Delta F_5=F_1h+C'U_5,          Delta G_5=G_1h+E'U_5.   (33)
```

Common retained values force `(U_5)_ret=0`.  Admissibility and
`Delta F_5'=F_1'h+3U_5'` then give

```text
Lambda(U_5')=-(1/3)Lambda(F_1'h).                       (34)
```

Substitution in the `m=5` version of `(21)` gives

```text
Lambda(T_5(Delta F_5,Delta G_5))
       =(5/3)Lambda(F_1'h)=0,                           (35)
```

because `h` has constant retained value, so
`Lambda(F_1'h)=h_ret Lambda(F_1')=0`.  Two compensations differ by `(15)`
and are invisible by `(23)`.  The normal sidecar of `F_5,G_5` first occurs
in `J_6`, so changing their representatives cannot affect `J_5`.  Hence

```text
boxed: the J_5 scalar is independent of all compatible F_4,F_5
       restriction choices and their actual representatives,             (36)
```

provided the earlier frozen actual data through `F_3,G_3` are fixed.
This does not cover changing `F_2,F_3` or their normal sidecars.

## 6. Strict boundary and reproduction

The theorem itself proves finite stagewise existence, not a controlled or
global series.  In particular it does not

- compute or prove vanishing of the now well-defined `J_5` scalar;
- freeze explicit actual target expressions for `F_4,G_4,F_5,G_5` or bound
  their degrees;
- prove a coherent all-order continuation or convergence/algebraization;
- construct a positive global pair on the Russell cylinder;
- construct a polynomial Keller map or prove its noninjectivity; or
- prove a counterexample to `JC(2)`, which remains open.

Downstream THM-4043 applies `w`-shifts of THM-3683's universal identity,
proves the `J_5` gate vanishes, and reaches `J_6`.  That later theorem does
not retroactively add a `J_5` computation to the present certificate.
THM-4046 then reaches `J_7` and proves that `J_8` cannot be cleared.

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_j3_lift_rigidity_thm4039.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_j3_lift_rigidity_thm4039.py
```

Normal and optimized exact runs byte-match the stored LF transcript.  The
script contains no Python `assert` statements.  **QED.**
