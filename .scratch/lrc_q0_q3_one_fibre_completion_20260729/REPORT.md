# q0-to-q3 one-fibre completion: exact selector and provenance theorem

**Status:** `FINITE-EXACT / INDEPENDENTLY HOSTILE-AUDITED; PROMOTED AS
THM-2880`, 2026-07-29.  This is not a row exclusion or an LRC(14)
conclusion.  It executes the final experiment requested by THM-2877.

## Verdict

The positive address completion exists exactly, and the entire
THM-2868/2874 coefficient atlas can be recovered after one particular
noncanonical source rechart.  Nevertheless there is no
named-right-origin-preserving, source-induced, physically typed completion
along this route.

There are three distinct outcomes.

1. **Pull the target selector back.**  The canonical q3 two-origin selector
   pulls back through `g` to `((12,5),(11,5))`, whose q0 occupancy is
   `(1,0)`.  It reproduces exactly the nonzero `P C_m` bank on all 26
   lawful samples, including `U=chi_3`, `V=1`, the Kummer ratio, and the
   THM-2874 Hermitian support `{0,3,10}`.  This changes the marked
   right-endpoint selector provenance: it is not the named q0 pair
   `((0,0),(12,0))`.  The physical source atom and scalar `P` are unchanged.
2. **Push the source selector forward.**  The named q0 pair maps to
   `((1,8),(0,8))`.  Both target charts are present and have the same
   coefficient, so their signed difference is zero on every one of the 26
   samples.
3. **Use the missing fibre.**  The missing target fibre is
   `F={9} x B3`, while
   `g^(-1)(F)={8} x B0` is wholly absent from q0.  Consequently every
   source-induced full-current term on `F` is zero, independently of how
   `F` is polarized.  Positive and address-`chi_3` endpoint atlases become
   nonzero only if one freely duplicates the source scalar `P`; that is a
   new source copy, not transport by `g`.

The physical gate is sharper still.  The q0 root ancestry has counts
`(966606,28534)`, whereas q3 has `(0,28535)`: its `U` factor is empty.
The `QA/QAB` semantic columns are zero in both q0 and q3 rows.  Thus the
one-fibre completion neither supplies a physical ancestry action nor
approaches the q11/q7 macro-truth attachment.

## 1. Exact finite universe

The companion works on:

```text
endpoint address plane:       F_13^2, 169 points;
source support:               81 points;
q3 target support:            90 points;
missing horizontal fibre:      9 points;
physical intervals:            4
  (q0 source/target and q3 source/target);
lawful raw Prony samples:      26;
formal frequency sections:     13;
certified endpoint fields:      2;
literal ancestry factors:       U and V;
semantic columns:              origin00, origin12, QA, QAB.
```

The relevant physical intervals are

```text
I_0^src = (142004992589460, 142005019034340),
I_0^tgt = (142005424522500, 142005450967380),

I_3^src = (210736584398340, 210736610843220),
I_3^tgt = (210737016331380, 210737042776260).
```

The q3 target has 56 full native-factor cells, and the source and target
weighted atoms both have weight `27581135604`.  All four displayed
intervals lie in macro `E3`.  Therefore the stopping point below is not
absence of the allocation carrier, a native-factor failure, or an
`E3/not-E3` crossing.

## 2. Address support and the origin obstruction

Write

```text
A0={0,1,2,3,4,5,6,7,12},
B0={0,1,2,3,4,5,8,9,10},

A3={0,1,2,3,4,5,6,7,8,9},
B3={0,3,4,5,8,9,10,11,12}.
```

For

```text
g(a,b)=(a+1,b+8),
```

the canonical endpoint evaluator gives the disjoint decomposition

```text
A3 x B3
 =g(A0 x B0) disjoint-union F,

F={9} x B3,                         |F|=9,              (1)

g^(-1)(F)={8} x B0.
```

This is the exact positive one-fibre completion.  It is not a permutation
of occupied supports: its new target fibre has no occupied source
preimage.

Let

```text
O=((0,0),(12,0))
```

be the named pair.  The four occupancy rows are

```text
q0 at O:                  (1,1),
q3 at O:                  (1,0),
q3 at g(O):               (1,1),
q0 at g^(-1)(O):          (1,0),

g(O)=((1,8),(0,8)),
g^(-1)(O)=((12,5),(11,5)).                              (2)
```

The fibre `F` meets neither `O` nor `g(O)`.  Thus adding any coefficient
supported on `F` cannot change either named-selector value.

Equation (2) also gives a map-independent hostile: no positive support
injection fixing the named pair can exist, because `(12,0)` is occupied
at q0 and absent at q3.  The failure is not an unfortunate choice of the
particular affine witness.

## 3. The coefficient bank needs no rechart

The q0-to-q3 target displacement is three allocation units.  Exact
endpoint exponent arithmetic gives

```text
12 R_dil (3 unit)=0 mod N,
26 R_dil (3 unit)=0 mod N.                              (3)
```

Hence both endpoint weights and both Prony nodes are literally identical
at q0 and q3.  In particular

```text
C_m(q0)=C_m(q3) !=0                                    (4)
```

for every raw sample in

```text
(1,2), (43,44), (85,86), (127,128), (170,171),
(211,212), (253,254), (295,296), (339,340),
(379,380), (421,422), (463,464), (505,506).
```

The q0 chart `(12,5)` restricts to the same physical source atom as the
canonical chart, so its x-sweep coefficient is the same pinned common
nonzero `P`.  The chart `(11,5)` is empty.  Consequently pulling the target
selector back gives

```text
P[C_m(q0,(12,5))-C_m(q0,(11,5))]
 =P C_m
```

on all 26 samples, exactly the THM-2868 q3 current bank.

The variable-offset Prony splits replay without alteration:

```text
U_(r+1)=omega^3 U_r,
V_(r+1)=V_r,
t_r=U_r/V_r,
```

and the normalized Hermitian edge has thirteen distinct nonzero values and
Fourier support

```text
{0,3,10}.                                                (5)
```

By contrast, pushing the q0 named selector forward gives

```text
C_m(q3,(1,8))-C_m(q3,(0,8))=C_m-C_m=0                  (6)
```

termwise on the full bank.

This is a right-selector provenance obstruction, not a source-origin or
Prony obstruction.  Equation (4) shows that no cyclotomic automorphism,
Galois rechart, section shift, or node transport is needed.

## 4. Arbitrary missing-fibre polarization is endpoint-only

Let `p_x` denote the source coefficient on address `x`: it is `P` on the
q0 support and zero off that support.  Let `w_y` be arbitrary scalar
weights on `F`.  Since every `g^(-1)y` with `y in F` lies in the absent
set `{8} x B0`,

```text
sum_(y in F) w_y p_(g^(-1)y) C_m(y)=0                  (7)
```

for every `m`.  This is pointwise before summation, so (7) holds for every
positive, signed, cyclotomic, or character polarization.

Two useful endpoint-only controls are nonzero:

```text
positive fibre:
  sum_(y in F) C_m(y)=9 C_m;

address-character-three fibre:
  sum_((9,b) in F) omega^(-3b) C_m
   =beta_3 C_m,            beta_3 !=0.                  (8)
```

Attaching a free copy of `P` to (8) carries the full Prony/Kummer/Hermitian
atlas by the fixed scalar `beta_3`.  But this replaces the zero
`p_(g^(-1)y)` in (7) by `P` by decree.  It duplicates source provenance
and is not induced by the address injection or the physical source
support.

Thus (8) is a clean positive control: the fibre is coefficient-capable.
Equation (7) identifies exactly why it is not current-capable.

## 5. Address `chi_3` and coefficient `chi_3` are separate

For the normalized address Fourier moment

```text
M_3(S)=sum_((a,b) in S) omega^(-3b),
```

translation by `(+1,+8)` gives

```text
M_3(gS0)=omega^2 M_3(S0),

M_3(gS0)=9 beta_3,
M_3(F)=beta_3,
M_3(E_q3)=10 beta_3.                                   (9)
```

Thus `g` preserves the vertical character-three **line** but changes its
normalization by `omega^2`.  Multiplication by `omega^(-2)` removes this
phase.  Since `omega` belongs to `F=Q(zeta_91)`, this is a common
`F`-rational scalar gauge.  It scales `U` and `V` together, preserves their
`chi_3/trivial` laws, leaves `t=U/V` and (5) unchanged, and cannot turn
either zero in (6)--(7) into a nonzero value.

No identification of the address Fourier axis with THM-2868's multiplier
frequency axis is used.  Equation (9) is an address-support statement;
equations (4)--(5) are physical coefficient statements.

## 6. Physical and semantic stopping gate

The literal q0 and q3 ancestry banks are

```text
q0: (|U|,|V|)
    =(966606,28534),
    digest 15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd;

q3: (|U|,|V|)
    =(0,28535),
    digest bdf95bb641f8d0ab0d8ebe764fecda782b9a7fe4995031c7f23e2f02868de931.
                                                                    (10)
```

Source and target frames agree within each q-row, but q3 has no `U`
labels at all.  The promoted THM-2859 strict forest audit independently
shows that the q3 pulled intervals are not same-labelled forest vertices.
Address reindexing acts on neither interval midpoint nor ancestry label, so
it cannot repair (10).

The THM-2847 coefficient-support frame has rows

```text
q0: (origin00,origin12,QA,QAB)=(1,1,0,0),
q3: (origin00,origin12,QA,QAB)=(1,0,0,0).               (11)
```

Therefore neither the pulled selector nor the missing fibre supplies a
`QA/QAB` transition.  Since both q0 and q3 remain in macro `E3`, this route
also never meets the complementary q7 block.

The route-specific first losses are now exact:

```text
preserve named source selector:
  coefficient cancellation (6);

pull back target selector:
  nonzero full coefficient atlas, but changed right-selector provenance,
  then empty q3 U ancestry (10);

use the missing fibre:
  source coefficient zero pointwise (7).
```

## 7. Cross-test with the directed u3 carry marker

The newly isolated factor-3/u3 marker does not change these conclusions.
On the complete target residue orbit, its state word is

```text
q=0,...,12:  SSSSSSSSSSSSD.
```

Thus u3 is dangerous only at q12 and its unique positive `D->S` unit exit
is

```text
q12 -> q0.                                               (12)
```

The present q0-to-q3 route is the nonwrapping path

```text
q0 -> q1 -> q2 -> q3
```

and has zero events in (12).  Furthermore every one of the 169 endpoint
representatives has `ell_u3=0`, so the event is independent of the address
map, selector, and missing fibre.

In particular, the address coordinate `12` in the pulled chart `(12,5)`
is not target residue `q12`.  Conflating those labels would create a false
connection.  The u3 exit is a genuine positive carry marker for wrapped
residue paths and may help base the THM-2851 triangle, but it supplies
neither the absent q3 `U` ancestry nor the zero `QA/QAB` columns in
(10)--(11).  It is orthogonal to all three routes tested here.

## 8. Theorem-ready formulation

> **q0/q3 one-fibre selector-provenance theorem.**  
> On the exact THM-2877 q0/q3 endpoint supports, the affine map
> `g(a,b)=(a+1,b+8)` has target complement `{9} x B3`.  Pullback of the
> canonical q3 right-endpoint selector reproduces the complete nonzero
> THM-2868/2874 coefficient atlas without coefficient rechart, but changes
> the q0 right selector to `((12,5),(11,5))`, while retaining the physical
> source atom and `P`.  Pushforward of the canonical q0
> selector cancels termwise.  Every source-induced coefficient on the
> missing fibre vanishes pointwise because its full preimage is q0-absent.
> The q3 `U` ancestry is empty and its `QA/QAB` columns are zero.  Hence
> this one-fibre route has no named-right-origin-preserving, source-induced
> physical completion, although both its pulled-selector and free-fibre
> coefficient shadows retain the exact projective/Hermitian atlas.

The conclusion is scoped to the transporter induced by `g` and its
one-fibre completion.  A broader correspondence that freely duplicates
source mass is a different object and must supply its own physical,
ancestry, and positivity law.

## 9. Connection contract

```text
source:
  THM-2877 q0 endpoint rectangle and the THM-2868 common source scalar P;

target:
  q3 endpoint rectangle, named two-origin selector, and 26-sample
  Prony/Kummer/Hermitian atlas;

map:
  g(a,b)=(a+1,b+8), followed separately by selector pushforward,
  selector pullback, or the complement fibre F;

preserved:
  address support inclusion, vertical character-three line, both endpoint
  nodes, all 26 raw q0/q3 coefficients, U=chi3, V=trivial, projective
  ratio, Hermitian support {0,3,10}, native q3 carrier, and macro E3;

destroyed / missing:
  simultaneous named-right-origin equivariance, source provenance on F,
  q3 U ancestry, same-labelled forest vertex, QA/QAB transition,
  q7 complement support, and a physical action;

cheapest surviving test:
  do not revisit this affine fibre.  A new route must either give a
  source-backed positive correspondence to F with an explicit no-duplication
  law and nonempty U/QA/QAB typing, or use the independent u3 wrapped-carry
  marker to base a genuine q11->q7 ancestry transition before coupling it
  to the endpoint coefficient atlas.
```

## Reproduction

Run from the repository root:

```bash
python3 .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.py \
  > .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.out
python3 -O .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.py \
  > .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.opt.out
cmp .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.out \
    .scratch/lrc_q0_q3_one_fibre_completion_20260729/audit.opt.out
```

The companion contains no executable Python `assert` statements.  Normal
and optimized outputs are byte-identical.  LF-normalized SHA-256:

```text
script  f69ce49927cd5a089ec1571fd8d40c07f7037644536c4967e2a5c4fa939bc1bd
output  7543068b42551ef8b0cb0a2d2c16e4b73ef268c8084e8e26a41d8066cb79eded
```
