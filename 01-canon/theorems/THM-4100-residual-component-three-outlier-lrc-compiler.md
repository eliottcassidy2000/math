---
id: THM-4100
title: "Residual-component three-outlier LRC compiler"
status: >
  PROVED + PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED.
  If a closed antipodal-safe interval J has length L and b<c<d are three
  new speeds, then J survives whenever
  1/b+min(4/c,2/c+4/d)<=7L. A surviving rational arrangement endpoint
  gives an adaptive even clock N<=max(Q,14d) with E_N=R_N=empty. In
  particular every eleven-core {1,...,8,b,c,d}, 93<=b<c<d, of arbitrary
  outlier parity closes every dyadic two-odd-tail seam through all common
  dilations. This is a structured sufficient supplier, not arbitrary-core
  AP extraction, a converse, or LRC(14).
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4087-lrc-literal-open-two-comb-compiler
  - THM-4092-parity-weighted-antipodal-k-comb-density-compiler
related:
  - THM-1017-ap-core-bridge-reduction
  - THM-4081-lrc-antipodal-height12-obstruction-and-six-speed-floor
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
  - THM-4101-ap7-weight-seven-gap-four-second-moment-absorption
  - THM-4112-antipodal-component-ancestry-chain-and-scale-separated-lrc-families
  - THM-1094-exact-two-comb-component-theorem
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - MISTAKE-126
  - MISTAKE-274
script: 04-computation/lrc_residual_component_three_outlier_thm4100.py
output: 05-knowledge/results/lrc_residual_component_three_outlier_thm4100.out
independent_audit_script: 04-computation/lrc_residual_component_three_outlier_thm4100_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_residual_component_three_outlier_thm4100_independent_audit.out
script_sha256: e51c47f058a0d27ad231891db788f8bd45d6eafbe7758bbe5881c6fbca943b4f
output_sha256: 9b4ee3b390f4c53a2e04e6b002d6d2be9d541fca8c5778ce599a840f04ceda46
independent_audit_script_sha256: de03fe6b491487193e32a373687db3553ed30115e5f4ce6ba1b89f544dad8eba
independent_audit_output_sha256: f5ed4898094eb7191662fa3970c8a9ff05f2430fdadf937c6470efc2fa782837
independent_semantic_sha256: f596182653208fc0b34d330ee20b4404e79602d14692a0d67a86f0f9e4370796
hash_basis: raw LF bytes
audit: >
  PASS. A Fraction-only independent implementation reproduces both pair
  bounds, all 73,066 pair and 1,390,420 circular triple components, the AP8
  cutoff and clocks, the spread gain, both scalar hostiles, and the THM-4101
  scope separation. Its midpoint seam gauge and the primary left-end gauge
  select tied maximum components of the same length; clipping accounts for
  exactly one spurious extra piece per triple. Normal/-O streams byte-match.
---

# THM-4100 -- the residual-component three-outlier compiler

Put `delta=1/14`. For a finite set `D` of positive integer speeds, define

```text
G_D^+-={theta in R/Z:
          ||v theta||>=delta and ||v(theta+1/2)||>=delta
          for every v in D},                                  (1)

U_v={theta in R/Z:
          min(||v theta||,||v(theta+1/2)||)<delta}.             (2)
```

Thus `U_v` is the literal **open** antipodal danger comb and `G_D^+-` is
`(R/Z) minus union_(v in D) U_v`. THM-4092 identifies its exact parity
anatomy. With

```text
omega(v)=1  if v is even,
omega(v)=2  if v is odd,                                     (3)
```

the teeth of `U_v` have centres, length, period, and intervening closed
safe-gap length

```text
centres:        k/(omega(v)v),
tooth length:   a_v=1/(7v),
period:         p_v=1/(omega(v)v),
safe gap:       g_v=p_v-a_v
              =6/(7v)    if v is even,
              =5/(14v)   if v is odd.                        (4)
```

In particular every `v`-safe gap has length at least `5/(14v)`. The new
operation is to retain the **component envelope** of the faster residual,
rather than compressing all three combs to one parity-weighted mass total.

## 1. Two exact component bounds

Let `c<d`. THM-4087 proves the literal-open bound

```text
every component of U_c union U_d has length < 2/(7c).         (5)
```

There is a simultaneous direct span bound:

```text
every component of U_c union U_d has length
  < 1/(7c)+2/(7d).                                           (6)
```

To prove `(6)`, a `d`-tooth cannot meet two different `c`-teeth: its length
`1/(7d)` is smaller than the least gap `5/(14c)` between consecutive
`c`-teeth. Hence a component contains at most one `c`-tooth. A component
with no `c`-tooth is one `d`-tooth and already satisfies `(6)`. Otherwise,
every attached `d`-tooth meets the one `c`-tooth. Strict overlap makes the
leftmost protrusion strictly shorter than one `d`-tooth, and likewise on the
right. Adding the `c`-tooth length proves `(6)`. Touching open endpoints do
not join components; the touching point is weak-safe.

Combining `(5)` and `(6)`, put

```text
P(c,d)=min(2/(7c), 1/(7c)+2/(7d)).                         (7)
```

Every component of `U_c union U_d` has length strictly less than `P(c,d)`.
The second branch improves the first exactly when `d>2c`.

## 2. The three-comb residual-component theorem

Let `D` be any finite set of positive integers. Let `J` be a closed circular
interval with a chosen real lift `[alpha,alpha+L]`, where `L>0`, and assume

```text
J subset G_D^+-.                                            (8)
```

Let `b<c<d` be distinct positive integers outside `D`. Define

```text
T(b,c,d)=1/(7b)+2P(c,d)
        =1/7 * (1/b+min(4/c,2/c+4/d)).                     (9)
```

> **General interval theorem.** If
>
> ```text
> boxed: 1/b+min(4/c,2/c+4/d) <= 7L,                      (10)
> ```
>
> then
>
> ```text
> J minus (U_b union U_c union U_d) is nonempty.            (11)
> ```
>
> Consequently `D union {b,c,d}` has an antipodal-safe point.

### Proof of the component bound

Let `H=U_c union U_d`. By `(5)` and `b<c`, every component of `H` has
length

```text
<2/(7c)<2/(7b)=4/(14b)<5/(14b)<=g_b.                      (12)
```

It therefore cannot meet two consecutive `b`-teeth: doing so would make its
connected lift span the entire intervening `b`-safe gap. Thus every component
of `U_b union H` contains at most one `b`-tooth.

If it contains none, it is an `H`-component and is shorter than `P(c,d)`.
If it contains one, all of its `H`-components meet that tooth. The extreme
left and right protrusions are each shorter than `P(c,d)`. Since the tooth
itself has length `1/(7b)`, every component obeys

```text
boxed:
component_length(U_b union U_c union U_d)
  <1/(7b)+2P(c,d)
  =1/(7b)+min(4/(7c),2/(7c)+4/(7d)).                       (13)
```

If the open union in `(11)` covered the connected closed interval `J`, then
`J` would lie in one of its open components. That component would have
length at least `L`, contradicting `(10)` and the strict inequality `(13)`.
This proves `(11)`. Equality is allowed in `(10)` because `(13)` is strict.
**QED.**

The map and its indispensable sidecar are now explicit:

```text
source:       one closed antipodal-safe interval J
operation:    group U_c union U_d, then reinsert U_b
preserved:    literal weak two-phase safety
lost by mass: which residual component meets which b-tooth
sidecar:      ordered component envelope P(c,d) and parity gap g_b
target:       a surviving arrangement endpoint and adaptive owner clock.   (14)
```

## 3. Rational survivor and adaptive even clock

Assume the two endpoints of the chosen lift of `J` admit rational
presentations with positive even denominators at most `Q`. The nonempty
closed set in `(11)` is a finite union of intervals and points. Choose an
endpoint of one of its components. It is either an endpoint of `J` or a
tooth endpoint

```text
(14k +/- 1)/(14v)       if v is even,
( 7k +/- 1)/(14v)       if v is odd,                       (15)
```

for `v in {b,c,d}`. Hence the survivor has an even presentation

```text
theta=r/N,        N<=max(Q,14d).                            (16)
```

For completeness, retain the inherited labelled clock objects. Write
`|x|_N=min(x mod N,N-(x mod N))` and put, for `C=D union {b,c,d}`,

```text
A_N(C)={0<=s<N: 14|vs|_N>=N for every v in C},              (17)

E_N(C)={odd z mod 2N:
          7|zs|_N<N for every s in A_N(C)}.                 (18)
```

THM-2066's owner-reversal relation satisfies

```text
R_N(C) subset E_N(C)^2.                                    (19)
```

Both labels `r` and `r+N/2 mod N` lie in `A_N(C)`. For every odd class `z`,

```text
|z(r+N/2)|_N=N/2-|zr|_N.                                  (20)
```

The two distances cannot both be strictly below `N/7`, since their sum is
`N/2>2N/7`. Eligibility in `(18)` requires both strict inequalities.
Therefore

```text
boxed: E_N(C)=R_N(C)=empty.                                (21)
```

When `|C|=11`, THM-2061/2066/2072 turn `(21)` into the dyadic seam
conclusion. For every positive integer `q` and every two distinct positive
odd integers `x,y`,

```text
2qC union {x,y}                                             (22)
```

is `1/14`-lonely. The common-dilation transport is forward and exact: use
`theta/q`; an odd `q` preserves the old half phase, while an even `q`
collapses it to the old zero phase. Both phases were safe. The clock becomes
`qN`. This is not an obstruction equivalence under even dilation.

## 4. AP8 plus three arbitrary-parity outliers from height 93

THM-4092 proves the exact interval

```text
J_8=[11/49,13/56] subset G_{ {1,...,8} }^+-,
|J_8|=3/392,             Q_8=98.                            (23)
```

Thus `(10)` becomes

```text
1/b+min(4/c,2/c+4/d) <= 3/56.                              (24)
```

If `93<=b<c<d`, then `c>=b+1` and

```text
1/b+min(4/c,2/c+4/d)
 <=1/b+4/c
 <=1/b+4/(b+1)
 <=1/93+4/94
 =233/4371
 <3/56,                                                     (25)

3/56-233/4371=65/244776.                                  (26)
```

No parity hypothesis on `b,c,d` occurs. Therefore every eleven-core

```text
boxed: C_(b,c,d)={1,...,8,b,c,d},       93<=b<c<d,          (27)
```

has an antipodal-safe point and an even adaptive clock

```text
N<=14d,       E_N=R_N=empty.                               (28)
```

The uniform threshold `93` is least for criterion `(24)`: at the next
ordered triple `(92,93,94)`, its left side exceeds `3/56` by

```text
37/119784.                                                 (29)
```

This is a boundary of the sufficient criterion, not a hostile to global
antipodal feasibility. The sharpened `d`-dependent branch is material. For

```text
(b,c,d)=(90,94,189),                                      (30)
```

the old envelope fails but the new one succeeds:

```text
1/90+4/94=227/4230 > 3/56,
1/90+2/94+4/189=4757/88830 < 3/56,
3/56-4757/88830=1/50760.                                  (31)
```

The exact witness is again `theta=11/49` at clock `N=98`.

### Exact comparison with THM-4092

THM-4092's AP8 density thresholds depend on the number `o` of odd outliers.
The present component theorem gives one threshold for every parity pattern:

| odd outliers `o` | THM-4092 threshold | THM-4100 threshold | stronger row |
|---:|---:|---:|---|
| 0 | `85` | `93` | THM-4092 |
| 1 | `106` | `93` | THM-4100 |
| 2 | `150` | `93` | THM-4100 |
| 3 | `281` | `93` | THM-4100 |

The mechanisms are complementary. THM-4092 retains parity weight and a
reciprocal boundary tax, but forgets component ancestry. THM-4100 retains
the ordered residual component and thereby crosses the expensive odd-parity
rows. It does not supersede THM-4092's all-even threshold `85`, its exact
nonuniform reciprocal criterion, or its AP5--AP7 families with four through
six outliers. It extends THM-4087 in a different direction: AP9 plus two
outliers becomes AP8 plus three.

THM-4101 is the adjacent but disjoint boundary advance. It absorbs **four**
outliers over AP7 at total parity weight seven by retaining a selected odd
gap-four pair, the split `theta=1/8`, and a second-moment overlap credit; its
uniform height is `264` and its declared parity/resonance hypotheses are
load-bearing. THM-4100 instead absorbs **three** arbitrary-parity outliers
over AP8 by retaining literal component ancestry, with uniform height `93`
and no pair resonance. Neither theorem contains the other. Their genuine
common lesson is that THM-4092's scalar first moment is not the last useful
carrier: THM-4100 repairs it by component localization, while THM-4101
repairs its weight-seven boundary by labelled overlap.

## 5. The first scalar failures and their missing sidecars

### Same weight, opposite antipodal feasibility

Let

```text
O=[12] minus {10},        S=[12] minus {4}.                  (32)
```

Both are primitive eleven-cores of height `12`, with six odd and five even
speeds, hence identical THM-4092 total weight

```text
sum_v omega(v)=17.                                         (33)
```

Nevertheless `O` contains THM-4081's obstruction

```text
D*={1,3,4,5,7,8,9,11,12},                                 (34)
```

whereas `S` is antipodal-safe at

```text
theta=18/77.                                               (35)
```

This is minimum-height exact: THM-4081 proves there is no obstruction of
height at most `11` and classifies height `12`. It is not a globally
minimum-cardinality claim, because antipodal obstruction cardinalities seven
and eight remain open.

The first failed implication is

```text
parity total + sorted order  -/->  antipodal feasibility.    (36)
```

Under the exact coordinate `t=2theta`, odd `v` contributes multiplier `v`
at threshold `1/7`, while even `v` contributes multiplier `v/2` at threshold
`1/14`. Speed `10` is only a weaker shadow of the already present odd speed
`5`; speed `4` supplies the genuinely new multiplier `2`. The missing
sidecar is transformed-multiplier support plus its new/shadow incidence with
the residual safe components.

### Same Cover14 scalars, opposite AP deletion shape

Compare

```text
A=[12] union {182},
C=([13] minus {3}) union {182}.                            (37)
```

Both are primitive thirteen-speed Cover14 families, have maximum `182`, six
odd and seven even speeds, and total weight `19`. Only `A` has maximum-
deletion core `[12]`. A complete exact tent-wall/pair-intersection sweep gives

```text
M(A)=14/183 at theta=14/183,       owners (1,182),
M(C)=2/17   at theta=6/17,         owners (6,11).             (38)
```

Thus

```text
Cover14 + maximum + parity total + sorted order
  -/-> AP maximum-deletion core.                             (39)
```

The missing sidecar is the normalized maximum-deletion support/gap word.
This collision does **not** refute THM-1017's open extraction target:
`C` has `M(C)>1/13` and is outside that target's deep residual. It diagnoses
the first scalar loss before the full hypothesis; AP extraction under
`M<1/13` remains open.

## 6. THM-4088 typing firewall

Sorting any distinct real labels produces a transitive order tournament.
Accordingly, `O` and `S` have the same rank-order fingerprint with `55`
edges, while `A` and `C` have the same one with `78` edges. THM-4088 proves
that such a bare order tournament is blind not only to these metric/support
coordinates but even to rational, algebraic-irrational, and transcendental
type. It is used here only as a **procedural schedule** for `b<c<d`.

No tournament invariant enters the proof of `(13)`. The load-bearing data
are the numerical component length, the parity safe-gap, and the endpoint
owner. Likewise THM-4088's strict-margin theorem says that every positive
witness interval contains times of every arithmetic type. Arithmetic type
cannot select the survivor; the adaptive rational endpoint is chosen only to
compile the owner clock.

The controlled-forgetting ledger is:

| source | map | preserved | destroyed | required sidecar |
|---|---|---|---|---|
| `U_c union U_d` | one scalar mass | total occupied length | component span/ancestry | `P(c,d)` |
| sorted outliers | transitive order tournament | order only | gaps, parity, support, owners | reciprocal speeds and `g_b` |
| survivor component | chosen endpoint | weak safety | interior arithmetic type | rational endpoint owner/clock |
| maximum-deletion core | `(W,max,Cover14)` | coarse scalars | AP support word | normalized sorted gaps |

## 7. Exact referee and correction of the seam census

Reproduce with

```bash
python3 -B 04-computation/lrc_residual_component_three_outlier_thm4100.py
python3 -B -O 04-computation/lrc_residual_component_three_outlier_thm4100.py
python3 -B 04-computation/lrc_residual_component_three_outlier_thm4100_independent_audit.py
python3 -B -O 04-computation/lrc_residual_component_three_outlier_thm4100_independent_audit.py
```

The assertion-independent `Fraction` referee uses only `require` gates,
literal open teeth, strict-overlap merging, and exact integer owner
inequalities. Its universes and results are:

```text
two-comb rows                              1,225
literal circular two-comb components      73,066
maximum pair-min-bound ratio               72/73
  row                                      (16,41)
  component                                [139/574,74/287]

three-comb rows                            19,600
literal circular three-comb components    1,390,420
sharpened d-span branch rows               4,600
maximum three-comb bound ratio             1050/1081
  row                                      (1,42,47)
  component                                [279/658,379/658]

AP8 rows, 93<=b<c<d<=120                  3,276
odd owner classes                          2,227,526
failures                                   0
maximum adaptive clock                     1,680
arrangement source sides                   left 1,771; right 1,505.  (40)
```

The boundary row `(93,94,95)` has

```text
theta=11/49=22/98,       half label=71,
|A_98|=34,               E_98=empty,
length-budget slack=65/1713432.                            (41)
```

An initial scratch report counted `1,410,020` three-comb components by
clipping a seam-crossing circular component into two lifted pieces. That is
exactly one extra piece for each of the `19,600` triples. The canonical
circular census in `(40)` selects one full lifted representative by its left
endpoint and is

```text
1,410,020-19,600=1,390,420.                               (42)
```

This correction changes no inequality, witness, or theorem conclusion. It is
recorded in the frozen output rather than silently replacing the exploratory
number. The independent path chooses the unique full lift whose midpoint lies
in `[0,1)`. At the maximizing row `(1,42,47)` it selects the tied seam
component `[-25/329,25/329]`, while the primary left-end gauge selects the
component in `(40)`. Both have length `50/329`; this is a representative tie,
not a count or bound discrepancy. The theorem itself is carried by the
analytic literal-open proof in Sections 1--3. Normal and optimized outputs are
byte-identical on both implementations.

Frozen raw-byte hashes:

```text
script  e51c47f058a0d27ad231891db788f8bd45d6eafbe7758bbe5881c6fbca943b4f
output  9b4ee3b390f4c53a2e04e6b002d6d2be9d541fca8c5778ce599a840f04ceda46
independent script de03fe6b491487193e32a373687db3553ed30115e5f4ce6ba1b89f544dad8eba
independent output f5ed4898094eb7191662fa3970c8a9ff05f2430fdadf937c6470efc2fa782837
```

## 8. Scope

THM-4100 is an arbitrary-`D` **conditional** interval compiler and a new
infinite AP8 three-outlier seam family. It does not prove that an arbitrary
eleven-core contains the interval `(8)`, does not extract AP structure from a
compact Cover14 residual, and does not promote parity weight or a transitive
order tournament into a structural invariant. It gives no converse to
`(10)`, no claim that `93` is a global AP8 threshold, and no result for an
unstructured fourth outlier.

THM-1017's AP extraction target and the residual-to-AP supplier in the live
proof map remain open. Consequently **LRC(14) remains OPEN.**
