---
id: THM-4112
title: "Antipodal component-ancestry chain and scale-separated LRC families"
status: >
  PROVED + PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. A parity-aware one-tooth/two-protrusion lemma
  iterates to arbitrary depth under explicit parent-gap gates. Ratio-two
  chains of depth at most six and every finite chain with adjacent ratio at
  least 12/5 give exact scale-separated suppliers. This yields AP7
  four-outlier eleven-cores, a primitive non-AP core plus six outliers, and
  AP8 plus five outliers; the latter two are direct thirteen-speed LRC(14)
  families. Arbitrary cores and LRC(14) remain open.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4092-parity-weighted-antipodal-k-comb-density-compiler
  - THM-4100-residual-component-three-outlier-lrc-compiler
related:
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-1233-all-prefix-component-span-compactification
  - THM-4098-weight-seven-antipodal-scale-escape-and-missing-parity-rows
  - THM-4101-ap7-weight-seven-gap-four-second-moment-absorption
  - THM-4109-ap7-weight-seven-gap-atlas-and-sharp-pair-overlap
  - THM-4110-sparse-reciprocal-phase-graph-saturation-and-ap13-torsion-tariff
  - THM-4117-physical-eleven-plus-two-primitive-support-obstruction
  - THM-4119-infinite-supplier-free-eleven-plus-two-residue-family
  - THM-4121-sharp-projective-tail-multiplier-residue-compiler
  - MISTAKE-169
  - MISTAKE-274
script: 04-computation/lrc_antipodal_component_ancestry_chain_thm4112.py
output: 05-knowledge/results/lrc_antipodal_component_ancestry_chain_thm4112.out
independent_audit_script: 04-computation/lrc_antipodal_component_ancestry_chain_thm4112_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_antipodal_component_ancestry_chain_thm4112_independent_audit.out
script_sha256: a036113d9918a81390b5d803c0e2c1a1b6e62c9168311e1f3599490e9e990195
output_sha256: cc7da6625a80008cd7644ad6a211f3bd0c3c119ff95ead7a51d22dbc30f76e3e
semantic_sha256: 9834738b13362e78c30f577bd4dc1219221c9c1d8413a18102dc1b56d65cf093
independent_audit_script_sha256: 9a084c2ef41ad54850a3ef0a444f1c35eaa2ffe017c3c35843d9398b7b6a719a
independent_audit_output_sha256: 4ccf4f419156ac2215f427639b8d7de7b138211d34c17043aa5969ce87dd54f3
independent_semantic_sha256: a716c2a1656d979a087674ebd3e92fb3a691e02a87f4bc867092d6b1dffcec6e
hash_basis: raw LF bytes
primary_audit: >
  PASS. Fraction-exact literal-open components, the exact one-comb base,
  symbolic ancestry through depth twelve, 1,474 numeric algebra rows, 1,932
  AP7 adaptive-root rows, the geometric-depth constants, both hostiles, the
  non-AP interval, and both direct thirteen-speed boundary families replay
  identically under normal and optimized Python and match the frozen output.
independent_audit: >
  ACCEPT. A clean-room Fraction-only wall-cell and endpoint-connectivity
  engine imports no primary code and independently reproduces the exact-span
  base, ancestry recursion, both scale regimes, adaptive AP7 split,
  parity-free suppliers, both phases of the non-AP interval, endpoint clocks,
  direct thirteen-speed rows, and both ancestry hostiles. Normal, optimized,
  and frozen outputs byte-match; the smallest theorem failure is none.
---

# THM-4112 -- antipodal component ancestry and scale separation

**PROVED + PROVED RELATIVE TO THM-2061/2066/2072 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

This theorem isolates the operation implicit in THM-4100 and iterates it:
a slow open tooth has only two sides, so a faster residual satisfying the
slow parent-gap gate can contribute at most one component envelope on each
side. The resulting coefficients are `1,2,4,...`. Unlike the blind suffix
estimate in THM-1233, every level retains its actual parent gap; unlike the
scalar density compiler in THM-4092, it retains component ancestry.

The theorem has three consequences that should be kept separate.

1. It gives an arbitrary-depth conditional component recursion.
2. It absorbs any finite number of sufficiently scale-separated outliers,
   and gives sharp finite-depth ratio-two rows.
3. It gives explicit eleven- and thirteen-speed `1/14`-lonely families.

It does **not** prove LRC(14) for an arbitrary thirteen-speed set, extract one
of the displayed cores from an arbitrary residual, or turn failure of one of
its sufficient inequalities into a cover.

## 1. Literal antipodal teeth and the one-parent lemma

Put `delta=1/14`. For a positive integer speed `v`, define the literal open
antipodal danger comb

```text
U_v={theta in R/Z:
       min(||v theta||,||v(theta+1/2)||)<delta}.             (1)
```

As in THM-4092 and THM-4100, put

```text
omega(v)=1 if v is even, and 2 if v is odd.                 (2)
```

The teeth of `U_v` have length, period, and intervening closed safe-gap
length

```text
a_v=1/(7v),
p_v=1/(omega(v)v),
g_v=p_v-a_v
   =6/(7v)   if v is even,
   =5/(14v)  if v is odd.                                  (3)
```

For a proper connected open subset of the circle, its **span** is the length
of any connected real lift. Endpoint contact does not join two literal open
components.

> **Lemma 1 (one tooth, two protrusions).** Let `H` be an open subset of the
> circle and let every connected component of `H` have span at most `B`,
> where
>
> ```text
> 0<B<=g_v.                                                 (4)
> ```
>
> Then every component of `U_v union H` has span strictly less than
>
> ```text
> a_v+2B.                                                   (5)
> ```

### Proof

An `H`-component cannot meet two consecutive `v`-teeth. Any connected lift
that meets the interiors of both has span strictly greater than the complete
intervening gap `g_v`, contradicting `(4)`. Consequently a component of
`U_v union H` contains at most one `v`-tooth: otherwise the portion crossing
one intervening safe gap belongs to one `H`-component meeting both teeth.

A union component containing no `v`-tooth is an `H`-component and has span
at most `B<a_v+2B`. A component containing one tooth consists of that tooth
and `H`-components that meet it. The leftmost protrusion is strictly shorter
than `B`, as is the rightmost protrusion; strict overlap with the open tooth
is what makes both inequalities strict. Their spans plus the tooth length
give `(5)`. Touching endpoints remain safe separators. **QED.**

This is a local component statement. It does not assert that an arbitrary
period window contains a complete safe gap; that false move is MISTAKE-169.
It also does not replace literal connected overlap by closed or almost-
everywhere overlap; that is the failure recorded in MISTAKE-274.

## 2. The arbitrary-depth ancestry chain

Let `v_1,...,v_r` be distinct positive integer speeds. Define envelopes
backwards by

```text
h_r=a_(v_r),
h_i=a_(v_i)+2h_(i+1),                 1<=i<r.              (6)
```

Equivalently,

```text
h_i=sum_(j=i)^r 2^(j-i)/(7v_j).                           (7)
```

> **Theorem 2 (component-ancestry chain).** Assume
>
> ```text
> h_(i+1)<=g_(v_i),                    1<=i<r.              (8)
> ```
>
> Then the components of `U_(v_i) union ... union U_(v_r)`
> have span strictly less than `h_i` for every `i<r`; the final one-comb
> components have span exactly `h_r`.

The distinction at the recursion base is load-bearing: a single open tooth
has span exactly `a_(v_r)`, not strictly less. Starting from that weak base,
Lemma 1 makes every reinserted level strict. This proves Theorem 2 by reverse
induction.

Now let `D` be a finite core, assume `{v_1,...,v_r}` is disjoint from `D`,
and let

```text
J=[alpha,alpha+L] subset
G_D^+-=(R/Z) minus union_(d in D) U_d                   (9)
```

be a closed lifted antipodal-safe interval of length `L>0`. If `(8)` holds
and

```text
h_1<=L,                                                    (10)
```

then

```text
J minus union_(i=1)^r U_(v_i) is nonempty.                 (11)
```

Indeed, if the open union covered connected `J`, then `J` would lie in one
of its open components. For `r>=2` that component has span strictly below
`h_1<=L`. For `r=1`, an open interval of span `h_1` still cannot contain a
closed interval of length at least `h_1`. Thus equality in `(10)` is legal.

The exact connection ledger is

```text
source:       a closed antipodal-safe core interval J
operation:    reinsert literal danger combs from fast to slow
preserved:    weak safety at both theta and theta+1/2
destroyed:    identities and multiplicities of attached faster components
sidecar:      every suffix envelope h_i and parent gap g_(v_i)
target:       a surviving rational arrangement endpoint and owner clock.  (12)
```

## 3. Finite doubling and arbitrary geometric depth

Suppose the outliers are increasing and

```text
v_(i+1)>=2v_i.                                             (13)
```

For a chain of length `r`, `(7)` gives

```text
h_i<=(r-i+1)/(7v_i),
h_(i+1)<=(r-i)/(14v_i).                                   (14)
```

Therefore every parent gap in `(8)` is automatic when `r<=6`, because

```text
h_(i+1)<=5/(14v_i)<=g_(v_i).                              (15)
```

We obtain the finite-depth compiler

```text
r<=6,  v_(i+1)>=2v_i,  r/(7v_1)<=L
    => D union {v_1,...,v_r} is antipodal-safe.            (16)
```

There is also an arbitrary-depth form. Let `lambda>2` and suppose

```text
v_(i+1)>=lambda v_i.                                      (17)
```

The infinite geometric majorants are

```text
h_(i+1)<=1/[7v_i(lambda-2)],
h_i<=lambda/[7v_i(lambda-2)].                             (18)
```

Thus, for every finite depth `r`,

```text
lambda>=12/5,
v_1>=lambda/[7L(lambda-2)]
    => D union {v_1,...,v_r} is antipodal-safe.            (19)
```

At `lambda=12/5`, the tail bound in `(18)` equals the smallest possible
odd parent gap `5/(14v_i)` and the top envelope is at most `6/(7v_1)`.
Equality is allowed. Formula `(19)` is uniform in `r`; it is not an
asymptotic statement and does not hide a finite search.

For the proved AP intervals of THM-4092 and the non-AP interval constructed
below, `(16)` gives the exact threshold ladder

| core `D` | safe interval length `L` | doubling outliers `r` | least `v_1` from `(16)` |
|---|---:|---:|---:|
| `{1,...,7}` | `9/490` | `4` | `32` |
| `{1,...,8}` | `3/392` | `5` | `94` |
| `{3,4,5,6,8,10,12}` | `3/56` | `6` | `16` |

At arbitrary depth and ratio `12/5`, the corresponding first-speed
thresholds are respectively `47`, `112`, and `16`. At ratio `3`, AP8 can
absorb arbitrarily many outliers from `v_1=56`.

## 4. A four-outlier sharpening

The bottom two levels can use THM-4100's sharper pair envelope. For
`x<y<z`, put

```text
P(y,z)=min(2/(7y),1/(7y)+2/(7z)),
Q(x,y,z)=1/(7x)+2P(y,z)
        =1/7[1/x+min(4/y,2/y+4/z)].                        (20)
```

The pair bound gives component span below `P(y,z)`, and

```text
P(y,z)<=2/(7y)<5/(14x)<=g_x.                              (21)
```

Lemma 1 therefore gives component span below `Q(x,y,z)` for the three-comb
residual. Let `s` be a fourth speed, not necessarily the smallest. If

```text
Q(x,y,z)<=g_s,                                             (22)
```

then all components of the four-comb union have span below

```text
a_s+2Q(x,y,z)
=1/7[1/s+2/x+min(8/y,4/y+8/z)].                            (23)
```

Consequently a core interval of length `L` survives these four outliers if

```text
1/x+min(4/y,2/y+4/z)
 <= 6/s       when s is even,
 <= 5/(2s)   when s is odd,                               (24)

1/s+2/x+min(8/y,4/y+8/z)<=7L.                             (25)
```

The first line is a genuine ancestry gate, not a cosmetic hypothesis.

### AP7 plus four outliers

THM-4092 gives

```text
J_7=[4/35,13/98],             |J_7|=9/490.                (26)
```

Let `84<=a<b<c<d`. If `a` or `b` is even, the four outliers satisfy
`(24)--(25)` with an adaptive root:

- if `a` is even, use root `s=a` and tail `(b,c,d)`;
- if `a` is odd, `b` is even, and `b<2a`, use root `s=b` and tail `(a,c,d)`;
- if `a` is odd, `b` is even, and `b>=2a`, use root `s=a` and tail `(b,c,d)`.

In the first case the worst **seven-scaled** survival margin is at `a=84`:

```text
9/70-[1/84+2/85+8/86]=1/8772.                             (27)
```

In the second case the seven-scaled left side is maximized at
`a=85,b=86,c=87`:

```text
9/70-[1/86+2/85+8/87]=650/445179.                         (28)
```

The parent inequality follows from

```text
1/a+4/(b+1)<=6/b,       a<b<2a,                           (29)
```

whose cleared numerator is `b(2a-b-1)+6a>0`. In the final case,

```text
1/b+4/c<5/b<=5/(2a),
1/a+2/b+8/c<1/a+10/b<=6/a<9/70.                           (30)
```

There is also a parity-free, less locally uniform family:

```text
a>=47,       b>=2a,       b<c<d.                          (31)
```

Use root `a`. The parent bound is at most `5/b<=5/(2a)`, and the survival
left side is at most `6/a<=6/47<9/70`. Hence every row in `(31)` survives
`J_7`, irrespective of parity.

The exact controls include

```text
(85,86,91,101),  outlier parity weight 7,
theta=81/707,     even clock 1414;                         (32)

(47,95,97,99),   four odd outliers, parity weight 8,
theta=4/35,       even clock 70.                           (33)
```

The first row has none of THM-4109's selected odd gaps `4,8,12`. Thus it
belongs neither to THM-4092's scalar weight-six row nor to THM-4109's
selected-pair weight-seven atlas.

## 5. A primitive non-AP seven-speed core

Put

```text
D_0={3,4,5,6,8,10,12},
J_0=[1/42,13/168],
|J_0|=3/56.                                                (34)
```

This core is primitive and is not an arithmetic progression. Exact affine
norm ranges on `J_0` are

| speed | range of `||v theta||` on `J_0` |
|---:|---:|
| `3` | `[1/14,13/56]` |
| `4` | `[2/21,13/42]` |
| `5` | `[5/42,65/168]` |
| `6` | `[1/7,13/28]` |
| `8` | `[4/21,1/2]` |
| `10` | `[19/84,1/2]` |
| `12` | `[1/14,1/2]` |

For even speeds the half-shift repeats the same phase. For the odd speeds
`3,5`, the displayed upper endpoints are below `3/7`, so their half-shift
norms are also at least `1/14`. Hence

```text
J_0 subset G_(D_0)^+-.                                     (35)
```

The four-outlier criterion gives, without parity assumptions,

```text
a>=16,       b>=2a,       b<c<d
    => D_0 union {a,b,c,d} is antipodal-safe.               (36)
```

The all-odd row

```text
(15,31,33,35)                                              (37)
```

crosses the scalar ceiling one step earlier: its ancestry slack is exactly
`19/95480`, and `theta=1/42` is safe. The adjacent row `(28,29,30,31)` has
criterion slack `89/170520`; `(26,27,28,29)` misses the same sufficient
criterion by `457/137592`, which is not a cover claim.

## 6. Direct thirteen-speed LRC families

Apply `(16)` at its two load-bearing finite-depth boundaries.

### Non-AP core plus six doubling outliers

For

```text
(v_1,...,v_6)=(16,32,64,128,256,512),                     (38)
```

the ancestry levels are

```text
(h_1,...,h_6)
=(3/56,5/224,1/112,3/896,1/896,1/3584).                   (39)
```

Thus `h_1=|J_0|`. Literal enumeration gives `416` circular danger
components, actual maximum span `1/112`, and the exact survivor

```text
theta=43/1792.                                             (40)
```

More generally, any six outliers with `v_1>=16` and
`v_(i+1)>=2v_i` give a thirteen-speed `1/14`-lonely set with core `D_0`.

### AP8 plus five doubling outliers

THM-4092 gives

```text
J_8=[11/49,13/56],             |J_8|=3/392.                (41)
```

For

```text
(v_1,...,v_5)=(94,188,376,752,1504),                      (42)
```

the levels are

```text
(5/658,1/329,3/2632,1/2632,1/10528),                      (43)
```

and

```text
|J_8|-h_1=1/18424.                                        (44)
```

Literal enumeration gives `1316` circular components, actual maximum span
`1/658`, and survivor `theta=11/49`. More generally, any five outliers with
`v_1>=94` and adjacent ratios at least two give a thirteen-speed
`1/14`-lonely set with core `{1,...,8}`.

These are direct thirteen-speed suppliers; they do not need the dyadic seam
compiler.

## 7. Rational endpoints and eleven-core seams

Suppose the endpoints of `J` have even presentations with denominators at
most `Q_0`. The survivor in `(11)` is a finite union of closed intervals
and points. Choose one of its endpoints. It is either an endpoint of `J`
or a danger-tooth endpoint

```text
(14k plus-or-minus 1)/(14v)   if v is even,
( 7k plus-or-minus 1)/(14v)   if v is odd.                 (45)
```

It therefore has an even presentation `theta=r/N` with

```text
N<=max(Q_0,14 max_i v_i).                                  (46)
```

Both `r` and `r+N/2` are safe owner labels. The complement identity

```text
|z(r+N/2)|_N=N/2-|zr|_N                                   (47)
```

for odd `z` makes THM-4100's endpoint argument apply verbatim:

```text
E_N=R_N=empty.                                             (48)
```

When the resulting core `C` has eleven speeds, THM-2061/2066/2072 turn
`(48)` into the dyadic two-tail family

```text
2qC union {x,y},                                           (49)
```

for every positive integer `q` and every two distinct positive odd speeds
`x,y`. In particular the AP7 four-outlier families above generate explicit
thirteen-speed LRC(14) families through all common dilations.

## 8. Parent-gap hostiles

The parent-gap hypothesis `(8)` cannot be deleted. For

```text
(85,87,89,91),                                             (50)
```

rooting at `85` gives

```text
Q(87,89,91)=437/54201 > g_85=1/238.                        (51)
```

The literal four-comb union has component

```text
(69/1274,46/637),          span=23/1274,                   (52)
```

which meets three `85`-teeth. It exceeds the false no-gap envelope

```text
a_85+2Q=11719/658155                                      (53)
```

by exactly

```text
207559/838489470.                                          (54)
```

Thus the first failed implication is `small residual component => at most
one parent tooth`; the missing coordinate is the parent safe gap.

There is also a sharp obstruction to merely choosing another root. For

```text
(43,45,47,49),                                             (55)
```

the other three combs have a component meeting two consecutive root teeth
for **each** of the four possible roots. Nevertheless `theta=4/35` is safe
for AP7 plus all four speeds. This is a method obstruction, not an
infeasibility witness.

## 9. Exact referee and scope

The primary Fraction-only referee freezes:

- the symbolic `1,2,4,...` recursion through depth `12`;
- every AP7 algebraic boundary and the exact rows `(32)--(33)`;
- the four-parent hostile `(50)` and all-root hostile `(55)`;
- the full affine norm ranges proving `(35)`;
- the D0 four-outlier crossings;
- all `416` and `1316` literal circular components in `(38)` and `(42)`;
- the exact witnesses, clocks, and equality/positive margins.

Run

```text
python3 -B 04-computation/lrc_antipodal_component_ancestry_chain_thm4112.py
python3 -B -O 04-computation/lrc_antipodal_component_ancestry_chain_thm4112.py
python3 -B 04-computation/lrc_antipodal_component_ancestry_chain_thm4112_independent_audit.py
python3 -B -O 04-computation/lrc_antipodal_component_ancestry_chain_thm4112_independent_audit.py
```

The separate AP8 threshold-25 finite census is deliberately absent. No
statement about that census is part of THM-4112. The geometric series in
`(18)--(19)` is a proved analytic majorant, not an extrapolation from the
depth-12 symbolic control.

The independent referee reproduces the open-component semantics, parent-gap
gates, non-AP interval, finite-depth equality boundary, endpoint-clock scope,
and exact hostiles with a separate wall-cell engine. The two implementations
and their frozen outputs agree on every theorem-bearing row. **QED.**
