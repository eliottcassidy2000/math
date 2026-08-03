---
id: THM-3223
title: "Fourth/fifth resonant PRS primitive walls, Pell content, and pivot resurrection"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The fourth and fifth fraction-free reciprocal pivots of the universal
  degree-d resonant factorial pair factor into inherited walls and two new
  primitive irreducibles W13,W20 of degrees 13,20.  Both new walls are 1
  modulo 2, so every fixed offset still has only a finite explicit exceptional
  prime set through row five.  At offset two, p=43 kills the fourth pivot but
  the fifth pivot returns as a unit by an exact whole-row clutching identity;
  the following row then vanishes identically, so the resurrection is sharp
  and lasts exactly one row.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  Two independent hostile audits rederived the factor walls, Rabin
  certificates, fixed-offset quantifiers, Pell/gauge recurrences, p=41 and
  p=43 reductions, and the whole-row clutch with its following terminal
  collapse.  Fresh ordinary and optimized replays byte-match the stored
  transcript and the declared LF-normalized hashes.
depends_on:
  - THM-3217-universal-resonant-degree-prs-wall-atlas-and-fixed-offset-exception-set
related:
  - THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return
  - THM-3214-two-jet-pseudo-division-locality-and-catalan-sharpness
script: 04-computation/factorial_fourth_fifth_prs_primitive_walls_thm3223.py
output: 05-knowledge/results/factorial_fourth_fifth_prs_primitive_walls_thm3223.out
script_sha256: 2afe51afe4719922b739e4aa0ea43fc285a130864f6db1d1ba281430d4e93fdc
output_sha256: 3938603f595f4fff35357ecc51aadb266fb2ea2e5b52e5f8fd8105765cff20e2
hash_basis: LF-normalized bytes
---

# THM-3223 -- fourth/fifth resonant PRS primitive walls, Pell content, and pivot resurrection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3217 identifies the first three degree-parametric reciprocal rows and
their fixed-offset exception divisor.  Two questions remained deliberately
open there: whether the atlas continues without a new structural operation,
and whether a selected pivot wall can be crossed by a later row.  Both have
exact positive answers through row five.

The calculation also exposes its internal clock.  Orders along the universal
`d=1` content divisor are Pell numbers, while a common normalization gauge
follows sums of consecutive Pell numbers.  This is forced by the cubic
homogeneity `E(cf,dg)=c^2dE(f,g)`, not a fit to numerical data.

## 1. Two further fraction-free rows

Use THM-3217's degree-`d` reciprocal pair `a,b`, normalized by `(2d-4)!`,
and its fraction-free operator

```text
P_1(f,g)=f_0g_1-f_1g_0,

E(f,g)=z^(-2){f_0^2g-[f_0g_0+P_1(f,g)z]f}.              (1)
```

Index the rows by

```text
R_(-1)=b,        R_0=a,        R_(k+1)=E(R_k,R_(k-1)),

r=R_1,           s=R_2,        t=R_3,
w=R_4=E(t,s),    x=R_5=E(w,t).                           (2)
```

Put `rho_k=[z^0]R_k`.  By THM-3214, `rho_4` uses precisely the initial
eight-jet and `rho_5` the initial ten-jet.  We work with `d>=12`, so every
one of the eleven required top coefficients exists.

Retain THM-3217's universal wall polynomials `H,J,U`.  Define

```text
W13(d)=
 140928614400d^13-4450810855424d^12+2187898619166720d^11
 -119479225810944000d^10+2401299899413954560d^9
 -26074123979998887936d^8+177312698542783856640d^7
 -809582018591761817600d^6+2562896130982456504320d^5
 -5668624022456998123776d^4+8630773322305700966400d^3
 -8649353474475612069600d^2+5150611778276800836000d
 -1384085705619177185625,                                (3)
```

and

```text
W20(d)=
 6734508720128000d^20-9560308579093708800d^19
 -104590212636279308288d^18+77325440060392923463680d^17
 -10990353099918947941089280d^16
 +585398859813499740975267840d^15
 -16234517220500771292093349888d^14
 +280269110820744266239872860160d^13
 -3319321835643659813075309035520d^12
 +28555759263867535936704384860160d^11
 -184855141179174118241589312618496d^10
 +920280708402038646121095278100480d^9
 -3567288992936033291861681133649920d^8
 +10820035211743102835629862352322560d^7
 -25630500379637406272477883155103744d^6
 +46975014922828624287489671780597760d^5
 -65372851083689767775383071292512000d^4
 +66843975116236841444114398164960000d^3
 -47400207255249718865601589483500000d^2
 +20842506093742030140968421821925000d
 -4284107145365556700983792046921875.                    (4)
```

## 2. Exact fourth and fifth pivots

Let

```text
D4=2^28
 (2d-19)(2d-17)^2(2d-15)^3(2d-13)^4
 (2d-11)^7(2d-9)^10(2d-7)^17(2d-5)^24,                 (5)

D5=2^69
 (2d-23)(2d-21)^2(2d-19)^3(2d-17)^4
 (2d-15)^7(2d-13)^10(2d-11)^17(2d-9)^24
 (2d-7)^41(2d-5)^58.                                    (6)
```

Then the two new pivots factor exactly as

```text
rho_4=w_0=-(d-1)^12 H(d)^4 J(d)^2 W13(d)/D4,             (7)

rho_5=x_0=-(d-1)^29 H(d)^10 J(d)^4 U(d)^2 W20(d)/D5.    (8)
```

### Proof

THM-3217 gives every reversed input coefficient as the rational function

```text
binom(d-2+epsilon,j)
 sum_(ell=0)^j binom(j,ell)d^(j-ell)(-1)^ell
 (2d-4+2epsilon-2j+ell)!/(2d-4)!,                       (9)
```

for `epsilon=0,1`.  Substitute the first eleven entries into five iterations
of `(1)` in the exact field `Q(d)`.  Reducing numerator and denominator to
primitive coprime polynomials gives `(7),(8)`.  This calculation never
inverts a row pivot, so the formulas remain valid after reduction on every
displayed wall.

An independent path starts from the raw multinomial coefficients of the
moment polynomials before reciprocal normalization and checks both pivots at
each integer degree `12<=d<=20`.  The 18 direct controls agree exactly with
`(7),(8)`.

## 3. The new factors are genuine primitive walls

Both `W13,W20` are primitive in `Z[d]`.  Rabin's finite-field criterion gives

```text
W13 irreducible modulo 47,       degree 13,
W20 irreducible modulo 29,       degree 20.               (10)
```

The reductions retain their degrees.  For `W13`, the exact check is
`X^(47^13)=X` modulo `W13` and the degree-one proper-subfield gcd is one.
For `W20`, it is `X^(29^20)=X` together with the degree-ten and degree-four
proper-subfield gcds.  Gauss's lemma proves that both polynomials are
irreducible over `Q`.

More simply but importantly,

```text
W13(d)=W20(d)=1                              (mod 2).     (11)
```

Every nonconstant coefficient in `(3),(4)` is even and each constant is
odd.  Thus neither new polynomial vanishes at an integer.  The two factors
are new irreducible chart divisors, not disguised products of `H,J,U`, and
they create no infinite fixed-offset degeneracy.

## 4. Pell content and normalization gauge

The exact orders of the first five pivots along `d=1` are

```text
ord_(d-1)(rho_1,...,rho_5)=1,2,5,12,29.                  (12)
```

These are `P_1,...,P_5` for the Pell recurrence

```text
P_0=0,             P_1=1,             P_(k+1)=2P_k+P_(k-1). (13)
```

There is no hidden additional `(d-1)` factor: the complementary wall
polynomials and odd-linear denominators in the displayed formulas are all
nonzero at `d=1`.

For a common rescaling of the two initial reciprocal rows by `u`, let `e_k`
be the exponent on `R_k`.  The exact homogeneity of `(1)` gives

```text
e_(-1)=e_0=1,              e_(k+1)=2e_k+e_(k-1),

(e_(-1),e_0,e_1,...,e_5)=1,1,3,7,17,41,99.              (14)
```

Thus the row gauge is the companion Pell clock.  Equations `(12)--(14)` do
not assert an unproved all-depth factor formula; they identify the exact
mechanism and verified content through row five.

## 5. Every fixed offset still has a finite row-five atlas

Let `Xi_s` be THM-3217's nonzero integer controlling its eight coordinates,
and define

```text
Xi_s^[5]=Xi_s W13(s)W20(s)
          (2s-19)(2s-21)(2s-23).                         (15)
```

For every integer `s>=2`, equation `(11)` and parity of the three new linear
factors give

```text
Xi_s^[5] !=0.                                             (16)
```

Fix `s>=2`, put `d=p+s`, and suppose

```text
p>max(11,2s-4),                    p does not divide Xi_s^[5]. (17)
```

Every factor in `(7),(8)` reduces modulo `p` to its value at `s`, and every
denominator is a unit.  THM-3217 handles the first eight coordinates; hence

```text
r_0,r_1,P_1(r,a),s_0,s_1,P_1(s,r),t_0,t_1,rho_4,rho_5
```

are all `p`-adic units.  The common normalization change from `(2d-4)!` to
`(2p)!` has exponent `e_k` from `(14)` and is itself a unit under `(17)`.
Therefore every fixed offset has only the explicit finite exception set
through row five.

This remains independent of THM-3148's height-zero resultant `delta_s`.
The present theorem concerns selected reciprocal top-jet pivots, not residual
common roots.

## 6. A pivot can die and return one row later

The new walls are arithmetically realized already at offset `s=2`.  Normalize
the four consecutive coordinates

```text
(t_0,w_0,w_1,x_0)=(rho_3,rho_4,[z^1]R_4,rho_5).          (18)
```

Exact reduction gives

```text
p=41:       (t_0,w_0,w_1,x_0)=(6,15,2,0),
p=43:       (t_0,w_0,w_1,x_0)=(14,0,24,23).              (19)
```

At `p=41`, `W20(2)=0` while `W13(2)=27`; the fifth selected pivot dies.
At `p=43`, `W13(2)=0` while `W20(2)=36`; the fourth pivot dies but the fifth
returns as a unit.  Both primes miss every numerator factor in `Xi_2` and
all displayed denominators, so these are isolated new walls.

The return at `p=43` is structural and stronger than a constant-coefficient
calculation.  For arbitrary series `f,g`, if `f_0=0`, then

```text
E(f,g)=f_1 g_0 z^(-1)f.                                 (20)
```

Indeed `P_1(f,g)=-f_1g_0` on that wall, so the first two terms in `(1)`
reduce exactly to the right side of `(20)`.  In particular

```text
[z^0]E(f,g)=f_1^2g_0.                                   (21)
```

Taking `(f,g)=(w,t)` in `(20)` gives, with `c=w_1t_0`,

```text
x=E(w,t)=c z^(-1)w,
x_0=w_1^2t_0=24^2*14=23                       (mod 43). (22)
```

Thus the whole shifted row, not only its first coordinate, crosses the wall.
The clutch is nevertheless terminal one step later.  When `c` is a unit,
`w=(z/c)x`, and direct substitution into `(1)` gives

```text
R_6=E(x,w)=E(x,(z/c)x)=0.                                (23)
```

If `ord_z(f)>=2`, then `f_0=f_1=0` and `(20)` already gives `E(f,g)=0`.
Hence a simple anchored zero is the sharp boundary: death of one pivot is
not death of the reciprocal state, the neighbor jet carries the entire row
back once, and exact shifted proportionality kills the following row.

## 7. Consequence and boundary

The proved connection is

```text
source:       the universal degree-d resonant factorial reciprocal pair;
operation:    two further fraction-free E iterates;
target:       primitive row-four/row-five divisors W13,W20;
preserved:    exact rational-function and fixed-offset p-adic chart data;
destroyed:    the rest of each row and any automatic Newton-slope choice;
sidecar:      neighboring coefficient w_1 on a rho_4 wall.                (24)
```

Equations `(7),(8)` extend the atlas, while `(20)--(23)` give the exact
one-row clutch and terminal boundary across one genuine wall.  They do
**not** prove that row five separates
every offset, that the primitive-wall pattern continues indefinitely, or
that one selected pivot controls the full exterior state.  No arbitrary
radial-coefficient `SFC`, new `NC(2)`, new `GMC(2)`, or `LRC(14)` consequence
is asserted.

The next all-depth question is now precise: identify the primitive factor of
`rho_k` after removing inherited content, and prove that its integer-offset
specialization is nonzero.  The Pell gauge is known; primitive-wall selection
is the remaining arithmetic holotopy.

## 8. Exact evidence

Run

```text
python 04-computation/factorial_fourth_fifth_prs_primitive_walls_thm3223.py
python -O 04-computation/factorial_fourth_fifth_prs_primitive_walls_thm3223.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins promoted THM-3217; constructs eleven reciprocal top coefficients and
five `E` rows in `Q(d)`; proves `(7),(8)` exactly; checks 18 raw-moment
coordinates by an independent expansion; verifies the two Rabin
irreducibility certificates and parity identities; checks the Pell and gauge
recurrences; checks 24 exact whole-row clutch, delayed-collapse, and
double-zero controls; and reproduces every residue in `(19)--(23)`.  It uses no
floating point, interpolation, randomness, or optimization-sensitive
`assert`.

QED.
