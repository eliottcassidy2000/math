---
id: THM-3217
title: "Universal resonant-degree PRS wall atlas and fixed-offset exception set"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the degree-d resonant factorial pair, the first three fraction-free
  reciprocal remainders have eight explicit pivot, neighbor, and connection
  coordinates.  Five universal wall polynomials H,J,K,U,V
  specialize at d=p+6 to the maintained offset-six H,J,K,H8,H9 walls.  For
  every fixed offset s>=2, a nonzero explicit integer Xi_s contains every
  prime at which any displayed coordinate can fail to be a p-adic unit after
  d=p+s.  This is a three-row chart theorem, not arbitrary-depth separation.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return
  - THM-3214-two-jet-pseudo-division-locality-and-catalan-sharpness
related:
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
script: 04-computation/factorial_universal_resonant_prs_wall_atlas_thm3217.py
output: 05-knowledge/results/factorial_universal_resonant_prs_wall_atlas_thm3217.out
script_sha256: ad83e262b691ae7d52b978f0251e5e56ccacf25ca4e6b304b28fffbdda5eb9e0
output_sha256: 7189a25d1d94ef07c2841d60f57588db600b7cafc3efaf9a719288495e594d41
hash_basis: LF-normalized bytes
---

# THM-3217 -- universal resonant-degree PRS wall atlas and fixed-offset exception set

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The offset-six scalars `H,J,K,H8,H9` in THM-3176 and THM-3192 look at first
like artifacts of one large Euclidean calculation.  They are not.  They are
the restriction to `d=p+6` of universal divisors on the degree line.  This
theorem computes that degree-parametric atlas and then pulls it back along
every fixed-offset section `d=p+s`.

The result is deliberately local in PRS depth.  It controls the first three
fraction-free reciprocal remainders, exactly the seven-jet budget visible to
THM-3214.  It does not assert that three rows separate every offset.

## 1. The degree-parametric reciprocal pair

Let

```text
L(t^k)=k!,

M_n^[d](v)=L((d-t+v t^2)^n).                              (1)
```

For `d>=9`, take the resonant pair

```text
A_d=M_(d-2)^[d],                    B_d=M_(d-1)^[d].       (2)
```

Reverse their coefficients from the top and use the same normalization:

```text
a_j=[v^(d-2-j)]A_d/(2d-4)!,
b_j=[v^(d-1-j)]B_d/(2d-4)!,

a(z)=sum_(j>=0) a_j z^j,            b(z)=sum_(j>=0) b_j z^j. (3)
```

Only the first eight entries are used.  Direct multinomial expansion gives,
for `epsilon in {0,1}` and `0<=j<=7`,

```text
[v^(d-2+epsilon-j)]M_(d-2+epsilon)^[d]
 = binom(d-2+epsilon,j)
   sum_(ell=0)^j binom(j,ell)d^(j-ell)(-1)^ell
                    (2d-4+2epsilon-2j+ell)!.              (4)
```

Thus every input jet in `(3)` is an explicit rational function of `d`.

Use THM-3214's fraction-free two-jet operator

```text
P_1(f,g)=f_0g_1-f_1g_0,

E(f,g)=z^(-2){f_0^2g-[f_0g_0+P_1(f,g)z]f}.               (5)
```

Define

```text
r=E(a,b),                 s=E(r,a),                 t=E(s,r). (6)
```

No pivot is inverted in `(5),(6)`, so every identity below remains an
identity on its own wall.

## 2. Eight universal coordinates

Put

```text
H(d)=24d-35,
G(d)=24d^2-259d+315,
I(d)=28d-45,                                                (7)

J(d)=256d^4-33792d^3+187360d^2-348768d+218295,

J_+(d)=256d^5-37120d^4+650464d^3-2875072d^2
       +4835439d-2837835,                                  (8)

K(d)=5120d^5-963840d^4+6841088d^3-18014016d^2
     +20587392d-8513505.                                   (9)
```

The two higher wall polynomials are

```text
U(d)=327680d^8+75694080d^7-5035040768d^6+61753251840d^5
     -341070457600d^4+1023103349760d^3-1737933934176d^2
     +1581643823760d-602628451425,                         (10)

V(d)=327680d^9+70123520d^8-6372368384d^7+152519243776d^6
     -1447233207040d^5+7075663957760d^4-19705993431648d^3
     +31775324380272d^2-27782756946945d+10244683674225.    (11)
```

The first remainder has

```text
r_0= -(d-1)H
      /[2(2d-7)(2d-5)^2],

r_1= (d-1)G
      /[4(2d-9)(2d-7)(2d-5)^2],

P_1(r,a)=d(d-1)I
      /[(2d-9)(2d-7)(2d-5)^3].                            (12)
```

The second has

```text
s_0=(d-1)^2 J
     /[16(2d-11)(2d-9)^2(2d-7)^3(2d-5)^4],

s_1=-(d-1)^2 J_+
     /[32(2d-13)(2d-11)(2d-9)^2(2d-7)^3(2d-5)^4],

P_1(s,r)=d(d-1)^3 K
     /[16(2d-13)(2d-11)(2d-9)^3(2d-7)^4(2d-5)^6].         (13)
```

Finally the first two coefficients of the third remainder are

```text
t_0=-(d-1)^5 H^2 U
     /[2048(2d-15)(2d-13)^2(2d-11)^3(2d-9)^4
             (2d-7)^7(2d-5)^10],

t_1=(d-1)^5 H^2 V
     /[4096(2d-17)(2d-15)(2d-13)^2(2d-11)^3(2d-9)^4
             (2d-7)^7(2d-5)^10].                          (14)
```

### Proof of `(12)--(14)`

Insert `(4)` into `(5)`.  One application consumes input jets through order
two, the next through order four, and the displayed third-row one-jet through
order seven, exactly as THM-3214 predicts.  Arithmetic in `Q(d)` gives
`(12)--(14)` after cancellation.  No specialization in `d`, interpolation,
floating point, or assumed nonzero pivot is used.  The exact companion carries
out this calculation in the rational-function field `Q(d)`.

## 3. Offset six is a literal section

Substitution `d=p+6` gives

```text
H(p+6)=24p+109,

J(p+6)=256p^4-27648p^3-365600p^2-1528800p-2096649,

K(p+6)=5120p^5-810240p^4-14447872p^3-92004672p^2
       -256323456p-265142241.                              (15)
```

Likewise `U(p+6)` and `V(p+6)` are exactly THM-3176's `H_8(p)` and
`H_9(p)`, coefficient for coefficient.  Therefore the maintained offset-six
walls are not merely analogous to `(7)--(11)`:

```text
(H,J,K,H_8,H_9)_(offset six)
       =(H,J,K,U,V)|_(d=p+6).                              (16)
```

There is one normalization change.  For a general fixed offset `s`, put
`d=p+s`.  Relative to a common `(2p)!` normalization, `(3)` has been divided
by the common factor

```text
u_(p,s)=(2d-4)!/(2p)!=prod_(i=1)^(2s-4)(2p+i).             (17)
```

The scaling law `E(xf,yg)=x^2yE(f,g)` gives common scaling exponents

```text
a,b,r,s,t:             1,1,3,7,17,
P_1(r,a),P_1(s,r):             4,10.                        (18)
```

If `p>2s-4`, then `u_(p,s)` is a `p`-adic unit.  Hence `(12)--(14)` detect
the same unit and zero charts as the common `(2p)!` normalization.  At
`s=6`, this recovers THM-3192's `H,J,K` chart ideals and THM-3176's two
third-row endpoint sidecars.

## 4. A finite exceptional-prime set for every fixed offset

Let

```text
D(s)=prod_(c in {5,7,9,11,13,15,17})(2s-c),

Xi_s=2s(s-1)D(s)H(s)G(s)I(s)J(s)J_+(s)K(s)U(s)V(s).       (19)
```

Then

```text
Xi_s !=0                         for every integer s>=2.    (20)
```

Indeed, the factors in `D(s)` are odd while `2s` is even.  The only roots of
the linear factors `H,I` are `35/24,45/28`.  The remaining six polynomials
have no residue root in the following finite fields:

```text
G mod 11,       J mod 13,       J_+ mod 37,
K mod 23,       U mod 17,       V mod 19.                  (21)
```

Their leading coefficients are nonzero in the displayed fields.  An integer
root would reduce to a root in every one of them, proving `(20)`.

Now fix `s>=2` and let `p` be prime with

```text
p>max(7,2s-4),                         p does not divide Xi_s. (22)
```

Because `W(p+s)=W(s) (mod p)` for every integral polynomial `W`, all
numerators and denominators in `(12)--(14)` are `p`-adic units.  Equation
`(17)` is also a unit.  Thus all eight displayed coordinates

```text
r_0,r_1,P_1(r,a),s_0,s_1,P_1(s,r),t_0,t_1                (23)
```

are `p`-adic units.  Consequently, for each fixed offset, every failure of
this displayed eight-coordinate three-row chart lies among the finitely many
prime divisors of the explicit nonzero integer `Xi_s` (together with the
bounded primes excluded by `(22)`).

This exception integer is independent of THM-3148's `delta_s`.  The latter
classifies common roots on the height-zero Frobenius residual; `Xi_s`
classifies selected top-jet PRS charts.  Neither condition implies the other.

## 5. The wall hierarchy is real

At `s=6`, five distinct primes isolate the five primary universal walls:

```text
p=109          hits H only among H,J,K,U,V,
p=232961       hits J only,
p=2767         hits K only,
p=1067703961   hits U only,
p=52511        hits V only.                                (24)
```

These are exact congruence checks at `d=p+6`; every listed number is prime.
Thus no one of the five primary wall factors can simply be deleted from the
atlas.  A wall in `(24)` kills its selected coordinate, not necessarily the
whole reciprocal exterior state or the final Newton separation.  THM-3192's
chart-cancellation hostile remains the governing boundary.

Geometrically, `(7)--(14)` form an arithmetic divisor arrangement over the
degree line.  A fixed-offset family is a section `d=p+s`; `(19)` is the
section's finite intersection divisor.  This is the useful holotopy here:
the complicated prime-by-prime formulas are shadows of one degree-parametric
object, while different charts can still exchange responsibility on a wall.

## 6. Consequence and boundary

For every fixed `s`, the first three fraction-free reciprocal rows are
uniformly regular outside a computable finite set of primes.  This extends
the structural meaning of the offset-six `H,J,K,H8,H9` calculation to all
fixed offsets without assuming the still-open unbounded resultant statement
`delta_s!=0`.

What it does **not** prove is equally important:

- three rows need not separate the Newton slopes for arbitrary `s`;
- if `p|Xi_s`, another Pluecker chart may survive and the theorem makes no
  rank-death claim;
- if `s` grows with `p`, `Xi_s` is no longer a fixed exception integer;
- no arbitrary radial-coefficient `SFC`, new `NC(2)`, or new `GMC(2)` claim
  follows from these selected top jets.

The next arithmetic target is an atlas-selection theorem across the full
`floor(s/2)` jet budget, not another offset-six expansion.

## 7. Exact evidence

Run

```text
python 04-computation/factorial_universal_resonant_prs_wall_atlas_thm3217.py
python -O 04-computation/factorial_universal_resonant_prs_wall_atlas_thm3217.py
```

and compare LF-normalized bytes with the declared output.  The companion
pins both proved dependencies; works in the exact rational-function field
`Q(d)`; reconstructs eight reciprocal PRS coordinates from eight input jets;
checks all formulas `(12)--(14)`; independently reconstructs 128 displayed
coordinates directly from raw moment coefficients at degrees `9` through
`24`; verifies the five offset-six
specializations; proves the modular no-integer-root certificates in `(21)`;
checks `Xi_s` on the independent control range `2<=s<=64`; verifies the
normalization exponents; and confirms all five isolated wall primes in `(24)`.
It contains no floating point and no optimization-sensitive `assert`.

An independent hostile audit rederived the top multinomial formula, the
orientation of `E`, all eight coordinates at five fresh degrees, the literal
`d=p+6` specializations, all six modular no-root certificates, and the five
isolating primes.  It also caught, before promotion, a stale dependency hash
and the omitted `P_1(r,a)` scaling exponent; the immutable candidate repairs
both.  Fresh normal and `-O` runs byte-match the stored transcript and the
declared LF-normalized hashes.

QED.
