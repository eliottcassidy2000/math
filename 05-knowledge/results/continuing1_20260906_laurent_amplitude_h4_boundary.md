# The simplest amplitude-interpolation recipe stops at the next actual support

**Status: FINITE-EXACT stopping boundary.** This does not change the frozen
[positive cubic amplitude repair](continuing1_20260906_laurent_real_power_boundary.md).
It refutes only the proposed coefficient positivity of the unique minimal-
degree amplitude interpolant. It does not exclude higher-degree positive
weights or combinations of other windows.

Take the next genuine A2B3 support **(-27,1,15)**, with x=1,h=4,g=14.
All three exponents are one modulo14, and the count triple(1,12,1)
witnesses its first support return. Keep the complete actual first/doubled
rows and alpha-completed hit response W_raw from the preceding source map.
Set T(s)=sQ_raw(-s), J_0(s)=sW_raw(-s). The monic first equation is

    p(s)=s^4-60s^3+84s^2-10s+1/11.

Let r_1<...<r_4 be its positive roots and z_i=sqrt(r_i)>0. There is a
unique polynomial A(z) of degree at most three satisfying

    A(z_i)=T(r_i)/J_0(r_i).

The exact certificate proves

    -344095/1000000 <= [z^2]A <= -344094/1000000 <0.       (1)

Thus the straightforward extension of the cubic coefficient-positive
amplitude interpolation already fails at this next actual support. All
four actual doubled responses and hit responses remain negative; (1) is
a failure of this particular positive-coefficient representation recipe,
not of actual noncancellation or of rootwise positivity of their ratio.

The proof uses ordinary Lagrange interpolation. The quadratic coefficient is

    sum_i [T(r_i)/J_0(r_i)]
       *[-sum_(j!=i) z_j]/[product_(j!=i)(z_i-z_j)].       (2)

The following disjoint rational intervals have opposite endpoint signs
of p and exhaust its degree. Integer-square-root enclosures give bounds
for every z_i; all divisions in (2) stay away from zero.

| Root | Left endpoint | Right endpoint |
|---|---:|---:|
| r_1 | 8419/849544 | 11993/1210189 |
| r_2 | 259526/2155711 | 291249/2419213 |
| r_3 | 11376199/8744207 | 1124341/864214 |
| r_4 | 120947883/2065060 | 91906768/1569213 |

The [standalone source](../../04-computation/continuing1_20260906_laurent_amplitude_h4_boundary.py)
constructs J_0 by full beta-path convolution and the actual alpha binomial
factor, including the negative beta index. It independently enumerates the
literal charge fibres at masses14 and28 to verify the first and doubled
rows. Rational interval arithmetic in (2) then proves (1).

```
python -B 04-computation/continuing1_20260906_laurent_amplitude_h4_boundary.py
python -B -O 04-computation/continuing1_20260906_laurent_amplitude_h4_boundary.py
```

The two runs pass41 always-active gates with byte-identical LF
[output](continuing1_20260906_laurent_amplitude_h4_boundary.out).
No numerical root or interpolation solver is used in the certificate.

    source 67038e1b3bbe6870e0e4d46998de47823a47bb0fca9ef06b468aaec3ee5358d3
    output 5587d5476994f1dee7b7448e7755ec3604fd426aff93ebae38cada3987cfc545

An exploratory h=4,...,8 calculation motivated this test; no signs beyond
this exact h=4 witness are promoted. The remaining operational question is
which additional phase weights or mixed-carrier generators can retain a
provable sign while crossing the appropriate actual-response cone wall.
