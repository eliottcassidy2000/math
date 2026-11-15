# Many unit shifts, but a small real exposure at increasing width

The [complete width-by-shift recurrence](ALL_WIDTH_TRANSPORTS.md)
counts all same-core global transports at width 12k. Here are
two consequences, one positive and one an obstruction to using
this grammar as the full target selector.

## Exact local multiplicity of each fixed correction at large width

Normalize a transport by assigning its free ternary units independently
and uniformly to the two tapes. Let D be its signed loss divided
by two. Its distribution is the bounded independent sum from the
displacement polynomial: 2k−1 bits at step −16, 2k−1 at step +25,
k−1 at step −34, and a deterministic constant 8(k−2). The mean
is −7/2 and the variance is (2918k−2037)/4. In particular for every
fixed even correction delta (including +2), the number of exact
matched suffix pairs at width 12k with loss delta is
```
(1+o(1)) * [4*32^(k−1)] / sqrt(2*pi * (1459/2) k).
```
This counts actual distinct address charts at one common width,
not independent representations of the same absolute target.

For a direct proof, Fourier invert the integer-valued D. After
centering, its characteristic function is the product of cosines
with arguments 16t/2, 25t/2, 34t/2, raised to the three numbers of
bits above. Because gcd(16,25,34)=1, their absolute values can
all be 1 simultaneously in [-pi,pi] only at zero. Away from a
fixed neighborhood of zero the product therefore decays
geometrically; inside that neighborhood it is bounded by
exp(−c k t²). For |t|<=k^(-2/5), the log cosine expansion gives
−Var(D)t²/2 with uniform error O(k t^4)=o(1). The phase for a
fixed requested half-step is asymptotically 1 after scaling
t by 1/sqrt(k). Integrating the Gaussian proves the displayed
formula. One could also read exact coefficients from the
polynomial without appealing to this limit.

## Allowing every matched chart still has an exponential width tariff

Take the union over all global matched suffixes at widths at least
12K, even allowing every odd integer between the two endpoint counts
as if it were already a witness. For a chart at width 12k, let
A=4096^k and y=A z+V. Every V is in the upper half of its block.
Thus the two-address formula gives
`F(x,y)>=c A(x+1)(z+1)` with one positive absolute c. Its
companion differs by O(k), by the exact displacement polynomial.
If their interval meets targets at most X, the upper endpoint
is at most X+O(k), so k=O(log X) and the pair of outer addresses
satisfies a divisor bound O(X/A). Each pair offers at most
O(k+1) odd integers in its interval.

There are exactly 4*32^(k−1) charts at width 12k, including
all signs. Therefore the union counted with multiplicity is bounded by
```
O( X log X * sum_(k>=K) (k+1) 128^(−k) ).
```
In particular the high-width bank above
K=(1+epsilon) log_128 log X has density zero through X.
This statement already allows overlapping suffix choices,
nonunit returns, both signs, and even unproved interpolation
between endpoints. Any finite composition of forward global
moves of this kind is itself a global move with summed width
and loss. Each intermediate departure is O(k) as well,
since the sum of the absolute step bounds is O(k). The
same estimate therefore also counts all such intermediate endpoints.

This is substantially more restrictive than the individual
width census: there are exponentially many charts for a
unit step, yet the set of address pairs admitting an expensive
repair shrinks even faster. Increasing the fixed-reader width
alone is not a substitute for choosing the scale or revising
the higher-order tournament interface.
