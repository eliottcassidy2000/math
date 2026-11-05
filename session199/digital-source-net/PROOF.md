# Canonical physical sources form a perturbed binary digital net

Session 199 continuation. **PROVED:** an exact bit-reversal product formula,
including its positive zero-run correction, and an all-order joint discrepancy
bound for actual one-apex sources. **OPEN:** exact two-source scalar coverage.
This is not a replacement for the Session198 forty-template lemma.

## 1. Inherited physical source and scope

Let g(1)=1, g(2h)=g(h), and g(2h-1)=g(h)+h-(h mod 2) for h>=2.
The reviewed Session198 two-index packet constructs, for every U>=2, the actual
source tau(U)=(g(U),U mod 2,U,1). It is the source obtained from the edge
(1,0,2,1) by reading the binary digits of U-1 after the leading one, with zero
mapped to Z and one mapped to O. The physical sidecar is therefore explicit.

The earlier independent-source family is exactly
F(J^i tau(U),J^j tau(V)), with i,j in {0,1} and the ten C2/C3 forms F.
It was already isolated in session198-source-pair-construction/two-index-coverage/
PROOF.md. Its open claim is coverage of all odd targets above 691. We do not
rename that existing family or infer its universal coverage from more samples.

A relevant older lesson is the exact occupancy/quantization distinction in
01-canon/theorems/THM-2048-fiber-quantization-tax-for-lrc-discrepancy.md.
The present result transports neither an LRC carrier nor a cyclic action. It
instead identifies a different exact carrier: the binary bit-reversal net,
with the physical source deviation retained as an explicit positive correction.
The older theorem is motivation, not a proof dependency. The companion
THM-729-density-second-moment-is-autocorrelation-riemann-discrepancy.md also
separates a distribution estimate from its stronger covering objective; its
empirical second-moment estimate is not used here.

## 2. Exact product and zero-run correction

Put n=U-1>=1, k=bit_length(n), B=2^k, and

    R=rev_k(n), C=B-1-R.

The reversal uses exactly k digits. A zero run is a maximal interval of zero
bits i,...,i+ell-1, counting positions from the least significant bit. Since
bit k-1 is one, every such run has a one immediately above it. Put
m_i=n mod 2^i (in particular m_0=0).

**Theorem 1.** The actual source Hamilton count satisfies

    g(U) = n C/B + delta(n),
    delta(n) = 1 + sum_zero_runs [1-(m_i/2^i)(1-2^(-ell))].       (1)

If r(n) is the number of zero runs, then

    1 <= delta(n) <= 1+r(n) <= 1+floor(k/2).                    (2)

The first equality delta=1 holds exactly when there are no zero bits.

**Proof.** Let f(n)=g(n+1). The inherited binary recurrence is

    f(0)=1, f(2x+1)=f(x), f(2x)=f(x)+x+(x mod 2).

Writing epsilon_i for bit i of n and stripping its least significant bits,

    f(n)=1+sum_(epsilon_i=0)
         [floor(n/2^(i+1))+epsilon_(i+1)],

where epsilon_k=0. The top bit is one and contributes no term. On the other
hand C/B=sum_(epsilon_i=0)2^(-i-1). Subtract nC/B. In each zero run, the
sum of the epsilon_(i+1) terms is exactly one, at its upper end. The sum of
the fractional parts n/2^(i+j+1), j=0,...,ell-1, is

    (m_i/2^i)(1-2^(-ell)).

This proves (1). Each bracket lies in (0,1], and every run uses at least one
zero and its immediately higher one, so 2r<=k. This proves (2), including
its equality boundary. QED.

The integer form suitable for exact verification is

    B g(U)-n C = B(1+r)
                 -sum_zero_runs m_i [B/2^i-B/2^(i+ell)].       (3)

No floating-point approximation is required.

## 3. A genuine joint physical discrepancy theorem

For fixed k>=1 let B=2^k and let n run from B/2 through B-1. These are B/2
actual sources, with U=n+1. Normalize each by

    x_n=n/B,              y_n=g(n+1)/n.

Use the following anchored counting function, for 1/2<=x<=1 and 0<=y<=1:

    A_k(x,y)=#{n: B/2<=n<B, x_n<x, y_n<y}.

**Theorem 2.** Uniformly over all such real x,y,

    |A_k(x,y)-B(x-1/2)y| <= 2k+4.                            (4)

Consequently, for every axis-aligned half-open rectangle
R=[a,b) x [c,d) in [1/2,1] x [0,1],

    |#{physical normalized sources in R}-B area(R)| <= 8k+16. (5)

In particular any such rectangle of area greater than (8k+16)/B contains
an actual source, with an O/Z construction word of length k-1. Both height
and the normalized Hamilton count are controlled jointly.

**Proof.** First replace y_n by z_n=C/B. On the full range 0<=n<B,
(n/B,(B-1-rev_k(n))/B) is a reflected binary bit-reversal net.
For a dyadic prefix interval in the first coordinate, fixing its h leading
bits leaves k-h free bits. The second coordinates form a uniformly spaced
grid of 2^(k-h) points, with one fixed shift. Its count below any y differs
from 2^(k-h)y by at most one.

The integer prefix [0,ceil(Bx)) is a disjoint union of at most k such dyadic
prefix intervals (the endpoint x=1 is a single interval). Thus its count
below y differs from Bxy by at most k+1, allowing one for the nonintegral
endpoint. Subtracting [0,B/2), which is itself one dyadic interval, proves
that the reflected-net anchored count differs from B(x-1/2)y by at most k+2.

By (1)-(2), throughout the physical shell,

    0 <= y_n-z_n <= (k+2)/B = epsilon.

Therefore the physical count is between the reflected-net counts at y and
max(0,y-epsilon). Their expected counts differ by at most
B(x-1/2)epsilon <= (k+2)/2. Combining this with the previous k+2 error gives
an absolute error at most 3(k+2)/2, and hence the simpler bound (4).
Four anchored counts give (5). Strict area excess ensures a positive integer
count. Every counted index is an actual source from Section1. QED.

## 4. Exact-target boundary

The new map is from a binary digit word to an actual source and a reflected
bit-reversal point. It preserves the source recipe, height, and Hamilton count
through the exact correction (1); omitting that correction destroys exact
scalar equality. The discrepancy statement is joint, not a marginal histogram.

For independent source choices, (4)-(5) give ordinary geometric control on
both copies of the physical source data. They do not count points on an exact
integer level of one of the forty bilinear forms. The required target strip
shrinks with the source height; its expected occupancy can be smaller than the
proved discrepancy error. This is the remaining arithmetic issue, not a hidden
physical-membership assumption. No exact coverage or new cofinite threshold is
asserted.
