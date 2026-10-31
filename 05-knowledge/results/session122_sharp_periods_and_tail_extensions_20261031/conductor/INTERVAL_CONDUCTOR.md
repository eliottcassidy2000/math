# Exact interval discrepancy and an improved many-tail conductor

Status: PROVED, with companion independent analytic audit and exact controls. This is a local interval theorem, not unrestricted Lonely Runner.

Let I=[b,b+W], W>0, and suppose every core speed is safe throughout I at radius r. Attached speeds v_i are positive real numbers; integer speeds give the ceiling versions below. Safety means ||v_i t||>=r, with endpoints retained. Arbitrary fixed phases may also be added to each attached speed.

## 1. Exact one-speed envelope

Put a=2r, 0<a<1. For speed v and x=vW=n+d, n=floor(x), 0<=d<1, the bad set ||vt+phi||<r has measure at most

    [a n + min(d,a)]/v.                                    (1)

The n full periods contribute exactly an/v. An interval of residual length d on the unit circle meets a single circular bad arc of length a in at most min(d,a). Some phase attains this measure, so (1) is the exact worst-phase envelope. Open endpoints change no measure.

## 2. General conductor for every k with 2kr<1

For k>=2 and ka<1 define

    N = ceil(a(k-1)/(1-ka)),    c = kaN.                    (2)

If every v_i W>=c, all k attachments have a common safe time in I. In the integer-speed case it suffices that min(v_i)>=ceil(c/W). For k=1 the exact elementary constant is c=a.

Proof of the envelope gate: let M(x)=[a floor(x)+min(frac(x),a)]/x. On each cell [n,n+1], M increases up to n+a and decreases thereafter (the n=0 initial part is constant1). For n>=1 its peak is a(n+1)/(n+a), decreasing in n, and is <=1/k exactly when n>=a(k-1)/(1-ka). All peaks from N onward are therefore admissible. On the previous cell the decreasing tail reaches1/k exactly at x=kaN. The inequalities N-1+a<c<=N follow from N-1<a(k-1)/(1-ka)<=N and ka<1 (the right inequality is immediate). Thus M(x)<=1/k for all x>=c.

Each of the k open bad sets consequently has measure <=W/k. They cannot cover I. Indeed, if they covered I, each would be a proper relatively open subset, since its measure is less than W. Connectedness of I then forces two of these sets to overlap in a nonempty relatively open subset, which has positive measure. A finite covering would have sum of measures strictly greater than W, contradicting the bound. This proves existence even at equality; it does not assert positive measure of the successful set.

For k=1, a closed interval at least as long as the largest open bad gap a/v meets the closed safe set.

## 3. Small-radius optimal interval constant

If 0<r<=1/(4k-2), then N=1, and the constant simplifies to

    c=2kr.                                                (3)

The equality radius is allowed. Relative to session121's sufficient min speed >=ceil(1/W), the new thresholds at r=1/14 for k=1,...,6 are

    1/7, 2/7, 3/7, 4/7, 10/7, 30/7,

each divided by W before taking the ceiling. The fifth and sixth entries use (2), not (3).

The constant2kr is optimal as a universal interval-only bound, even for unshifted positive integer speeds. This is not a sharpness statement for a fixed core's entire safe set.

To prove optimality, fix c0<2kr and choose positive rational s with c0<s<2kr. On a unit interval choose z_i=(2i-1)/(2k), i=1,...,k, and phases phi_i congruent to -s z_i modulo1. Every x in[0,1] has a nearest z_i with |x-z_i|<=1/(2k), hence the corresponding shifted speed satisfies ||s x+phi_i||<=s/(2k)<r.

Choose a sufficiently large positive rational b such that 1/b<r-s/(2k). Select rational delta_i in[0,1) with s b+delta_i congruent to phi_i modulo1 and let v_i=s+delta_i/b. At t=b+x the phase differs from s(x-z_i) modulo1 by delta_i x/b, of absolute value <1/b. Thus at least one speed remains strictly bad at every point of[b,b+1], while min(v_i)>c0. Clear the common denominator D of all v_i: the speeds Dv_i are positive integers and the interval [(b)/D,(b+1)/D] has the same dimensionless products v_i W and the same complete failure. This disproves every smaller universal constant c0. Repeated speeds are allowed in the theorem; the construction need not assume distinctness.

## 4. What the reduction preserves and loses

The source is the complete closed core-safe interval; the target is the exact worst-phase measure envelope. The map retains width, speed, radius, and open/closed endpoint conventions but discards mutual phase correlations and the rest of the core's safe set. Connectedness restores equality-case existence. The sharpness construction concerns arbitrary intervals and does not realize them as the entire safe set of a prescribed core.

The c in(2) is the minimal cutoff for this particular individual-measure condition M(x)<=1/k on every x>=c. For N>1 it is not claimed to be the optimal simultaneous-attachment constant. No claim of literature priority is made. The earlier session121 conductor and the independently derived partial-period conductor remain valid, but weaker in the small-radius range.
