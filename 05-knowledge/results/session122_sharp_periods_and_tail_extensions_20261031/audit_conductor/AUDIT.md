# Independent audit: the optimal interval-only clock constant

**PASS, analytic theorem and independent finite rational controls.**
This audit imports no producer module. It does not assert unrestricted LRC,
optimality for a specified core's entire safe set, or literature priority.

## Statement audited

Let k>=2, 0<r<=1/(4k-2), and let I be a closed interval of width W>0.
For any positive real speeds v_1,...,v_k and any real phase shifts theta_i,
if every v_i W>=2kr, there exists t in I such that

    ||v_i t+theta_i||>=r  for all i.

For positive integer speeds and zero shifts this gives the advertised
cutoff ceil(2kr/W). If I lies inside a core's safe set, the witness retains
all the original core constraints. The conclusion is nonemptiness, not
positive measure: boundary examples below have only isolated safe points.

## Independent proof

Put a=2r. On an interval of length x in the phase variable, a period-one bad
set consisting of an open arc of length a has measure at most

    a floor(x)+min(x-floor(x),a).

Indeed a block of each whole period contributes exactly a. The residual
arc has length d=x-floor(x), and its intersection with the bad arc has
measure at most min(d,a). The bound is attained by a suitable starting
phase, so it is an exact phase-independent envelope.

Divide by x and write M(x) for this fraction. For 0<x<1, M(x)=min(1,a/x).
On each interval n<=x<=n+1 with n>=1, the fraction increases until x=n+a
and decreases thereafter. Its maximum there is a(n+1)/(n+a). These
maxima decrease with n since a<1. Thus for x>=ka,

    M(x)<=max(1/k,2a/(1+a))=1/k,

where the equality uses (2k-1)a<=1, exactly r<=1/(4k-2).
The first term bounds the initial n=0 interval. All subsequent intervals
are bounded by the n=1 peak. This also covers every endpoint x=n and x=n+a.

Consequently each bad set B_i inside I has measure at most W/k<W.
Their sum of measures is at most W. Suppose that these relatively open
sets covered I. Each is proper, because its measure is less than W.
Connectedness of I implies that two members of this finite open cover
intersect: otherwise its nonempty members partition I into disjoint open
sets, a separation. Their intersection is a nonempty relatively open
subset of a nondegenerate interval and therefore has positive measure,
even if it contains an endpoint. Hence the sum of the bad-set measures
would strictly exceed the measure W of their union. Contradiction.
The retained witness lies in the complement of every open bad set,
which is why equality at a safe endpoint is permitted.

For k=1 the separate elementary bound is vW>=2r: no open bad interval of
length 2r/v can contain the whole closed I. No k>=2 connected-cover step
is needed in that case.

## General envelope conductor through the full strict-density range

**PASS.** For k>=2 and 0<r<1/(2k), put a=2r and

    T=a(k-1)/(1-ka),  N=ceil(T),  C=kaN.

Then v_i W>=C for every speed still guarantees a common safe time,
with arbitrary real phase shifts. This does not claim C is optimal
outside the preceding low-radius range.

The peak at x=n+a has value a(n+1)/(n+a), which is at most 1/k
exactly when n>=T. The last possible offending peak is n=N-1.
Since N-1<T<=N and ka<1,

    N-1+a < kaN=C < N.

Thus C lies on the decreasing branch after that peak, where the envelope
is aN/x and equals 1/k at x=C. It stays at most 1/k until x=N, and
every subsequent peak has index n>=N>=T. This proves the bound on all
x>=C without skipping an intervening peak. The same relative-open-cover
argument handles equality. If T is itself integral, the next peak merely
attains 1/k, which is harmless.

At radius 1/14, this gives C=10/7 for five tails and C=30/7 for six tails.
The strict density assumption is retained. No claim for k=7 at that
radius follows from the formula, whose denominator would be zero.

## Sharpness of the coefficient, with distinct unshifted integer speeds

The constant 2kr cannot be replaced by any c<2kr in a theorem about an
arbitrary supplied interval and the minimum attached speed alone.
This statement is distinct from any optimal-radius or optimal-core theorem.

Choose rational s with max(c,0)<s<2kr. Set

    d_j=(2j+1)s/(2k),  j=0,...,k-1,
    m=r-s/(2k)>0.

Since 2kr<=k/(2k-1)<1, all d_j are strictly between zero and one.
Choose an integer D with every D d_j integral. Put a_j=D(1-d_j),
choose an integer N>s/m, and define

    U=ND,  v_j=U+a_j,  I=[1/D, 1/D+s/U].

These v_j are distinct positive integers. Every t in I is t=1/D+x/U
for some x in [0,s]. The nearest d_j is at distance at most s/(2k),
while modulo one

    v_j t = x-d_j+(a_j/U)x.

The error term has magnitude at most s/N<m because 0<a_j<D.
Therefore for that j, ||v_j t||<r. Every t in I is bad for at least
one speed, so there is no common safe time. Yet

    min_j(v_j) width(I)>s>c.

Thus any smaller universal coefficient fails, even with distinct integer
speeds and zero phase shifts. The construction does not claim I is the
full safe set of a prescribed nonempty core. The parent's equivalent
rational-speed construction followed by common-denominator rescaling is
also sound, with the same strict margin condition.

## Comparisons and source check

The immediate predecessor `session121_near_frontier_closures_20261031/lrc/LRC.md`
T1 used the sufficient cutoff ceil(1/W), at the same radius range.
`precision_lifts_ramsey_deficits_20261030/LRC.md` L2 and
`clique_palette_joint_ports_20261030/LRC.md` L4 used that cutoff for pairs.
Those source statements were inspected at commit
`d345f993475dd7b5272ac6c3b4ab693899515a38`.
The independently proposed partial-period constant
2kr/(1-2r(k-1)) is also valid but larger than 2kr. It is not needed here.
A limited source check is not an exhaustive priority search.

At r=1/14 the new constants for k=1,2,3,4 are respectively
1/7,2/7,3/7,4/7. The k>=2 theorem applies only through k=4 at this radius.

## Reproduction and controls

    python3 -B check.py > normal.json
    python3 -B -O check.py > optimized.json
    cmp normal.json optimized.json

The script uses rational arithmetic and explicit exceptions, without producer
imports or removable assertions. It checks the envelope at all its defining
breakpoints across k=2,...,12 and three radii per k; literal bad-set measures
against the envelope; actual full-comb intersections at the exact conductor
cutoff with independent phase shifts; point-only boundary examples; the full strict-density-range envelope at 56 additional parameter pairs;
actual shifted intersections at the corresponding general cutoff; and
seven explicit unshifted distinct-integer-speed sharpness counterexamples.
The two hostile scalar controls show why the radius condition and clock
cutoff cannot simply be erased from this proof. Normal and optimized
receipts are byte-identical.
