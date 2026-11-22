# Mod-3 pressure excludes almost all prime heights with order-six parts

A shared-gap chain may use an arbitrary ordinal assembly of strong
components of order at most six in each gap, with arbitrary singleton
padding and arbitrary zero runs. Its number of distinct prime heights
through \(X\) is \(O(X^s)\) with \(s=65535/65536\); in particular it
misses a relative density-one set of primes. This extends the
[order-five obstruction](RESET-PRESSURE.md) and assumes the same
[physical shared-gap transfer](SHARED-CONTINUANTS.md). It is not
an obstruction for all tournaments or all larger pieces.

## Physical reductions and real growth

A nonempty gap of profile \((h,Q)\), ordinary and two-slot counts,
applied to a state of height \(F\) updates it to \((Q(F+U),0,hF)\).
Since \(Q\ge2h\), a positive boundary has \(r=G/(D+G)\in[0,1/3]\).
Both profile reads multiply under ordinal assembly, so any factor
\(\gcd(h,Q)>1\) persists after the gap and through the whole suffix.
For a prime output the gap would already have produced a height
strictly greater than that factor. The same argument excludes adjacent
positive gaps by the exact reset factorization. Terminal zeros can be
normalized as in [RESET-PRESSURE](RESET-PRESSURE.md): \(ZZ\) is scalar
equivalent to a transitive width-one gap, and a lone final \(Z\) following
a nonempty module doubles that module's \(Q\) via padding. On the initial state a leading \(ZZ\) equals the
transitive width-one positive state, while a leading \(Z\) followed
by any \((h,Q)\) gap equals the initial \((h,2Q)\) gap state. These
padding rules retain the same palette. Thus count a first positive gap followed by successive nonempty
zero runs and positive gaps.

Exhaustive subset/last-vertex DP over all labelled strong tournaments
through order six gives these coprime nontrivial profiles:
\[
(5,32),(11,100),(17,224),(19,256),(23,272),(23,296),(25,304),
(29,320),(31,328),(37,352),(37,376),(41,416),(43,424).
\]
The transitive lane starts with \((1,2)\). All singleton padding is
included by multiplying \(Q\) by \(2^k\), \(k\ge0\). In bounding
composite gaps it is safe to overcount ordered lists of nontrivial
pieces, even if the same resulting profile is counted repeatedly or
has inadmissible content.

For \(j\ge1\), let \(a,b\) be the sums and \(c,d\) the middle entries
of the first and third columns of \(Z^j\). For input ratio \(r\) put
\(f=a+(b-a)r\), \(u=c+(d-c)r\). The growth across the next zero run
and gap is \(g=h f+Q(f+u)\), and the new ratio is \(r'=h f/g\).

## Keep the projective state modulo three

For a prime target greater than three, the positive pair \((D,G)\)
never vanishes modulo three since that common factor persists.
Retain its projective type \(t\), with representatives
\((0,1),(1,0),(1,1),(1,2)\). The next type is the projectivization
of \((Q(F_*+U_*),hF_*)\) reduced modulo three, using the actual
integer columns of \(Z^j\). Reject a zero pair. Use positive weights
\[
 W_t(r)=(A_t+B_t r)/1000,\quad
 A=(1000,1362,1374,1227),\qquad B=(116,194,195,186);
\]
then \(1\le W_t(r)\le 1.44\).

The exact integer replay proves the following pointwise bounds for
the finite box \(1\le j\le48,\ 0\le k\le13\):
* The sum of \(g^{-1}W_{t'}(r')\) over surviving one-piece transitions
  is less than \(0.966 W_t(r)\) for every input type and ratio.
* Without the type filter or \(W\), the corresponding sum of
  \(g^{-1}\) over nontrivial one-piece transitions is less than \(0.158\).

To certify the continuum, use the 48 cells with endpoints \(i/144\).
On each cell lower bound \(g\) at its smaller endpoint value and
upper bound \(r'\) at its larger endpoint value. The latter is valid
because a linear-fractional function with positive denominator is
monotone or constant. Round each term upwards to a billionth by
integer division and compare against the lower endpoint value of
\(W_t\). The replay performs all 192 comparisons and exhausts all
labelled order-six profiles independently.

## Infinite tails and composite gaps

Because \(Z^3\ge2I\), \(g\ge Q_i 2^k 2^{\lfloor j/3\rfloor}\).
Thus the omitted unweighted sum is at most
\[
 \sum_i Q_i^{-1}\left(10/65536+5/8192\right).
\]
The base reciprocal sum is below \(0.077\) without the transitive lane
and \(0.577\) with it. In particular the weighted one-piece tail
costs less than \(0.001 W_t(r)\), and the full unweighted nontrivial
one-piece pressure is less than \(0.159\).

An additional nontrivial piece of profile \((h_i,Q_i)\) in a gap
multiplies \(g\) by at least \((h_i+2Q_i)/3\), since at least two
thirds of the current growth comes from the \(Q(f+u)\) term and
\(Q_i\ge h_i\). The sum \(\beta\) of reciprocal multipliers is less
than \(0.109\), again by exact rounding. Starting from the nontrivial
one-piece bound, all lists with at least two pieces cost at most
\[
 1.44\cdot0.159\cdot 0.109/(1-0.109)\ W_t(r).
\]
There is only one padding exponent for the whole composite gap; the
estimate assigns it to the first piece. The total weighted transfer
is therefore strictly below \(0.996 W_t(r)\) at exponent one.

## An explicit exponent below one

Take \(s=65535/65536\) and \(c=1001/1000\). The finite rectangle has
\(g<2^{80}\) by the replay, and \(c^{65536}>2^{80}\), so each term's
increase from exponent one is at most \(c\). The same multiplier
bounds the increase of \(Q_i^{-1}\) and of the reciprocal factors
in \(\beta\), since those bases also lie below \(2^{80}\).
Use \(r_0=c/2>2^{-s}\). For the infinite tail replace
\(10/65536+5/8192\) by the geometric bound
\[
 \frac{2r_0^{16}+3r_0^{17}/(1-r_0)}{1-r_0}
 +\left(2+\frac{3r_0}{1-r_0}\right)\frac{r_0^{14}}{1-r_0}.
\]
The first term covers \(j\ge49\) at all paddings, and the second
covers \(k\ge14\) at all zero lengths. With these replacements
the same exact rational comparisons give total pressure less than
\(0.996\) at this explicit exponent as well.

Count only profile presentations, not topologically different
placements that have the same profile. The first positive gap has
finite inverse-\(s\) weight: \(Q\) multiplies, the sum of inverse
\(s\)-powers of the finitely many nontrivial base \(Q_i\) is below
one, and singleton padding is geometric. Iterating the contraction
gives a finite sum of inverse \(s\)-powers over all surviving profile
words. It bounds distinct prime outputs through \(X\) by \(O(X^s)\).

To replay the catalogue and both rational pressure certificates,
run python order-six-check.py --profiles.
