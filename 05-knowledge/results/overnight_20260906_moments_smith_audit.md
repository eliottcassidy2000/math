# Independent audit of the first odd-prime mixed two-jet cluster

**Verdict: PASS, analytic proof audit.** Moments lane, 2026-09-06.
Source: [the wildcard mixed-cluster report](overnight_20260906_smith_mixed_cluster.md).
Scope: arbitrary lifts and duplicate depth, every determinantal rank,
the competing tariff boundary, and the consecutive-node corollary.
This report claims no additional computational census.

The argument applies to an odd prime `p`, geometric scale `e>=1`, and
`p+1` lifted nodes with every residue present and exactly one duplicate
pair of difference valuation `d>=1`. No canonical residue representatives
are needed in the proof.

The potentially delicate row division is valid integrally: for a
derivative row `D(u)=(q u^(q-1))_q`,

```text
(D(u+delta)-D(u))/delta = (q(q-1)u^(q-2))_q mod p.
```

The unit part of `delta` introduces no additional loss. Modulo `p` this
is the ordinary second-derivative functional at the repeated residue.
The polynomial `X^(p+1)/(p+1)-X^2/2` is integral over `Z_(p)` and separates
this divided row from the old derivative-value bank: its first derivative
is `X^p-X` modulo `p`, while its second derivative is `-1`.

For the consecutive column set below degree `h`, the representative
derivative bank costs precisely one power of `p` when `p<h<2p` and none
at `h=2p`. The all-node bank costs exactly `d+1` for `p+2<=h<2p` and
`d` at `h=2p`. The displayed factorial-Vandermonde minors prove exact
attainment, and the modular ranks prove the lower bounds. Saturating the
bank cannot discard its old reduction: every old row is an integral
combination of saturated rows. Because the complete value/first-jet
observer has rank `h` modulo `p`, the saturated bank can be extended by
the required number of actual value rows. Undoing its row transformation
then gives an original minor of exactly the claimed valuation.

I checked all alternatives in the weighted middle-rank lower bound.
With `L=h(h-1)/2-(p+1)`, the possible costs are bounded below by

```text
eL+d+1                 all derivatives, first columns;
eL+e+d                 all derivatives, later columns;
eL+e+1                 p derivatives, first columns;
eL+2e                  every remaining possibility.
```

Each is at least `eL+1+min(e,d)`. Both competing minima have explicit
witnesses. This includes `d=e`; the equality is a tie of two lawful
certificates, with no continuity or genericity argument required.

At `h=2p+1`, deleting one row leaves at least one complete duplicated
pair of the same Hasse order. Its difference is divisible by `delta`,
independently of the columns. The minimum geometric scale is
`e(2p^2-1)`. The proposed witness contains values and first derivatives
at every residue and, after division, a second derivative at the repeated
residue. Since `2` is a unit, its kernel would have multiplicity three
there and two elsewhere, requiring degree at least `2p+1`. Thus its
determinant is a unit after exactly one `delta` division. This proves the
penultimate determinantal divisor without a hidden residue-lift condition.

The full determinant has valuation `2ep(p+1)+4d`, agreeing with the
pairwise-difference product. Consecutive differences give the displayed
partition. For `p=3`, the complete list is exactly

```text
(0,0,e,2e+1,3e+min(e,d),5e-1,
 6e+d-min(e,d),7e+3d).
```

The intervening exponent range is empty, while the middle determinantal
range still has rank five. The entries remain sorted. At `e=1`, the two
inherited entries `4,4` coincide; at `d=e`, the old `2p`-entry partition
is preserved exactly. If `d<e`, only the old entry at index `p+2`
decreases, by `e-d`, before the two terminal entries are appended.

The `p=2` exclusion is necessary at both uses of second derivatives:
the denominator `2` in the separating polynomial is not a unit, and an
ordinary second derivative does not supply a second Hasse condition.
The actual `(0,2,4)` partition consequently differs at rank five, as
the source correctly records.

Finally, for consecutive nodes through `n=p(p+1)`, the global split into
residue classes uses unit resultants and is a valid integral CRT. Each
class has at most `p+1` nodes. A class of size `p+1` is, after translation,
`p*(0,...,p)`, so it has exactly one internal duplicate of depth one;
the mixed theorem applies with `e=d=1`. The stated whole next-band
extension therefore follows for every odd prime. It does not justify
splitting a single residue class into an unchanged old Smith prefix.

No mathematical repair was found. The arbitrary-lift quantifier, exact
minor witnesses, boundary tariff, and retained CRT hypotheses all match
the claimed conclusions.
