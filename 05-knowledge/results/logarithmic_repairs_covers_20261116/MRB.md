# MRB tail scale and the cost of an absolute repair

For the [MRB constant](https://marvinrayburns.com/ThegeometryV11.pdf),
put a_n = exp(log(n)/n)−1. Its alternating tail starting at n
has absolute value S_n = a_n−a_(n+1)+a_(n+2)−... for n≥5.
The sequence is decreasing and convex there: t=log(x)/x is
decreasing and t''=(2log(x)−3)/x³>0; exponentiation preserves
these two properties. Algebra yields
2S_n−a_n = (a_n−a_(n+1))−(a_(n+1)−a_(n+2))+...
The differences are nonnegative and decreasing, hence
a_n ≤ 2S_n ≤ a_n+(a_n−a_(n+1)). In particular

$$
2n S_n = \log n+O((\log n)^2/n).
$$

That is a concrete way to get a logarithmic *absolute* debit
out of a shrinking analytic relative error. In contrast the
[triangle budget](../triangular_backbone_chanpak_20261116/CORE.md)
forces any positive induced repair from count H to cost at least
2 ceil(log_3 H), and the [factored-gap theorem](../triangular_backbone_chanpak_20261116/FACTORED_GAPS.md)
may force much more. Thus the literal budget 2H S_H could not pay
for positive induced growth at sufficiently large H: its leading
coefficient of log H is one, versus 2/log(3). This compares scales;
the alternating series does not produce a tournament move.

A stronger absolute selector test also uses the previous MRB
pressure theorem. Suppose the set of anchor counts A has
Σ_(a∈A) a^(−s)<∞ with 0<s<1. Allow **all** odd targets within
R(X) of an anchor for targets in [X,2X]. If R(X)=O(X^(1−s)),
then the repair covers only o(X) targets, even before testing
whether the repairs are physically available. For sufficiently
large X the radius is below X/2. Every eligible anchor is then
in [X/2,5X/2]; their number is o(X^s), by bounding it by
(5X/2)^s times the remaining series tail. Multiplying by
2R(X)+1 gives o(X). Thus the previous 18/23 and 37/43 pressure
certificates require radii growing faster than X^(5/23) and
X^(6/43), respectively, for full coverage in those banks.
In particular a logarithmic signed repair, even if perfect
locally, cannot cover the holes of those sparse anchor banks.
One needs more anchors, a larger movable radius, or a joint
construction that is not a bounded thickening of the old bank.
