# Exact inverse and sparsity for a one-site unit repair

In the one-apex transitive-spine family fix an arbitrary binary prefix U
and suffix V and flip the intervening single mark. Write u=bin(U) and
v=revcomp(V); either integer may be zero. The marked seam identity in
[the graceful session](../chanpak_graceful_ternary_20261115/GRACEFUL.md)
gives, exactly,

    H(U0V) = 1+2uv+p(u)+p(v)+u+(u mod 2)
    H(U1V) = 1+2uv+p(u)+p(v)+v+(v mod 2).

Indeed use the seam identity with first-address integer 2u or 2u+1
and second address v; the odd case has seam parity cost v mod 2.
This is a compact exact contrasting analogue of the earlier two-site
gap selector: here the gap constrains the two source addresses
against one another instead of fixing a single absolute address.

## Exact one-square-root target inverse per parity

A prescribed positive gap 2g forces ceil(v/2)−ceil(u/2)=g after
interchanging sides if necessary. Put

    u=2a−epsilon, v=2(a+g)−eta, epsilon,eta in {0,1}.

For a specified lower endpoint the quadratic baseline is

    B(a)=1+2(2a−epsilon)(2(a+g)−eta)+2a.

For the higher endpoint add 2g. The residual p(u)+p(v) is nonnegative
and at most 4a+2g, by the exact binary recurrence. Consecutive
integer baselines differ by 16a+8g−4(epsilon+eta)+10, strictly
more than this entire uncertainty. Thus for each parity pair the
greatest a with B(a)≤N is the unique possible candidate. Take its
integer square root from the quadratic formula, check a≥0 and
nonnegative addresses, and test the exact checksum. This is a
complete target-first inverse with at most four parity candidates,
at every prescribed gap. Both words are physically reconstructible;
leading address zeros are neutral pads, not new choices.

## Architecture capacity

Across **all** spines, all positions and all contexts, only
O(G(sqrt(X)+1)) distinct endpoint values ≤X can be incident to a
nonzero one-site gap of magnitude ≤2G. Indeed symmetry lets us
take 1≤g≤G, and for all a≥1 we have u≥a, v≥a+g, so the inequality
2uv≤X forces a≤sqrt(X/2); there are only four parity pairs.
The a=0 boundary has only four possibilities per g as well.
For each fixed positive g the count is in fact Θ(sqrt(X)), because
one parity lane has quadratic growth and disjoint uncertainty bands.
In particular **unit one-site repairs occur on a zero-density
set of absolute counts**, no matter where the marked site is placed
or how the spine is grown. More generally, allowing every gap up
to o(sqrt(X)) still gives zero density in this lane. This applies
to the one-site member of the 721–727 [interleaved escape](NONLINEAR_ESCAPE.md),
not to unit switches that change several marks or to arbitrary tournament
constructions.
