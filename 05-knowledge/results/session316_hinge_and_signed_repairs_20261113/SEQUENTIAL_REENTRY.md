# Sequential reentry after exposing a fixed set of internal ports

The [neutral cyclic-five bank](NEUTRAL_BANK.md) supports more than a
one-time equality. Fix the total weight $k$ of a $d$-bit mirror word.
At some chosen set $J$ of left positions, retain the values of the
bits and the literal exterior probe or local enlargement used there.
In comparisons between words these augmented blocks, and any augmented
mirror blocks, must be identical at the marked positions. All other
sites still use the original cyclic-five chart. Then every word in
this marked fiber has the same full $Q$ and the same pointed apex
profile, even after exposing any number of such local positions.

Indeed for an unprobed individual site $i$, its only changed term
is the split around the common apex in that site. With current block
polynomials $Q_j(z)$, the exact change is
$$
 2z^2\prod_{j\ne i} Q_j(z).
$$
All unprobed sites have the same $Q_i=Q_R$, so their gains have the
same weight in the whole polynomial. In a fixed-weight comparison
the signed gains cancel. The three noninterior apex statuses still
agree one by one, and the interior status follows by subtracting
from the full polynomial, just as before any probe. Exterior extension
through the retained apex remains legal by the same ordered-run gluing
argument. A fixed local modification may be made at another chosen
position, then comparisons can continue within each fiber with that
new bit held fixed.

With $|J|=s$ and $r$ selected bits among the revealed positions,
the number of marked realizers in the fiber is
$\binom{d-s}{k-r}$, with the usual zero convention. Mark the probes
so deleting them recovers the unprobed mirror word, which was rigid and
intrinsically decoded from its apex and ordered components.
Thus at any nontrivial residual fiber neutral moves can be reused:
move one selected bit to one unselected site while leaving all
revealed positions fixed. The one-probe formula in the companion
note proves a genuine signed readout from the original fiber; this
extension proves persistent neutrality of the unrevealed fiber, not
that every later pair of probe outputs has the same magnitude or can
be forced to hit an arbitrary integer. Adding probes itself generally
changes the count; the equality is between completions of a fixed
revealed pattern.

## A read order that provably keeps a nonzero signed separation

Use the positive local vertex-0 probe in the companion note and read
left sites from right to left, never probing the right half. Compare
two words agreeing on all already-probed sites, differing at the next
left site $k$, and compensating at unprobed sites to keep the weight
fixed. Normalize the output count change by $\prod_{i\ne k}h_i$
with the current block counts. In the one-apex formula the varying
terms are:
* $\Delta s$ at an initial boundary, or $(e_{k-1}/h_{k-1})\Delta s$;
* $\Delta e(s_{k+1}/h_{k+1})$;
* $\Delta S\sum_{i>k}(E_i/h_i)\prod_{k<j<i}(q_j/h_j)$, and
  $\Delta E\sum_{i<k}(S_i/h_i)\prod_{i<j<k}(q_j/h_j)$;
* the compensation $-2h'_k/15$, with $h'_k=33$.

Here $s,S,e,E$ are as in the companion note, and the local changes
are $+1,+1,-1,-5$, respectively. The compensation weight is unchanged:
all differing remote sites are unprobed, each with count 15, and
their signed number is still minus one. Earlier sites in the negative
overlap sum are unprobed. Already-probed later sites have
$q/h=360/33>8$, $E/h\ge226/33>63/15$, and
$s/h\le9/33<6/15$. Thus the positive overlap can only increase
and the negative adjacent term can only shrink in magnitude relative
to the unprobed comparison. That comparison had strictly positive
change at every site. The next two possible bit values therefore
give distinct exact counts whenever both are allowed, and the
unrevealed completions are still neutral within each branch.
This is a repeatable read and reentry rule, not arbitrary-target
coverage: the two exact branch counts are dictated by the marked
configuration.
