# A factored star seam closes the defect-14 ray

**Proved:** all counts $38\cdot2^j+1$, $j\ge0$, occur as
Hamiltonian path counts of tournaments. The large-height
construction is a literal marked seam whose slope is a product
and whose intercept is a signed two-term expression. This finishes
the entire previously unresolved defect-14 layer, not Chan–Pak.

## Two independently chosen contact profiles

For a tournament $D$ and a subset $S$ of tagged vertices define
$h=H(D)$ and $k=H(D+u)$, where $u$ beats exactly $S$ in $D$.
For a nonempty subset $U\subseteq V(D)$, take each path through
$D[U]$ starting in $S$, and weight it by $H(D-U)$. Sum over
even $|U|$ for $e$, and odd $|U|$ for $o$. Take $H(\varnothing)=1$.
These four integers $(h,k,e,o)$ are all computed in one finite
module, independently of the eventual growing width.

Use two such profiles $A,B$. Orient $B\to X\to A^{*}$ and
$B\to A^{*}$, where $X$ is a transitive block of width $m\ge1$
and $A^{*}$ is the converse of $A$. A single blocker $u$ loses
to $X$, beats precisely the tagged vertices of $B$, and loses to
precisely the tagged vertices of $A^{*}$. Its remaining contacts
are forced. Thus every return between visits to $X$ passes through
$u$. This is exactly the star configuration from the
[matching criterion](MARKED_VISIT_MATCHING.md), with the blocker
on the out-neighbor side (converse covers the other case).

## The exact product and intercept

By the [cycle lift](SPECIALIZATION_AND_CYCLE_LIFT.md), the count is
$U(2)-V(2)+2^{m-1}(V(2)+W(2))$. The three terms classify respectively
no block use, an odd outside return, and an even outside return.

The outside path from the out-neighbors to the in-neighbors of $X$
has form $P_{A^{*}}\,u\,P_B$. The first piece is optional and when
nonempty ends at a tagged vertex; the second is nonempty and starts
at a tagged vertex. Reverse the first piece to read its counts
in $A$. Such an outside path has odd order exactly when the two
pieces have the same parity. Hence
$$
V(2)=2\big((e_A+h_A)e_B+o_Ao_B\big),\quad
W(2)=2\big((e_A+h_A)o_B+o_Ae_B\big).
$$
For $U$, packs either avoid $u$ entirely, have its cycle in just
one arm, or take a cycle through $u$ that visits both arms.
In the last case there is a nonempty piece in each arm with the
same parity. Thus
$$
U(2)=h_Ah_B+(k_A-h_A)h_B+(k_B-h_B)h_A
       +2(e_Ae_B+o_Ao_B).
$$
The mixed terms cancel literally, giving
$$
\boxed{H(X_m)=h_A(k_B-2e_B)+h_B(k_A-h_A)
       +(h_A+e_A+o_A)(e_B+o_B)2^m.}
$$
The factorization does not say that arbitrary integers in these
slots can be chosen independently. Each four-tuple must come from
one literal tagged tournament. The separability makes it possible
to look for the two four-tuples independently before gluing.

## The outstanding ray

A small pair is given by the following outgoing masks, indexed
from zero, and tag masks:

| Module | Rows | Tag mask | $(h,k,e,o)$ |
|---|---|---:|---|
| $A$ | $0,21,25,3,9$ | 2 | $(5,19,7,7)$ |
| $B$ | $2,4,9,19,7,31$ | 56 | $(9,101,63,65)$ |

The intercept is $5(101-126)+9(19-5)=-125+126=1$.
The coefficient is $19\cdot128=2432=38\cdot2^6$.
Use width $m=j-6$ for all $j\ge7$. The accompanying replay also
constructs literal tournaments for $j=0,\ldots,6$, so this proof
of the entire 38-ray does not rely on a prior bulk coverage cap.
A second pair with coefficient $9728$ is included for comparison.

## All defect layers through 2055

The same factorization yields a bank of 1,421 independently checkable
pairs of actual tagged modules. One pair covers each needed odd part
of a dyadic slope. Both arms have order at most eight. The
[bank](star_seam_bank.tsv) records the vertex masks, not just proposed
integer profiles. The replay recomputes its four-tuples with path
subset DP, verifies the two factors, and checks each intercept to be
exactly one. No completeness claim about the source search is used.

Here is why the checked bank applies to infinitely many targets in
each defect layer. In the canonical binary encoding, an odd target
other than an all-one word consists of a prefix ending in 0 followed
by a nonempty run of ones. Only the first 1 in this final run changes
the defect: any further 1 keeps it fixed and sends the marked count
$t$ to $2t-1$. Hence each such prefix gives one whole family
$1+C\,2^j$, where $C$ is found by appending a single 1.

For a target of defect $D$, let its last zero run have length $b$,
starting at marked count $u$. Its contribution to the defect is
$u(2^b-1)-1$, or $u(2^b-1)$ if it is the initial run. All such
contributions are nonnegative, by the
[canonical recurrence](LOW_DEFECT_DESCENT.md). Therefore the marked
count at the end of that run is $u2^b\le2(D+1)$. After $k$ letters,
the marked count is at least $2^k+1$, so for $D\le2055$ the prefix
ends within 12 letters. The replay considers every binary prefix
ending in zero of length at most 12, keeping exactly those whose
defect after the first final 1 is at most 2055. There are 1,421.

For each prefix the verified coefficient $c$ has the same odd part
as $C$. The star construction reads $1+c\,2^m$ at every positive
width $m$, so it covers the original ray at every sufficiently large
height. Any earlier value in that ray is at most $1+c$.
Across the bank this is at most 6,572,033. The repository's
[corrected coverage certificate](https://github.com/eliottcassidy2000/math/blob/main/session203/coverage-recertification/README.md)
already covers the odds through $2^{32}-1$ except 7 and 21, so it
covers all these earlier terms. All-one words have defect zero and
already give their own ordinary count. Consequently **any further
missing odd count must have canonical defect at least 2056**.
The explicit 38-ray construction above does not depend on this
inherited finite coverage step.

The next open task is to select a pair directly from the prescribed
odd part $q$, and control how many factors of two its coefficient
spends. The equation for the intercept is an exact physical
Bézout-type constraint and the coefficient is a product:
$f_A f_B=q\,2^v$. A numerical factor pair or a rational solution
to the intercept equation does not prove that both tagged modules
exist. This bank is a constructive all-height result on finitely
many slopes, not a proof of the unbounded selector.

Run python3 star_factor_replay.py next to star_seam_bank.tsv.
