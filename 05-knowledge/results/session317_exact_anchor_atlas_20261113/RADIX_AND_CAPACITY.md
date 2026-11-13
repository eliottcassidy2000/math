# Mixed-power digits can cover the integers; a fixed independent tournament bank cannot

These are two separate theorems. The arithmetic one isolates a
sufficient interval-gluing rule. The tournament one shows that simply
making one fixed endpoint-matched palette larger will not meet it.

## An exact mixed-power decoder with no normalization step

Let $q>h\ge2$ be coprime integers. Given $D\ge(q+1)(h-1)$, **every
nonnegative integer** belongs to some set
$$
 S_n=\left\{\sum_{i=0}^n d_iq^i h^{n-i}:0\le d_i\le D\right\}.
$$
This is an unnormalized integer assertion, not the ordinary rational
base expansion of $h^n$ times an input. The mixed powers are the
same as in the older
[rational-base lonely-runner atom tree](../../../01-canon/theorems/THM-3848-rational-base-prefix-atom-tree-and-lonely-runner-separation.md),
but the proof here keeps both an integer residue and a real overlap
interval instead of a single allowed sheet.

Define $A_0=0,\ B_0=D$ and
$$
 A_n=hA_{n-1}+(h-1)q^n,\qquad
 B_n=hB_{n-1}+(D-h+1)q^n.
$$
Inductively $[A_n,B_n]\cap\mathbb Z\subseteq S_n$.
In $hS_{n-1}+q^n[0,D]$, a desired residue modulo $h$ selects
one residue of the last digit. Its least allowed digit is at most
$h-1$ and its largest at least $D-h+1$. Successive such digits
are $h$ apart, so their previous-level intervals, after dividing
by $h$ within the residue, are shifted by $q^n$.
They overlap or touch as integer intervals if
$B_{n-1}-A_{n-1}\ge q^n-1$. This holds at $n=1$ and propagates:
$D\ge q+h-1$ implies
$$
 B_n-A_n\ge h(q^n-1)+(D-2h+2)q^n\ge q^{n+1}-1.
$$
Thus any prescribed integer in the displayed interval has an exact
digit selected by the residue and the interval bounds, then descends.

It remains to connect levels. Substituting
$A_n=(h-1)q(q^n-h^n)/(q-h)$ and
$B_n=D h^n+(D-h+1)q(q^n-h^n)/(q-h)$ shows that at
$D=(q+1)(h-1)$ one has $B_n-A_{n+1}=(h-1)h^n\ge0$.
Larger $D$ only helps. The intervals begin at 0 and tend to infinity.
This proves complete absolute coverage in the arithmetic model.

## The capacity inequality for the physical one-sided construction

Take a finite bank of types. Type $t$ has fixed
$(h_t,q_t,s_t,S_t)$ and a finite menu of actual blocks/cuts varying
only in $(e_t,E_t,b_t)$, in the notation of the
[common-apex additive law](../session313_zeckendorf_target_selection_20261113/ADDITIVE_PORTS.md).
Assume nontrivial cuts, so $S_t,E_t>0$, and let $M_t$ be the number
of distinct $E_t$ values. Any sequence of types is allowed. Then
$$
 \sum_t \frac{M_t}{q_t}<1
 \quad\Longrightarrow\quad
 \#\{H\le X\text{ in the bank}\}=O(X^\alpha)
 \text{ for some }\alpha<1.
$$

For a fixed ordered type sequence, fix the word of $E$ responses.
All long overlap terms are then fixed, including the factors $S$ and
$q$. The remaining sum, divided by $\prod h_i$, consists of the
boundary terms $s_1/h_1,e_n/h_n$, the terms $b_i/h_i$, and
the terms $e_i s_{i+1}/(h_i h_{i+1})$. For a finite bank these lie
in a fixed rational lattice of bounded denominator and their total
range is $O(n)$. Thus at most $C(n+1)\prod_i M_{t_i}$ different
counts come from a type sequence. Conversely the first-to-last
overlap term alone gives $H\ge c\prod_i q_{t_i}$ for $n\ge2$,
with $c>0$ uniform over the bank. Choose $\alpha<1$ close enough
to 1 so that $\sum M_t/q_t^\alpha<1$. Sum the count bound
over $n$ and over all type sequences with $\prod q_{t_i}\le X/c$,
using
$\prod M_{t_i}\le (X/c)^\alpha\prod (M_{t_i}/q_{t_i}^\alpha)$.
The resulting weighted geometric sum converges, proving the claim.

In particular **every single fixed one-sided type is sparse**:
$E$ counts a subset of ordered two-slot covers in which the marked
slot is nonempty, so $1\le E\le q-h$ and $M\le q-h<q$.
Correlated split/end offsets do not rescue it; they contribute only
the polynomial factor just accounted for. Zero-cut degeneracies are
excluded from the stated theorem. With $S$ fixed, $E=0$ would force
the cut to be full and $S=q-h$, which forces every cut in that
type to be full. If $S=0$ it forces all cuts empty.

Thus the old modular selectors were not merely missing a clever
choice of digits: an entire single-type mechanism fails on real
counting capacity. To translate the positive arithmetic decoder,
one must use genuine changing type, retained interactions, or a
different scale. The inequality is a rejection test only; violating
it does not imply coverage or valid endpoint overlap.
