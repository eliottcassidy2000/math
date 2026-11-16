# Every even gap by one apex-arc flip, with an outer bank still free

Use the notation of [SHIELD](SHIELD.md). Write rev(s) for ordinary order
reversal of a tape and rc(s) for reversal followed by complement
z→3−z. When a tape consists of common letters 0 and 3, its binary
address reads those letters as 0 and 1 respectively, left to right.

**Theorem.** For every even integer d, there is a central two-apex
spine patch of length O(log(2+|d|)) whose mutual-apex flip has signed
gap d. After choosing d one may wrap the patch by an arbitrary
common-letter word of length m starting in 3 and by its reverse
complement on the other side. All 2^(m−1) literal exterior choices
retain the signed gap exactly and give strong tournaments. No other
arc needs to change. The lower and upper counts are not assigned
absolute targets by this theorem.

## Two arithmetic charts

Choose a common-letter word s of length k≥2, of binary address T,
and set P=s c rev(s), where the central core c is either 12 or 0123.
Both cores are invariant under reversal followed by swapping the
two rows, as is P. In particular this implies

    R_a(P)−L_a(P)=R_b(P)−L_b(P),                        (1)

because the two binary addresses of P are reversals of one another
and R complements the reversed row.

For c write e=L_b(c)−L_a(c), and d_0 for the difference b−a
at its first site. Since c is reverse-swapped, its last difference
is −d_0 and R_a(c)−R_b(c)=−e. Put C=2^k−1 and ℓ=T mod2.
The common prefix s has L=(T,T), and rev(s) has R=(C−T,C−T).
Both outer blocks have zero intrinsic gap, and their weighted
cross determinant vanishes. The two weighted cross determinants
with c sum to e(C−2T). The two local adjacent corrections from
(SHIELD, eq. 3) sum to d_0(1−2ℓ). Therefore

    Δ(P)=Δ(c)+e(C−2T)+d_0(1−2ℓ).                     (2)

Here are the two evaluations (direct use of the contraction count):

| Core c | e | d_0 | Δ(c) | Δ(P) |
|---|---:|---:|---:|---|
| 12 | 1 | 1 | 0 | 2^k − 2(T+ℓ) |
| 0123 | 2 | 0 | 0 | 2^(k+1) − 2 − 4T |

In the first row T+ℓ takes every even integer from 0 to 2^k,
so the gaps are exactly the multiples of four in [−2^k,2^k].
In the second, T runs from 0 to 2^k−1, giving every integer
congruent to two modulo four in [−2^(k+1)+2,2^(k+1)−2].
Choose k with 2^k≥|d|+2 and invert the appropriate displayed
formula. This realizes every signed even d with k logarithmic.

## Preserve the gap while the outside moves

Take any common-letter word u of length m≥1 starting in 3.
The larger tape is W=u P rc(u). The outer blocks have zero
intrinsic gap and their direct weighted cross determinant is
zero. The other two weighted cross determinants sum to

    L_b(u)[R_a(P)−L_a(P)]−L_a(u)[R_b(P)−L_b(P)] = 0

by (1) and L_a(u)=L_b(u). The first and last letters of P are
common, so the local adjacent corrections at both new boundaries
vanish. Thus Δ(W)=Δ(P)=d at all external addresses. No terminal condition on u is needed: the boundaries of P
are common, so either common terminal letter of u works.

The outer first site is 3 and the final site is 0. Both apices
beat the first spine site, and the last site beats both of them.
Every spine site is connected along that spine, so both apex-arc
orientations are strong. The first and last spine sites also
give an odd cycle through one apex for each choice of an even
number of other spine vertices. Thus for a tape of length n,

    1+2^(n−2) ≤ H ≤ 2·3^n,

the upper bound coming from the two three-run colorings in the
[census proof](CENSUS.md). This establishes unbounded genuinely
strong anchors even for each fixed d. The count scale can grow
by adding exterior bits, but those bit choices do not now solve
the independent equation requiring H at a named absolute target.
