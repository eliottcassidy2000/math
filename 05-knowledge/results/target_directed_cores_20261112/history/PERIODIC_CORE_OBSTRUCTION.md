# A ternary ready ladder cannot use a periodic funnel extension rule

2026-11-12. **PROVED**, in the exact operation language below. This
strengthens the finite-bank obstruction by allowing unbounded core orders.
It does not exclude aperiodic extensions, rewriting or deleting parts of a
core, changing the outer clock model, or target-directed core replacement.

## 1. Operation and theorem

Use the inherited one-start funnel construction and six-coordinate state
from [TARGET_RETURN](../../berggren_target_return_20261111/reset/TARGET_RETURN.md):

    (h,e,s,p,f,t)=(P(1),E(1),S(2),P(2),E(2),S(3)).

P(m) counts ordered m-slot spanning path covers, allowing empty slots.
E(m), S(m) prescribe respectively the ending or starting root in the first,
nonempty slot. A legal funnel tail w is empty or ends in1. Its new root
beats the old root, loses to every other old-core vertex, and beats exactly
the zero-bit tail vertices. The old core dominates the transitive tail.
Append means applying this construction to the whole current rooted core,
with the newborn as the new core root.

At each stage surround the core with an outer-tail ready chart
Q(C,L,b,R;k) as in [FAMILIES F1–F5](../../multicore_metagraph_20261113/FAMILIES.md).
The growing TT_k module is the outer clock. Write

    H(Q)=A(C,L,b,R) 2^k+B(C,L,b,R).

**Theorem.** There is no eventual period of nonempty core-extension blocks
and outer contexts that realizes H↦3H+2 at every step for every k≥1.
More generally, no fixed finite nonempty block of legal funnel extensions,
repeated on the core, can multiply the dyadic slope A by 3^q each time,
for any positive integer q, while returning to the same outer context.

The core may gain arbitrarily many vertices. The extension block may
contain several funnel operations, with arbitrary fixed tail words. Empty
tail words are allowed. What is fixed is the repeated operation block and
the outer context at its boundaries, not the core itself. A periodic
collection of outer contexts reduces to one boundary context after a full
period. The clock k is retained; it is not reindexed.

The same conclusion holds for a single ready-count ladder at a fixed
outer tail word, without asking for a family in k; see the corollary in
section4. Thus the obstruction is not caused solely by the dyadic-family
ansatz.

## 2. The response block and its exact 3-adic cost

For a tail word w of length l, the inherited coefficients are

    L_m=(m−1)^l+sum_(i:w_i=1) m^i(m−1)^(l−i−1),
    D_m=(m−1)^l+sum_(i:w_i=0) (m−1)^i m^(l−i−1).

B_m is the inherited nonnegative six-state tail-word count. Its state is
(last,started), initially(0,0). At bit b there are transitions: retain the
state with weight m−1; replace last by1 if b=1 or2 if b=0 with weight1;
and set started=1 with weight1 if it is already1 or b=0. Sum final states
with last≠2. This convention is fixed in the attached checker.

The last three coordinates evolve independently of the first three:

    [p']   [2B_2  −2D_2  2L_3] [p]
    [f'] = [ L_2    −1      0 ] [f] = M_w [p,f,t]^T.             (1)
    [t']   [ D_3     0     3^l] [t]

Let r(w) be the number of initial1 bits, including r=0 for the empty
word. Then

    v_3(det M_w)=r(w).                                         (2)

Proof. Direct expansion gives

    det M_w=2[(D_2L_2−B_2)3^l+L_3D_3].                         (3)

If w is nonempty and ends in1, D_3≡2^l modulo3. If r<l, sum the initial
run in L_3 to obtain

    L_3=3^r 2^(l−r)+terms divisible by3^(r+1).

The zero at position r supplies no term. Thus v_3(L_3)=r, and the first
summand in(3) is divisible by3^l, of strictly greater valuation. If r=l,
the word is all ones, L_2=B_2=2^l, D_2=1, L_3=3^l, D_3=2^l;
then det M_w=2^(l+1)3^l. The empty word has determinant2. This proves(2),
including nonvanishing of every determinant.

## 3. The same initial run has a strictly larger physical growth cost

For every actual nonempty rooted tournament,

    E(2)≤P(2)/2,       S(3)≥P(2)/2.                            (4)

The first inequality follows by swapping the two slots: covers ending at
the root in slot1 and in slot2 are disjoint equally large subsets of all
two-slot covers. For the second, split the slot containing the root just
before it. Put the root-starting suffix in slot1, its prefix in slot2, and
the other old slot in slot3. Retain one bit specifying which old slot held
the root. This injects all two-slot covers into two copies of the S(3)
covers. Empty prefixes and empty other slots are permitted.

For a nonempty tail ending in1, B_2≥2^(l−1): on its first l−1 bits use
only the retain/replace-last transitions, then replace last by1 on the
final bit. These are 2^(l−1) accepted transition histories. Also

    D_2≤1+2+4+...+2^(l−1)=2^l−1.

Hence 2B_2−D_2≥1. The same inequality holds for the empty word. Using(1)
and(4),

    p'≥(L_3+2B_2−D_2)p≥(L_3+1)p>3^r p.                       (5)

The last comparison uses L_3≥3^r, directly from its nonnegative expansion.
Thus every leading1 supplies exactly one factor3 to the determinant's
3-adic valuation, but the positive path-cover response grows by **strictly
more** than that factor. The strict excess is the useful lost coordinate.

## 4. Proof of the theorem

Let a repeated nonempty extension block have product matrix M, and put
R=sum r(w) over its words. Then det M is nonzero with valuation R. By(5)
there is a fixed real number c>3^R such that

    p_j≥c^j p_0                                                   (6)

along its actual repeated core orbit.

For the fixed outer context, FAMILIES F4–F5 gives

    A=alpha_L(S(2)+P(2))+alpha_B H,
    alpha_L>0,       alpha_B≥0.

Both H and S(2) are at most P(2), so

    alpha_L p_j≤A_j≤(2alpha_L+alpha_B)p_j.                       (7)

Suppose A_j=A_0 lambda^j with lambda=3^q. Equations(7) bound p_j both
above and below by positive constants times lambda^j. Its generating
function is the rational function

    sum_(j≥0) p_j z^j = e_1^T (I−zM)^(-1)[p_0,f_0,t_0]^T.

The positive lower bound forces this function to diverge as real
z approaches1/lambda from below; the upper bound gives convergence at
every smaller nonnegative z. Therefore det(I−M/lambda)=0, so lambda
is an eigenvalue of M. Since lambda is an integer, division of the monic
integer characteristic polynomial by X−lambda shows that lambda divides
det M. Consequently q≤R.

But(6) now grows strictly faster than lambda^j, contradicting the upper
bound from(7). This proves the theorem. Notice that no bounded core bank,
finite census, or generic spectral assumption occurs in the argument.

**Fixed-target corollary.** Let a fixed legal outer word w complete each
successive core to its ready source of count H_j. A repeated nonempty
funnel extension block cannot satisfy

    H_j+1=3^(qj)(H_0+1),       q≥1.

Indeed its exact count formula is

    H_j+1=L_2(S(2)+P(2))+B_1 H−D_1 E(1)−E(2).

Here L_2≥1, B_1≥1 and D_1 is0 or1. Since E(1)≤H and E(2)≤P(2)/2,

    P(2)/2≤H_j+1≤(2L_2+B_1)P(2).

The same rational-series/determinant/growth argument applies verbatim.
Periodic outer words reduce to one word after a full period. This does
not cover a growing or aperiodically selected outer word, an interior
rewrite, or a change of port language.

## 5. What this says about the 462 family

The inherited all-k chain is

    462·2^k−5 → 1386·2^k−13 → 4158·2^k−37.

Its core addresses, written from the rooted singleton outward, are

    (1),       (11,1),       (1,11,empty).

Their six-states are respectively

    (3,1,3,12,3,7),
    (15,5,15,216,57,115),
    (27,12,30,564,171,327).

These are actual **rewrites** of the inner address, not a repeated
extension of the preceding core. The second and third addresses even
have the same total number3 of initial1 bits in their constituent words,
but order those words differently. That is precisely a direction not
excluded by the theorem. A prospective infinite completion should retain
the noncommuting response block, not merely the scalar target or total
number of1 bits.

The shifted target H+1 is multiplied by3 at each ready step. This matches
the old −1-anchor collision algebra only at its affine target coordinate;
the determinant and positive cover-growth data are additional tournament
sidecars. Scalar Collatz collisions cannot be inserted as physical core
rewrites without checking those responses.

**Remaining positive task.** Find an aperiodic extension compiler, or an
exact address-rewriting rule with changing interior core, that supplies
the next slope/intercept pair. The result above does not rule out either.
It explains why a newly found finite chain should not be completed by
guessing a periodic append-only tail rule.
