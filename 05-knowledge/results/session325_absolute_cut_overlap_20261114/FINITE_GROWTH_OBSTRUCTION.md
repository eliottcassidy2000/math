# Finite reservoir banks cannot chase the growing common interface

The two primitive [selective chart collisions](ABSOLUTE_OVERLAP.md) match at a fixed pair of source cores for every common reservoir. Appending the same source block already destroys overlap over one common padded reservoir. Here arbitrary *different* source lengths and *different* reservoir cores and padding are allowed. For a fixed finite reservoir bank the resulting overlap can occur at only finitely many pairs of source growth lengths.

## Precise statement

Let X fill the γ=2 chart and Y fill the γ=6 chart. Suppose their H counts equal h>0. Write their optional slot counts as (q,r) and (q',r'). Assume

    3q r'/(q' r) is not an integral power of 2/3.

Append TT_z to X and TT_w to Y, with z,w≥0 independent. On each side choose any nonempty fixed core from a finite reservoir bank, and append independently any number of singleton SCCs to it. Then only finitely many distinct pairs (z,w) admit any exact common target. There can still be infinitely many targets at each of the finitely many pairs, notably (0,0) in both primitive examples.

For the disclosed profiles with h=571 and h=401, respectively, the arithmetic certificates for the hypothesis are the reduced fractions 26139456/13062793 and 67815/33716. After removing all factors 2 and 3, the residual numerator and denominator differ (15127 and 13062793; 7535 and 8429). Thus both primitive overlaps provably have the finite-growth obstruction.

## Proof

Fix reservoir cores with endpoint intercepts d=h H(V) and e=h H(W) and odd slot parts a=oddpart(Q_V(2)), b=oddpart(Q_W(2)). Absorb their fixed two-adic valuations into the lower bounds on exponents x,y. A common target is exactly

    a·2^x[-3h+6q·2^z+2r·3^z] − b·2^y[-3h+6q'·2^w+6r'·3^w] = d−e.

Both bracketed quantities are odd. If infinitely many growth pairs work, pass to one fixed pair of reservoir cores, and then to a sequence on which at least one growth length tends to infinity. The absolute relation expanded on the left has six terms, plus the constant term if nonzero. Use the [finite normalized block theorem](EXPONENTIAL_OVERLAP.md) and pass to one fixed minimal zero-sum partition with fixed ratios in each block. The term containing the unbounded factor 3^z or 3^w has unbounded three-adic valuation. The only other monomial that can have unbounded three-adic valuation is the opposite 3-term. Thus both growth lengths tend to infinity and those two terms share a block. Fixed ratios imply both z−w and x−y are constant on the subsequence.

If d≠e, the difference of two odd-times-dyadic terms forces min(x,y)≤v₂(|d−e|) (strict equality when x≠y), so fixing x−y bounds both padding exponents. Pass to one pair. If d=e, divide out the common power of two instead. In either case the infinite sequence now gives a fixed equation

    C_0 + C_2·2^z + C_3·3^z = 0.

An unbounded increasing sequence of z forces all coefficients to vanish, by dividing successively by the largest remaining base power. The two nonconstant coefficients vanish only if, with t=w−z,

    a·2^(x−y)/b = (q'/q)·2^t = (3r'/r)·3^t.

Rearranging gives 3q r'/(q' r)=(2/3)^t, contradicting the hypothesis. This proves finiteness. The argument does not produce the last exceptional pair, and it does not apply if new strong reservoir cores are allowed without a finite bank. That is the exact structural escape it points to.

For the special restriction w=z there is an effective finite test: if intercepts differ the normalized exponents have min fixed or bounded by v₂(|d−e|), and both bracketed quantities are bounded between positive constant multiples of 3^z. Thus the larger normalized exponent is also bounded by a fixed computable constant in the off-diagonal case. If intercepts agree, oddness forces x=y. For each of the resulting bounded exponent differences, an equation C_0+C_2·2^z+C_3·3^z=0 has a computable root bound by dominance: beyond a finite threshold the leading term exceeds the sum of the other magnitudes.
