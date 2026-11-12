# How the family grammars intersect

**PROVED elementary corollary of the reviewed births and wildcard grammars.**

Let D be the ordinal/C3 hereditary class in `wildcard/PROOF.md`. Let T(b) be a marked apex over a transitive spine with binary word b, using `births/PROOF.md`'s orientation. Let r(b) be the number of adjacent10 transitions.

Then

    T(b) belongs to D  iff  r(b)<=1.                         (B1)

If r=0 the tournament is transitive. If r=1 its normal form has a single C3 quotient and transitive children, with ordinal prefix and suffix, hence belongs to D. If r>=2 its induced-prime inventory contains W5. Every induced subobject of D lies in D, whereas a prime five-vertex tournament cannot be made by ordinal/C3 substitution (any nontrivial expansion leaves a proper module). Hence T(b) is outside D. The apex marking does not change membership in D.

This computes the complete intersection, not merely the first hostile. The earliest escape is a five-vertex W5, but the same criterion decides the intersection at arbitrary order.

For every m>=0 the number of length-m binary words with exactly r descents10 is

    binom(m+1,2r+1),                                        (B2)

with out-of-range binomials zero. To see this, adjoin initial0 and final1. The padded word has2r+1 change points among its m+1 adjacent positions, and any odd-sized choice of change points reconstructs exactly one such padded word. Therefore the number of apex-marked length-m words lying in D is

    (m+1)+binom(m+1,3),

and their ordinary generating function by total tournament size is

    z/(1-z)^2 + z^3/(1-z)^4.                               (B3)

These are MARKED isomorphism counts: forgetting the apex can merge words. The proportion of the2^m marked one-apex classes lying in D tends to zero, although both source families possess complete concise grammars. Neither complete family is representative of every large tournament.

The batch grammar supplies the next hierarchy: tournaments with a specified deletion set of k vertices leaving a transitive spine are encoded losslessly by a tournament on those k vertices and a word over2^k signatures. The k=1 slice is exactly the feedback-word family; k growing gives an exhaustive marked language, but classification of its prime inventories and unmarked quotients is not uniformly solved by the k=1 formula. This identifies a concrete direction for further structural understanding, rather than an unsupported extrapolation.
