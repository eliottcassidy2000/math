# Literal ×2 and ×3 are synchronized in the fixed-step port

The [six-vertex physical route](PORT_ROUTES.md) reads the optional
slot profile of its module D as

    4Q₃+9Q₂−3H  →  4Q₃+9Q₂−H  →  4Q₃+9Q₂+H.

Fix the modulus M=2h by fixing H(D)=h. There is a complete effective
absolute inverse for this lane, and an all-size obstruction to
treating its ×2 and ×3 as separate clocks.

Every nontrivial strong component of D has at least three paths
and order at most 1+3 log_5(h) by
[Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
There are at most log_3(h) nontrivial components, since H
multiplies along the ordered condensation. Enumerate all their
finite lists with H-product h. Each such list contributes a pair
(B,C) of products of optional slot counts Q₂ and Q₃. A singleton
component leaves H unchanged and multiplies Q_j by j. Indeed a
cover slot in an ordinal sum concatenates its restrictions, so
Q_j is multiplicative. Thus all possible module profiles are exactly

    H=h,     (Q₂,Q₃)=(B·2^s,C·3^s),      s≥0,

over the finite catalogue. For the empty core require s≥1 so D is
nonempty. This is the joint version of the
[fixed-H two-slot decoder](../chanpak_graceful_ternary_20261115/BERGGREN.md):
the **same** exponent s is mandatory in both coordinates.

To hit a prescribed middle target N, solve for each core list

    N+h=4C·3^s+9B·2^s.

The right side is strictly increasing, so exact bisection finds the
only possible s for that list; test the equality before reporting
a route. If a list succeeds, the actual distinguished module and
both bracketing witnesses are reconstructed by its ordinal sum and
the two literal arc flips. For fixed h and N≤X there are only
O_h(log X) middle counts, even allowing every D. Their residue
class can mix while the actual heights remain exponentially
separated.

This is a literal ×2×3 appearance with neither of the independence
hypotheses of a full ×2×3 action: it moves on one synchronized
semigroup ray in profile space. Choosing Q₂'s s and Q₃'s s
independently would insert a nonexistent physical module. The
failure is not a numerical prime-order coincidence; it is the
constraint that both clocks come from the same singleton budget.
