# Both S320 returns enter a strictly ordered cover cone

The two [S320 common-context returns](../session320_second_moments_and_native_contexts_20261113/PRIME_RETURN.md) equalize Hamilton counts but do not remain incomparable. Let X,Y denote the output order used there (the source called X produces X here). For the order-13 return at 532,623 and the order-15 prime return at 4,591,771, the entire ordered-slot differences are respectively

$$
Q_Y(z)-Q_X(z)=2z(z^2-1)(z^8+38z^6+340z^4+708z^2+240),
$$

$$
Q_Y(z)-Q_X(z)=2z(z^2-1)(z^{10}+48z^8+828z^6+4752z^4+7016z^2+720).
$$

These identities are independently checked by two exact recurrences in the companion. First compute Hamilton paths on every subset by last vertex, then obtain ordered nonempty covers by subset convolution. Separately compute every directed cycle count on each subset by rooting at its smallest vertex. The odd-cycle packing recurrence either leaves the smallest vertex uncovered, with weight z, or picks a directed odd cycle through it with weight 2z and recurses on the complement. The two computations agree at every integer evaluation from 0 to the order, sufficient for polynomial identity. The literal masks and both quotients are embedded in the replay.

## Four visits cannot rescue either return, and neither can more visits

For either difference write the right-hand side as z times g(z−1). The displayed factorization gives a polynomial g with nonnegative coefficients, zero constant term, and positive leading term of degree n−3 for output order n. In fact, if the parentheses equal R(z²), then g(u)=2u(u+2)R((u+1)²). Writing F_j for ordered covers by j nonempty paths, inversion of the choice of nonempty slots gives

$$
F_j(Y)-F_j(X)=j\Delta^{j-1}g(0).
$$

The forward difference of u^d at 0 counts onto maps to j−1 labelled bins. Consequently **every** F_j difference is strictly positive for 2≤j≤n−2, and is zero at j=1,n−1,n. In particular the four-visit differences are 26,045,760 and 808,292,160. The sign is not a numerical extrapolation from the low slots.

For any common substitution quotient with one marked vertex, each ordered output path cover projects to quotient walks. The other factors are fixed, and the counts are nonnegative combinations of the marked factor's F_j. Thus this one-sided order is preserved in all output cover coordinates. A common continuation preserves equality of Hamilton counts **if and only if** the marked quotient vertex belongs to a singleton strongly connected component:

* If its component is singleton, the vertex cannot recur in any directed walk, so only the equal F_1 coordinate contributes.
* Otherwise it lies on a directed triangle. Begin with that triangle's walk returning to the marked vertex and insert all unused quotient labels into the walk. This gives a two-visit itinerary with positive coefficient: each other factor fills its single run by a Hamilton path. Thus the strict F_2 coordinate already separates the counts.

Every common strong quotient on at least three vertices therefore breaks equality. The inequality then persists through all further common homogeneous substitutions because the one-run Hamilton contribution has positive coefficient. This strengthens the inherited three-visit obstruction to an exact all-visit classification. A selective contact, unequal context, or different source pair remains possible.

## Packing interpretation

For odd equal-order outputs, expand Q in powers z^(n−2i) using the positive odd-cycle packing formula. If a_i denotes the signed difference of the weighted packing counts at deficit i, equality at z=1 means their sum is zero. Dividing Q_Y−Q_X by z(z²−1) produces, in descending powers of z², the successive cumulative sums a_0, a_0+a_1, … of these differences (the last, zero sum is omitted). The positive parentheses above thus certify a literal ordering of weighted packing prefixes, stronger than merely comparing Q at a handful of positive integers. The resulting positive cone is a useful rejection test when searching for a second exact return: a candidate must leave it before a common cyclic quotient can balance opposing signs.
