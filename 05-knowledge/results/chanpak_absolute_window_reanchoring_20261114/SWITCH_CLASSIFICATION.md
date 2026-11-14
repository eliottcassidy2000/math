# Exact contact-switch independence: two-site switches and the carry exception

For an ordered transitive spacer of m vertices between fixed marked
cores, the only dependence of its apex count that is nonlinear in
the contact bits is a positive common multiple of g(z). See
[the spacer calculation](INTERLEAVED_CONTACT_RETURN.md).
Index contacts from least significant bit. The weighted-inversion
form of g gives a mixed second difference for disjoint switch
vectors X,Y (entries -1,0,1):

    -h_left h_right sum_(i in supp X,j in supp Y) X_i Y_j W(|i-j|),
    W(d)=2^(d-1)+[d=1].                                 (1)

Each changed bit switches between its binary endpoints. No disjoint
switch interferes with its admissibility. On a simultaneous cube of
such switches, (1) is the only interaction: the response is
quadratic, so it is additive iff every pair has zero mixed term.
This is the precise transfer from the old knot-curvature note:
a mixed difference detects interaction even when both marginal
increments are well behaved. No knot equality is used.

## Separated supports: the precise exception

Suppose supp X lies strictly below supp Y. Let i0 be the last active
position of X, j0 the first of Y, and d=j0-i0. Put
A=sum X_i 2^(i0-i), B=sum Y_j 2^(j-j0); both are odd nonzero integers.
The generic mixed sum is 2^(d-1) A B. The only extra term occurs at
d=1 and equals X_i0 Y_j0. Thus cancellation is impossible for d>=2,
and for d=1 requires |A|=|B|=1 and opposite signs for AB and the
boundary product.

If such a signed binary sum is on a nonsingleton support and has
absolute value one, its sign is opposite to its lowest digit:
the excess from the higher powers is nonzero and even, so cannot
have left the signed integer equal to the lowest digit. Moreover
its active places must be consecutive and the leading digit must
oppose all the lower ones, since only the sum of *all* smaller powers
can reduce the largest power to absolute value one. This is exactly
one binary carry. A singleton, in contrast, has the sign of its
only digit. Therefore separated supports can cancel **exactly**
when they touch, one flips a singleton, and the other is a full
signed carry. This corrects the preliminary overstatement that
every separated pair interacts. Example: positions 0,1 switch
(+1,-1) and position 2 switches +1. Since W(1)=W(2)=2, sum (1)
is zero. Two separated switches that each affect at least two
bits cannot be independent.

## Complete two-site classification

A two-site switch, up to direction, has entries (+1,sigma) at its
left and right sites, with sigma=+1 (symmetric) or -1
(antisymmetric). Pairwise disjoint two-site switches with zero
mixed term must nest.

- Separated pairs are excluded by the preceding theorem.
- For a crossing a<c<b<d, put L=c-a, G=b-c, R=d-b. The four
  distances are L,G,R,L+G+R. When any gap is at least two,
  W(L+G+R) is at least four times the largest of the other three,
  so cancellation fails. If all three gaps are one, the sum of
  three signed 2s cannot cancel a signed 4.

For nested sites a<c<d<b, write L=c-a, R=b-d, G=d-c.
If both switches are symmetric, every term is positive.
If both are antisymmetric, the mixed sum is
W(L)+W(R)-W(L+G)-W(R+G), zero iff L=R=G=1.
If signs are opposite, the equation is equality at L,R for either
D_G(n)=W(n+G)-W(n) (outer symmetric) or
U_G(n)=W(n+G)+W(n) (outer antisymmetric).
Both are strictly increasing for n>=1: for n>=2 the terms scale
by two, and the first increase follows directly from W(1)=W(2)=2
and W(3)=4. Thus the opposite-sign pairs work iff L=R, i.e. they
share a center.

Consequently there are at most three mutually independent,
disjoint two-site switches, at every spacer length. At most one
can be symmetric; two antisymmetric switches must occupy the
same four consecutive sites, and no third antisymmetric switch
can share that requirement with either. In a three-switch cube
the third switch is necessarily symmetric, surrounds the four
sites, and has their center. This classifies the support of the
six-site physical return, without claiming optimality for switches
with more than two sites. In particular, the difference of two
selected corners can be tiny despite the individual increments
being widely separated: exact additivity is not the same as a
full interval of offsets.
