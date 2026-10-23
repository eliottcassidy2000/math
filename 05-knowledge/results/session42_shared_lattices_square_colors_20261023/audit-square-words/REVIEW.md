# Ordered square-label reconstruction: independent audit

Status: **PASS**, with exact carrier restrictions below. This packet audits the
Session42 root discovery independently; it makes no literature-priority claim.

## 1. The all-size reconstruction theorem

Let p=(p_0,...,p_(N-1)) be a permutation of {1,...,N}. Retain the **ordered**
adjacent sums s_i=p_(i-1)+p_i, i=1,...,N-1. Define c_0=0 and
c_i=s_i-c_(i-1). Any realization has

    p_i=(-1)^i x+c_i,     x=p_0.

Consequently an actual vertex value at any specified position recovers the
entire path. Without that marked value:

* For odd N, sum_i p_i=N(N+1)/2 has coefficient one on x, so there is at most
  one realization.
* For even N, the fixed sum of squares N(N+1)(2N+1)/6 gives a quadratic in x
  with leading coefficient N. Therefore there are at most two realizations.

The moment proof works for any prescribed finite set of distinct real vertex
values, with its own fixed first and second moments. The next classification
uses the special interval {1,...,N}.

## 2. Exactly which interval paths admit a twin

Suppose two realizations p,q have q_0=p_0+d, d>0. Then
q_i-p_i=(-1)^i d. Write A for the values in p's even positions and B for its
odd-position values. Both A disjoint-union B and (A+d) disjoint-union (B-d)
are {1,...,N}.

Using the indicator polynomials A(z),B(z), and P(z)=z+...+z^N gives

    A+B=P,     z^d A+z^(-d) B=P.

Multiplying the second relation by z^d, subtracting the first, and cancelling
z^d-1 in the integral polynomial domain yields

    (1+z^d) A(z)=P(z).

If zeta is a primitive 2d-th root of unity, then zeta^d=-1 and zeta is not one.
Thus P(zeta)=0 forces zeta^N=1, proving 2d divides N. The quotient is now exact:

    A(z)=(z+...+z^d) sum_(j=0)^(N/(2d)-1) z^(2dj).

Hence A is precisely the union of alternating blocks of length d, beginning
with {1,...,d}. B is its complement. This is necessary, not merely one
convenient construction.

Conversely, any permutation alternating those two block sets has a twin:
add d at A values and subtract d at B values. This preserves every adjacent
sum and permutes {1,...,N}. If the original starts in B, reverse the signs of
the first displacement. The involution exchanges the two realizations.

The proof shows why a missing endpoint/position owner is significant even
when every edge's exact arithmetic label has been retained.

## 3. Exact count and divisor generating function

Put N=2m. For each d dividing m, the permutations starting in A_d and
alternating A_d,B_d number (m!)^2. Each is the lower-start member of exactly
one twin pair. Distinct d cannot count the same pair, since d is the difference
of its initial values; nor can one path participate in different pairs,
because Section 1 permits at most two realizations.

Thus, writing tau(m) for the divisor count:

    double-fibre ordered-sum words = tau(m) (m!)^2;
    distinct ordered-sum words = (2m)! - tau(m) (m!)^2;
    fraction of permutations with a twin = 2 tau(m)/binom(2m,m).

For odd N all N! permutations have distinct sum words. The normalized
counting series has the exact formal Lambert-series identity

    sum_(m>=1) [double-fibre count/(m!)^2] z^m
      = sum_(m>=1) tau(m) z^m
      = sum_(d>=1) z^d/(1-z^d).

The probability has the asymptotic

    2 tau(m)/binom(2m,m) ~ 2 tau(m) sqrt(pi m)/4^m.

The divisor factor is retained; it is not an asymptotically constant
coefficient. This follows directly from the central-binomial Stirling
asymptotic. In particular the ambiguity becomes exponentially rare among
all interval permutations, despite the exact infinite families.

### Reversal and palindromic clocks

For N>1, a palindromic adjacent-sum word makes p and its reversal distinct
realizations. Odd N is therefore impossible. For N=2m the reversal must be
the unique twin tau_d(p). For each divisor d of m, choose an ordering of the
m pairs of tau_d (m! choices), and choose the initial block side (2 choices).
The first m positions select one member from each ordered pair, with sides
forced to alternate. The last m positions are tau_d applied to the reversed
first half. This yields exactly all palindromic-clock permutations for d:
2m! permutations, or m! sum words after dividing by the two reversed paths.
The twin displacement makes different d disjoint.

Thus the number of palindromic words is tau(m)m!, and Burnside for word
reversal gives

    number of sum-word classes modulo reversal
      = ((2m)! - tau(m)(m!)^2 + tau(m)m!)/2.

For odd N>1 the count is N!/2. At N=1 there is one empty sum word, fixed by
reversal, so that initial case is separate. These formulas are independently
checked in the complete permutation census through N=8. They are not counts
of square-restricted Hamiltonian paths.

These counts concern all interval permutations. Restricting sums to squares
selects Hamiltonian paths in a smaller graph and does not preserve the count.

## 4. The first square-sum ambiguity is exactly order 36

Let Q_(N,d) retain precisely those edges of the square-sum graph whose
endpoints belong to opposite alternating d-block sets. By Section 2, a
square-sum twin exists if and only if some Q_(N,d), d dividing N/2, has a
Hamiltonian path.

All odd orders are unique by Section 1. For every even N<36 and every divisor
d of N/2, the independent checker constructs Q_(N,d) and supplies one of
three elementary impossibility certificates:

1. an isolated vertex;
2. more than two degree-one vertices;
3. a proper connected component.

There are exactly 52 such (N,d) pairs. Their complete certificates are stored
in receipt.json. No Hamiltonian search is needed for this finite exclusion.

At N=36,d=6 the root-supplied path is

    3,33,16,9,27,22,14,11,25,24,1,35,29,20,5,31,18,7,
    2,23,26,10,6,19,30,34,15,21,4,12,13,36,28,8,17,32.

Applying the block involution produces the independently checked twin

    9,27,22,3,33,16,20,5,31,18,7,29,35,14,11,25,24,1,
    8,17,32,4,12,13,36,28,21,15,10,6,19,30,34,2,23,26.

Both are permutations, all adjacent sums are squares, and their ordered sums
agree. The second is not the reversal of the first. A direct reconstruction
from each possible initial value 1,...,36 finds exactly these two solutions.
Together with the 52 obstruction certificates, this proves the minimality
claim, rather than merely reporting the first example found by a search.

## 5. Hostile scope checks and reproduction

This is not a statement about a multiset of sums or unpositioned edge labels.
It is also not true for partial paths: on {1,...,9}, the single sum 10 has
eight distinct two-vertex realizations. Full vertex-set coverage is essential
to the two moment constraints.

The executable independently groups all permutations through N=8 by their
sum words, reconstructs all fibres, checks their size and block form, and
checks the divisor count. It then verifies all 52 finite obstruction
certificates and the order-36 witness. Every condition raises an exception
rather than relying on Python assertions.

    python3 check.py
    python3 -O check.py

The two outputs are byte-identical. The packet does not enumerate all
Hamiltonian paths of Q_36, classify all square-sum graphs, prove a Berggren
island map, or infer tournament metagraph equivalence from this path code.
