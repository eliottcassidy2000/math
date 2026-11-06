# A complete digit-and-gcd inverse for positive unit readouts

Status: **PROVED, independently audited in the stated domain**. This is an exact arithmetic
inverse and a complexity bound, not a proof of universal occupancy. The
twelve-template all-target lemma remains **OPEN**. Session 201.

## 1. Physical readouts and the scalar domain

Use the reviewed function g(1)=1, g(2h)=g(h), and
g(2h-1)=g(h)+h-(h mod2). The odd recurrence is used for h>=2;
its extension at h=1 is the true scalar identity1=1. For all h>=1,

    1<=g(h)<=h,                 g(h) is odd.                 (1)

These follow inductively; in particular g(2h)<=h handles the even-parent
case in the upper bound for an odd argument.

Fix integers alpha,beta,gamma,delta satisfying

    alpha>=1 odd,   beta>=alpha,
    -alpha<=gamma<=0,   alpha-beta<=delta<=0.                (2)

We seek all integer v>=2 such that

    n=alpha g(v)+beta v+gamma(v mod2)+delta.                 (3)

This domain includes every physical unit-source readout. For an actual
source s=(D,E,U,d), the coefficients of K(s,tau(v)) are

    alpha=d, beta=U-d, gamma=E-d,
    delta=D-E-U+d=-f.

The reviewed nonnegative coordinates imply U>=2d, 0<=E<=d and
0<=f<=U-2d, giving (2). The source operations preserve odd D,d and
U=E modulo2 from the seed, so alpha is odd. For a central template,
adjunction replaces s by the actual source M_Wdagger tau(u).

The theorem below is numerical on the entire domain (2); it does not
declare every coefficient tuple in (2) to be a physical source.

## 2. Binary peeling preserves a fixed gcd

Split the current index h by its parity. For h=2j, (3) becomes a row with

    (alpha,beta,gamma,delta) -> (alpha,2beta,0,delta).         (4)

For h=2j-1 it becomes

    (alpha,alpha+2beta,-alpha,delta+gamma-beta).              (5)

Both rows still satisfy (2). The first coefficient alpha stays fixed,
beta at least doubles, and

    d0=gcd(alpha,beta)                                      (6)

is unchanged, since alpha is odd. At depth k,
beta_k=2^k beta_0 modulo alpha. Original indices are tracked exactly as

    v=s h-r,  s>=1, 0<=r<s.

The even branch sends (s,r) to(2s,r); the odd branch to(2s,s+r).
These are disjoint parity branches, so they lose or duplicate no original
index. An internal h=1 is permitted as scalar arithmetic, including the
tautological odd identity. It is not used as a physical Z operation on a
singleton. Every returned v>=2 is compiled afterward by its canonical
nontrivial unit-source word.

## 3. A short complete interval and at most two leaf candidates

For every h>=1, equations(1)-(2) give

    beta(h-1)+alpha <= alpha g(h)+beta h+gamma(h mod2)+delta
                    <= (alpha+beta)h.                     (7)

Therefore a solution lies in the exact necessary interval

    ceil(n/(alpha+beta)) <= h <= 1+floor((n-alpha)/beta).     (8)

Intersect this with h>=1 and the requested original-index range under
v=s h-r. A target n<alpha has no solution. Let L,H be the endpoints of
the un-intersected interval. Its width satisfies

    H-L < 1+n alpha/beta^2.                                (9)

For a parity p in{0,1}, write h=2t+p. A solution necessarily satisfies

    2beta t = n-delta-(beta+gamma)p  (mod alpha).            (10)

If d0 does not divide the right side, that parity is impossible. Otherwise
divide by d0 and invert modulo m=alpha/d0. This gives one t progression
modulo m, or all integers t if m=1. Thus possible h of this parity are
spaced by2alpha/d0. The original exact equality, including g(h), must be
checked afterward.

Stop binary splitting when

    beta^2 >= n d0.                                       (11)

At such a leaf, (9) gives H-L<1+alpha/d0<=2alpha/d0. Hence the interval
contains at most one congruence candidate of each parity, and at most two
g evaluations decide it completely. This does not accept a congruence
candidate without evaluating g.

The gcd divisibility filter can also be used before the stopping depth;
an impossible parity branch is discarded immediately. No empirical
height cutoff or probabilistic primality assumption enters this inverse.

## 4. An explicit search-node bound

Let K be the least nonnegative integer for which

    4^K beta_0^2 >= n d0.                                 (12)

Every branch stops by depth K, so the complete binary search visits at
most2^(K+1)-1 nodes and evaluates g at most2^(K+1) times. In particular
the number of nodes and g evaluations is

    O(1+sqrt(n d0)/beta_0).                               (13)

This is a bound for one fixed positive readout, not an existence theorem.
The modular inverse can be prepared once per depth: since
beta_k/d0=2^k beta_0/d0 modulo m, its inverse is the initial inverse times
2^(-k) modulo the odd integer m. All scalar g evaluations terminate in
O(log n) recurrence steps in the target-bounded applications below.

For the fixed finite central-word family, use its complete product-bounded
domain2<=u<=v and restore the full transpose-closed matrix menu. For each
u,W apply the preceding procedure to the right index. The physical row
satisfies beta>=u-1: for W empty this is equality; for nonempty central W,
the final J in Wdagger has beta equal to its predecessor's U>=u.
Also d0<=alpha<=beta. Thus the per-row node bound is at most

    O(1+sqrt(n/(u-1))).

There are a fixed number of matrices and u<=sqrt(2n). Summing yields a
complete target membership algorithm with

    O(n^(3/4) polylog n) bit operations.                   (14)

The polylog factor accounts for exact integer arithmetic, modular inversion
and recurrence evaluation. An average-gcd improvement is not asserted.
This improves the naive scalar/progression census without altering its
target universe. It produces actual index certificates or a complete
nonmembership result for the selected finite template family.

## 5. Verification, targeted use, and limits

The author checker compares complete results with direct index enumeration
on104 rows throughout the numerical cone (2), at all targets1 through201.
It checks every visited row, the invariant gcd, exact original-index
reconstruction, at most one leaf candidate per parity, and the full binary
node budget. It also compares complete template answers with the older
progression inverse on named positive and hostile targets.

The larger runs are targeted prime-constellation probes: n,3n-2,5n-4 are
selected prime together, excluding the two exact rank-one divisor templates
identified in the separate factorization packet. These are individual
constructive tests, not a larger interval census. Their returned indices
are checked directly by (3) and the physical graph compiler.

Finite successful probes do not prove occupancy, and a first-witness run
does not enumerate all representations or certify optimal indices. Only
an exhausted complete run can certify nonmembership in the finite menu.
The original tournament spectrum and the twelve-form occupancy assertion
remain OPEN.
