# The all-row depth initial module and its two-step staircase

**PROVED IN THIS NOTE.** This is an all-row homogeneous-interface theorem,
not an all-row inhomogeneous solution or termination result. It upgrades
the finite homogeneous depth computations and the conditional weight
staircase in the initial Session254 note.

## 1. Exact initial module at every depth

Let K be a field, x,t independent indeterminates, and set

    u=x^2t, p=t(1+x^2t), y=xtp.

For an integer d>=0 define the K-vector subspace

    P_d = sum_(a,b>=0, a+b<=d) x^a u^b K[p,y]  ⊂ K[x,t].

For n>=0 let I_(d,n) be its pure top-row coefficient space:

    I_(d,n)={ [t^n]f : f in P_d and [t^i]f=0 for all i<n }.

For every d>=0 and n>=d,

    I_(d,n) = K[x]_(<= floor((n+3d)/2)).                   (1)

This describes exact leading coefficients of genuine elements of P_d,
not merely an upper bound from a finite matrix. The theorem is valid in
any characteristic. Its later application to the normalized Jacobian
chart is over characteristic zero.

### Proof: the two-point polynomial bound

Invert x temporarily and substitute t=u/x^2. A generator becomes

    x^a u^b p^c y^e
      = x^(a-2c-3e) u^(b+c+2e)(1+u)^(c+e).               (2)

Decompose a finite linear combination into Laurent x-weight sectors.
In the sector of weight -h, every generator obeys

    2c+3e=h+a.

Its polynomial in u has degree

    b+2c+3e=b+h+a <= d+h.                                (3)

Suppose x^r occurs with nonzero coefficient in [t^n]f for an f having
all earlier rows zero. This term lies in the weight -h sector with
h=2n-r. The polynomial F_h(u) in that sector has u-adic order exactly n.
Indeed different Laurent x-weights cannot cancel, and a term u^i in
that sector corresponds to x^(2i-h)t^i. Equation (3) implies
h>=n-d>=0.

For every generator in this sector,

    c+e >= (2c+3e)/3 = (h+a)/3 >= h/3.

Consequently (1+u)^ceil(h/3) divides F_h. The points 0 and -1 are distinct
in every field. The orders at these two points therefore satisfy

    n+ceil(h/3) <= deg F_h <= d+h.

In particular 3n<=3d+2h. Since r=2n-h, this yields

    2r<=n+3d.

Every nonzero leading coefficient thus has the degree bound in (1).
Cancellation of earlier rows cannot bypass it: that only increases
the zero order at u=0 while retaining the zero at u=-1 and degree cap.

### Proof: every allowed monomial is attained

A generator with valuation exactly n is obtained by choosing

    0<=b<=d, 0<=a<=d-b,
    0<=e<=floor((n-b)/2), c=n-b-2e.

Its leading coefficient is exactly x^(a+2b+e). For each fixed b these
exponents fill the integer interval

    [2b, d+b+floor((n-b)/2)].                             (4)

All these intervals are nonempty because n>=d>=b. Consecutive intervals
overlap or touch: the upper endpoint for b-1 is

    d+b-1+floor((n-b+1)/2) >= 2b-1.

Their upper endpoints are nondecreasing. Their union is therefore

    [0,2d+floor((n-d)/2)]
       = [0,floor((n+3d)/2)].

For each exponent choose the corresponding single generator; its lower
rows vanish and its coefficient at n is exactly that monomial. Summing
these witnesses realizes every polynomial in (1). QED.

**Mechanism.** This is a two-point vanishing budget, not a matrix-rank
extrapolation. The retained factor (1+u) records the second point. Dropping
it would leave only deg F<=d+h and the much weaker bound r<=n+d.

## 2. All-row terminal tangent space

Return to a characteristic-zero field and the inherited normalized curve

    A0=1+x^2/4, C0=-3x/4-x^3/8.

The homogeneous top-row tangent at row n is

    (t^n theta A0', t^n theta C0'), degree(theta)<=n.

It preserves projected depth precisely when

    theta A0' in I_(2,n), theta C0' in I_(3,n).

For every n>=3, (1) makes these conditions exactly

    deg(theta)+1 <= floor((n+6)/2),
    deg(theta)+2 <= floor((n+9)/2).

For a nonzero theta the leading terms cannot cancel: A0'=x/2 and
C0'=-3/4-3x^2/8 have nonzero leading coefficients. The first inequality
is the stricter one when n is odd; they coincide when n is even. Thus

    ker D_n = K[x]_(<= floor(n/2)+2),                      (5)
    dim ker D_n = floor(n/2)+3,
    rank D_n = ceil(n/2)-2.

The cap degree(theta)<=n is automatic at this bound for n>=3.
These formulas prove the row15/16 computations uniformly. They also
exhibit the first failure of the proposed linear continuation n-6:

    n=17: kernel degree<=10, dimension11, rank7;
          theta=x^11 is forbidden, though 11=17-6.

This hostile needs only the A-side: theta A0'=x^12/2, whereas
I_(2,17) has maximum degree floor(23/2)=11. Hence the failure is proved
without a large-rank computation; depth_semigroup.py independently
checks the entire row17 selector as a finite control.

As before, (5) describes the direction space of a *nonempty* affine
terminal fiber. It does not establish compatibility of any given
inhomogeneous depth debt at any row.

## 3. The bracket operator remains universal

The inherited early prefix is

    A1=4/3+2x^2, C1=-4x-3x^3/2, q=-(x^2+6)/2.

For m>=3, varying the row m-1 tangent by theta changes the quantities
in the source-normal bracket identity G_m=T_m-q B_m/m by

    delta B_m=m(-x+3x^3/4)theta+theta',
    delta T_m=(2x-7x^3/4-3x^5/8)theta.

These follow by isolating the row m-1 factors in the determinant and
cubic expressions. There are no terms quadratic in that top-row
variation, since 2(m-1)>m for m>=3. Subtraction gives

    delta G_m=-x theta+(x^2+6)theta'/(2m)=L_m(theta).       (6)

Thus no inhomogeneous higher-row data enters this homogeneous map.

## 4. All-row repair quotient and unconditional homogeneous weight staircase

For m>=5, put n=m-1 and

    d_m=floor((m-1)/2)+2.

By (5), this is the precise terminal tangent degree bound for every
compatible source-normal prefix with the prescribed early rows.
Apply the exact Student sequence in PROOF.md §3 with d=d_m and N=m.
The complete row-m bracket quotient has dimension

    m-d_m=floor(m/2)-1.                                  (7)

It consists of one truncated Student moment and

    floor(m/2)-2

rigid high-degree coefficients, at degrees d_m+2,...,m. The full
prefix-preserving valuation-m source packet consists of

    p^(m-2b)y^b, 0<=b<=floor(m/2).

Its leading terms lie below degree d_m+2, so its quotient image is
exactly rank one. The odd b are invisible and every even b is visible.
Consequently the **least visible source weight at every m>=5** is

    w_min(m)=2m-2 floor(m/4).                             (8)

For every k>=2 the successive rows 4k-1 and 4k both have weight6k;
their least-weight visible payers are p^3 y^(2k-2) and y^(2k).
The repeated weight24 of rows15/16 is the k=4 instance.

The word “unconditional” here applies to the homogeneous depth and
response formulas, not existence of prefixes or lifts. Equations
(7)–(8) do **not** say any actual source satisfies the rigid debts,
that later depth debts vanish, or that infinitely many successive
scalar payments terminate in a polynomial Hamiltonian. Indeed the
number of untouched bracket debts grows roughly as half the row.
This gives an all-height obstruction budget for a repair strategy,
not a proof of a large conjecture.

## 5. Reproducibility and source boundaries

The definitions and normalization are the source-pinned ones in
SOURCES.md and PROOF.md. No new external literature is invoked.

    python3 05-knowledge/results/session254_four_frontiers_20261111/jacobian/depth_semigroup.py

The script checks explicit monomial attainment for every d=0,...,8 and
n=d,...,40, verifies the universal operator identities symbolically,
and independently builds the complete row17 depth matrices and selector.
The all-parameter statements follow from the proof above, not the finite
controls. No inherited inhomogeneous source reconstruction is performed.
