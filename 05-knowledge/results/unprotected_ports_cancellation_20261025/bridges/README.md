# Periodic integration retains one scalar; deeper cancellation can retain more

2026-10-25. **PROVED:** a complete cyclic parametrization of the physical
three-direction common kernel; an explicit wraparound dimension formula
for physical primes p>=19; and the characteristic-dependent rank and Smith
form of powers of the cyclic difference. **FINITE-EXACT:** independent
physical incidence and matrix controls. These results do not identify a
kernel vector with a binary error pattern, a detour, or a section.

Concurrent provenance: the arriving Session62 packet gives a stronger
arbitrary-support theorem and a wider interval range. The overlap and the
distinct difference-depth/Hasse refinements here are recorded in
[INCOMING.md](../INCOMING.md); this note does not supersede that stronger
support classification.

## 1. Inherited mechanism and the two different primes

The source is the [cube-plus-bridge common-kernel theorem](../../trijection_hourglass_ports_20261025/incidence/common_kernel.md).
It proves the no-wrap formula and leaves cyclic closure open. The exact
objects are points in the *physical* plane F_p^2, where p>3 is prime:

    T_L={(-a-3b,j-9b): a,b in {0,1}, 0<=j<L},
    2<=L<=p.

The three directions are u=(1,0), v=(0,1), w=(3,9). Arrays take values in
an arbitrary coefficient field k, whose characteristic need not equal p.
Let J be the space of arrays with zero sums on every actual physical line
in those three directions. The four columns have distinct x-coordinates
for p>3, and their L points are distinct for L<=p.

The canonical hostile is a periodic integration constant discarded by the
finite-support proof. The corrected near miss is averaging by 1/p when
the coefficient field has characteristic p. The least-used sidecar is one
evaluation of the integrated quantity, rather than its total sum. The
concept board is physical incidence / periodic difference / support holes /
invariant scalar / nilpotent depth / exact continuation observables.

## 2. Complete periodic parametrization for every physical p>3

All subscripts in this section are in Z/p. Extend each column by zero
outside I={0,...,L-1}. As in the source, write

    A=f_(0,0), B=f_(0,1), C=f_(1,0), D=f_(1,1),
    h=A, P=A+B, R=A+C, S=B+D.

Let (t f)_j=f_(j-1), and K=t^3+t^6+t^9. The physical u-line equations
give S=-t^9 R. Put r=R. The w-line equations give

    (1-t^3)P = -t^3(1-t^9)r,
    (1-t^3)(P+Kr)=0.

Because addition by three is one cycle on Z/p, the kernel of 1-t^3 is
the constant sequences over every coefficient field. Thus there is a
unique scalar lambda such that

    P=-Kr+lambda*1,        R=r,        S=-t^9 r.       (W1)

The complete parameter conditions are

    supp h subset I,             sum h=0;
    supp r subset E=I intersect (I-9),    sum r=0;
    (Kr)_j=lambda for every j outside I;
    p*lambda=0 in k.                                  (W2)

Here I-9 means {j:j+9 belongs to I}, with cyclic arithmetic. Recover the
array by

    A=h, B=P-h, C=r-h, D=S-P+h.                       (W3)

**Theorem W.** Conditions (W2) and formulas (W1),(W3) give a linear
bijection onto J for every p>3 and 2<=L<=p, in every coefficient field.

**Proof.** The two displayed physical line equations hold by direct line
labels, including their modular aliases. The support conditions say that
R,S,P and h vanish outside I, hence so do all four columns. Conversely,
the support of R and S forces exactly supp r subset E, and support of P
is exactly the outside-I condition in (W2). The v-line equations are the
four column sums. Since sum P=-3 sum r+p lambda and sum S=-sum r, they
are equivalent to sum h=0, sum r=0 and p lambda=0. This calculation uses
no division by three or p. Finally h=A, r=A+C, and lambda is recovered
from P+Kr at any one coordinate, proving injectivity. QED.

This is also a complete dimension algorithm for the small physical primes
5,7,11,13,17: take the rank of the explicitly given parameter equations
(W2). The simplified table below is not asserted for those primes, where
the support holes can overlap differently.

## 3. Closed wraparound formula for p>=19

The source already proves, without a characteristic restriction,

    dim J=L-1+max(0,L-10),        2<=L<=p-9.           (W4)

For the remaining lengths put d=p-L, so 0<=d<=8. Then

    dim J=2L-11+c_d+epsilon_d,
    c_d=max(0,9-2d,6-d),
    epsilon_d=1 if char(k)=p and d<=6, and 0 otherwise. (W5)

Thus the values of c_d, in order d=0,...,8, are

    9,7,5,3,2,1,0,0,0.

Together (W4),(W5) cover every 2<=L<=p for p>=19. At full length this
gives dim J=2p-2+[char(k)=p]. The additional characteristic-p mode first
appears at L=p-6, not at the first modular alias.

### Structural proof of the boundary table

Write the omitted positions as H={-d,...,-1} using representatives near
the end of the period. The r-support excludes H and H-9. These two holes
are disjoint in the stated range. The remaining coordinates divide into

* a free interval of length f=p-d-9=L-9;
* a local interval {-9,...,-d-1} of length g=9-d.

All outside-I equations reduce to

    r_(j-3)+r_(j-6)=lambda,        -d<=j<=-1,          (W6)

because r_(j-9)=0. In this equation coordinates outside the local interval
are zero. There are no other local constraints. For d<=3 these are d
independent, disjoint two-variable equations. For d=4 they reduce to

    r_-7=lambda, r_-6+r_-9=lambda, r_-5+r_-8=lambda.

For d=5 they reduce to

    r_-8=lambda, r_-7=lambda, r_-6+r_-9=lambda.

For d=6 they reduce to r_-9=r_-8=r_-7=lambda. At d=7 or 8 an equation
has both left-hand coordinates zero, forcing lambda=0, and the remaining
local coordinates are then all zero. These reductions use coefficient
one and repeated equations only, so they remain valid in characteristics
two and three.

For lambda=0 the local free dimension is exactly c_d. When d<=6, every
lambda admits a local solution, adding one parameter precisely when
p lambda=0 does not already force lambda=0. The free interval is nonempty
for p>=19 and d<=8, so sum r=0 removes exactly one independent free
coordinate, without restricting these local choices. Its contribution is
f-1. Adding L-1 h parameters gives (W5).

### An explicit characteristic-p hostile

Take p=19,L=13,d=6, h=0, lambda=1, and

    r_0=-3, r_10=r_11=r_12=1, all other r_j=0.

This has sum r=0 and satisfies every support and outside-I condition.
The columns obtained from (W1),(W3) have integer sums

    0,19,0,-19.

Every cyclic u- and w-line sum is already zero over the integers. The
remaining v-line sums vanish after reduction modulo 19, producing a
common-kernel mode absent in other characteristics. Setting lambda=0
would discard it. The checker verifies that a lambda mode exists exactly
in the cases specified by (W5).

The inherited three-space weighted-coupling lemma gives dim ker M=2 dim J
when all three coupling weights are nonzero. For the physical weights
(-3,-9,1), this applies in the physical coefficient field F_p and more
generally in characteristic different from three. The J formula remains
valid in characteristic three, but those particular coupling weights
vanish and the factor-two conclusion must not be imported there.

## 4. One scalar repairs first differences, not arbitrary depth

The quantity Q=P+Kr is translation invariant. Its value at one chosen
coordinate is exactly lambda, in every characteristic. Thus

    Q -> ((1-t^3)Q, Q_0)

is injective. This recovers the missing coordinate without an average or
division by p. The scalar is a k-valued coordinate, not a Boolean toggle.

Now put Delta=1-t^3 on all p cyclic coordinates. For every m>=1,

    dim ker Delta^m = 1              if char(k) != p;
    dim ker Delta^m = min(m,p)       if char(k) = p.    (W7)

**Proof.** The cyclic module is k[t]/(t^p-1). In characteristic different
from p, t^p-1 is squarefree; its only common root with 1-t^3 is 1, because
gcd(p,3)=1. Raising the latter polynomial to a power therefore retains a
one-dimensional kernel. In characteristic p, t^p-1=(t-1)^p and
1-t^3=(1-t)(1+t+t^2), whose second factor is a unit in this local ring
because its value at t=1 is 3!=0. Powers of Delta have exactly the same
kernel dimensions as powers of t-1. QED.

Consequently one scalar is insufficient for depth m>=2 in characteristic
p. For m<=p, values at m successive positions along the step-three cycle
are sufficient: the order-m finite-difference recurrence propagates them
uniquely. At least m independent linear sidecar coordinates are needed,
because the kernel has dimension m. At m=p the difference is identically
zero, so all p coordinates are hidden. In other characteristics, one
coordinate remains sufficient at every depth.

The integer operator makes this boundary particularly explicit. For
1<=m<=p, the nonzero Smith invariant factors of Delta^m are

    1 repeated p-m times, p repeated m-1 times,       (W8)

and there is one zero factor. To prove this, the rational kernel and left
kernel are both the primitive constant vector. All maximal minors
therefore have the same absolute value. The product of nonzero complex
eigenvalues of Delta is p, so for Delta^m it is p^m. Its adjugate is a
constant matrix and its trace is that eigenvalue product, giving each
maximal minor absolute value p^(m-1). This is the product of its nonzero
Smith factors. By (W7), reduction modulo p has rank p-m; hence exactly
m-1 of those factors are divisible by p. Their product is p^(m-1), forcing
each to equal p and every other factor to equal one. This includes m=p;
it is not a claim about the different higher-depth Smith pattern m>p.

## 5. Exact connections and boundaries

There is an exact older bridge to the July 24 theorem
[THM-2201, cyclic root-fibre Hasse jet](../../../../01-canon/theorems/THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier.md).
Encode a cyclic vector Q as Phi(Q)=sum Q_j*g^j. In characteristic p,
put epsilon=g-1. The constant vector Q=lambda*1 becomes

    Phi(Q)=lambda*(1+g+...+g^(p-1))
          =lambda*epsilon^(p-1).

This is exactly the full-fibre mode isolated in that theorem at p=13.
Its augmentation is zero, but its highest Hasse coordinate is lambda.
The Pascal coefficients J_i(Q)=sum_j Q_j*binom(j,i) retain it. More
generally ker Delta^m consists of the top m Hasse coordinates when m<=p;
these provide an alternative to the m evaluation coordinates in section 4.
Thus the map is an explicit change of basis preserving cyclic translation
and its cancellation kernel. It forgets none of those linear data when the
full jet is retained. The rooted owners, admissible event word and physical
geometry from the older application remain separate sidecars; the map
does not transport its dynamical conclusions to this strip.

| Source | Map | Preserved predicate | Lost coordinate and repair |
|---|---|---|---|
| Physical three-direction array | (h,r,lambda) from (W1)--(W3) | Every actual line sum and support condition | Dropping lambda loses a periodic mode; retain one evaluation |
| Cyclic sequence | Delta^m applied over the specified coefficient field | Its observed m-th difference | Retain a basis of ker Delta^m; one scalar only where (W7) permits |
| Integral cyclic operator | Smith factors and reduction modulo a prime | Exact rank change and torsion | Rational rank alone loses the p-primary factors |

The adjacent [finite-horizon clone-moment transfer](../counting/README.md)
has a different observable: aggregate future counting contributions. Its
rank is min(S,h) over Q, and an invisible signed bank at horizon two can
be visible at horizon three. The common principle is to measure the
kernel of the *actual observable* before deleting coordinates. Its moment
rank is not the dimension in (W7), and lambda is not its clone-cut toggle.

The [clone-cut involution](../../ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md)
preserves the inflated tournament while reversing a sign. That is an
equivariant pairing of combinatorial objects, not cyclic differentiation.
No assertion that one pairing certifies arbitrarily deep cancellation is
obtained from (W1). For involution differences themselves, the separate
[THM-2206, quadratic-deck characteristic audit](../../../../01-canon/theorems/THM-2206-quadratic-deck-augmentation-and-hamiltonian-characteristic-incompatibility.md)
gives (sigma-1)^m=(-2)^(m-1)(sigma-1), explaining why the order-two and
order-p cancellation operators must be distinguished.

For ports, the relevant transferable rule remains the exact
[shared-source relation](../../ports_regular_actions_hn_20261025/hn/RELATIONAL_LIFT.md):
an annihilated local statistic does not discharge an unrecorded boundary
condition. The physical line-sum equations are over k, whereas literal
tournament join and anti-reflection constraints are over F_2. A relation
between these carriers would need to preserve their actual source and
failure predicates. No such identification is asserted here.

## 6. Reproduction and scope

Run `python check.py > normal.json` and `python -O check.py > optimized.json`.
The script builds actual line partitions by physical cross products,
independently of the periodic parameter builder. Across 384 cases it
checks complete span, independence, and every line sum. The universe is
every 2<=L<=p for p=5,7,11,13,17,19,23, plus all nine wrap lengths for
p=29; coefficient primes are the distinct members of {2,3,5,p}.

It also checks 150 difference-power ranks, including m=p,p+1,p+2, and 23
integer cofactor determinants through m=p at p=5,7,11. Gates survive
optimized Python. These controls accompany the proofs above; they do not
extend the simplified support table to unproved small-prime cases.
