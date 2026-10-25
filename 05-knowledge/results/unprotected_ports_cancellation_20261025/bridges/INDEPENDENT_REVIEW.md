# Independent proof audit of cyclic wrap, coefficient characteristic, and Hasse jets

2026-10-25. Reviewer: `nonconsecutive_ports`, independently of the bridge
producer. **VERIFIED: PASS**, with no mathematical gap found in the stated
scope of [README.md](README.md). The unbounded conclusions follow from the
proofs audited below; the finite checks are corroboration. No producer edit
was required.

## Physical equations and the retained scalar

For the actual point (-a-3b,j-9b), a horizontal line has label j-9b and
a w=(3,9) line has cross-product label -3(j+3a). Dividing the latter label
by -3 is performed in the physical field F_p, where p>3; this does not
divide by three in the coefficient field. The four x-columns are distinct
and their points remain distinct for L<=p.

These labels give the stated cyclic equations

    S=-t^9 r,
    (1-t^3)P=-t^3(1-t^9)r.

The identity `(1-t^3)(t^3+t^6+t^9)=t^3(1-t^9)` is exact. Since step three
is a single p-cycle, the kernel of 1-t^3 consists precisely of constant
sequences over any coefficient field. Hence P=-Kr+lambda*1 is complete,
and lambda is unique. Recovering lambda by evaluation uses no average
and remains valid in characteristic p.

The support of r and t^9 r is exactly `I intersect (I-9)`. Support of P
is exactly the outside-I equation Kr=lambda. The four column sums reduce
to sum h=0, sum r=0, and p*lambda=0, because sum P=-3 sum r+p*lambda.
Conversely those conditions make every recovered column sum zero. Thus
W1--W3 are a genuine bijection, with no omitted support or sum condition.

## The wrap table

For p>=19 and d=p-L in 0..8, the holes H and H-9 are disjoint. Their
complement consists of the free interval of length L-9 and the local
interval {-9,...,-d-1}. Every outside-I equation is precisely

    r_(j-3)+r_(j-6)=lambda, -d<=j<=-1.

I independently checked the reductions for each d. At d=0..3 there are
d disjoint equations on 9-d local variables. At d=4 the three equations
are r_-7=lambda, r_-6+r_-9=lambda, r_-5+r_-8=lambda. At d=5 they are
r_-8=lambda, r_-7=lambda, r_-6+r_-9=lambda. At d=6 the three local
variables all equal lambda. At d=7,8 a zero left-hand side first forces
lambda=0, after which all local variables vanish.

The homogeneous local dimensions are therefore exactly
`9,7,5,3,2,1,0,0,0`. All reductions use coefficients one and duplicate
equations, so characteristics two and three create no hidden rank loss.
For d<=6 every lambda has a local solution; it survives p*lambda=0
exactly in characteristic p. The free interval is nonempty throughout
this range, so sum r=0 removes one independent free coordinate without
constraining the local solution. Adding the independent L-1 h dimensions
gives W5. The prior no-wrap proof supplies W4, and the two ranges cover
all L. No extension of the simplified table to physical p<19 is used.

The explicit p=19,L=13 hostile also checks directly: its r has sum zero,
Kr=1 on the six omitted positions, and the recovered column sums over
the integers are 0,19,0,-19. Thus it is a real extra characteristic-19
mode rather than a failure of physical incidence.

## Coefficient characteristic and higher differences

The arbitrary coefficient field may have characteristic zero, two, three,
p, or a different prime. The physical labels remain in F_p throughout.
For the particular coupling weights (-3,-9,1), all three are nonzero
exactly when the coefficient characteristic is not three. The document
correctly restricts the factor-two weighted-coupling consequence while
retaining the J parametrization in characteristic three.

For W7, when the coefficient characteristic is not p, the polynomial
t^p-1 is squarefree and its gcd with (1-t^3)^m has degree one. This proves
the one-dimensional kernel without assuming the coefficient field is
algebraically closed. In characteristic p, the cyclic algebra is
k[epsilon]/(epsilon^p), and 1-t^3 is epsilon times a unit, because p>3.
Its m-th power has kernel dimension min(m,p). Successive evaluations
along the step-three cycle determine a solution of the order-m recurrence
from m initial values. The statement that Delta^p is zero is understood
in this characteristic-p case; in other characteristics its rank remains
p-1, as W7 states.

For W8, the rational left and right kernels are both the primitive constant
vector, so the adjugate is a constant multiple of the all-ones matrix.
The product of the nonzero eigenvalues is p^m; taking its adjugate trace
gives each cofactor absolute value p^(m-1). This is the product of the
nonzero Smith factors. Reduction modulo p has rank p-m for 1<=m<=p, so
exactly m-1 factors are divisible by p. Their product forces those factors
to equal p, and every other nonzero factor to equal one. The stated range
m<=p is essential and is retained; no higher-depth Smith pattern is claimed.

## The July bridge is an actual change of basis

I read [THM-2201, cyclic root-fibre Hasse jet](../../../../01-canon/theorems/THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier.md),
including its owner and dynamical limitations. The map Phi(Q)=sum Q_j*g^j
and the Pascal coefficients J_i=sum Q_j*binom(j,i) give exactly the same
invertible basis change, now written for general physical p in coefficient
characteristic p. Translation corresponds to multiplication by g, and
step-three translation to g^3.

The full constant vector becomes lambda*epsilon^(p-1). Its augmentation
is zero while its highest Hasse coordinate is lambda. More generally
ker Delta^m is the ideal generated by epsilon^(p-m) for m<=p, so precisely
the top m jet coordinates remain free on that kernel. This verifies the
claimed alternative sidecar to m evaluations. It is a faithful linear
carrier, not a transport of the older theorem's owner rules, event word,
or LRC dynamical conclusion. Those limitations are stated correctly.

## Reproduction and limits

Both normal and optimized executions of `check.py` were independently
rerun and passed: 384 physical incidence cases, 150 difference-power rank
cases, 23 integer cofactor determinants, and four full-fibre Hasse tests.
The checker constructs actual physical line partitions by cross products
separately from its periodic parameter builder, then tests dimension,
independence, and every line sum. Its gates use explicit exceptions.

The promoted objects are linear arrays, their actual incidence constraints,
and characteristic-specific kernels. This audit does not identify a kernel
vector with a binary error pattern, a section, a tournament detour, or a
clone-cut toggle. It establishes the exact bridge stated in the source,
including the information that must be retained when coordinates are removed.
