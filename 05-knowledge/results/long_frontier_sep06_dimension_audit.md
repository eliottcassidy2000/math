# Independent audit of the sharp finite-length stability penalty

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT — PASS.**

Auditor: `three_ray_geometry`, separate from the proof/source author,
`certificate_audit`. I accept the full
[dimension asymptotic](long_frontier_sep06_dimension.md), including its
first-order equality conditions and the added first-moment consequence.
The inputs are **THM-4454**,
[sharp global signed-root duplication stability](../../01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md),
and the independently audited
[minimizer classification](long_frontier_sep06_minimizers.md).
No reserved namespace is used as a proved dependency.

For the actual normalized finite-list domain p1=p2=1, E>0, and total
length at most N, the theorem proves

    I_N=K3+C/N+o(1/N),
    C=-22-(16/3)sqrt(6)+10sqrt(2)+(40/3)sqrt(3)>0.

It does not claim finite-N attainment or identify an exact finite-N
optimizer. Zero padding, arbitrary tail signs and varying lengths are
retained throughout.

## 1. Ordered local geometry and uniformity

Let a>=b>=c be the three leading positive entries, ell their mean,
m the remaining square mass, and h=ell-c. Writing the deviations as
alpha,beta,-h gives alpha+beta=h and alpha>=beta>=-h. Therefore
alpha lies in [h/2,2h], and direct maximization gives

    v=alpha^2+beta^2+h^2<=6h^2.

This justifies the crucial conversion of every symmetric splitting
effect into O(h^2), while the ordered top-two sum is exactly 2ell+h.

I independently expanded the three leading moments, using
3ell^2+v=1-m. The tail estimates

    |sum_tail r_i^3|<=m^(3/2), sum_tail r_i^4<=m^2

are uniform in the number and signs of tail entries. They give the
displayed moment and distance expansions. Their energy and distance
denominators stay away from zero in a fixed neighborhood of the
three-atom limit. In particular,

    J=(3-4z)+(10z-4)m+O(m^(3/2)+h^2),
    R-K3=A m+B h+O(m^(3/2)+h^2).

The exact coefficients agree with both the displayed radicals and

    A=(10z-4-K3*sqrt(2)*z)/d0,
    B=K3*sqrt(2)/d0.

Both are strictly positive. Products mh and m^2 are absorbed by the
stated error uniformly. Shrinking the neighborhood then gives the
coercive lower bound (A/2)m+(B/2)h. No dimension-dependent remainder
constant has been hidden in that step.

The positive B term is linear in the ordered splitting h. It cannot
be removed by symmetry of the moments, since the denominator selects
only two of the three leading entries.

## 2. Exact upper family and splitting controls

I independently solved the normalization equation with n equal negative
tails. Its positive solution is exactly

    ell=[1+sqrt(n((n+2)-(n+3)v)/3)]/(n+3).

For v=0 it lies in (1/3,1/sqrt(3)) for every n>=1. The actual list
has three positive entries and n nonzero negative entries, positive
energy and total length n+3. Both moment normalizations hold exactly.

The explicit family is analytic in 1/n near zero. Its tail cubic and
quartic contributions are O(n^-2) and O(n^-3), improving the generic
local remainder sufficiently to prove

    R=K3+A*(sqrt(3)-1)^2/n+O(n^-2).

The identity C=A*(sqrt(3)-1)^2 is correct. Replacing n by total
length N=n+3 changes only the O(N^-2) term.

For each fixed lambda>0, both displayed ordered deviation patterns with
h=lambda/n are actual eligible lists for all sufficiently large n.
Their variances are respectively 3h^2/2 and 6h^2, so their common
first-order splitting penalty is B*lambda. Thus the exact controls
really expose the linear ordered effect and do not supply evidence
for a false quadratic-only approximation.

## 3. Infimum lower bound and equality conditions

The approximate-minimizer construction is legitimate without compact
attainment: the explicit upper family makes I_N finite for large N,
and THM-4454 bounds it below by K3. Choosing errors N^-2 makes every
such sequence tend to K3, so the classification applies.

The uniform local coercivity gives m+h=O(1/N). Multiplying the
Taylor error by N then gives O(N^-1/2+N^-1)=o(1). Padding to total
length N and applying signed Cauchy--Schwarz to all N-3 tails gives

    (N-3)m>=(3ell-1)^2.

This remains valid with positive or zero tails. Since ell tends to z,
the lower bound on Nm combines with A,B>0 to produce the sharp
coefficient C along every subsequence of approximate minimizers.

I also checked both directions of the first-order equality statement:

    N(R-K3)->C
    iff Nm->(sqrt(3)-1)^2 and Nh->0.

For necessity, Cauchy gives the asymptotic lower bound on Nm, and
the nonnegative B Nh term cannot compensate for excess radial mass.
For sufficiency, the ordered variance bound and the norm identity
force ell to z, bringing the list into the fixed Taylor neighborhood.

The exact Cauchy defect is m-S^2/(N-3), S=1-3ell. Under these
conditions it is o(1/N). Every zero or nonnegative padded tail then
contributes at least S^2/(N-3)^2 to this defect, so there are only
o(N) such coordinates. Both the used length and the negative-coordinate
count are therefore N-o(N).

The added first-moment corollary is valid:

    sum_tail |r_i-S/(N-3)|<=sqrt((N-3)*defect)=o(1).

The reference tail is eventually negative. Hence the positive tail
sum tends to zero, and its negative magnitude tends to sqrt(3)-1.
This conclusion concerns the sharp first finite-length order; it
does not contradict the divergent mixed-dust examples for plain
R->K3 minimization.

## 4. Source audit and independent replay

I read the complete
[source](../../04-computation/long_frontier_sep06_dimension.py).
The formal quotient identities, direct atom-gradient calculation,
independent moment expansion, exact normalization formula and symbolic
upper/split derivatives agree. The atom-gradient calculation uses the
explicit moment expression as a local extension; it does not declare
the positive-only limiting row to satisfy p1=1.

The rational interval code rounds every operation outward. Its square-root
brackets are checked by squaring, and division rejects intervals containing
zero. The twelve actual controls use n=10,100,1000,10000 with exactly
the three declared deviation patterns. They retain both moment laws,
all multiplicities and positive denominators. Their finite comparisons
are correctly separated from the analytic asymptotic proof.

Independent normal and optimized runs reproduce the frozen output
byte-for-byte. All **117 always-active exact gates PASS**. Independently
computed file hashes are:

```text
source SHA256 badb2489aeb8e4a27a0e5b61c203323dfb988be693488637056a1f9f3c19e667
output SHA256 d2e72609bb6b248bed6d721e75e0a2b0350316b046cd6bed5ad06ddae8973e30
semantic SHA256 15251e617e95cae1f478f4ec83c7ce9b884fde71c9521555a5c78f4c65f21c88
```

No mathematical correction, quantifier defect or source defect was found.
