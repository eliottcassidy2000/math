# A periodic grammar for the first defect not determined by four moments

**PROVED.** This classifies the sharp first **rational moment-interval ambiguity**
for every initial odd backbone with top b>=27. It does not assert that rounding
the interval fails at that same order, or that A000568 becomes unknown there.

Let B_b={3,5,...,b}, b odd, and a=b+2>=29. After the exact backbone is removed,
residual cycle lengths are odd numbers >=a. For a residual core let R be the
number of positive defects remaining after subtracting the minimum positive
defect, as in Session279. A degree-four scalar moment certificate is rationally
exact when every core has R<=4; an R>=5 core has a positive ambiguity by the
inherited finite-difference perturbation.

## The complete finite classification rule

Let P consist of all 53 tuples

    0<=x1<x2<x3<x4,  all xi even,  x1+x2+x3+x4<=28.

For x in P define its six reduced defects

    d_ij(a,x)=xi-gcd(a+xi,xj-xi), i<j.

The actual positive defects are a+d_ij. Let t(a) be the least sum xi among
tuples for which the minimum of these six numbers occurs exactly once. Then

    first rational-ambiguity order N_b = 4a+t(a).       (1)

In particular t(a) is always defined, belongs to
`{12,16,18,20,22,24,26,28}`, and is periodic on the odd integers with period6930.
This rule is an all-backbone theorem, not fitting a numerical table.

## Proof of completeness

Choose q in{0,2,4} so that 3 divides a+q. The four lengths

    a+q, a+q+2, a+q+4, a+q+6

have defects `(a+q-1,a+q-1,a+q-3,a+q+1,a+q+1,a+q+3)`.
Thus their minimum is unique and they give R=5 at total mass
4a+4q+12<=4a+28. Therefore a first obstruction exists by that order.

No five-cycle residual core can occur by then, since 5a>4a+28 for a>=29.
With at most three cycles there are at most three edges. With four cycles
and a repeated length there are at most five positive edges; extracting their
minimum deletes at least one, leaving R<=4. Consequently the first obstruction
has four **distinct** lengths. Subtracting a from each gives exactly a tuple
in P. Their maximum pair difference is22<a, so no length divides another;
all six defects are positive. R=5 holds exactly when the minimum is unique.
This proves (1), including its completeness and the initial cutoff a>=29.

At the first core all five nonzero shifted Bernoulli parameters lie strictly
between0 and1. Its active-count law has positive mass at every node0,...,5.
The fifth finite-difference signed perturbation changes its zero atom while
preserving its first four factorial moments. Holding every other core fixed, this perturbation remains feasible in the
aggregated R=5 moment cell. Its positive weight prevents cancellation. Therefore
the lower and upper rational moment bounds are strictly distinct. At every smaller n all R<=4,
so inversion is exact. Backbone and singleton padding do not alter R.

## Period and exact finite certificate

Every gcd in the rule has second argument one of2,4,...,22. Because a+xi is
odd, only the odd part of that difference matters. Their least common multiple
is3^2*5*7*11=3465. Shifting a by6930 preserves oddness and every gcd, hence
preserves the full defect ordering and t(a). The algorithm evaluates all3465
odd residue classes in one period and all53 patterns, and checks the shifted
period as a control. The classification itself follows from the symbolic gcd
rule; the histogram merely summarizes it:

| t | number of odd residue classes mod6930 |
|---:|---:|
|12|1155|
|16|462|
|18|462|
|20|792|
|22|99|
|24|99|
|26|132|
|28|264|

Examples: B through31 first fails at144 with core(33,35,37,39);
B through35 is exact through167 and first fails at168 with(39,41,43,45).
Thus order150 can be computed by the degree-four interval alone with B through35.
B through49 is exact through215 (a=51 gives t=12, N=216).

## What the classification says about growth

The next obstruction is determined by finitely many small difference divisors,
not by the absolute size of the cycle lengths. The family translates with a
periodic arithmetic correction. This is an example of a genuine finite grammar
for an infinite set of certificate failure types. It classifies neither all
tournaments nor all symmetry groups; it identifies exactly the information
required to keep this one counting certificate rationally exact as size grows.

## The exact adaptive schedule

For n>=0 let b_*(n) be the least odd b>=27 for which n<N_b. This is a
certificate-optimal initial backbone **within this b>=27, degree-four family**,
not a globally optimal counting algorithm. Since20<=N_b-4b<=36, start at the
least odd b>=27 strictly greater than(n-36)/4, then test successive odd b.
The first successful value is b_*(n). By b>(n-20)/4 success is guaranteed, so
at most three successive odd candidates are inspected: the starting b is
strictly above(n-36)/4, so b+4 is strictly above(n-20)/4. The monotonicity of N_b follows
also directly: enlarging the backbone removes permitted residual cores.

Moreover N_(b+6930)=N_b+27720. It follows that

    b_*(n+27720)=b_*(n)+6930, for every n>=144.

Indeed N_27=144; for n>=144 the translated target excludes every b<=27+6930,
and all later candidates translate back into the same admissible regime.
This is a finite periodic controller for arbitrarily large exact-certification
backbones, not a claim that computing their moment packets takes constant time.
Examples: n150 selects b33 (first ambiguity156); n200 selects b45 (210);
n250 selects b57 (260). At n1,000,000 it selects b249995 (first ambiguity
1,000,008); this last statement concerns the proved schedule only, not a
computed millionth tournament count.
