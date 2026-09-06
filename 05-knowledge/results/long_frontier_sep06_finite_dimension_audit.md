# Independent audit of the finite-dimension stopping obstruction

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT replay.**
The [proof note](long_frontier_sep06_finite_dimension.md),
[source](../../04-computation/long_frontier_sep06_finite_dimension.py), and
[frozen output](long_frontier_sep06_finite_dimension.out) pass without a
mathematical correction. The proposed inequality `N(R-K3)>=C` for every
eligible length and every N is refuted. The sharp asymptotic theorem is
unchanged, and the corresponding question restricted to N>=6 remains open.

## 1. Actual domain and exact witnesses

The lists `(256,16,1,-16)/257` and `(675,30,1,-15,-15)/676` have respective
lengths four and five, no zero coordinates, numerator sum equal to the
denominator, and numerator square sum equal to its square. Thus p1=p2=1
exactly. Their first two displayed coordinates are the ordered two largest
positive entries. At least two coordinates are nonzero, so p4<1 and the
energy is positive. Their distance from the positive two-atom target is
also strictly positive.

An independent calculation formed the actual polynomial
`G(s)=product_i(1+r_i*s)` by Fraction convolution, then evaluated
`D=-[s^4]G(s)^2` and `E=sum_(i<j)r_i^2*r_j^2`. It recovered

```text
length four: D/E = 66305/65793,
length five: D/E = 444131/456301.
```

Substitution into the actual ordered-distance quotient, followed by exact
radical comparisons, independently verified the stated intervals

```text
21419/10000 < 4(R-K3) < 21420/10000 < C,
21649/10000 < 5(R-K3) < 21650/10000 < C.
```

The comparisons with the equal-three family also pass. These are actual
finite eligible lists, so they refute both exact-length and length-at-most
versions of the proposed all-N bound. No limiting zero-energy list is
substituted into an undefined quotient.

## 2. Why the rational curve produces the obstruction

For N>=4 and integer k>=1, set `T=(N-3)(N-2)/2`. The displayed curve's
first moment and second moment identities follow respectively from
cancellation of `(N-3)k` and from `(N-3)^2+(N-3)=2T`. All coordinates
are nonzero. The first two displayed entries are largest among the three
positive entries because

```text
T*k^2 / ((N-3)*k) = (N-2)*k/2 >= 1,
(N-3)*k >= 1.
```

For each fixed N it converges to one positive unit atom; the already-proved
one-atom boundary limit therefore gives `R -> L=sqrt(2)-2/3`. The exact
comparisons `5(L-K3)<C<6(L-K3)` were independently recomputed. In
particular the strict limiting defect is negative for N=4 and N=5, so
every sufficiently large integer k in either of these two curves is a
counterexample. This inference does not assert monotonicity of the curve.
The source's first-violating-parameter statements concern only the two
explicit finite integer heads, not all rational lists.

The equal-three formula is the positive solution of
`3a^2+(3a-1)^2/(N-3)=1`; its negative coordinates are
`-(3a-1)/(N-3)`. This verifies its normalization independently of the
interval checks. All bounded comparisons concern the correct positive
two-root distance. The N=4,...,7 comparisons refute finite-N optimality of
that particular family; the N=8,...,20 comparisons do not prove an optimizer
or a global transition.

## 3. Source and replay scope

The source uses exact rational interval arithmetic. The integer-square-root
construction encloses each positive radical, every reciprocal rejects an
interval containing zero, and each strict comparison separates interval
endpoints. Repeated interval multiplication may enlarge the interval but
cannot invalidate its enclosure. Both the exact moment formula for J and
the ordering of the two positive coordinates are correct. The checks are
explicit exceptions and remain active with Python optimization enabled.

The complete declared universe is k=1,...,16 at N=4, k=1,...,15 at N=5,
and precisely two named controls at each N=4,...,20. Independent normal
and optimized replays agree byte-for-byte with the saved output, passing
**445 always-active gates**. No additional census or numerical inference
was used in this audit.

```text
source SHA256 bc3ba159919527741e8b7439823759f13575916549ecdb666b9667f9fb5054f1
output SHA256 490e74d700b33019a726f75e78b3672aff13b1e3ecc4737cfe22f83583d9d024
semantic SHA256 1dc1ef2b4c92cf87562a81c1248f567b290475d797d15a5569e83f27b422cf6c
```

The obstruction concerns an extrapolation of the asymptotic coefficient
to all finite dimensions. It neither retracts the sharp asymptotic formula
nor settles the surviving N>=6 inequality or any exact finite optimizer.
