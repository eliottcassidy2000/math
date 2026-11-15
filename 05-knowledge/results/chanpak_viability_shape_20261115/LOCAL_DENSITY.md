# The unrestricted two-tape bank has trivial local factors and a divisor-law real factor

Use the physical fixed carrier from
[absolute steering](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md):
```
F(x,y)=B+Qxy+Tx+Sy+h[p(x)+p(y)],   x,y nonnegative integers,
(B,Q,T,S,h)=(2855,5760,4034,3914,45).
```
Here p is the exact even binary surcharge. Presentations below are
**ordered address pairs**, not distinct counts or isomorphism
classes. In the spirit of [Vaughan's account of Hardy–Littlewood
local and real factors](https://assets.cambridge.org/97805215/73474/excerpt/9780521573474_excerpt.pdf),
the result identifies both factors for this first moment; it is
not a pointwise circle-method asymptotic.

## Main theorem

For any fixed M and any odd residue a modulo 2M, the number of
presentations with F<=X and F≡a mod 2M is, uniformly in a,
```
 (1/M){(X/Q) log X + K X} + O_M(X^(1−epsilon_M)),
 epsilon_M>0,
 K = c(S)+c(T)−(1+log Q)/Q,
 c(b) = [log Gamma(b/Q)−log Gamma((b+h)/Q)]/h.
```
Without any residue restriction the remainder can be taken
O(sqrt(X) log X). In particular the fixed-modulus factor is 1 after
normalizing for odd counts, even at primes dividing 45 or Q.
The full real first moment has logarithmically growing mean
multiplicity across target integers. It does not imply positivity
at any specified large integer.

## Interval mixing for the exact checksum

Write c0=p/2. Modulo (2M,M) the two scalar digit transitions on
(u,c0) are u→2u+bit, with c0 unchanged for bit 1 and increased
by (u+(u mod 2))/2 for bit 0. They give a Markov operator for
uniform independent bits. Every state has exactly two incoming
preimages with multiplicity: the new u determines the bit, then
the doubling equation on u has two preimages, and c0 is forced.
Thus the uniform law is stationary. The
[fixed-length prefix theorem](../chanpak_absolute_steering_20261115/PREFIX_PROGRAMMING.md)
gives a positive common transition length for all state pairs.
A block of that length dominates a fixed positive fraction
of the uniform law, hence iteration contracts the deviation
geometrically. On any aligned interval of 2^k addresses the
relative discrepancy is O_M(2^(-eta_M k)), uniformly in the
high prefix. Decomposing an integer interval into at most
two dyadic blocks per scale gives absolute discrepancy
O_M((V+1)^(1−eta_M)) on intervals of length V (reduce eta_M
if necessary so 0<eta_M<=1/2).

## A vertical lattice estimate, uniform in the other address

Fix x, let A=Qx+S, Y=X−B−Tx−h p(x), and W=Y/A.
For large W the number of y with A y+h p(y)<=Y has the estimate
```
 (Y/h) log(1+h/A)+O(1+sqrt(Y)/A+log(W+2)/A).
```
The constant is uniform over the rows, since A>=S>0 is fixed
away from zero and h is fixed. To see it, partition y=qu+b at a
power of two q (at least 2). From the exact
binary recurrence,
p(y)=y(1−rho(y))+O(log(y+2)), where rho is binary reverse fraction.
In the class b mod q, replace rho(y) by rev(b)/q incurring error
O(1/q). Thus the smooth row cutoff uses slope
A_b=A+h(1−rev(b)/q), and the two membership tests can differ
only in a y-interval of width O((W/q+log(W+2))/A).
Both y ranges are bounded by W. Summing the boundary discrepancies
costs O(q+W/(Aq)+log(W+2)/A). The Riemann sum Y/q sum_b 1/A_b
has error O(W/(Aq)), since the derivative of the integrand is
O(1/A²). Take q comparable to max(2,sqrt(Y)/A).

For the row with the additional congruence, use the same cutoffs
but q comparable to sqrt(W+1) (and at least 2).
At each fixed b, the exact binary split expresses
c0(qu+b) as c0(u) plus a function of u mod 2M (and b).
Thus interval mixing applies to the actual congruence. Let
d=gcd(h,M). After division of F−B by 2, the linear coefficient
A/2=2880x+1957 is a unit modulo d. The function depending on b
must still be retained for congruences modulo M, but all terms
coming from h c0 vanish modulo d. Consequently exactly 1/d of
u residues modulo 2M can satisfy the condition modulo d,
because q is a unit mod the odd d. For each such u residue,
a proportion d/M of c0(u) gives the specified full modulus.
This gives exactly weight 1/M in each low-word class.
The mixing error summed over classes is O_M((W+1)^(1−eta_M/2)),
since each u interval has length at most W/q and there are q
classes. The same proof applies to horizontal rows using
T/2=2017, also a unit modulo d.

## Hyperbola decomposition and the second term

Let L=floor(sqrt(X/Q)). Every witness is in the vertical strip
x<=L or the horizontal strip y<=L because F>=Qxy. Their intersection
contains all but O(1) integer pairs of the square [0,L]^2.
Indeed, for a fixed sufficiently large constant b0,
if one address is <=L−b0, the loss Q b0 L in the product dominates
all other terms, which are at most B+(T+S+2h)L. Thus only a fixed
corner square can be absent. The full square contributes
X/Q+O(sqrt X), or X/(QM)+O_M(X^(1−eta_M/2)) for one compatible
residue by interval mixing (condition on one address).

In each vertical row x<=L, replacing Y by X in the main term
costs O(1), since B+Tx+h p(x)=O(x+1) and its logarithmic slope
factor is O(1/(x+1)). The sum of the unconditional row errors
is O(sqrt(X) log X); the sum of the modular row errors is
O_M(X^(1−eta_M/4)). Do the same for horizontal rows.
For b=S,T the finite main term is evaluated by the product
identity for Gamma:
```
 sum_{n=0}^L (1/h) log(1+h/(Qn+b))
 = (1/Q) log L + c(b) + O(1/L).
```
Equivalently subtract the harmonic sum termwise and use the
absolutely convergent constant; the Gamma identity follows
by repeated Gamma(z+1)=z Gamma(z) and the fixed-shift ratio
asymptotic. Subtract the square contribution, and substitute
log L=(log X−log Q)/2+O(X^(−1/2)). This proves the asserted
coefficient and error.

## How to use the result, and what it forbids claiming

Each prime-power local density in this presentation problem
is exactly the uniform compatible density in the limit; all
finite products agree with that normalization because the
same conclusion holds at every composite M. There is no
Euler-factor obstruction left to fix by repeating the same
residue arguments. The real first moment is much richer than
that of the nearly balanced windows: its logarithm comes
from the sum over size scales of the shorter address, not
from a single tight tangent neighborhood.

For example the main-term contribution of x in [A,2A] with
A tending to infinity and A<=L/2 is (log 2) X/Q up to lower
order, regardless of the scale A. Thus a meaningful major
arc needs to combine many genuinely different aspect ratios.
A pointwise coverage claim would still need to bound the
remaining Fourier contribution or construct a deterministic
representative. A large mean cannot rule out collisions or
sparse exceptional targets, and the physical brackets do not
supply an exact coefficient.

### The shape of the real factor

The ratio is not stationary on a bounded sector. For each fixed
0<theta<=1/2, the vertical row estimate through x<=X^theta gives
(theta/Q)X log X+O(X), also with factor 1/M in each fixed odd class.
For 1/2<theta<1, use the transpose and the fact that
log((x+1)(y+1))/log X tends in probability to 1: the number of
pairs with product below X^(1−epsilon) is O(X^(1−epsilon)log X),
and every pair has product O(X). Thus log(x+1)/log X converges
to the uniform law on [0,1] over the presentations, and
log((x+1)/(y+1))/log X to uniform on [−1,1].
This convergence also holds after a fixed residue condition.

For a fixed ratio bound 1/R<=x/y<=R with R>1 and x,y positive,
there are only (X/Q)log R+O_R(sqrt X) presentations through X.
The region xy<=X/Q in that angular sector has the asserted area
in coordinates u=xy, v=x/y, with Jacobian 1/(2v). The
boundary length is O_R(sqrt X), giving the stated lattice
error; the actual linear and checksum corrections change the
height by at most O_R(sqrt X). Consequently even R=10^6
permits fewer than 0.5% of the odd counts in this carrier.
The thin tangent band |x−y|<=C sqrt(x+y) has only O_C(X^(3/4))
ordered pairs, by summing O_C(sqrt s) choices at sums
s=x+y<=O_C(sqrt X). Its exact support has density zero despite
the earlier all-target bracketing result. These bounds are on
restricted presentations, not on the wider tournament spectrum.

There is even a crude explicit range where the full bank
cannot cover: for X>=Q, the standard divisor estimate and
the two axes give at most
`
(X/Q)(1+log(X/Q)) + X/S + X/T + 1
`
presentations. For X<=exp(2800) this is below 0.49X
(and in particular below the number of odd targets when X is
large enough). This is a very long finite shortage of the bank,
not proof of infinitely many holes. See the separate
[strong-word countermodel](MARKED_WORD_COUNTERMODEL.md) for an
infinite failure of the unweighted local/mean inference.
