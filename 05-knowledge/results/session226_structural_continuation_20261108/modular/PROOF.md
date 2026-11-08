# Two-J valuation total positivity and sharp universal modular floor

Status: PROVED mathematical continuation; [independent audit PASS](../modular_audit/AUDIT.md).

## Scope and dependency

The sole mathematical dependency is the exact fixed-length two-J law and its
explicit kernels in `session225_joint_modular_law_recovery_20261108/PROOF.md`
and `KERNELS.md`. We use their letters A=K00, B=K01=K10, C=K11 and weights
(9,6,1)/16. For odd q and residue r, the limiting density is

L_q(r) = (9 product A_i(r) + 6 product B_i(r) + product C_i(r))/16.

This packet proves properties of that limiting probability law L_q(r)/q.
It gives no height-truncated law, integer coverage, cutoff, finite-bank
completeness, or Chan–Pak conjecture proof. The inherited O_q(1/L) convergence
is only for fixed q; none of our uniform density estimates upgrades it.

## 1. Strict monotone likelihood ratios at every odd prime power

Fix p>=3 and k>=1. Put alpha=1-1/p, d=1-1/p^2,
b=1+1/p-1/p^2, c=1+2/p. For valuation v<k put t=v+1.
The inherited formulas are

A(v)=c-p^(-t)(c+t d), B(v)=b-p^(-t)d;
A(k)=c-p^(-k-1)(2+k alpha), B(k)=b-alpha p^(-k-1).

C(0)=1+p^(-3); C(v)=1-p^(-2)+p^(-3) for 1<=v<=k.
All values are positive. The ratios B/C and A/B are strictly increasing in
v on {0,...,k}.

Proof. B strictly increases. C drops at the first step and is then constant,
so B/C strictly increases. At a nonterminal step v to v+1 (v<=k-2),

Delta B=p^(-t) alpha d,
Delta A/Delta B=t+1+p/(p^2-1)>2.

At the terminal step k-1 to k,

Delta B=p^(-k) alpha,
Delta A/Delta B=k+1/alpha>2.

Meanwhile A(v)<1+2/p<=5/3 and B(v)>=1-1/p^2+p^(-3)>8/9,
so A(v)/B(v)<15/8<2. Therefore
B Delta A-A Delta B>0 at every step, proving the assertion. In particular
A/C=(A/B)(B/C) also strictly increases.

The conditional probability of valuation v in branch H is
p_H(v)=n_v H(v)/p^k, where n_v=(p-1)p^(k-v-1) for v<k and n_k=1.
The common positive n_v factor cancels from every likelihood ratio.
Thus, in the latent order C<B<A, the valuation laws are strictly ordered
by monotone likelihood ratio, not merely by their zero-residue probabilities.

## 2. Strict multivariate total positivity

Let q=product p_i^k_i and V_i=min(v_p_i(F),k_i) in the limiting two-J law.
Write f(v) for the joint probability mass on product_i {0,...,k_i}.
For all vectors x,y,

f(x meet y) f(x join y) >= f(x) f(y).

The inequality is strict whenever x,y are incomparable in coordinatewise
order. Thus f is MTP2 (multivariate totally positive of order two), strictly
so on every nontrivial incomparable pair. For comparable pairs equality is
an identity. This statement is empty of strict content when q has only one
prime factor.

Elementary proof. Write f(v)=sum_H w_H product_i p_H,i(v_i), with
w_C=1/16,w_B=6/16,w_A=9/16. Expand both products in the claimed inequality.
Same-branch terms cancel exactly. Pair two distinct branches H<J. Set
r_i(s)=p_J,i(s)/p_H,i(s), a strictly increasing positive function. Factor out
the common positive quantity w_H w_J product_i p_H,i(x_i)p_H,i(y_i).
Let I={i:x_i>y_i}, Jset={i:x_i<y_i}, and E the remaining indices. Put
R_I^+=product_I r_i(max(x_i,y_i)), R_I^-=product_I r_i(min(x_i,y_i)),
and similarly R_Jset^+,R_Jset^-. The residual difference is

(product_E r_i(x_i)) (R_I^+-R_I^-)(R_Jset^+-R_Jset^-).

This is nonnegative and is strictly positive if both I and Jset are nonempty.
Summing over latent branch pairs proves the theorem without invoking a
marginalization theorem for MTP2.

## 3. Association, conditioning, and posterior reinforcement

Any two coordinatewise nondecreasing real functions of the valuation vector
have nonnegative covariance. This remains true after specifying arbitrary
values of any subset of coordinates. In particular, any two nonconstant
nondecreasing functions of distinct remaining coordinates have strictly
positive covariance. All conditioning events mentioned here have positive
probability.

Proof. Under each branch, coordinates are independent. Product measures are
associated: induct on the number of coordinates, using the covariance
decomposition and the one-coordinate identity

Cov(a(X),b(X))=1/2 E[(a(X)-a(X'))(b(X)-b(X'))]>=0.

Increasing conditional means remain increasing when integrating independent
coordinates. The likelihood-ratio order from section 1 implies stochastic
order, by the same identity applied to an increasing function and the
increasing likelihood ratio. Therefore, for increasing f,g, their conditional
means given H increase in the order C<B<A. Decomposing covariance over H
proves association. If f,g are nonconstant increasing functions of separate
coordinates, both branch means strictly increase; their covariance over the
three positive latent weights is strictly positive. Fixing other coordinates
only replaces w_H by positive weights proportional to w_H product p_H,i(v_i),
leaving the same proof intact.

Increasing any observed valuation multiplies each higher/lower branch odds
by a factor strictly greater than one. Consequently it stochastically raises
every unobserved valuation. This gives posterior reinforcement, not just the
previously known unconditional positive correlation of divisibility events.
No assertion is made after arbitrary event selection or height conditioning.

## 4. The exact residue minimizers

For every odd q>1, L_q(r) is strictly increasing in each truncated prime
valuation, holding all others fixed. Its minimum is attained exactly on
units modulo q and its maximum exactly at r=0 modulo q.

Proof. For all valuation increases beyond the first, A and B strictly increase
and C is constant. At the first step, write a,b,c for the products of A,B,C
over the other prime powers. Then

16 Delta L=9 a Delta A+6 b Delta B-c/p^2.

The inherited lower bound A>(1-p^(-2))^2 gives

a>(product_(odd primes)(1-p^(-2)))^2=(8/pi^2)^2>1/2.

Also C<=1+p^(-3), whence
c<=product_(odd primes)(1+p^(-3))
 <= product_(n>=2)(1+n^(-3)) <= exp(sum_(n>=2)n^(-3)) < exp(1/4)<2.

Here sum_(n>=2)n^(-3)<1/8+integral_2^infinity x^(-3)dx=1/4.
For k=1 the first Delta A=(2p-1)/p^2>1/p^2. For k>=2,
p Delta A=2-1/p-3/p^2+2/p^3>1, so again Delta A>1/p^2.
Therefore 16 Delta L>(9a-c)/p^2>0. This proves strict coordinate increase.
(The elementary inequality pi^2<10 suffices for the displayed a bound.)

## 5. Sharp universal floor as an Euler-product constant

Define the absolutely convergent positive products over all odd primes

A_* = product_p (1-2/p^2+1/p^3),
B_* = product_p (1-1/p^2+1/p^3),
C_* = product_p (1+1/p^3),
and C_*sharp=(9 A_*+6 B_*+C_*)/16.

Then

inf_{q odd positive, r mod q} L_q(r)=C_*sharp,

with strict L_q(r)>C_*sharp for every finite q. This is a sharp constant,
not the coarser 36/pi^4+7/(2pi^2) bound inherited from Session225.

Proof. By section 4 it suffices to take units r. Their density depends only
on the prime support of q. If a,b,c are the three products for a finite
support, adjoining a new prime p changes 16 L by

p^(-3)[-9a(2p-1)-6b(p-1)+c],

which is negative, since a>1/2,c<2,p>=3,b>0. The finite-support values
therefore strictly decrease under inclusion. Absolute convergence follows
from the summability of p^(-2). Exhausting all odd primes proves the infimum
formula; adjoining even one omitted prime proves strictness at finite q.
The q=1 value is 1 and is also strictly above the infimum.

For numerical evaluation, let C(N) use the primes <=N. Then

0<C(N)-C_*sharp<3/(2N), for N>=3.

Indeed product_tail(1-a_p)>=1-sum_tail a_p, with a_p=2/p^2-1/p^3
for A and a_p=1/p^2-1/p^3 for B. Their deficits are bounded by 2/N
and 1/N using sum_(n>N)n^(-2)<=1/N. The positive C-tail only reduces
the total deficit. The coefficients give (9*2+6)/(16N)=3/(2N).
Floating point output is descriptive only; the rational exact tests and
analytic error bound do not depend on floating point calculations.

## Inheritance / connection / loss ledger

Closest proved mechanism: the shared three-branch CRT mixture and explicit
prime-power kernels. Canonical hostile: replacing the mixture of products
by a product of mixtures fails already at q=15. Corrected near miss: branch C
is not increasing in valuation, so coordinate monotonicity does not follow
termwise. Least-used sidecar: the ordering of all branch likelihood ratios,
rather than just their probabilities at zero.

Concept board: (1) shared CRT latent label, (2) prime valuations and likelihood
ratios, (3) total positivity / conditional dependence, (4) extremal residue
and Euler-product scale. The likelihood-ratio proof turns (1) into (3);
the quantitative first-step check is separately needed for (4).

Map: full normalized O/Z/J fixed-length model -> inherited limiting scalar
law -> truncated prime valuations. Preserved predicates: modular divisibility,
coordinate association, likelihood-ratio reinforcement, density extrema.
Lost information: unit classes, SCC labels, integer size, tree shape, word
history. Restoration sidecars for further transfers are unavailable here.
There is no bridge to height conditioning or coverage.

Reproduce with `python3 check.py`. Finite checks validate implementations and
hostile examples, not the all-modulus proof. No literature-priority claim.

### Certified numerical enclosure

`check.py` additionally evaluates each finite product using integer interval
arithmetic at scale 10^30, rounding the lower endpoint downward and upper
endpoint upward after every factor. Combining this with the proved tail
bound at N=10^6 certifies

0.771866451842 < C_*sharp < 0.771867951843.

These endpoints are rounded outward to twelve decimal places. This enclosure
is rigorous and separate from the descriptive floating point estimate.
