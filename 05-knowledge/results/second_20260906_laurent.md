# Endpoint 33: the actual carried doubled row is negative at every first root

**Status: PROVED in the stated scope; complete exact certificate and analytic
proof independently audited.** No theorem ID is
reserved. This closes the next six-first-channel family in the existing
`A=2,B=3` progression, with an unbounded opposite endpoint. It also records
two FINITE-EXACT failures of attempted sign certificates on actual normalized
rows. General trinomial two-rung separation remains **OPEN**.

## 1. Statement and exact scope

Let `g>=17` be an integer with `gcd(g,33)=1`, and let
`alpha,beta,gamma` be arbitrary nonzero complex numbers. For

```text
f(u)=alpha*u^(-33)+beta*u^(2g-33)+gamma*u^(3g-33),
```

the first nonzero constant-term moment occurs at exactly `g` or `2g`.
Both alternatives occur for every allowed g. Thus

```text
first nonzero mass <=2g<3g=(negative endpoint)+(positive endpoint).
```

There are six channels in the first support-return fiber and twelve in the
doubled fiber, including its lower carry. There are exactly five distinct
canceling coefficient phases at mass g, and every one has a nonzero moment
at mass `2g`. The sign below applies after the displayed coefficient
normalization; it is not a real sign assertion about a complex raw moment.

The first-mass statement uses the gcd condition. The polynomial certificate
itself does not use it. For example at `g=18`, the support `(-33,3,21)`
already has first support return six, so eighteen is not its first
admissible mass.
No claim about separation at every later pair of masses is made.

## 2. Inheritance, novelty and the retained board

The closest proved mechanisms are
[THM-4436, complete factorial-row roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md),
the [endpoint-27 carried certificate](overnight7_20260906_laurent_quartic_carry.md),
and the [all-h characteristic degree interface](overnight7_20260906_laurent_midpoint_transport.md),
Section 5. The new [minimal derivative cone](creative_20260906_laurent_bridge.md)
was the attempted entry mechanism. The positive final certificate uses the
actual multiplication operator in the first-root quotient, after that
entry failed.

The canonical hostile is now the
[combined-pencil free-scale obstruction](combined_pencil_empty_core_morning_sep06.md):
Euler coupling, all-weight compatibility and the same ordinary carrier do
not retain the crossing-edge normalization. The corrected near miss is
therefore substituting a flexible midpoint model for the actual carried
row. The least-used sidecar in this session is the **literal coefficient
of the lower carry**, retained before quotient reduction.

The board is: actual carried counts; alpha-completed square; minimal
coupled observers; positive phase multipliers; and quotient characteristic
coefficients. The actual cone failure below caused the transition to the
last object; the phase-multiplier scale probe prevented an unsupported
shortcut back. Searches for endpoint33, width33, the literal support and
six-first-channel closure in current theorem and result files found no
existing proof of this family. Endpoint15 and endpoint27 are explicitly
inherited, not counted again. No external-priority claim is made.

## 3. Exact first and doubled coefficient maps

Put

```text
x=g-16>=1,       tau=alpha*gamma^2/beta^3,
X=alpha^x beta^15 gamma,       K=(2g)_22>0.
```

A return channel of total mass g satisfies
`33*g=2g*n_beta+3g*n_gamma`, hence
`2*n_beta+3*n_gamma=33`. Its nonnegative solutions are exactly

```text
(n_alpha,n_beta,n_gamma)=(x+j,15-3j,1+2j),   j=0,...,5.
```

At mass `2g`, every channel is

```text
(2x+e,30-3e,2+2e),                         e=-1,...,10.    (1)
```

The lower channel e=-1 is valid because `2x-1>=1`.
These are complete fibers, not a selected semigroup subfamily.

Define the monic first polynomial and normalized doubled Laurent row by

```text
p_x(t)=sum_(j=0)^5 [11!*(x+5)_(5-j)/((15-3j)!(1+2j)!)] t^j,

q_x(t)=sum_(e=-1)^10 [(2x+10)_(10-e)/((30-3e)!(2+2e)!)] t^e. (2)
```

Here `(a)_r` denotes a falling factorial and `(a)_0=1`. Literal
multinomial expansion, including (1), gives

```text
CT(f^g)    = X * binom(g,11) * p_x(tau),
CT(f^(2g)) = X^2 * K * q_x(tau).                          (3)
```

The row in THM-4436 with `A=2,B=3,h=5,r=0,z=1` is exactly a
positive multiple of `p_x`. It has five distinct strictly negative roots
for every integer x>=1. In particular `p_x(0)>0`, so the inverse t in
q_x is valid in every first-root quotient. The five phases are attained
by taking `beta=gamma=1` and alpha equal to a selected root. Taking
all three coefficients equal to one gives the first-detection alternative g.

For an arbitrary positive return mass m, the charge equation modulo g
implies `g | 33m`. The gcd hypothesis therefore gives `g | m`. The
six actual channels above show g is feasible; no positive earlier mass is.
Off the five first roots the moment at g is nonzero. It remains to prove
`q_x(rho)<0` at each of them.

## 4. A polynomial multiplication operator with no hidden denominator

Initially work over `Q(x)[t]/(p_x(t))`, where t is invertible. The
remainder constructed below lies in `Q[x,t]`; inversion is not asserted
over the unspecialized ring `Q[x][t]/(p_x(t))`. Its constant
coefficient is

```text
p_0(x)=(11!/15!)*(x+5)_5.
```

The only prospective denominator in reducing q_x is its inverse carry.
It cancels as the exact polynomial identity

```text
q_(-1)(x)/p_0(x)
 = [64*15!/(33!*11!)] * x * product_(j=0)^4(2x+2j+1).       (4)
```

To verify (4), split `(2x+10)_11` into its six even and five odd
factors. The five factors `(x+1),...,(x+5)` cancel p_0; the remaining
even factor is `2x`.

Reduce the nonnegative powers of q_x by ordinary monic division, and use

```text
t^-1 = -sum_(j=1)^5 p_j(x)t^(j-1)/p_0(x)
```

for the carry. The result is

```text
R_x(t)=sum_(j=0)^4 R_j(x)t^j in Q[x,t],
q_x(t)=R_x(t) mod p_x(t),          deg_x R_j<=10-j.          (5)
```

The degree bound is structural. Give x and t weight one. The coefficient
of t^e in q_x for e>=0 has degree `10-e`; each monic relation coefficient
p_j has degree `5-j`. Monic reduction therefore preserves total weight
at most ten. For the inverse term, (4) has degree six and p_(j+1) has
degree `4-j`, giving the same bound in (5).

Let `T_x` be the matrix of multiplication by R_x in the basis
`1,t,...,t^4`, and write

```text
C_x(z)=det(zI-T_x)=z^5+c_1(x)z^4+...+c_5(x).                (6)
```

Its entry in row i, column j has degree at most `10+j-i`.
In a k-by-k principal determinant the row and column index sums cancel;
therefore `deg_x c_k<=10k`. This supplies the independent finite identity
audit bound, rather than inferring a theorem from a short parameter scan.

## 5. The exact positive certificate and its consequence

The [producer](../../04-computation/second_20260906_laurent.py)
constructs (2)--(6) over exact rational polynomial rings and verifies
Cayley-Hamilton as a second symbolic consequence identity. For
`k=1,...,5`, it obtains

```text
c_k(x)=d_k sum_(j=0)^(10k) a_(k,j) (x-1)^j,
d_k>0,             every a_(k,j) is a strictly positive integer. (7)
```

There are 155 positive coefficients in all. The **complete**, exact
d_k and a_(k,j) are frozen in five `CHAR_CERTIFICATE` JSON lines in the
[matching output](second_20260906_laurent.out); each line gives k, degree,
positive content and the integer coefficients in ascending order.
This finite certificate is part of the proof. It establishes strict
positivity of all five c_k for every real x>=1, not merely sampled x.

Consequently C_x(z)>0 for every real z>=0. At any real root rho of p_x,
evaluation of the multiplication operator, or Cayley-Hamilton followed
by evaluation, gives

```text
C_x(q_x(rho))=0.
```

The real number q_x(rho) therefore cannot be nonnegative. Thus
`q_x(rho)<0`. For integral x>=1, THM-4436 places all five roots on
this real negative axis, so (3) proves the stated noncancellation.
The polynomial proof does not need root simplicity for its real-root
implication; simplicity supplies the exact count of five attained phases.

The [independent reconstruction](second_20260906_laurent_audit.md)
starts from literal count fibers at x=1,...,51, constructs its own quotient
multiplication matrices, computes characteristic coefficients by principal
minors, and reconstructs the shifted polynomials by Newton interpolation.
The degree bounds above make this a complete identity certificate; it
checks the entire 155-coefficient positivity object. It imports no
producer. This differs from the producer's symbolic matrix powers and
trace recurrence.

## 6. Why the attempted derivative-cone entry is insufficient

This boundary is on an **actual**, fully normalized carried support,
already covered by the inherited endpoint15 theorem. Take h=2,x=1,g=8,
with support `(-15,1,9)` and monic first row `p=t^2+10t+1`. The
alpha-completed common carrier is

```text
H=(1+u)^8 (u^6+7t*u^4+10t^2*u^2+t^3),       k=5.
```

The carrier has degree fourteen. At either negative root of p it is
real-rooted: the cubic `y^3-7y^2+10y-1` has one positive root in each
of `(0,1),(1,2),(5,6)`, so `1+7t+10t^2+t^3` has three negative
roots, and the quadratic pullback in H has only real roots.

Reduce the five derivative-square generators and the actual normalized
row `t^2 Q_raw` modulo p. Independently divide each nonzero vector by
its positive integer gcd. In constant/linear coordinates the results are

| Generator | Primitive coordinate vector |
|---|---|
| r=0 | `(43382112,429439087)` |
| r=1 | `(115957212,1147859323)` |
| r=2 | `(211037892,2089062211)` |
| r=3 | `(3652956,36160583)` |
| r=4 | `(785924,7779857)` |
| Actual target | `(29132922,288386281)` |

For the linear functional `L(v)=det(v_0,v)`, its five ray values are

```text
0, 2458284732, 11149685028, 372200124, 142706596,
```

whereas its actual-target value is `-9483716742`. Hence the target
lies outside the nonnegative cone generated by all these rays.
The [minimal-cone theorem](creative_20260906_laurent_bridge.md) proves
that these five rays generate every permitted two-sided window, so
adding higher coupled windows does not repair this obstruction. Quotient
reduction is allowed here: even a polynomial identity valid only at the
first roots cannot produce the desired constant nonnegative combination.
The actual row is nevertheless negative at both first roots by inherited
endpoint15; cone membership was only sufficient.

A second tempting repair also fails at scale. For h=5, set y=t^2 and
seek the unique normalized identity

```text
q_x(t) (1+A_1 y+A_2 y^2)
 = W_x(t) (B_0+B_1 y+B_2 y^2) mod p_x(t),                    (8)
```

where W_x is the actual alpha-completed beta-hit response, retaining its
full binomial coefficients. Explicitly, with `g=x+16`,

```text
W_x(t)=sum_(e=-1)^10 binom(2g,2+2e) t^e
         *sum_(a+b=10-e) binom(x+5+2a,3a) binom(x+5+2b,3b),
```

where a,b are nonnegative and the binomials use their ordinary integer
values. At x=10000 the exact five-by-five solve gives
coefficient signs, in order `(B_0,B_1,B_2,A_1,A_2)`,

```text
(+,-,-,-,-).
```

The first implication to fail is extrapolating positive multiplier
coefficients from a short x bank. The retained exact hostile exposes the
missing large-x coordinate. The checker verifies the exact identity (8)
and this sign witness using fixed-length integer binomial sums. This
refutes that particular positive-coefficient multiplier certificate,
not the actual return theorem, not all phase multipliers, and not any
possible inequality valid only at the first roots.

## 7. Verification, files and next boundary

The producer's literal bank contains x=1,...,24, including nonprimitive
rows for coefficient identities; fifteen are primitive first-return
instances. Every first count is checked against its charge and mass,
every doubled coefficient against the literal multinomial, and the lower
carry is retained. The g=18 early-return hostile, actual derivative-cone
separator, and large-x multiplier hostile are explicit controls.

```bash
python3 -B 04-computation/second_20260906_laurent.py
python3 -B -O 04-computation/second_20260906_laurent.py
```

The producer passes 397 explicit gates. The independent audit supplies
the separate interpolation certificate and a comparison of all coefficients.
Normal and optimized producer outputs agree byte for byte. Raw LF hashes:

```text
source 02be5390867e0e7fb9e1f8610d7ef6ce1e822e0cf92e84bcaf1815fd46588407
output 9d5ece94d9174f75cecf06a62954181547b8a1d57c9c16628143cf6db201bc93
semantic f26269010352004ee03320c86de7b2b789fb23d967538ec4a679fb357e089ace
```

The connection preserves the complete actual coefficient row, not only
support or real-rootedness. Passing to the first-root quotient forgets
values away from the first cancellation locus; the polynomial (2), its
carry identity (4), and the coefficient monomial in (3) restore the
necessary interpretation. Characteristic positivity tests the actual
signed response at every first root. This is the concrete consumer of
the quotient representation.

The next endpoint in this progression is 39, with seven first channels.
No certificate for it was computed here. The remaining general obligation
is an all-h positivity proof or another uniform same-root mechanism;
the two explicit observer failures above constrain that search.
