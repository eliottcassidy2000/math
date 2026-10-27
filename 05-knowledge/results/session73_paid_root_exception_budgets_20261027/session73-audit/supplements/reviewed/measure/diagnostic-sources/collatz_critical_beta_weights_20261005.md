# Parameter-averaged critical weights and the remaining positivity problem

2026-10-05. **PROVED, scoped:** a critical beta mixture, finite total mass,
an exact mass conservation identity, polynomial comparison with every fixed
parameter, and finite rational approximation bounds. **FINITE-EXACT:** the
implementation controls. **OPEN:** strict positivity at every integer, hence
universal Collatz. These are research results, not a canon promotion or a
claim of literature priority.

## 1. Inheritance and the change of target

The anchor is a positive summable incoming flow on every positive odd
nonroot source. The niche is universal parameter averaging of route codes;
the wildcard is the distinction between numerical approximation and support
discovery. The board is **rooted support / inverse branching / refuel depth /
parameter mixture / finite certificate / positive lower bound**.

The closest proved mechanism is [C8 of the critical-flow note](collatz_three_bits_critical_flow_20261005.md),
integrated from `origin/main` commit `69f6b3903`. It already supplies a
canonical nonnegative, summable, computable-real flow, positive exactly on
the rooted component. [The refuel-bill audit](collatz_refuel_bill_code_20261005.md)
identifies its maximality and corrects interpretations of its logarithm.
Those are inherited results, not new conclusions of this note.

The companion [Green-operator refinement](collatz_green_weight_20261005.md)
improves their total-mass bound by using the root's actual inverse column.
The [discounted mixtures](collatz_mixed_refuel_weights_20261005.md) provide a
different useful tradeoff: stronger mass bounds and strict incoming payment.
Here the new mixture retains the critical, undiscounted endpoint.

The hostile example is a zero-weight unrooted component: excellent global
approximation need not reveal a positive weight at any particular missing
source. The corrected near miss is interpreting contraction of the
*backward weight solver* as termination of the *forward integer orbit*.
The least-used sidecar is the full certified address beneath a compressed
pair of route statistics.

The earlier measure route has external precedent. **CITED:** Assani's
*Collatz map as a non-singular transformation*, Theorem 4.1, characterizes
Collatz using a finite measure equivalent to counting measure, a power
bounded composition operator, and the specified conservative part.
Equivalence to counting measure retains positivity at every integer.
Our killed odd-map inequalities and beta formulas below are a separate,
explicit construction; we do not import positivity from that theorem.
[Primary text](https://arxiv.org/html/2208.11675).

## 2. The critical family and its sharper root estimate

Write \(U(n)=\operatorname{oddpart}(3n+1)\), \(S(n)=4n+1\).
Every positive odd integer is uniquely \(n=S^j(b)\), where
\(v_2(3b+1)\in\{1,2\}\). At a nonroot base set

\[
 U(b)=S^{k(b)}G(b).
\]

On a rooted base, let \(T\) be the number of \(G\)-edges to its first
hit of 1 and \(K\) the sum of their sibling-removal depths. For
\(n=S^j(b)\), put \(M=K+j\). These are base-edge counts, not an
assertion that the actual odd Collatz path has exactly \(T\) steps.
Keep the actual source and full address when compiling a ROOT certificate.

For \(0<r<1\), the inherited critical flow is

\[
 f_r(n)=(1-r)^T r^M
 \quad\text{on rooted sources},\qquad f_r(n)=0\text{ otherwise}.       \tag{1}
\]

Thus \(f_r(1)=1\), and \(f_r(S^j1)=r^j\). Define the killed operator
on \(V=\{3,5,7,\ldots\}\) by

\[
 (\mathcal K f)(y)=\sum_{n\in V:U(n)=y}f(n),\qquad y\in V.
\]

Every nonroot target not divisible by three has one primitive inverse
base and all its siblings. Summing that geometric fibre gives

\[
 \mathcal K f_r(y)=f_r(y)\quad(3\nmid y),\qquad
 \mathcal K f_r(y)=0\quad(3\mid y).                           \tag{2}
\]

This includes unrooted targets: all their predecessors are also unrooted.
There is no payment equation at the killed root.

Put \(D=1+r+r^2\). The base backward operator has exact norm
\(\kappa=(1+r)/D<1\), but its root column has the smaller sum
\(a_0=r(1+r^2)/D\). The root self-loop has been omitted. Therefore
the Green refinement gives

\[
 \sum_b g_r(b)\le1+\frac{a_0}{1-\kappa}
 =1+r+r^{-1}=D/r,
 \qquad
 \sum_n f_r(n)\le\frac{D}{r(1-r)}.                           \tag{3}
\]

This is unconditional mass control. It does not supply a nonzero value
at a previously unresolved base.

## 3. A parameter-free weight with total mass at most eleven

Average (1) using the probability density \(6r(1-r)\) on \((0,1)\):

\[
 E(n)=\int_0^1 6r(1-r)f_r(n)\,dr.
\]

**B1.** On every rooted source this is the explicit rational number

\[
 \boxed{E(n)=6B(M+2,T+2)
       =\frac{6(M+1)!(T+1)!}{(M+T+3)!}.}                    \tag{4}
\]

It is zero off the rooted component. Tonelli's theorem and (3) yield

\[
 \boxed{\sum_{n\text{ odd}>0}E(n)
 \le6\int_0^1(1+r+r^2)\,dr=11.}                           \tag{5}
\]

The factors in the density cancel both endpoint costs in (3). An
unweighted parameter average would already make the root ray divergent:
\(\int_0^1r^jdr=1/(j+1)\). Choosing an integrable parameter law is
part of the construction, not an optional normalization.

Linearity and nonnegativity preserve (2). In particular
\(\mathcal K E\le E\). The mixture is positive at exactly the same
sources as every \(f_r\); averaging creates no new positive support.

**B2 — exact boundary mass.** The root ray has weights and total

\[
 E(S^j1)=\frac6{(j+2)(j+3)},\qquad
 \sum_{j\ge0}E(S^j1)=3.                                   \tag{6}
\]

Consequently exactly two units of nonroot mass flow into the killed root.
Summing (2) over all nonroot targets, justified by (5), gives

\[
 \boxed{\sum_{n>0\text{ odd}:3\mid n}E(n)=2.}              \tag{7}
\]

Indeed \(\sum_V\mathcal K E=\sum_V E-2\), whereas (2) makes the
left side \(\sum_{y\in V,3\nmid y}E(y)\). Equation (7) does not
assert that all multiples of three are rooted. It describes the total
weight on the rooted multiples of three, with zero assigned to any others.
Nor does it determine the total mass in (5) exactly.

## 4. Averaging pays for choosing a good parameter

**B3.** For every rooted source and every \(0<r<1\),

\[
 E(n)\ge
 \frac{6(M+1)(T+1)}{(M+T+1)(M+T+2)(M+T+3)}f_r(n)
 \ge\frac6{(M+T+2)(M+T+3)}f_r(n).                         \tag{8}
\]

The first inequality is the elementary binomial bound
\(\binom{M+T}{M}r^M(1-r)^T\le1\); the second uses
\((M+1)(T+1)\ge M+T+1\). It also holds trivially off the support.
Thus one averaged rule loses at most a quadratic factor against the best
fixed parameter separately for each certified route. The earlier fixed
choice imposed exponential costs in both counters.

The source of this connection is the family of geometric address
probabilities; the target is one incoming subinvariant flow. The map is
positive integration. It preserves summability and payment, and discards
the chosen parameter. The pair \((T,M)\) suffices to evaluate the
weight, but discards the order and integer realizability of an address.
The full guarded ROOT certificate is the necessary sidecar. The cheapest
decisive test is the exact inverse-fibre identity in (2), including its
infinite remainder, not a comparison of typical trajectories.

This is a concrete universal-parameter coding move. It is not a claim
that a universal prefix machine or a randomness theorem certifies every
integer's route.

## 5. Numerical approximation terminates without a Collatz assumption

**B4 — pointwise interval.** Decode \(n=S^j(b)\) and follow at most
\(N\) base edges. A root hit gives (4) exactly. If no root is seen,
let \(K_N\) be the sibling depths already removed. Any eventual root
would need at least one more base edge. Hence

\[
 0\le E(n)\le6B(j+K_N+2,N+3)
 \le\frac6{(N+3)(N+4)}.                                  \tag{9}
\]

The same interval contains zero for an unrooted source. Its effective
width tends to zero uniformly. A zero lower endpoint therefore means
**unresolved at this cutoff**, not **proved nonconvergent**. The program
uses 27 as a hostile control: cutoff zero returns such an interval, while
the later ROOT certificate gives \(E(27)=1/21296186450>0\).

**B5 — finite support in total variation.** Construct the inverse tree
from ROOT, keeping only base depth \(T<N\), each inverse sibling depth
\(k\le B\), and outer sibling depth \(j\le B\). This is a finite
set of certified sources. Give each its full weight (4), and give other
sources zero. Denote the resulting submeasure by \(E_{N,B}\).

For \(N,B\ge1\) and rational \(0<\epsilon<1\), a computable bound is

\[
\begin{aligned}
 \|E-E_{N,B}\|_1\le{}&12\epsilon+
 12(1-\epsilon^2/3)^{N-1}\\
 &+6\left(\frac1B+\frac2{B+1}+\frac4{B+2}
                  +\frac3{B+3}+\frac2{B+4}\right).          \tag{10}
\end{aligned}
\]

To prove it, split the omitted sources into depth, inverse-width and
outer-sibling tails; overlaps can only improve this upper bound.

* The depth tail of (1) is at most
  \(a_0\kappa^{N-1}/[(1-\kappa)(1-r)]\). After mixing it is at most
  \(6\int_0^1(1+r^2)\kappa^{N-1}dr\).
  Since \(\kappa\le1-r^2/3\), splitting at \(r=\epsilon\)
  gives the first line of (10).
* If \(A_B\) keeps inverse depths through \(B\), then
  \(\|A-A_B\|_1\le r^{B+1}\). The resolvent identity and (3)
  give \(\|g_r-g_{r,B}\|_1\le D^2r^{B-2}\).
  After sibling extension and mixing this is bounded by
  \(6\int_0^1D^2r^{B-1}dr\), whose five coefficients are
  \(1,2,3,2,1\).
* Omitting \(j>B\) costs at most
  \(6\int_0^1D r^{B+1}dr
    =6[(B+2)^{-1}+(B+3)^{-1}+(B+4)^{-1}]\).

First choose \(\epsilon\), then \(N\), then \(B\). This proves
effective approximation by finite rational measures. The finite measures
are approximations; truncation is not claimed to preserve every incoming
inequality. A small error in E-mass cannot exclude a missing source whose
E-mass is zero. For coverage, the inherited independent positive atomic
reference measure must still assign that missing integer a positive atom.

## 6. What remains useful to prove next

The global analytic tasks of defining the weights, summing them, and
computing them to precision have solutions. The exact remaining condition
is

\[
                         E(n)>0\quad\text{for every positive odd }n. \tag{11}
\]

By C3 of the inherited critical-flow note, a strictly positive summable
critical flow forces every orbit to root: an infinite distinct orbit
contradicts summability; an extra incoming predecessor rules out any
remaining nonroot cycle. For this particular E, (11) is also immediately
equivalent to ROOT support by construction. Neither argument proves (11).

The [quantified lookahead obstruction](collatz_finite_lookahead_weight_obstruction_20261005.md)
explains why a finite observer alone cannot supply the missing lower
bounds. It constructs arbitrarily large positive negative-cycle shadows
with identical retained data at expanding endpoints. A local anchor-depth
correction pays whole guards, but its entry refuel bill is explicit.

The productive next test is therefore a positive lower bound on a named
unbounded source family across **its refuel boundary**, retaining the
actual integer and anchor precision. The parameter mixture eliminates
the need to guess one globally efficient geometric parameter. It leaves
that family-support obligation visible and unchanged.

## 7. Exact checks

[Program](../../04-computation/experiments/collatz_critical_beta_weights_20261005.py)
and [saved output](collatz_critical_beta_weights_20261005.out):

```text
python 04-computation/experiments/collatz_critical_beta_weights_20261005.py
python -O 04-computation/experiments/collatz_critical_beta_weights_20261005.py
```

There are 20,224 explicit checks. The universe includes all
\(0\le T,M\le24\) for independent polynomial integration, five rational
parameter controls, five exact inverse-fibre truncations with their
remainders, one hundred root-ray terms with its infinite tail, and 1,536
backward-generated sources at \(T<5\) and individual depths at most five.
Their addresses are independently recovered by forward iteration. Root
handling, zero-support components, unresolved 27, malformed inputs, and
coarse tail controls are hostiles. No finite check is used to establish
(11); the unbounded statements have proofs above.
