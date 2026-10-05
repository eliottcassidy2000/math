# A computable summable Green weight, with positivity as the remaining Collatz obligation

2026-10-05. **PROVED:** the exact operator norm, the unconditional rooted
Green construction and its effective error bounds, its maximality, the
characterization of its positive support, and the finite inverse-kernel
compiler. **FINITE-EXACT:** the explicit controls below.
**OPEN:** positivity at every primitive base, equivalently universal
positive Collatz. A convergent numerical weight approximation does not
prove that every coordinate is nonzero.

Artifacts: [program](../../04-computation/experiments/collatz_green_weight_20261005.py)
and [output](collatz_green_weight_20261005.out). All arithmetic in the
program is integral or rational and remains checked under optimized Python.

## 1. Inheritance and the strengthened target

The closest mechanism is C3--C4 of
[three-bit sibling flow](collatz_three_bit_sibling_flow_20261005.md).
It reduces a positive summable incoming flow to individual inequalities
on sibling bases. Its conditional construction used already terminating
paths to allocate positive weights. The earlier
[effective prefix mass](collatz_effective_prefix_mass_20261005.md), P6,
gives the underlying discounted-flow equivalence.

Incoming C8 of [three-bits critical flow](collatz_three_bits_critical_flow_20261005.md)
and the corrected [refuel-bill code](collatz_refuel_bill_code_20261005.md),
read at origin/main commit 69f6b3903, already prove the unconditional
rooted product construction, critical summability, computable real
coordinates, maximality among summable feasible weights, and support
equal to the root basin. The construction below independently recovers
those results through a positive operator; it does not claim them as new.

The additional refinements here are the exact root column, its improved
critical total-mass bound, and a finite sibling-depth truncation with a
certified global error in the summable-sequence norm. In particular,
at critical rho=1 the bound improves from
\((1+r+r^2)/r^2\) to \(1+r+1/r\). Strict positivity everywhere
remains distinct from constructing and computing the nonnegative weight.

The inherited hostile is an arbitrary long positive climb: weights
bounded above and below by fixed multiples of the source-height prior
cannot pay every edge. The new hostile is a contractive operator whose
root cannot be reached from another component. The least-used sidecar is
the root boundary condition in the operator, together with the forbidden
ternary phase in each inverse sibling fibre. The board is **source /
sibling base / incoming fibre / rooted operator / precision / nonvanishing**.

## 2. The base operator and its exact norm

Let \(U(n)=(3n+1)/2^{v_2(3n+1)}\) on positive odd integers and
\(S(n)=4n+1\). Let

\[
 {\cal B}=\{b>0:b\text{ odd},\ v_2(3b+1)\in\{1,2\}\}.
\]

Every positive odd integer has a unique decomposition \(S^k(b)\),
\(b\in{\cal B}\). A total decoder repeatedly applies \((n-1)/4\) while
\(n\equiv5\pmod8\). For \(b\ne1\), define \(G(b)\) and \(k(b)\) by

\[
 U(b)=S^{k(b)}G(b),\qquad G(b)\in{\cal B}.
\]

Fix rational \(0<r<1,\ 0<\rho\le1\), put \(c_0=\rho(1-r)\), and set
\(d_b=c_0r^{k(b)}\). On sequences indexed by \({\cal B}\), define

\[
 (Af)(1)=0,\qquad (Af)(b)=d_b f(G(b))\quad(b\ne1).       \tag{1}
\]

The root row is zero. It is not the actual root self-loop with a discount
silently attached.

For a target \(y>0\) odd, a predecessor exists exactly when \(3\nmid y\).
Its unique minimal sibling base is

\[
 b_0(y)=
 \begin{cases}
 (2y-1)/3,&y\equiv5\pmod6,\\
 (4y-1)/3,&y\equiv1\pmod6.
 \end{cases}                                         \tag{2}
\]

Thus for a fixed parent \(c\in{\cal B}\), possible children are exactly
\(b_0(S^k(c))\), omitting targets divisible by 3 and the child 1.
The children for different \(k\) are distinct: their actual U-targets differ.
Since \(S^k(c)\equiv c+k\pmod3\), precisely one depth class modulo 3 is
forbidden.

Write \(j(c)\in\{0,1,2\}\) for \(-c\bmod3\). The exact column sum is

\[
 C(c)=\sum_{b:G(b)=c}d_b
 =\rho\left(1-\frac{r^{j(c)}}{1+r+r^2}\right)
   -{\bf1}_{c=1}\rho(1-r).                            \tag{3}
\]

The subtraction at 1 removes depth zero, which would be the root source.
For nonroot parents the largest value occurs at \(c\equiv1\pmod3\);
the base 7 realizes it. Positivity and Tonelli's theorem therefore give

\[
 \boxed{\ \|A\|_{\ell^1\to\ell^1}
     =\kappa=\frac{\rho(1+r)}{1+r+r^2}<\rho\le1.\ }       \tag{4}
\]

Indeed, the norm of a positive matrix on \(\ell^1\) is the supremum of its
column sums: the upper bound follows by summing absolute values, and
a unit vector at a maximizing column attains it. The row bound is also

\[
 \|A\|_{\ell^\infty\to\ell^\infty}\le c_0<1.
\]

The root column is the sharper starting bound

\[
 a_0=\|A\delta_1\|_1
     =\frac{\rho r(1+r^2)}{1+r+r^2}.                   \tag{5}
\]

At \(r=1/16,\rho=1/2\), these constants are
\(c_0=15/32,\ \kappa=136/273,\ a_0=257/8736\).
The improvement over \(\rho\) comes from the missing ternary phase, not
from a numerical spectral experiment.

## 3. The unconditional Green construction

Define

\[
 H=\sum_{j\ge0}A^j\delta_1.                            \tag{6}
\]

This series converges in \(\ell^1\), without any Collatz hypothesis.
It is nonnegative and satisfies

\[
 H=\delta_1+AH,\qquad H(1)=1.
\]

The geometric norm estimates give

\[
 \|H\|_1\le1+\frac{a_0}{1-\kappa},\qquad
 \left\|H-\sum_{j=0}^N A^j\delta_1\right\|_1
 \le\frac{a_0\kappa^N}{1-\kappa}.                     \tag{7}
\]

For the default parameters the total-mass upper bound is \(4641/4384\).
At the critical endpoint \(\rho=1\), the exact root column gives the
particularly useful improvement

\[
 \boxed{\ \|H\|_1\le1+\frac{a_0}{1-\kappa}
       =1+r+\frac1r.\ }                              \tag{7a}
\]

This is \(273/16\) for \(r=1/16\) and \(21/4\) for \(r=1/4\).
The inherited general-column estimate is
\(1/(1-\kappa)=(1+r+r^2)/r^2\), respectively 273 and 21.
The improvement uses the smaller first-step mass from the killed root,
followed by the worst-column norm for all later steps. It is a rigorous
upper bound, not a claim that the Green vector attains that bound.

All these are weight norms, not the remaining probability under either
source prior.

There is at most one nonzero term of (6) at a fixed base. If the first
G-hit of root occurs after \(\tau\) steps, then

\[
 H(b)=\prod_{i=0}^{\tau-1}d_{G^i(b)}>0.                \tag{8}
\]

If there is no such finite hit, every term at b is zero, so \(H(b)=0\).
Root row deletion is what prevents padding a first hit by extra root loops.

The U-orbit and G-orbit have the same root-membership predicate. If a
G-step has parent c, then \(U(b)=S^k(c)\) and
\(U(S^k(c))=U(c)\); a finite G-root path therefore transports an actual
root proof. Conversely, along a rooted U-orbit, its first U-step decreases
the odd root rank, and stripping sibling depth preserves that rank or
reduces it further at base 1. Hence its G-path terminates too. Consequently

\[
 \boxed{\ H(b)>0\quad\Longleftrightarrow\quad b
         \text{ has an actual finite ROOT route}.\ }  \tag{9}
\]

Equations (6)--(9) construct H unconditionally. They do not establish the
left side for every b.

### Maximality and the original flow

Let g be any bounded nonnegative feasible weight with \(g(1)=1\) and
\(g(b)\le d_b g(G(b))\). Iterating to a finite root gives \(g(b)\le H(b)\).
On an orbit avoiding root, iteration gives
\(g(b)\le c_0^j\|g\|_\infty\) for every j, hence \(g(b)=0\).
Therefore H is the pointwise greatest bounded feasible weight with root
value 1; in particular it is greatest among the summable ones.

Extend H to all odd sources by \(f(S^k b)=r^kH(b)\). Its total mass is
\(\|H\|_1/(1-r)\). For the incoming operator with source and target root
deleted, every nonroot 3-unit target y has

\[
 {\cal K}f(y)=\rho f(y).
\]

At a target divisible by 3 there are no predecessors, so
\({\cal K}f(y)=0\le\rho f(y)\). The root value \(f(1)=1\) is bookkeeping;
no inequality is imposed at that target. If one separately counts its
nonroot predecessors, their total weight is \(r/(1-r)\).

Thus every desired flow inequality already holds. At rho=1 this is a
critical incoming inequality rather than a strictly discounted one;
the base operator still contracts because its missing ternary phase
gives kappa<1 and each base-edge factor is at most 1-r<1.
The outstanding requirement is positive weight at every input, exactly (9).

## 4. An always-terminating precision algorithm

To approximate a coordinate H(b), trace at most N actual G-edges. If
root is hit, return the finite product (8), exactly. Otherwise return
the interval

\[
 [0,c_0^{N+1}].                                      \tag{10}
\]

If a later hit exists it needs at least N+1 edges and its product is at
most \(c_0^{N+1}\); if it never hits, the true value is zero. Thus (10)
is sound in both cases. For every rational error tolerance, N can be
chosen effectively before starting the trace. This proves uniform
computability of H as a real-valued function without using an unknown
stopping time.

A zero approximation does not certify a zero limit. The base 27 is an
explicit control: depth zero returns a zero lower approximation, while
its checked G-root depth is 40 and its exact weight is positive. It can
be much smaller than a routine error tolerance.

Nor does contraction force positive support in general. On the abstract
graph with root 1 isolated and \(G(n)=n+1\) for \(n\ge2\), take \(d=1/2\).
Then \(\|A\|_1=1/2\) but the Green vector is simply \(\delta_1\), zero on
every nonroot vertex. This is a hostile to that inference, not a Collatz
counterexample.

If one proves any positive explicit lower bound \(L(b)\le H(b)\), then
choosing N with \(c_0^{N+1}<L(b)\) forces a root hit within N edges.
A lower bound valid at every base would therefore discharge the global
obligation. The precision theorem supplies an upper error bound; it
does not supply such a positive lower bound.

## 5. Finite-support approximation with an exact infinite error bound

Computability in \(\ell^1\) also needs control of the infinitely many
children, not just a bounded path depth. Let \(A_K\) retain sibling
depths \(0\le k\le K\). Each retained column has finitely many children,
and the omitted geometric tail gives

\[
 0\le A_K\le A,\qquad
 \|A-A_K\|_1\le\varepsilon_K=\rho r^{K+1}.             \tag{11}
\]

The resolvent identity, or the telescoping identity for powers, yields

\[
 \left\|\sum_{j\ge0}A^j\delta_1
       -\sum_{j\ge0}A_K^j\delta_1\right\|_1
 \le\frac{\varepsilon_K}{(1-\kappa)^2}.               \tag{12}
\]

For example,
\(\|A^j-A_K^j\|_1\le j\kappa^{j-1}\varepsilon_K\);
summing this elementary geometric derivative proves (12).

The finite rational kernel

\[
 H_{N,K}=\sum_{j=0}^N A_K^j\delta_1
\]

therefore satisfies

\[
 0\le H_{N,K}\le H,\qquad
 \|H-H_{N,K}\|_1
 \le\frac{a_0\kappa^N}{1-\kappa}
     +\frac{\rho r^{K+1}}{(1-\kappa)^2}.              \tag{13}
\]

Choose N and K to make each summand at most half the requested tolerance.
This is a terminating finite-support precision compiler. It discovers
bases by the exact inverse formula (2), not by assuming future forward
searches terminate. Each retained base has its checked parent equation,
and every retained parent is also retained. The supported weights obey
the exact edge equality; weights at all other bases are zero and satisfy
the inequality. Every sibling of a retained base is therefore rooted,
with the total source mass counted once using the unique base.

For the default parameters, tolerance 1/100 uses 16 bases at N=4,K=2;
tolerance 1/1000 uses 128 bases at N=7,K=2. These are certified Green-weight
errors. In particular the latter kernel still omits base 27. Fast
approximation of H is not a bound on unresolved input probability under
the gamma or height source priors.

## 6. Cancelling the sibling depth in a prior-relative coordinate

The earlier source prior satisfies
\(\nu(S^k n)=16^{-k}\nu(n)\). At r=1/16, write H(b)=\(\nu(b)h(b)\).
For a nonroot base the equality becomes

\[
 h(b)=a_b h(G(b)),\qquad
 a_b=\rho\frac{15}{16}\frac{\nu(U(b))}{\nu(b)}.         \tag{14}
\]

The unbounded k cancels from this scalar multiplier. Since a base's
actual valuation is 1 or 2, the target's binary height differs by at most
one. At \(\rho=1/2\),

\[
 a_b\in\{15/128,\ 15/32,\ 15/8\}.
\]

The last value exceeds one. The target G(b) still requires the full
integer and its exact sibling-depth decoder; a three-letter multiplier
alphabet is not a finite-state proof of coverage.

Because H(1)=1 and \(\nu(1)=2/3\), this factored coordinate has
h(1)=3/2. If root value 1 is desired in the factored coordinate, use
\((2/3)H/\nu\); the nonroot multipliers are unchanged. This distinction
prevents changing a boundary normalization while leaving its recurrence
apparently intact.

## 7. Reproduction and scope

    python -B 04-computation/experiments/collatz_green_weight_20261005.py
    python -B -O 04-computation/experiments/collatz_green_weight_20261005.py

Declared universes:

* all 192 primitive bases below 512, at depth cutoffs 0 through 8, for
  nine rational parameter pairs, including critical rho=1; direct children plus exact geometric
  tails agree with every column formula and the attained norm;
* three critical root-column bounds and finite critical inverse kernels;
* 1,728 pointwise precision intervals checked against explicitly grounded
  longer controls, with the actual path products retained;
* five finite kernel choices (N,K)=(0,0),(1,1),(2,2),(4,3),(5,5),
  independently checked by their forward parent decoder, plus larger
  nested controls;
* two finite global-precision compilations, a large exact tolerance control,
  and fifteen malformed-input or
  killed-root controls;
* the abstract isolated-root example showing why contraction alone is
  insufficient for positive support.

The result provides a particular weight everywhere as a computable real,
and proves all its inequalities and summability. Its positivity set is
exactly the already-defined root basin. The remaining research target is
a source-sensitive nonvanishing bound or another argument that this
explicit Green weight has no zero coordinates.
