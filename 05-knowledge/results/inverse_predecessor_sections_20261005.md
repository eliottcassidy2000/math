# Small inverse predecessors, exact incoming mass, and an adaptive repair

2026-10-05. **PROVED, scoped:** sharp predecessor-size bounds, a complete
three-class inverse palette, an increasing unit inverse ray, failure of a
uniform size-local flow comparison, and an adaptive finite-head repair.
**FINITE-EXACT:** the declared controls. **OPEN:** positivity at every
integer and universal Collatz. Inverse branching transports a common-future
predicate; it does not provide a new grounded endpoint.

Artifacts: [program](../../04-computation/experiments/inverse_predecessor_sections_20261005.py)
and [output](inverse_predecessor_sections_20261005.out).

## 1. Inheritance and the exact question

The closest proved mechanism is the full three-channel inverse table in
[inverse rays and ternary addresses, section 2](inverse_ray_ternary_addresses_20261004.md).
It inherits the older [reverse-tree pieces](collatz_mod6_20260917_reverse_tree_pieces.md).
The least multiple-of-three predecessor and its section were already proved
in [fusion helpers, Proposition S1](collatz_fusion_helpers_20261005.md).
Incoming commit bde61bf16,
[adaptive mixture flow, P3 and P5](collatz_adaptive_mixture_flow_20261005.md),
adds the exact infinite-fibre payment and the injection measure on rooted
multiples of three. Those constructions are inherited here.

The phrase “every odd 3-unit has a predecessor below 22 times itself”
conceals a useful qualification: the inherited predecessor can be chosen
**divisible by three**. Unrestricted predecessors admit a much smaller bound.

The canonical hostile is the removed self-loop at 1. The corrected near
miss is inferring a uniform positive weight ratio from a bounded size ratio.
The least-used sidecar is the completed route's sibling-depth counter.
The board is **target / inverse residue class / size / rooted support /
incoming mass / adaptive cutoff**. No identification with the primes
2, 3, 11 or an automorphic object is used: 22 here is a convenient integer
above the arithmetic supremum 64/3.

## 2. The complete small-predecessor palette

Write \(U(x)=\operatorname{oddpart}(3x+1)\), \(S(x)=4x+1\).
A positive odd target n has an odd predecessor iff \(3\nmid n\).
For such n, all predecessors are

\[
 b_j(n)=\frac{2^{a(n)+2j}n-1}{3},\qquad j\ge0,\qquad
 a(n)=\begin{cases}2,&n\equiv1\pmod3,\\1,&n\equiv2\pmod3.\end{cases} \tag{1}
\]

Their exact valuations are \(a(n)+2j\). They are strictly increasing with j
and satisfy \(b_{j+1}=S(b_j)\).

For each desired source class \(c=0,1,2\pmod3\), there is a unique exponent
\(\kappa_c(n)\in\{1,\ldots,6\}\) with \(2^\kappa n\equiv1+3c\pmod9\).
This follows from the order-six cycle of 2 modulo 9:

| n mod 9 | 1 | 2 | 4 | 5 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| source 0 mod 3 | 6 | 5 | 4 | 1 | 2 | 3 |
| source 1 mod 3 | 2 | 1 | 6 | 3 | 4 | 5 |
| source 2 mod 3 | 4 | 3 | 2 | 5 | 6 | 1 |

All sources in that class use exponents \(\kappa_c+6t\), \(t\ge0\).
This is the inherited palette, not a new three-state completeness theorem
for iterated inverse dynamics.

For \(n>1\), the first three predecessors in (1) are **exactly** all the
predecessors below \(22n\), and contain one of each source class modulo 3.
Indeed their exponents are at most 6, so each is below \(64n/3<22n\).
The fourth exponent is at least 7, and
\((128n-1)/3>22n\). Multiplication by 4 cycles the three source classes.

### Root boundary

At target 1, (1) begins \(1,5,21,\ldots\). A strict first-hit certificate
omits the self-loop 1. The smallest unrestricted or unit predecessor is
then 5; the smallest three-divisible predecessor is 21. In the prescribed
source class \(1\pmod3\), removing 1 moves the first source to 85
(exponent 8). Thus a claimed uniform exponent-six bound for every
prescribed class must retain that root exception.

## 3. Sharp size bounds and an infinite bounded-growth inverse ray

For n>1, the least unrestricted predecessor is below \(4n/3\).
The least predecessor coprime to three has the exponent table

\[
 (2,1,2,3,4,1)
 \quad\text{on target classes }(1,2,4,5,7,8)\pmod9.
\]

Thus its size is below \(16n/3\). The least predecessor divisible by
three, given by the first palette row, is below \(64n/3\).
Both latter bounds also hold at target 1 after removing its self-loop.

These constants are sharp suprema for the stated sections:
targets \(n=18t+1\) give the unrestricted ratio
\(4/3-1/(3n)\) and the three-divisible ratio \(64/3-1/(3n)\);
targets \(n=18t+7\) give the least unit ratio \(16/3-1/(3n)\).
Here t tends through positive integers. No equality with the supremum
occurs at a finite target.

The corresponding **size-induction boundary** is exact for n>1:

* a smaller unrestricted predecessor exists iff \(n\equiv2\pmod3\);
* a smaller unit predecessor exists iff \(n\equiv2\) or \(8\pmod9\);
* a smaller three-divisible predecessor exists iff \(n\equiv5\pmod9\).

For exponent a=1 the inverse is \((2n-1)/3<n\); for every a>=2 it
is at least \((4n-1)/3>n\). The tables identify exactly when exponent
one belongs to the required class. Since all later predecessors are
larger, these are impossibility statements for that predecessor type,
not merely failures of a particular selection rule. The other five unit
target classes force expansion in the three-divisible section. Its
factor-22 bound cannot by itself provide ordinary-size induction.

There is also a canonical **strictly increasing** unit predecessor R(n).
For n>1 choose exponents

\[
 (2,3,2,3,4,5)
 \quad\text{on the same six target classes},                    \tag{2}
\]

and set \(R(1)=5\). These are the least exponents that give a unit source
larger than n. Directly,

\[
 n<R(n)<\frac{32}{3}n,\qquad 3\nmid R(n),\qquad U(R(n))=n.        \tag{3}
\]

The constant 32/3 is sharp along \(n=18t+17\). Iteration therefore produces
an infinite increasing sequence of distinct unit predecessors with

\[
 n<R^d(n)<(32/3)^d n \qquad(d\ge1).                            \tag{4}
\]

Every finite segment is an actual route back to the supplied n.
It is rooted iff n is rooted. Even an infinite collection of inverse
ancestors retains that same unresolved terminal obligation. This is a
constructive inverse ray, not a proof that any supplied n reaches 1.

The inherited ternary-address notes explain the changing phase: computing
a child's residue consumes an additional ternary digit of its target.
For example, 7 and 25 agree modulo 9, but their least unit predecessors
37 and 133 have different classes modulo 9. Full source/carry data, not
just a palette label, must accompany composition.

## 4. Fixed size caps lose all relative mass on an explicit rooted family

Use the incoming mixture W. A **completed** nonroot target has counters
\((L,K)\) and weight

\[
 w(L,K)=\frac{2K!(L+1)!}{(L+K+2)!}.
\]

Its j-th predecessor in (1) has counters \((L+1,K+j)\).
P3 of the incoming note proves the full row identity. Its finite head gives

\[
 \frac{\sum_{j=0}^{J-1}W(b_j(n))}{W(n)}
 =1-R_J(L,K),\qquad
 R_J(L,K)=\frac{w(L,K+J)}{w(L,K)}
 =\prod_{i=0}^{L+1}\frac{K+1+i}{K+J+1+i}.                    \tag{5}
\]

The equality concerns nonroot targets, so it never inserts the forbidden
root loop into a first-hit route.

Take \(n_K=S^K(1)\), with \(K\ge1\) and \(K\equiv0\) or 1 modulo 3.
Then \(3\nmid n_K\), \(U(n_K)=1\), and its counters are \((0,K)\).
It is an explicit certified family. Formula (5) becomes

\[
 R_J(0,K)=\frac{(K+1)(K+2)}{(K+J+1)(K+J+2)}.
\]

For every fixed J, the captured fraction tends to zero as K grows.
In particular the entire inverse palette below \(22n_K\) has fraction

\[
 \boxed{1-\frac{(K+1)(K+2)}{(K+4)(K+5)}\longrightarrow0.}     \tag{6}
\]

The same obstruction holds for **every fixed size cap C>0**. Choose an
integer J with \(2\cdot4^J\ge3C+1\). If \(b_j(n)<Cn\), (1) implies

\[
 2^{a(n)}4^j<3C+1/n\le3C+1,
\]

so j<J. Thus the C-sized head is contained in a fixed finite sibling head,
and its relative mass tends to zero on the same rooted family.

Consequently no constant c>0 can make “the total W-mass of predecessors
below Cn is at least cW(n)” hold for all rooted unit targets. This is an
exact obstruction to this size-local flow comparison, not to all possible
weights or all positivity arguments. The first failed implication is
**bounded arithmetic size ratio implies uniformly controlled flow ratio**.

### A nearby three-divisible source can have tiny weight

For \(K=9t+1\), \(n_K\equiv5\pmod9\). Its least predecessor is

\[
 z_t=\frac{2n_K-1}{3}
     =\frac{2\cdot4^{K+1}-5}{9},\qquad
 3\mid z_t,\quad U(z_t)=n_K,\quad U(n_K)=1.
\]

Despite \(z_t<n_K\) and \(z_t/n_K\to2/3\),

\[
 \frac{W(z_t)}{W(n_K)}=\frac2{K+3}\longrightarrow0,\qquad
 \lambda(z_t)=W(z_t)=\frac4{(K+1)(K+2)(K+3)}.                \tag{7}
\]

These distinct, two-odd-step-certified injection atoms already show that
the inherited probability measure lambda has infinite positive ordinary
size moments: \(\sum n^p\lambda(n)=\infty\) for every p>0.
It also has infinite squared logarithmic size moment. Along this family,
size is exponential in t, atom mass is proportional to \(t^{-3}\),
and squared log size times mass is proportional to \(t^{-1}\).
This is compatible with the inherited finite expected **odd-step** root
time; size and route length are different observables.

## 5. Adaptive incoming capture, with its exact cost

There is a positive repair when the completed counters are retained.
Given rational \(0<\delta<1\), choose the first J siblings until the exact
tail (5) is at most delta. Monotonicity makes the least such J accessible
by exact integer/rational binary search.

A convenient explicit sufficient bound is

\[
 \boxed{J=
 \left\lceil(\delta^{-1}-1)\frac{K+L+2}{L+2}\right\rceil.}    \tag{8}
\]

Indeed every factor in (5) is at most
\((1+J/(K+L+2))^{-1}\), so Bernoulli's inequality gives

\[
 R_J(L,K)\le
 \left(1+\frac{J}{K+L+2}\right)^{-(L+2)}
 \le\left(1+\frac{(L+2)J}{K+L+2}\right)^{-1}\le\delta.
\]

The actual size cost remains explicit:

\[
 b_{J-1}(n)<\frac{4^J}{3}n.                               \tag{9}
\]

On the root-ray targets, the exact minimal count satisfies

\[
 \frac{J_{\min}(0,K,\delta)}K\longrightarrow
                 \delta^{-1/2}-1.                       \tag{10}
\]

This follows directly from the quadratic expression for \(R_J(0,K)\);
rounding changes J by at most one. Thus fixed captured mass can require
exponentially large inverse size factors in the retained counter K.
The type of the input matters: (8) takes certified route counters, not an
independently known bound on an arbitrary integer's unknown route.

| Source | Map | Preserved fact | Lost or missing information |
|---|---|---|---|
| Unit target n | Least inverse section | Actual one-step common future | Whether n is rooted |
| Completed target and (L,K) | Adaptive first-J inverse head | Chosen fraction of incoming W-mass | Source order unless the certificates stay attached |
| Fixed C-sized inverse palette | Numerical size restriction | Concrete integer predecessors | No uniform positive fraction of W |
| Rooted injection atoms z_t | Size-moment test | Explicit two-step certificates | No conclusion about unseen injection atoms |

A positive lower bound on any one of these predecessors would transport
to n through the incoming inequality. The existence and size bound alone
do not supply that positive lower bound.

## 6. Reproduction and declared universes

    python -B 04-computation/experiments/inverse_predecessor_sections_20261005.py
    python -B -O 04-computation/experiments/inverse_predecessor_sections_20261005.py

The program checks all positive odd 3-unit targets below 4096 against
independent direct exponent searches; all three residue classes; sharpness
identities; the exact factor-22 palette; 12 increasing unit rays through
16 inverse steps; and all formal pairs \(0\le L\le12,\ 0\le K\le48\).
The formal array controls compare factorial weights against product tails,
and test exact optimal cutoffs for delta=1/2,1/4,1/10.
They do not assert that every formal pair is an arithmetic ROOT word.
Six explicit root-ray targets check actual inverse rows and size costs.
Sixteen exact injection atoms check (7). Malformed exact-type and
root-self-loop controls are included. No orbit census is used to infer
universal positivity.
