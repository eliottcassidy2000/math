# Run coordinates, vertex coordinates, and certified pole correction

**Session277. PROVED IN NOTES; exact computation reported separately.**
These results concern transitive inflations of a fixed tournament quotient.
They do not compute A000568, prove Chan–Pak, or classify all tournaments.

## Inheritance and the new junction

The quotient-walk/separator law and full path-cover response are inherited
from [Session268](../../session268_composable_tournaments_20261111/SYNTHESIS.md).
The canonical modular tree, infinite prime alphabet and finite-alphabet
limitations are inherited from [Session271](../../session271_predictive_tournaments_20261111/structures/ATLAS.md).
We do not claim a new substitution identity. We specialize it to transitive
children, expose its Stirling coordinate change, and prove a recurrence-class
dichotomy and a certified absolute-error repair for the resulting infinite
families. These are scoped refinements of those mechanisms.

A tangential source is the repo's moment-transform thread: THM-1600,
`01-canon/theorems/THM-1600-laplace-layer-exact-at-degree-one-and-span-two-eliminated.md`,
uses truncated exponentials to replace an asymptotic moment argument with an
exact sequence. Its usable idea here is the **choice of moment measure**, not
any GMC claim. Our measure is a discrete geometric distribution; no GMC or
Jacobian conjecture is assumed. The distinction between finite recurrences
before and after a change of coordinate is the connection to retain.

## 1. The exact multivariate compiler

Let Q have q named vertices and adjacency matrix A (zero diagonal). Let
T_Q(r)=Q[TT_(r_1),...,TT_(r_q)], allowing r_i=0 by deleting that module.
Set H(empty)=1. Define

    F_Q(x) = sum_(r_i>=0) H(T_Q(r)) product_i x_i^r_i/r_i!,
    u_i = exp(x_i)-1, U=diag(u_i).

Then, as a formal power series,

    F_Q(x) = 1 + 1^T U (I-AU)^(-1) 1.                 (1)

**Proof.** Collapse each maximal consecutive visit to one module in a
Hamiltonian path into a run. The resulting run word is a directed walk in Q.
If it visits module i in k_i runs, distributing its r_i totally ordered
vertices into k_i nonempty runs has k_i! S(r_i,k_i) possibilities: assign each
vertex to a run and put the assigned vertices in increasing order within that
run. Runs are not required to be intervals of the module's total order.
The EGF of those assignments is u_i^k_i. Summing walk words proves (1).
Every path has exactly one run word and assignment, so there is no overcount.
For fixed r, k_i<=r_i makes coefficient extraction finite. This yields an
exact physical constructor/evaluator, not merely a formal analogy.

The full vector r and quotient orientation are sidecars that must be retained.
Setting all x_i equal forgets the placement of sizes; it is useful for an
aggregate but not for the Chan–Pak target selector.

## 2. A sharp recurrence-class dichotomy

Define the integer aggregate

    B_n(Q) = sum_(sum r_i=n) n!/product r_i! * H(T_Q(r)).

This is a multinomial-weighted path total over named-module assignments, NOT
an ordinary isomorphism-class count and NOT the set of individual H values.
Let W_0=1 and W_k=1^T A^(k-1)1 for k>=1. Then

    B_n = sum_(k=0)^n W_k k! S(n,k),
    sum B_n t^n/n! = R_Q(exp(t)-1),
    R_Q(u) = 1+u 1^T(I-uA)^(-1)1.                    (2)

The run sequence W is always constant-coefficient recurrent by Cayley–Hamilton.
Nevertheless:

**Theorem.** B_n(Q) is P-recursive (satisfies a nonzero linear recurrence with
polynomial coefficients, eventually) if and only if Q is transitive.

If Q is transitive then A is nilpotent, R is polynomial, and (2) is a finite
linear combination of exp(jt). More precisely B_n=q^n: every inflation is
transitive and has exactly one Hamiltonian path.

If Q is not transitive then A contains a directed cycle. Its spectral radius
rho is positive. The power series sum W_k u^k has finite radius 1/rho and a
singularity at the positive real point 1/rho (nonnegative coefficients;
alternatively restrict to a maximal Perron SCC). Since it is rational, this is
a genuine pole. Composition with exp(t)-1 produces a pole at every

    t = log(1+1/rho) + 2*pi*i*j, j in Z.

Its derivative there is nonzero, so poles do not disappear. A D-finite analytic
function can have finite-plane singularities only at the finitely many zeros
of the leading polynomial of its differential equation. Hence this EGF is
not D-finite. A sequence is P-recursive iff its EGF is D-finite, proving the
claim. This includes eventual recurrences: changing finitely many coefficients
adds a polynomial EGF and cannot remove infinitely many poles.

**Practical meaning.** Failure to discover a polynomial-coefficient size
recurrence is sometimes structural, even for a three-vertex quotient with a
one-dimensional run recurrence. It does not mean that exact computation is
impossible. Rationality in exp(t), finite run recurrence, and nonholonomicity
in vertex size coexist.

## 3. Regular quotients: exact geometric moments, not a heuristic

For a regular quotient of order q=2d+1>=3, every vertex has d successors, so
W_k=q d^(k-1). Regardless of the isomorphism type of that quotient,

    F(t) = 1 + q(exp(t)-1)/(1-d(exp(t)-1)),
    B_n = q/d * sum_(k=1)^n d^k k! S(n,k), n>=1.     (3)

Thus every regular quotient of the same order has the same aggregate. This
is an explicit loss of structure: an exact recurrence for the aggregate does
not distinguish those quotients or their individual inflation target sets.

Put tau=log(1+1/d) and c=q/(d(d+1)). Expanding the denominator gives

    B_n = c sum_(m>=0) m^n exp(-tau*m)
        = c n! sum_(j in Z) (tau+2*pi*i*j)^(-n-1), n>=1.    (4)

The first equality is an absolutely convergent geometric moment. The second
follows by differentiating the classical partial-fraction expansion of
1/(1-exp(-tau)); after one derivative the pole sum is absolutely convergent.
It can also be proved by Poisson summation for x^n exp(-tau*x) on x>=0
(the value at zero is zero for n>=1). This is an exact infinite correction
formula, with no claim that a finite leading term already rounds correctly.

For K>=1 let P_(n,K) be the symmetric pole sum with |j|<=K, multiplied by c n!.
Then

    |B_n-P_(n,K)| <= 2c n! / (n (2*pi)^(n+1) K^n).    (5)

Indeed |tau+2*pi*i*j|>=2*pi*|j| and the tail of sum j^(-n-1) is at most
integral_K^infty x^(-n-1)dx. The errors are oscillatory: the two j-th poles
contribute 2c n! (tau^2+4*pi^2*j^2)^(-(n+1)/2)
cos((n+1) arctan(2*pi*j/tau)). There is no universal one-sided leading-pole
error. A tiny relative error still need not be less than half an integer.

**Certified recovery rule.** Compute a directed interval [L,U] for P_(n,K),
expand it by the bound in (5), and accept only if the expanded interval
contains exactly one integer. That integer is B_n. For every fixed n>=1,
taking K sufficiently large and arithmetic precision sufficiently high must
succeed. This is an all-order terminating certificate, not an assertion that
a particular finite table proves termination. Roughly linear K in n suffices
for fixed d asymptotically, by Stirling's formula applied to (5); the output
itself has Theta(n log n) digits. Our implementation chooses K adaptively and
retains all integer digits; if the analytic tail is below1/8 but the
interval is still ambiguous it doubles its extra precision and restarts.
Thus there is no fixed precision floor or arbitrary pole-count cutoff.
It does not claim constant-time generation.

The executable uses only Python's standard library. Log tau is enclosed by
the positive atanh series with a rational geometric tail. Pi is enclosed by
Machin's formula and alternating arctangent tails. Every subsequent real and
complex interval operation uses explicitly directed Decimal arithmetic.
The final certificate is checked against the independent finite Stirling sum
(3). The finite comparison is a control, not the basis of the certificate.

## 4. Chan–Pak: a finite target inverse for every fixed strong quotient

For a fixed nontrivial strong Q (q>=3) and positive module sizes r_i, every vertex of Q lies in
a directed triangle. (The out-neighbor/in-neighbor cut must have an edge from
the former to the latter; otherwise strong connectivity fails.) For each i,
T_Q(r) contains the induced tournament C3[TT_(r_i),1,1]. Its H is
2^(r_i)+1: there are three unsplit run words, plus the unique run word with
the large module split into two nonempty runs, contributing 2^(r_i)-2.

Hamiltonian path count is monotone under adding a vertex to a tournament:
each old directed path admits an insertion of the new vertex, and choosing
one deterministic legal slot gives an injection (deleting the new vertex
recovers the old path). Consequently

    H(T_Q(r)) >= 2^(max r_i)+1.                       (6)

The singleton quotient is excluded: its transitive inflations have H=1 at
every size and do not obey (6). For q>=3, targets c<=1 are rejected.

**Exact inverse theorem.** To decide whether a prescribed c is realized by a
positive transitive inflation of this fixed strong quotient, it suffices to
check the finite box

    1 <= r_i <= floor(log_2(c-1)).                    (7)

Use integer bit lengths rather than approximate logarithms in code. Formula
(1) supplies the exact count for each candidate, and the actual inflated
matrix supplies the witness. This is a complete target decision procedure
for this family, not just a search with an arbitrary size cap.

It also gives an explicit capacity obstruction: at most
floor(log_2(X-1))^q parameter vectors, hence at most that many distinct H<=X,
come from this quotient. A fixed finite bank of such fixed-order transitive
inflation quotients has density zero in the odd targets. Thus a universal
Chan–Pak proof via these ingredients must permit growing quotient complexity,
nontransitive children, or nonuniform interfaces. No reweighting of the exact
aggregate B_n can remove that missing ingredient. This makes the next task
precise rather than mistaking a recurrence for a target selector.

### A productive face: exact Mersenne factorization, with witnesses

The C3 face with one singleton module specializes (1) to

    H(C3[TT_a,TT_b,1]) = 2 + (2^a-1)(2^b-1), a,b>=1.  (8)

This is a derived explicit specialization, not a claim that cyclic blow-up
families or the underlying identity are new. To see it directly, extract one
copy of the singleton variable from the rational run series. The remaining
nonzero terms have run numbers (1,1), (2,1), (1,2), (2,2), with multiplicities
3,1,1,1. Their contribution is
3+(2^a-2)+(2^b-2)+(2^a-2)(2^b-2), which simplifies to (8).
Thus this face has a terminating exact arithmetic inverse: c-2 must factor
as two positive Mersenne numbers. No approximate count evaluation is needed.
For example 23=2+3*7 gives module sizes(2,3,1), and
47=2+3*15 gives(2,4,1). The finite box search classifies every transitive C3
inflation, not only this face: c=25,27,29,31 have no such inflation at *any*
size, whereas 45 has sizes(2,2,2). These are family exclusions, not exclusions
from the tournament spectrum. Actual witness matrices and exact exhaustive
box checks are included in `TARGETS.json`.

Ordinal sums multiply Hamiltonian counts, so actual matrices from this family
can immediately be composed into all finite products of these counts. An
implementation must retain factors as graph witnesses, not assume that an
integer factor of an attainable H is attainable (7 divides 49). This yields
constructive source material and a complete inverse for one soluble boundary,
while honestly identifying where a richer source interface is necessary.

The coefficient compiler composes with Session268's response calculus by
replacing u_i^k with each child's actual ordered k-path-cover count. Transitive
modules are the soluble Stirling boundary of that larger theory. The promising
constructive extension is to retain source interfaces in this replacement,
not to discard them at the aggregate stage.

## Scope and reproducibility

Run `python3 recurrence_bridge.py` in this directory. It checks the multivariate
compiler against independent Hamiltonian subset DP for all named quotients of
orders1–4 and module sizes0–2 with total order<=7; checks the C3 aggregate
through total order7; and produces directed pole certificates at n=10,20,30,
50,100,250 for regular orders3,5,7. The JSON includes exact integers, interval
endpoints, absolute analytic tail bounds and signed leading-pole error estimates.
Only the signed display estimates are rounded numerically; certificate intervals
are directed. A000568, individual H-value coverage, and all-class enumeration
are deliberately not inferred from this weighted aggregate.
